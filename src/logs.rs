//! Borehole logs: the one input format, parsed strictly.
//!
//! ```text
//! # comments and blank lines are ignored
//! units: topsoil, clay, sand, gravel, chalk
//! hole,x,y,collar,unit,from,to
//! BH01,120.0,340.0,52.30,topsoil,0.0,0.8
//! BH01,120.0,340.0,52.30,clay,0.8,6.2
//! ```
//!
//! `units` lists the stratigraphic column from top to bottom; the last unit is
//! the basement, which has no base. `from`/`to` are depths below the collar
//! (ground elevation at the hole). A unit missing from a hole between two units
//! it did log is pinched out there (thickness 0). The unit a hole stops in has
//! an unknown base: only that it lies deeper than the hole.
//!
//! Anything ambiguous is an error with a line number, never a guess: a typo in
//! a unit name, a gap or overlap between intervals, a unit logged twice or out
//! of order, a hole whose location changes between rows.

use std::collections::HashMap;
use std::fmt::Write as _;

pub const HEADER: &str = "hole,x,y,collar,unit,from,to";
/// Holes closer than this are treated as the same location.
pub const MIN_SEPARATION_M: f64 = 0.01;
const DEPTH_TOLERANCE_M: f64 = 1e-6;

#[derive(Debug, Clone, PartialEq)]
pub struct Interval {
    pub unit: usize,
    pub from: f64,
    pub to: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Hole {
    pub id: String,
    pub x: f64,
    pub y: f64,
    pub collar: f64,
    pub intervals: Vec<Interval>,
}

impl Hole {
    pub fn xy(&self) -> [f64; 2] {
        [self.x, self.y]
    }

    /// Total depth drilled.
    pub fn depth(&self) -> f64 {
        self.intervals.last().map_or(0.0, |i| i.to)
    }

    /// Index of the unit the hole stopped in.
    pub fn last_unit(&self) -> usize {
        self.intervals.last().map_or(0, |i| i.unit)
    }

    /// Thickness of `unit` if this hole determines it: logged, or pinched out
    /// above the unit the hole stopped in. `None` for the unit it stopped in
    /// and everything below.
    pub fn thickness(&self, unit: usize) -> Option<f64> {
        if unit >= self.last_unit() {
            return None;
        }
        Some(self.intervals.iter().find(|i| i.unit == unit).map_or(0.0, |i| i.to - i.from))
    }

    /// Elevation of the base of `unit`, if this hole determines it.
    pub fn base(&self, unit: usize) -> Option<f64> {
        if unit >= self.last_unit() {
            return None;
        }
        let depth = self.intervals.iter().take_while(|i| i.unit <= unit).last().map_or(0.0, |i| i.to);
        Some(self.collar - depth)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Site {
    pub units: Vec<String>,
    pub holes: Vec<Hole>,
}

impl Site {
    pub fn unit_index(&self, name: &str) -> Option<usize> {
        self.units.iter().position(|u| u == name)
    }

    pub fn without_hole(&self, index: usize) -> Site {
        let mut s = self.clone();
        s.holes.remove(index);
        s
    }

    pub fn to_csv(&self, comment: &str) -> String {
        let mut out = String::new();
        for line in comment.lines() {
            let _ = writeln!(out, "# {line}");
        }
        let _ = writeln!(out, "units: {}", self.units.join(", "));
        let _ = writeln!(out, "{HEADER}");
        for h in &self.holes {
            for i in &h.intervals {
                let _ = writeln!(
                    out,
                    "{},{:.1},{:.1},{:.2},{},{:.1},{:.1}",
                    h.id, h.x, h.y, h.collar, self.units[i.unit], i.from, i.to
                );
            }
        }
        out
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ParseError {
    pub line: usize,
    pub message: String,
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.line == 0 {
            write!(f, "{}", self.message)
        } else {
            write!(f, "line {}: {}", self.line, self.message)
        }
    }
}

impl std::error::Error for ParseError {}

fn err<T>(line: usize, message: impl Into<String>) -> Result<T, ParseError> {
    Err(ParseError { line, message: message.into() })
}

fn valid_name(s: &str) -> bool {
    !s.is_empty() && s.chars().all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-' || c == '.')
}

fn number(line: usize, field: &str, s: &str) -> Result<f64, ParseError> {
    match s.parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(v),
        _ => err(line, format!("{field} '{s}' is not a finite number")),
    }
}

pub fn parse(text: &str) -> Result<Site, ParseError> {
    let mut lines =
        text.lines().enumerate().map(|(i, l)| (i + 1, l.trim())).filter(|(_, l)| !l.is_empty() && !l.starts_with('#'));

    let Some((ln, first)) = lines.next() else {
        return err(0, "no data: expected a 'units:' line");
    };
    let Some(list) = first.strip_prefix("units:") else {
        return err(ln, "expected 'units: <top>, ..., <basement>' first");
    };
    let units: Vec<String> = list.split(',').map(|u| u.trim().to_string()).collect();
    for u in &units {
        if !valid_name(u) {
            return err(ln, format!("unit name '{u}' is not valid (letters, digits, _ - .)"));
        }
    }
    if units.len() < 2 {
        return err(ln, "need at least two units: one with a base, and the basement");
    }
    for (i, u) in units.iter().enumerate() {
        if units[..i].contains(u) {
            return err(ln, format!("unit '{u}' is listed twice"));
        }
    }
    let index: HashMap<&str, usize> = units.iter().enumerate().map(|(i, u)| (u.as_str(), i)).collect();

    match lines.next() {
        Some((_, h)) if h.replace(' ', "") == HEADER => {}
        Some((ln, h)) => return err(ln, format!("expected header '{HEADER}', got '{h}'")),
        None => return err(0, format!("missing header '{HEADER}'")),
    }

    let mut holes: Vec<Hole> = Vec::new();
    let mut seen: HashMap<String, usize> = HashMap::new(); // id -> first line
    for (ln, row) in lines {
        let f: Vec<&str> = row.split(',').map(str::trim).collect();
        if f.len() != 7 {
            return err(ln, format!("expected 7 fields ({HEADER}), got {}", f.len()));
        }
        let id = f[0];
        if !valid_name(id) {
            return err(ln, format!("hole id '{id}' is not valid (letters, digits, _ - .)"));
        }
        let (x, y, collar) = (number(ln, "x", f[1])?, number(ln, "y", f[2])?, number(ln, "collar", f[3])?);
        let Some(&unit) = index.get(f[4]) else {
            return err(ln, format!("unit '{}' is not in the units list ({})", f[4], units.join(", ")));
        };
        let (from, to) = (number(ln, "from", f[5])?, number(ln, "to", f[6])?);
        if from < 0.0 {
            return err(ln, format!("depth 'from' {from} is negative"));
        }
        if to <= from {
            return err(ln, format!("interval {from}..{to} is empty or inverted"));
        }

        let continues = holes.last().is_some_and(|h| h.id == id);
        if !continues {
            if let Some(first) = seen.get(id) {
                return err(
                    ln,
                    format!("hole {id} already ended (its rows start at line {first}); keep a hole's rows together"),
                );
            }
            if from.abs() > DEPTH_TOLERANCE_M {
                return err(ln, format!("hole {id} must start at depth 0, starts at {from}"));
            }
            seen.insert(id.to_string(), ln);
            holes.push(Hole { id: id.to_string(), x, y, collar, intervals: vec![Interval { unit, from, to }] });
            continue;
        }
        let h = holes.last_mut().expect("continues implies a hole");
        if (h.x, h.y, h.collar) != (x, y, collar) {
            return err(ln, format!("hole {id} changes location or collar between rows"));
        }
        let prev = h.intervals.last().expect("hole has an interval");
        if (from - prev.to).abs() > DEPTH_TOLERANCE_M {
            let what = if from > prev.to { "gap" } else { "overlap" };
            return err(ln, format!("hole {id}: {what} between {} and {from} m", prev.to));
        }
        if unit <= prev.unit {
            return err(
                ln,
                format!(
                    "hole {id}: {} below {} is out of stratigraphic order ({})",
                    units[unit],
                    units[prev.unit],
                    units.join(", ")
                ),
            );
        }
        h.intervals.push(Interval { unit, from: prev.to, to });
    }
    if holes.is_empty() {
        return err(0, "no holes");
    }
    for i in 0..holes.len() {
        for j in 0..i {
            let d = (holes[i].x - holes[j].x).hypot(holes[i].y - holes[j].y);
            if d < MIN_SEPARATION_M {
                return err(0, format!("holes {} and {} are at the same location", holes[j].id, holes[i].id));
            }
        }
    }
    Ok(Site { units, holes })
}

#[cfg(test)]
mod tests {
    use super::*;

    const GOOD: &str = "\
# a site
units: topsoil, clay, sand, chalk
hole,x,y,collar,unit,from,to
A,0,0,50,topsoil,0,1
A,0,0,50,clay,1,5
A,0,0,50,sand,5,9
A,0,0,50,chalk,9,10
B,100,0,48,topsoil,0,1
B,100,0,48,sand,1,4
";

    #[test]
    fn parses_and_derives() {
        let s = parse(GOOD).unwrap();
        assert_eq!(s.units, ["topsoil", "clay", "sand", "chalk"]);
        let (a, b) = (&s.holes[0], &s.holes[1]);
        assert_eq!(a.thickness(1), Some(4.0));
        assert_eq!(a.base(2), Some(41.0));
        assert_eq!(a.thickness(3), None, "basement has no base");
        // B skipped clay: pinched out, known. B stopped in sand: unknown.
        assert_eq!(b.thickness(1), Some(0.0));
        assert_eq!(b.base(1), Some(47.0));
        assert_eq!(b.thickness(2), None);
        assert_eq!(b.base(2), None);
    }

    #[test]
    fn round_trips_through_csv() {
        let s = parse(GOOD).unwrap();
        assert_eq!(parse(&s.to_csv("again")).unwrap(), s);
    }

    fn fails(text: &str, line: usize, needle: &str) {
        let e = parse(text).expect_err("should fail");
        assert_eq!(e.line, line, "{e}");
        assert!(e.message.contains(needle), "{e}");
    }

    #[test]
    fn rejects_what_it_cannot_be_sure_of() {
        let head = "units: a, b, c\nhole,x,y,collar,unit,from,to\n";
        fails("", 0, "no data");
        fails("hole,x,y\n", 1, "units:");
        fails("units: a\n", 1, "two units");
        fails("units: a, a\n", 1, "twice");
        fails("units: a, b c\n", 1, "not valid");
        fails("units: a, b\nhole,x,y\n", 2, "header");
        fails(&format!("{head}H,0,0,1,d,0,1\n"), 3, "not in the units list");
        fails(&format!("{head}H,0,0,1,a,0,1,9\n"), 3, "7 fields");
        fails(&format!("{head}H,0,nan,1,a,0,1\n"), 3, "finite");
        fails(&format!("{head}H,0,0,1,a,1,1\n"), 3, "empty");
        fails(&format!("{head}H,0,0,1,a,0.5,1\n"), 3, "start at depth 0");
        fails(&format!("{head}H,0,0,1,a,0,1\nH,0,0,1,b,1.5,2\n"), 4, "gap");
        fails(&format!("{head}H,0,0,1,a,0,1\nH,0,0,1,b,0.5,2\n"), 4, "overlap");
        fails(&format!("{head}H,0,0,1,b,0,1\nH,0,0,1,a,1,2\n"), 4, "stratigraphic order");
        fails(&format!("{head}H,0,0,1,a,0,1\nH,0,0,1,a,1,2\n"), 4, "stratigraphic order");
        fails(&format!("{head}H,0,0,1,a,0,1\nH,0,1,1,b,1,2\n"), 4, "location");
        fails(&format!("{head}H,0,0,1,a,0,1\nG,5,0,1,a,0,1\nH,0,0,1,b,1,2\n"), 5, "keep a hole's rows together");
        fails(&format!("{head}H,0,0,1,a,0,1\nG,0,0.001,1,a,0,1\n"), 0, "same location");
        fails(head, 0, "no holes");
    }
}
