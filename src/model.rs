//! The layered model: ground surface plus one thickness per unit.
//!
//! Kriging each contact elevation on its own lets layers cross, which is a
//! geological impossibility. So the model kriges the ground surface and each
//! unit's thickness (an isopach), clamps thicknesses at zero, and stacks them:
//! base of unit k = ground - sum of thicknesses down to k. Layers cannot cross
//! by construction, and a unit can pinch out.
//!
//! The variance of a base is the sum of the variances above it. That treats
//! the thicknesses as independent errors, which is an assumption, so the
//! intervals are checked, not trusted: `strata cv` against held-out holes and
//! `strata bench` against synthetic sites whose truth is known.

use crate::kriging::{Estimate, Kriging};
use crate::logs::Site;
use crate::variogram::{self, Variogram};

/// Fewer observations than this and a variogram fit is noise.
pub const MIN_OBSERVATIONS: usize = 6;
/// Two-sided 95% normal quantile.
pub const Z95: f64 = 1.959_963_985;
pub const DEFAULT_MAX_SIGMA: f64 = 1.5;
/// Bounds on the leave-one-out variance correction (see `variance_scale`).
pub const VARIANCE_SCALE_LIMITS: (f64, f64) = (0.25, 4.0);

/// One kriged variable and the holes behind it.
#[derive(Debug, Clone)]
pub struct Fitted {
    pub observations: usize,
    pub kriging: Kriging,
    /// Multiplier on the kriging variance, from the variable's own
    /// leave-one-out errors (see `variance_scale`). Clamped to
    /// `VARIANCE_SCALE_LIMITS`, so it can be below 1 as well as above.
    pub variance_scale: f64,
}

impl Fitted {
    pub fn variogram(&self) -> &Variogram {
        self.kriging.variogram()
    }

    pub fn predict(&self, p: [f64; 2]) -> Estimate {
        let e = self.kriging.predict(p);
        Estimate { value: e.value, variance: e.variance * self.variance_scale }
    }
}

/// Kriging variance takes the fitted variogram as exact, and a variogram
/// fitted to a few dozen holes is not. Rough surfaces fitted with a sill that
/// comes out low give intervals that are too narrow; smooth surfaces fitted
/// with an exponential or spherical model, which overstate short-range
/// variation, give intervals that are too wide. Either way the error is
/// measured, not assumed: krige each hole from the others under the same
/// variogram, and scale the variance by the mean squared standardised error.
/// The factor is clamped to [1/4, 4] so a handful of holes cannot make the
/// model wildly sure or wildly unsure.
fn variance_scale(points: &[[f64; 2]], values: &[f64], v: Variogram) -> Result<f64, String> {
    let mut sum = 0.0;
    for i in 0..points.len() {
        let (mut p, mut z) = (points.to_vec(), values.to_vec());
        p.remove(i);
        z.remove(i);
        let e = Kriging::new(p, z, v)?.predict(points[i]);
        if e.variance > 0.0 {
            sum += (e.value - values[i]).powi(2) / e.variance;
        }
    }
    Ok((sum / points.len() as f64).clamp(VARIANCE_SCALE_LIMITS.0, VARIANCE_SCALE_LIMITS.1))
}

#[derive(Debug, Clone)]
pub struct Model {
    pub units: Vec<String>,
    pub ground: Fitted,
    /// Thickness of every unit except the basement; `Err` says why a unit
    /// could not be modelled (and so no base at or below it can be given).
    pub thickness: Vec<Result<Fitted, String>>,
}

/// A surface: 0 is the ground, k >= 1 is the base of unit k - 1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Surface {
    Known(Estimate),
    /// The model has no basis for it (too few holes reached it).
    NoData,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Column {
    pub surfaces: Vec<Surface>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Answer {
    Known(Estimate),
    /// Estimated, but sigma exceeds the limit: the model will not vouch for it.
    TooUncertain(Estimate),
    NoData,
}

impl Column {
    /// The answer for a surface under a sigma limit. Once a surface is too
    /// uncertain, everything below it is too: errors only accumulate downward.
    pub fn answer(&self, surface: usize, max_sigma: f64) -> Answer {
        match self.surfaces[surface] {
            Surface::NoData => Answer::NoData,
            Surface::Known(e) if e.sigma() > max_sigma => Answer::TooUncertain(e),
            Surface::Known(e) => Answer::Known(e),
        }
    }

    /// How many surfaces from the top are answered under the limit.
    pub fn answered(&self, max_sigma: f64) -> usize {
        (0..self.surfaces.len()).take_while(|&s| matches!(self.answer(s, max_sigma), Answer::Known(_))).count()
    }
}

fn fit(points: Vec<[f64; 2]>, values: Vec<f64>) -> Result<Fitted, String> {
    let v = variogram::fit(&points, &values);
    let variance_scale = variance_scale(&points, &values, v)?;
    let observations = points.len();
    Ok(Fitted { observations, kriging: Kriging::new(points, values, v)?, variance_scale })
}

impl Model {
    pub fn build(site: &Site) -> Result<Model, String> {
        if site.holes.len() < MIN_OBSERVATIONS {
            return Err(format!("{} holes; the ground surface needs at least {MIN_OBSERVATIONS}", site.holes.len()));
        }
        let ground = fit(site.holes.iter().map(|h| h.xy()).collect(), site.holes.iter().map(|h| h.collar).collect())?;
        let mut thickness = Vec::new();
        for unit in 0..site.units.len() - 1 {
            let (pts, vals): (Vec<_>, Vec<_>) =
                site.holes.iter().filter_map(|h| h.thickness(unit).map(|t| (h.xy(), t))).unzip();
            thickness.push(if pts.len() < MIN_OBSERVATIONS {
                Err(format!("only {} holes got through {}; need {MIN_OBSERVATIONS}", pts.len(), site.units[unit]))
            } else {
                fit(pts, vals)
            });
        }
        Ok(Model { units: site.units.clone(), ground, thickness })
    }

    pub fn surface_count(&self) -> usize {
        self.units.len()
    }

    pub fn surface_name(&self, surface: usize) -> String {
        if surface == 0 {
            "ground surface".into()
        } else {
            format!("base of {}", self.units[surface - 1])
        }
    }

    /// Parses "ground" or a unit name (meaning its base).
    pub fn surface_by_name(&self, name: &str) -> Option<usize> {
        if name == "ground" {
            return Some(0);
        }
        self.units[..self.units.len() - 1].iter().position(|u| u == name).map(|i| i + 1)
    }

    pub fn column(&self, p: [f64; 2]) -> Column {
        let g = self.ground.predict(p);
        let mut surfaces = vec![Surface::Known(g)];
        let (mut elevation, mut variance, mut ok) = (g.value, g.variance, true);
        for t in &self.thickness {
            match (ok, t) {
                (true, Ok(f)) => {
                    let e = f.predict(p);
                    elevation -= e.value.max(0.0);
                    variance += e.variance;
                    surfaces.push(Surface::Known(Estimate { value: elevation, variance }));
                }
                _ => {
                    ok = false;
                    surfaces.push(Surface::NoData);
                }
            }
        }
        Column { surfaces }
    }

    /// Holes that stopped inside a unit say its base is deeper than the hole.
    /// The model never uses that for fitting; here it is checked: a hole whose
    /// bottom lies more than 1.96 sigma below the predicted base contradicts
    /// the model there.
    pub fn contradictions(&self, site: &Site) -> Vec<String> {
        let mut out = Vec::new();
        for h in &site.holes {
            let unit = h.last_unit();
            if unit + 1 >= self.units.len() {
                continue; // reached the basement
            }
            let bottom = h.collar - h.depth();
            if let Surface::Known(e) = self.column(h.xy()).surfaces[unit + 1] {
                if bottom < e.value - Z95 * e.sigma() {
                    out.push(format!(
                        "hole {} stopped in {} at {:.2} m, but the model puts the base of {} at {:.2} ± {:.2} m",
                        h.id,
                        self.units[unit],
                        bottom,
                        self.units[unit],
                        e.value,
                        Z95 * e.sigma()
                    ));
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::logs::parse;

    fn site() -> Site {
        let mut text = String::from("units: soil, clay, rock\nhole,x,y,collar,unit,from,to\n");
        for i in 0..8 {
            let (x, y) = ((i % 4) as f64 * 100.0, (i / 4) as f64 * 100.0);
            let c = 50.0 + i as f64 * 0.3;
            let clay = 4.0 + (i % 3) as f64;
            text += &format!("H{i},{x},{y},{c},soil,0,1\nH{i},{x},{y},{c},clay,1,{}\n", 1.0 + clay);
            text += &format!("H{i},{x},{y},{c},rock,{},{}\n", 1.0 + clay, 3.0 + clay);
        }
        parse(&text).unwrap()
    }

    #[test]
    fn honours_the_holes() {
        let s = site();
        let m = Model::build(&s).unwrap();
        for h in &s.holes {
            let c = m.column(h.xy());
            for k in 0..2 {
                let Surface::Known(e) = c.surfaces[k + 1] else { panic!() };
                assert!((e.value - h.base(k).unwrap()).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn too_few_holes_through_a_unit_means_no_data_below() {
        let mut s = site();
        // Stop most holes in the clay: its base is unknown at all but two.
        for h in s.holes.iter_mut().skip(2) {
            h.intervals.truncate(2);
        }
        let m = Model::build(&s).unwrap();
        assert!(m.thickness[0].is_ok());
        assert!(m.thickness[1].is_err());
        let c = m.column([150.0, 50.0]);
        assert!(matches!(c.surfaces[1], Surface::Known(_)));
        assert_eq!(c.surfaces[2], Surface::NoData);
        assert_eq!(c.answered(f64::INFINITY), 2);
    }

    #[test]
    fn refusal_cascades_downward() {
        let m = Model::build(&site()).unwrap();
        let c = m.column([5000.0, 5000.0]);
        let n = c.answered(0.0);
        assert_eq!(n, 0);
        assert!(matches!(c.answer(2, 0.0), Answer::TooUncertain(_)));
        assert_eq!(c.answered(f64::INFINITY), 3);
    }

    #[test]
    fn surface_names() {
        let m = Model::build(&site()).unwrap();
        assert_eq!(m.surface_by_name("ground"), Some(0));
        assert_eq!(m.surface_by_name("clay"), Some(2));
        assert_eq!(m.surface_by_name("rock"), None, "the basement has no base");
        assert_eq!(m.surface_name(2), "base of clay");
    }
}
