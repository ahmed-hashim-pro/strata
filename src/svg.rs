//! Pictures: a vertical cross-section and a plan map of where the model will
//! answer. Both draw only what the model vouches for; everything else is
//! hatched as refused, so the picture cannot look more certain than the model.

use std::fmt::Write as _;

use crate::logs::{Hole, Site};
use crate::model::{Column, Model, Surface, Z95};
use crate::synth::Truth;

/// Unit fills, in stratigraphic order (a fixed categorical order: adjacent
/// layers are adjacent pairs, which this order keeps distinguishable).
const UNIT_COLOURS: [&str; 8] =
    ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"];
/// Sequential ramp for sigma, light (sure) to dark (unsure).
const SIGMA_RAMP: [&str; 5] = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"];
const INK: &str = "#24292f";
const MUTED: &str = "#57606a";
const GRID: &str = "#d8dee4";
const FONT: &str = "-apple-system, 'Segoe UI', Helvetica, Arial, sans-serif";

fn n(v: f64) -> String {
    format!("{v:.1}")
}

fn unit_colour(u: usize) -> &'static str {
    UNIT_COLOURS[u % UNIT_COLOURS.len()]
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// A tick step of 1, 2 or 5 times a power of ten giving about `target` ticks.
fn nice_step(span: f64, target: f64) -> f64 {
    let raw = span / target;
    let p = 10f64.powf(raw.log10().floor());
    [1.0, 2.0, 5.0, 10.0].iter().map(|m| m * p).find(|s| *s >= raw).unwrap_or(10.0 * p)
}

fn hatch_def(out: &mut String) {
    let _ = writeln!(
        out,
        "<defs><pattern id=\"refused\" width=\"7\" height=\"7\" patternUnits=\"userSpaceOnUse\" patternTransform=\"rotate(45)\">\
<rect width=\"7\" height=\"7\" fill=\"#f6f8fa\"/><line x1=\"0\" y1=\"0\" x2=\"0\" y2=\"7\" stroke=\"#8c959f\" stroke-width=\"1.6\"/></pattern></defs>"
    );
}

fn header(out: &mut String, w: f64, h: f64, title: &str, subtitle: &str) {
    let _ = writeln!(
        out,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\" font-family=\"{FONT}\">",
        n(w),
        n(h),
        n(w),
        n(h)
    );
    hatch_def(out);
    let _ = writeln!(out, "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>");
    let _ = writeln!(
        out,
        "<text x=\"24\" y=\"30\" font-size=\"17\" font-weight=\"600\" fill=\"{INK}\">{}</text>",
        esc(title)
    );
    let _ = writeln!(out, "<text x=\"24\" y=\"50\" font-size=\"12.5\" fill=\"{MUTED}\">{}</text>", esc(subtitle));
}

/// Legend entries laid out left to right, wrapping at `right`.
struct Legend {
    x: f64,
    y: f64,
    left: f64,
    right: f64,
}

impl Legend {
    fn item(&mut self, out: &mut String, swatch: &str, label: &str) {
        let width = 30.0 + 6.6 * label.len() as f64 + 18.0;
        if self.x + width > self.right {
            self.x = self.left;
            self.y += 22.0;
        }
        let _ = writeln!(out, "<g transform=\"translate({},{})\">{swatch}</g>", n(self.x), n(self.y - 6.0));
        let _ = writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"12\" fill=\"{MUTED}\">{}</text>",
            n(self.x + 30.0),
            n(self.y + 4.0),
            esc(label)
        );
        self.x += width;
    }
}

pub struct SectionOptions<'a> {
    pub from: [f64; 2],
    pub to: [f64; 2],
    /// Holes within this distance of the line are projected onto it.
    pub buffer: f64,
    pub max_sigma: f64,
    /// Synthetic truth to overlay, if the site came from `strata synth`.
    pub truth: Option<&'a Truth>,
}

const STATIONS: usize = 240;

/// Consecutive index runs where `keep` holds.
fn runs(len: usize, keep: impl Fn(usize) -> bool) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut start = None;
    for i in 0..=len {
        let k = i < len && keep(i);
        match (k, start) {
            (true, None) => start = Some(i),
            (false, Some(s)) => {
                out.push((s, i - 1));
                start = None;
            }
            _ => {}
        }
    }
    out
}

fn est(c: &Column, s: usize) -> Option<(f64, f64)> {
    match c.surfaces[s] {
        Surface::Known(e) => Some((e.value, e.sigma())),
        Surface::NoData => None,
    }
}

pub fn section(site: &Site, model: &Model, o: &SectionOptions) -> String {
    let (ax, ay) = (o.from[0], o.from[1]);
    let (dx, dy) = (o.to[0] - ax, o.to[1] - ay);
    let len = dx.hypot(dy);
    let at = |t: f64| [ax + dx * t / len, ay + dy * t / len];
    let ts: Vec<f64> = (0..STATIONS).map(|i| len * i as f64 / (STATIONS - 1) as f64).collect();
    let cols: Vec<Column> = ts.iter().map(|&t| model.column(at(t))).collect();
    let answered: Vec<usize> = cols.iter().map(|c| c.answered(o.max_sigma)).collect();
    let surfaces = model.surface_count();
    let truth: Option<Vec<Vec<f64>>> = o.truth.map(|tr| ts.iter().map(|&t| tr.surfaces(at(t))).collect());

    // Holes near the line, by distance along it.
    let mut near: Vec<(f64, &Hole)> = site
        .holes
        .iter()
        .filter_map(|h| {
            let (hx, hy) = (h.x - ax, h.y - ay);
            let t = (hx * dx + hy * dy) / len;
            let off = (hx * dy - hy * dx).abs() / len;
            (off <= o.buffer && (0.0..=len).contains(&t)).then_some((t, h))
        })
        .collect();
    near.sort_by(|a, b| a.0.total_cmp(&b.0));

    // Vertical extent: every estimate with its interval, the holes, the truth.
    let (mut top, mut bottom) = (f64::NEG_INFINITY, f64::INFINITY);
    for c in &cols {
        for s in 0..surfaces {
            if let Some((v, sg)) = est(c, s) {
                top = top.max(v + Z95 * sg);
                bottom = bottom.min(v - Z95 * sg);
            }
        }
    }
    for (_, h) in &near {
        top = top.max(h.collar);
        bottom = bottom.min(h.collar - h.depth());
    }
    if let Some(tr) = &truth {
        for col in tr {
            top = top.max(col[0]);
            bottom = bottom.min(*col.last().unwrap());
        }
    }
    let pad = 0.06 * (top - bottom);
    top += pad + 2.0;
    bottom -= pad;

    let (left, right, plot_top, plot_h, w) = (64.0, 24.0, 96.0, 400.0, 980.0);
    let plot_w = w - left - right;
    let px = |t: f64| left + t / len * plot_w;
    let py = |z: f64| plot_top + (top - z) / (top - bottom) * plot_h;
    let ve = (plot_h / (top - bottom)) / (plot_w / len);

    let mut out = String::new();
    let h_total = plot_top + plot_h + 120.0;
    header(
        &mut out,
        w,
        h_total,
        &format!("Section A–B through {} holes", near.len()),
        &format!(
            "Kriged layers with 95% intervals. Hatched where the model refuses (σ > {} m). Holes within {} m of the line are projected onto it.",
            o.max_sigma, o.buffer
        ),
    );

    // Grid and axes.
    let zstep = nice_step(top - bottom, 6.0);
    let mut z = (bottom / zstep).ceil() * zstep;
    while z <= top {
        let y = py(z);
        let _ = writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{GRID}\" stroke-width=\"1\"/>",
            n(left),
            n(y),
            n(left + plot_w),
            n(y)
        );
        let _ = writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"11\" text-anchor=\"end\" fill=\"{MUTED}\">{z:.0} m</text>",
            n(left - 8.0),
            n(y + 4.0)
        );
        z += zstep;
    }
    let tstep = nice_step(len, 8.0);
    let mut t = 0.0;
    while t <= len + 1e-9 {
        let _ = writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"11\" text-anchor=\"middle\" fill=\"{MUTED}\">{t:.0}</text>",
            n(px(t)),
            n(plot_top + plot_h + 18.0)
        );
        t += tstep;
    }
    let _ = writeln!(
        out,
        "<text x=\"{}\" y=\"{}\" font-size=\"11\" text-anchor=\"end\" fill=\"{MUTED}\">distance along section (m), vertical exaggeration ×{:.0}</text>",
        n(left + plot_w),
        n(plot_top + plot_h + 36.0),
        ve
    );

    // Units, between their top and base, where both are answered.
    // Runs reach half a station past their ends so neighbouring runs meet.
    let half = 0.5 * plot_w / (STATIONS - 1) as f64;
    let poly = |out: &mut String,
                a: usize,
                b: usize,
                upper: &dyn Fn(usize) -> f64,
                lower: &dyn Fn(usize) -> f64,
                attrs: &str| {
        let (xa, xb) = ((px(ts[a]) - half).max(left), (px(ts[b]) + half).min(left + plot_w));
        let mut pts = format!("{},{} ", n(xa), n(py(upper(a))));
        for i in a..=b {
            let _ = write!(pts, "{},{} ", n(px(ts[i])), n(py(upper(i))));
        }
        let _ = write!(pts, "{},{} {},{} ", n(xb), n(py(upper(b))), n(xb), n(py(lower(b))));
        for i in (a..=b).rev() {
            let _ = write!(pts, "{},{} ", n(px(ts[i])), n(py(lower(i))));
        }
        let _ = write!(pts, "{},{}", n(xa), n(py(lower(a))));
        let _ = writeln!(out, "<polygon points=\"{}\" {attrs}/>", pts.trim_end());
    };
    let value = |i: usize, s: usize| est(&cols[i], s).map_or(bottom, |e| e.0);
    for u in 0..model.units.len() {
        let basement = u + 1 == surfaces;
        let keep = |i: usize| if basement { answered[i] == surfaces } else { answered[i] > u + 1 };
        for (a, b) in runs(STATIONS, keep) {
            let lower = |i: usize| if basement { bottom } else { value(i, u + 1) };
            poly(
                &mut out,
                a,
                b,
                &|i| value(i, u),
                &lower,
                &format!("fill=\"{}\" fill-opacity=\"0.42\"", unit_colour(u)),
            );
        }
    }
    // Refused: from the deepest answered surface (or the top) down.
    for (a, b) in runs(STATIONS, |i| answered[i] < surfaces) {
        let upper = |i: usize| if answered[i] == 0 { top } else { value(i, answered[i] - 1) };
        poly(&mut out, a, b, &upper, &|_| bottom, "fill=\"url(#refused)\"");
    }
    // Intervals and estimates for each answered surface.
    for s in 0..surfaces {
        for (a, b) in runs(STATIONS, |i| answered[i] > s) {
            let e = |i: usize| est(&cols[i], s).expect("answered implies estimated");
            poly(
                &mut out,
                a,
                b,
                &|i| e(i).0 + Z95 * e(i).1,
                &|i| e(i).0 - Z95 * e(i).1,
                &format!("fill=\"{INK}\" fill-opacity=\"0.13\""),
            );
            let mut d = String::new();
            for i in a..=b {
                let _ = write!(d, "{}{},{} ", if i == a { "M" } else { "L" }, n(px(ts[i])), n(py(e(i).0)));
            }
            let _ = writeln!(out, "<path d=\"{}\" fill=\"none\" stroke=\"{INK}\" stroke-width=\"1.3\"/>", d.trim_end());
        }
    }
    if let Some(tr) = &truth {
        for s in 0..surfaces {
            let mut d = String::new();
            for (i, col) in tr.iter().enumerate() {
                let _ = write!(d, "{}{},{} ", if i == 0 { "M" } else { "L" }, n(px(ts[i])), n(py(col[s])));
            }
            let _ = writeln!(
                out,
                "<path d=\"{}\" fill=\"none\" stroke=\"#cf222e\" stroke-width=\"1.4\" stroke-dasharray=\"5 3\"/>",
                d.trim_end()
            );
        }
    }

    // Holes as logged.
    for (t, h) in &near {
        let x = px(*t);
        for iv in &h.intervals {
            let _ = writeln!(
                out,
                "<rect x=\"{}\" y=\"{}\" width=\"9\" height=\"{}\" fill=\"{}\" stroke=\"{INK}\" stroke-width=\"0.8\"/>",
                n(x - 4.5),
                n(py(h.collar - iv.from)),
                n((py(h.collar - iv.to) - py(h.collar - iv.from)).max(0.5)),
                unit_colour(iv.unit)
            );
        }
        let _ = writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10.5\" text-anchor=\"middle\" fill=\"{INK}\">{}</text>",
            n(x),
            n(py(h.collar) - 6.0),
            esc(&h.id)
        );
    }
    let _ = writeln!(
        out,
        "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"none\" stroke=\"{INK}\" stroke-width=\"1\"/>",
        n(left),
        n(plot_top),
        n(plot_w),
        n(plot_h)
    );
    for (label, t, anchor) in [("A", 0.0, "start"), ("B", len, "end")] {
        let p = at(t);
        let _ = writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"12\" font-weight=\"600\" text-anchor=\"{anchor}\" fill=\"{INK}\">{label} ({:.0}, {:.0})</text>",
            n(px(t)),
            n(plot_top - 8.0),
            p[0],
            p[1]
        );
    }

    let mut lg = Legend { x: left, y: plot_top + plot_h + 64.0, left, right: w - right };
    for (u, name) in model.units.iter().enumerate() {
        lg.item(
            &mut out,
            &format!("<rect width=\"22\" height=\"12\" fill=\"{}\" fill-opacity=\"0.6\"/>", unit_colour(u)),
            name,
        );
    }
    lg.item(&mut out, &format!("<rect width=\"22\" height=\"12\" fill=\"{INK}\" fill-opacity=\"0.18\"/><line x1=\"0\" y1=\"6\" x2=\"22\" y2=\"6\" stroke=\"{INK}\" stroke-width=\"1.3\"/>"), "estimate and 95% interval");
    lg.item(
        &mut out,
        "<rect width=\"22\" height=\"12\" fill=\"url(#refused)\" stroke=\"#8c959f\" stroke-width=\"0.6\"/>",
        "refused",
    );
    if truth.is_some() {
        lg.item(&mut out, "<line x1=\"0\" y1=\"6\" x2=\"22\" y2=\"6\" stroke=\"#cf222e\" stroke-width=\"1.4\" stroke-dasharray=\"5 3\"/>", "true contacts (synthetic site)");
    }
    out.push_str("</svg>\n");
    out
}

pub struct MapOptions {
    /// Surface index (0 = ground, k = base of unit k - 1).
    pub surface: usize,
    /// Grid cell size in metres; default about 1/70 of the extent.
    pub cell: Option<f64>,
    pub max_sigma: f64,
    pub line: Option<([f64; 2], [f64; 2])>,
}

pub fn map(site: &Site, model: &Model, o: &MapOptions) -> String {
    let (mut x0, mut x1, mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for h in &site.holes {
        x0 = x0.min(h.x);
        x1 = x1.max(h.x);
        y0 = y0.min(h.y);
        y1 = y1.max(h.y);
    }
    let margin = 0.08 * (x1 - x0).max(y1 - y0).max(1.0);
    let (x0, x1, y0, y1) = (x0 - margin, x1 + margin, y0 - margin, y1 + margin);
    let cell = o.cell.unwrap_or((x1 - x0).max(y1 - y0) / 70.0);
    let (nx, ny) = (((x1 - x0) / cell).ceil() as usize, ((y1 - y0) / cell).ceil() as usize);

    let plot_w = 620.0;
    let scale = plot_w / (nx as f64 * cell);
    let plot_h = ny as f64 * cell * scale;
    let (left, plot_top) = (40.0, 76.0);
    let w = left + plot_w + 40.0;
    let px = |x: f64| left + (x - x0) * scale;
    let py = |y: f64| plot_top + (y0 + ny as f64 * cell - y) * scale;

    let name = model.surface_name(o.surface);
    let mut out = String::new();
    header(
        &mut out,
        w,
        plot_top + plot_h + 96.0,
        &format!("Where the model answers: {name}"),
        &format!(
            "Standard error of the {name} between the holes; hatched where σ > {} m and the model refuses.",
            o.max_sigma
        ),
    );

    let bin = |sigma: f64| (((sigma / o.max_sigma) * SIGMA_RAMP.len() as f64) as usize).min(SIGMA_RAMP.len() - 1);
    let (mut refused, mut total) = (0usize, 0usize);
    let _ = writeln!(out, "<g shape-rendering=\"crispEdges\">");
    for j in 0..ny {
        let yc = y0 + (j as f64 + 0.5) * cell;
        let fills: Vec<&str> = (0..nx)
            .map(|i| match model.column([x0 + (i as f64 + 0.5) * cell, yc]).surfaces[o.surface] {
                Surface::Known(e) if e.sigma() <= o.max_sigma => SIGMA_RAMP[bin(e.sigma())],
                _ => "url(#refused)",
            })
            .collect();
        total += nx;
        refused += fills.iter().filter(|f| f.starts_with("url")).count();
        // One rect per run of equal cells; rows overlap by a hair so no seams show.
        let mut i = 0;
        while i < nx {
            let k = (i..nx).find(|&k| fills[k] != fills[i]).unwrap_or(nx);
            let _ = writeln!(
                out,
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\"/>",
                n(px(x0 + i as f64 * cell)),
                n(py(yc + cell / 2.0)),
                n((k - i) as f64 * cell * scale + 0.4),
                n(cell * scale + 0.4),
                fills[i]
            );
            i = k;
        }
    }
    let _ = writeln!(out, "</g>");
    if let Some((a, b)) = o.line {
        let _ = writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{INK}\" stroke-width=\"2\" stroke-dasharray=\"7 4\"/>",
            n(px(a[0])),
            n(py(a[1])),
            n(px(b[0])),
            n(py(b[1]))
        );
        for (label, p) in [("A", a), ("B", b)] {
            let _ = writeln!(
                out,
                "<text x=\"{}\" y=\"{}\" font-size=\"13\" font-weight=\"700\" text-anchor=\"middle\" fill=\"{INK}\" stroke=\"#ffffff\" stroke-width=\"3\" paint-order=\"stroke\">{label}</text>",
                n(px(p[0])),
                n(py(p[1]) - 9.0)
            );
        }
    }
    let unit = o.surface.checked_sub(1);
    for h in &site.holes {
        let determines = match unit {
            None => true,
            Some(u) => h.base(u).is_some(),
        };
        let fill = if determines { INK } else { "#ffffff" };
        let _ = writeln!(
            out,
            "<circle cx=\"{}\" cy=\"{}\" r=\"4\" fill=\"{fill}\" stroke=\"{INK}\" stroke-width=\"1.5\"/>",
            n(px(h.x)),
            n(py(h.y))
        );
    }
    let _ = writeln!(
        out,
        "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"none\" stroke=\"{INK}\" stroke-width=\"1\"/>",
        n(left),
        n(plot_top),
        n(plot_w),
        n(plot_h)
    );
    let _ = writeln!(
        out,
        "<text x=\"{}\" y=\"{}\" font-size=\"11\" text-anchor=\"end\" fill=\"{MUTED}\">{:.0} m × {:.0} m, north up · refused over {:.0}% of the area</text>",
        n(left + plot_w),
        n(plot_top - 8.0),
        x1 - x0,
        y1 - y0,
        100.0 * refused as f64 / total as f64
    );

    // Colour key: sigma bins and refused on one row, the hole symbols below.
    let ky = plot_top + plot_h + 30.0;
    let step = o.max_sigma / SIGMA_RAMP.len() as f64;
    let kw = 52.0;
    let _ =
        writeln!(out, "<text x=\"{}\" y=\"{}\" font-size=\"12\" fill=\"{MUTED}\">σ (m)</text>", n(left), n(ky + 4.0));
    for (i, c) in SIGMA_RAMP.iter().enumerate() {
        let x = left + 44.0 + i as f64 * kw;
        let _ = writeln!(out, "<rect x=\"{}\" y=\"{}\" width=\"{kw}\" height=\"12\" fill=\"{c}\"/>", n(x), n(ky - 6.0));
    }
    for i in 0..=SIGMA_RAMP.len() {
        let _ = writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10.5\" text-anchor=\"middle\" fill=\"{MUTED}\">{}</text>",
            n(left + 44.0 + i as f64 * kw),
            n(ky + 20.0),
            format!("{:.2}", i as f64 * step).trim_end_matches('0').trim_end_matches('.')
        );
    }
    let mut lg =
        Legend { x: left + 44.0 + SIGMA_RAMP.len() as f64 * kw + 28.0, y: ky, left: left + 44.0, right: w - 20.0 };
    lg.item(
        &mut out,
        "<rect width=\"22\" height=\"12\" fill=\"url(#refused)\" stroke=\"#8c959f\" stroke-width=\"0.6\"/>",
        "refused",
    );
    lg.y += 34.0;
    lg.x = left;
    lg.item(&mut out, &format!("<circle cx=\"11\" cy=\"6\" r=\"4\" fill=\"{INK}\"/>"), "hole reached this surface");
    lg.item(
        &mut out,
        &format!("<circle cx=\"11\" cy=\"6\" r=\"4\" fill=\"#ffffff\" stroke=\"{INK}\" stroke-width=\"1.5\"/>"),
        "hole stopped above it",
    );
    out.push_str("</svg>\n");
    out
}
