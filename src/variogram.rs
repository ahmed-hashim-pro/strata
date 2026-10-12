//! Semivariograms: how fast a quantity stops resembling itself with distance.
//!
//! The variogram is the whole of the uncertainty story. Kriging is exact given
//! one; the stated intervals are only as honest as the fit. So the fit is
//! deterministic, weighted toward short lags (which drive the predictions
//! between nearby boreholes), and every fitted model is printed by
//! `strata validate` so a reader can disagree with it.

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Exponential,
    Spherical,
}

impl fmt::Display for Kind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Kind::Exponential => "exponential",
            Kind::Spherical => "spherical",
        })
    }
}

/// gamma(h) = nugget + sill * shape(h / range) for h > 0, and 0 at h = 0.
/// `sill` is the partial sill; `range` is the practical range (where the
/// model reaches 95% of its sill for the exponential, 100% for the spherical).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Variogram {
    pub kind: Kind,
    pub nugget: f64,
    pub sill: f64,
    pub range: f64,
}

impl Variogram {
    pub fn gamma(&self, h: f64) -> f64 {
        if h <= 0.0 {
            0.0
        } else {
            self.nugget + self.sill * shape(self.kind, h / self.range)
        }
    }

    pub fn total_sill(&self) -> f64 {
        self.nugget + self.sill
    }
}

impl fmt::Display for Variogram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}, nugget {:.3}, sill {:.3}, range {:.0} m", self.kind, self.nugget, self.sill, self.range)
    }
}

fn shape(kind: Kind, r: f64) -> f64 {
    match kind {
        Kind::Exponential => 1.0 - (-3.0 * r).exp(),
        Kind::Spherical => {
            if r >= 1.0 {
                1.0
            } else {
                1.5 * r - 0.5 * r * r * r
            }
        }
    }
}

/// One lag bin of the experimental variogram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bin {
    pub lag: f64,
    pub gamma: f64,
    pub pairs: usize,
}

pub const BINS: usize = 12;
const MIN_PAIRS: usize = 3;

/// Classical (Matheron) estimator over `bins` equal-width bins out to half
/// the largest pair distance. Bins with fewer than three pairs are dropped.
pub fn experimental(points: &[[f64; 2]], values: &[f64], bins: usize) -> Vec<Bin> {
    let max_lag = 0.5 * max_distance(points);
    if max_lag <= 0.0 {
        return Vec::new();
    }
    let width = max_lag / bins as f64;
    let mut sum_h = vec![0.0; bins];
    let mut sum_g = vec![0.0; bins];
    let mut count = vec![0usize; bins];
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            let h = dist(points[i], points[j]);
            if h > max_lag {
                continue;
            }
            let b = ((h / width) as usize).min(bins - 1);
            sum_h[b] += h;
            sum_g[b] += 0.5 * (values[i] - values[j]).powi(2);
            count[b] += 1;
        }
    }
    (0..bins)
        .filter(|&b| count[b] >= MIN_PAIRS)
        .map(|b| Bin { lag: sum_h[b] / count[b] as f64, gamma: sum_g[b] / count[b] as f64, pairs: count[b] })
        .collect()
}

/// Fits nugget, sill and range for both model kinds by weighted least squares
/// (weights N/h^2, as gstat's default) over a fixed geometric grid of ranges,
/// solving nugget and sill in closed form for each range. Deterministic: the
/// same data always gives the same model.
pub fn fit(points: &[[f64; 2]], values: &[f64]) -> Variogram {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0).max(1.0);
    let max_lag = 0.5 * max_distance(points);
    // A floor keeps the kriging system regular when every value is identical
    // (a unit that is absent everywhere): the model then says "almost no
    // variation", which is what the data says.
    let floor = 1e-6 * (1.0 + var);

    let bins = experimental(points, values, BINS);
    if bins.len() < 3 {
        return Variogram { kind: Kind::Exponential, nugget: 0.0, sill: var.max(floor), range: max_lag.max(1.0) };
    }

    let mut best: Option<(f64, Variogram)> = None;
    const STEPS: usize = 60;
    for kind in [Kind::Exponential, Kind::Spherical] {
        for s in 0..STEPS {
            let range = max_lag * 0.05 * 40f64.powf(s as f64 / (STEPS - 1) as f64);
            let (nugget, sill, sse) = fit_linear(&bins, |h| shape(kind, h / range));
            if best.as_ref().map_or(true, |(b, _)| sse < *b) {
                best = Some((sse, Variogram { kind, nugget, sill, range }));
            }
        }
    }
    let mut v = best.expect("grid is not empty").1;
    if v.total_sill() < floor {
        v.sill = floor - v.nugget;
    }
    v
}

/// gamma = nugget + sill * f(h), both non-negative, weighted by N/h^2.
fn fit_linear(bins: &[Bin], f: impl Fn(f64) -> f64) -> (f64, f64, f64) {
    let (mut sw, mut sf, mut sff, mut sg, mut sfg) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for b in bins {
        let w = b.pairs as f64 / (b.lag * b.lag);
        let fv = f(b.lag);
        sw += w;
        sf += w * fv;
        sff += w * fv * fv;
        sg += w * b.gamma;
        sfg += w * fv * b.gamma;
    }
    let det = sw * sff - sf * sf;
    let (mut c0, mut c1) = if det.abs() > 1e-12 * sw * sff {
        ((sg * sff - sf * sfg) / det, (sw * sfg - sf * sg) / det)
    } else {
        (sg / sw, 0.0)
    };
    if c0 < 0.0 {
        c0 = 0.0;
        c1 = if sff > 0.0 { sfg / sff } else { 0.0 };
    }
    if c1 < 0.0 {
        c1 = 0.0;
        c0 = sg / sw;
    }
    let sse = bins.iter().map(|b| b.pairs as f64 / (b.lag * b.lag) * (b.gamma - c0 - c1 * f(b.lag)).powi(2)).sum();
    (c0.max(0.0), c1.max(0.0), sse)
}

pub fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

fn max_distance(points: &[[f64; 2]]) -> f64 {
    let mut m = 0.0f64;
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            m = m.max(dist(points[i], points[j]));
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gamma_is_zero_at_zero_and_reaches_sill() {
        for kind in [Kind::Exponential, Kind::Spherical] {
            let v = Variogram { kind, nugget: 0.5, sill: 2.0, range: 100.0 };
            assert_eq!(v.gamma(0.0), 0.0);
            assert!(v.gamma(1e-9) >= 0.5);
            assert!(v.gamma(100.0) >= 0.5 + 2.0 * 0.95 - 1e-12);
            assert!((v.gamma(1e6) - 2.5).abs() < 1e-9);
        }
    }

    #[test]
    fn gamma_is_monotone() {
        for kind in [Kind::Exponential, Kind::Spherical] {
            let v = Variogram { kind, nugget: 0.1, sill: 1.0, range: 50.0 };
            let mut prev = 0.0;
            for i in 0..500 {
                let g = v.gamma(i as f64 * 0.5);
                assert!(g >= prev);
                prev = g;
            }
        }
    }

    #[test]
    fn experimental_on_a_line() {
        // z = x on a regular line: gamma(h) = h^2 / 2 exactly.
        let pts: Vec<[f64; 2]> = (0..41).map(|i| [i as f64, 0.0]).collect();
        let vals: Vec<f64> = (0..41).map(|i| i as f64).collect();
        for b in experimental(&pts, &vals, 10) {
            let expect: f64 = 0.5 * b.lag * b.lag;
            assert!((b.gamma - expect).abs() / expect < 0.05, "{b:?}");
        }
    }

    #[test]
    fn constant_data_fits_a_regular_model() {
        let pts: Vec<[f64; 2]> = (0..20).map(|i| [(i % 5) as f64 * 10.0, (i / 5) as f64 * 10.0]).collect();
        let v = fit(&pts, &[3.0; 20]);
        assert!(v.total_sill() > 0.0);
    }
}
