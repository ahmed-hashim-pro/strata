//! Synthetic sites with a known truth.
//!
//! Real borehole data never tells you what lies between the holes, so it can
//! only test the model where there is data (cross-validation). A synthetic
//! site can test it everywhere. Each surface is a smooth Gaussian random
//! field (Matérn 3/2 covariance), drawn with random Fourier features; boreholes
//! are sampled from it, logged to 0.1 m like a driller would, and some stop
//! before the basement. `bench` then compares the model's stated intervals
//! with the truth at points between the holes.

use crate::cv::SurfaceStats;
use crate::kriging::Estimate;
use crate::logs::{Hole, Interval, Site};
use crate::model::{Answer, Model};
use crate::rng::Rng;

/// The site is a square of this side, in metres.
pub const SIZE: f64 = 1000.0;
const FEATURES: usize = 400;
/// Minimum spacing between synthetic holes.
const SPACING_M: f64 = 30.0;

/// Covariance family of a synthetic field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Smoothness {
    /// Matérn 3/2: smooth, like most sedimentary contacts. The default.
    Smooth,
    /// Exponential: rough at every scale, like an erosion surface.
    Rough,
}

/// A stationary Gaussian random field with a Matérn 3/2 covariance,
/// sd^2 (1 + a) exp(-a) with a = sqrt(3) r / l, scaled so the correlation is
/// 0.05 at `range`. It is smooth, as sedimentary contacts are, and it is not
/// one of the variogram families the model fits, so the model is never handed
/// the true answer's shape. Drawn as a sum of random cosines whose
/// frequencies come from the covariance's spectral density: a bivariate
/// Student t with three degrees of freedom, scaled by 1/l. `Field::with`
/// can draw an exponential (rough) field instead.
#[derive(Debug, Clone)]
pub struct Field {
    mean: f64,
    amplitude: f64,
    waves: Vec<([f64; 2], f64)>,
}

/// (1 + a) exp(-a) = 0.05 at a = 4.7439.
const MATERN_RANGE_A: f64 = 4.743_864_9;

impl Field {
    pub fn new(rng: &mut Rng, mean: f64, sd: f64, range: f64) -> Field {
        Field::with(rng, Smoothness::Smooth, mean, sd, range)
    }

    /// The exponential covariance exp(-3 r / range) has as its spectral
    /// density a bivariate t with one degree of freedom; Matérn 3/2 one with
    /// three. A t with k degrees of freedom is a normal over sqrt(chi2_k / k).
    pub fn with(rng: &mut Rng, smoothness: Smoothness, mean: f64, sd: f64, range: f64) -> Field {
        let (dof, l) = match smoothness {
            Smoothness::Smooth => (3, 3f64.sqrt() * range / MATERN_RANGE_A),
            Smoothness::Rough => (1, range / 3.0),
        };
        let waves = (0..FEATURES)
            .map(|_| {
                let chi2: f64 = (0..dof).map(|_| rng.normal().powi(2)).sum::<f64>().max(1e-300);
                let s = 1.0 / (l * (chi2 / dof as f64).sqrt());
                let omega = [rng.normal() * s, rng.normal() * s];
                (omega, rng.range(0.0, std::f64::consts::TAU))
            })
            .collect();
        Field { mean, amplitude: sd * (2.0 / FEATURES as f64).sqrt(), waves }
    }

    pub fn at(&self, p: [f64; 2]) -> f64 {
        let s: f64 = self.waves.iter().map(|(w, b)| (w[0] * p[0] + w[1] * p[1] + b).cos()).sum();
        self.mean + self.amplitude * s
    }
}

/// The geology of one synthetic site.
#[derive(Debug, Clone)]
pub struct Truth {
    pub units: Vec<String>,
    ground: Field,
    thickness: Vec<Field>,
}

/// (name, mean thickness, sd, practical range). Clay is thick but variable
/// enough to pinch out; the basement (chalk) has no thickness.
const COLUMN: [(&str, f64, f64, f64); 4] =
    [("topsoil", 0.7, 0.25, 250.0), ("clay", 4.0, 2.2, 500.0), ("sand", 5.0, 1.2, 400.0), ("gravel", 3.0, 0.8, 400.0)];
const BASEMENT: &str = "chalk";

impl Truth {
    pub fn new(seed: u64) -> Truth {
        Truth::with(seed, Smoothness::Smooth)
    }

    pub fn with(seed: u64, smoothness: Smoothness) -> Truth {
        let mut rng = Rng::new(seed);
        let ground = Field::with(&mut rng, smoothness, 60.0, 2.0, 800.0);
        let thickness = COLUMN.iter().map(|&(_, m, sd, r)| Field::with(&mut rng, smoothness, m, sd, r)).collect();
        let mut units: Vec<String> = COLUMN.iter().map(|c| c.0.to_string()).collect();
        units.push(BASEMENT.into());
        Truth { units, ground, thickness }
    }

    /// Elevations of the ground and of each base, top down.
    pub fn surfaces(&self, p: [f64; 2]) -> Vec<f64> {
        let mut z = self.ground.at(p);
        let mut out = vec![z];
        for t in &self.thickness {
            z -= t.at(p).max(0.0);
            out.push(z);
        }
        out
    }

    /// Drills `count` holes. About one in five stops before the basement.
    pub fn site(&self, seed: u64, count: usize) -> Site {
        let mut rng = Rng::new(seed.wrapping_mul(0x2545_F491_4F6C_DD1D) ^ 0xB0E);
        let mut spots: Vec<[f64; 2]> = Vec::new();
        while spots.len() < count {
            let p = [round(rng.range(0.0, SIZE), 1), round(rng.range(0.0, SIZE), 1)];
            if spots.iter().all(|q| (p[0] - q[0]).hypot(p[1] - q[1]) >= SPACING_M) {
                spots.push(p);
            }
        }
        let holes = spots
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                let s = self.surfaces(p);
                let collar = round(s[0], 2);
                // Contact depths as logged: to 0.1 m, never decreasing.
                let mut depths: Vec<f64> = s[1..].iter().map(|z| round(collar - z, 1)).collect();
                for k in 1..depths.len() {
                    depths[k] = depths[k].max(depths[k - 1]);
                }
                let into_basement = round(rng.range(1.0, 3.0), 1);
                let mut end = round(depths.last().copied().unwrap_or(0.0) + into_basement, 1);
                if rng.uniform() < 0.2 {
                    end = round(rng.range(0.3, 0.95) * end, 1).max(0.1);
                }
                let mut intervals = Vec::new();
                let mut top = 0.0;
                for (unit, base) in depths.iter().copied().chain(std::iter::once(f64::INFINITY)).enumerate() {
                    let bottom = base.min(end);
                    if bottom > top {
                        intervals.push(Interval { unit, from: top, to: bottom });
                        top = bottom;
                    }
                    if top >= end {
                        break;
                    }
                }
                Hole { id: format!("BH{:02}", i + 1), x: p[0], y: p[1], collar, intervals }
            })
            .collect();
        Site { units: self.units.clone(), holes }
    }
}

fn round(v: f64, places: i32) -> f64 {
    let f = 10f64.powi(places);
    (v * f).round() / f
}

/// Comparison of stated intervals with the truth, per surface.
#[derive(Debug, Clone)]
pub struct Bench {
    pub surfaces: Vec<(String, SurfaceStats)>,
    pub sites: usize,
}

/// Builds a model for each of `sites` synthetic sites and scores it at
/// `points` random locations inside the site (100 m in from the edge, where
/// a real survey would be asked questions).
pub fn bench(
    sites: usize,
    holes: usize,
    points: usize,
    max_sigma: f64,
    smoothness: Smoothness,
) -> Result<Bench, String> {
    let truth0 = Truth::new(1);
    let names: Vec<String> = std::iter::once("ground surface".to_string())
        .chain(truth0.units[..truth0.units.len() - 1].iter().map(|u| format!("base of {u}")))
        .collect();
    let mut stats = vec![SurfaceStats::default(); names.len()];
    for seed in 1..=sites as u64 {
        let truth = Truth::with(seed, smoothness);
        let model = Model::build(&truth.site(seed, holes))?;
        let mut rng = Rng::new(seed ^ 0xBE4C);
        for _ in 0..points {
            let p = [rng.range(100.0, SIZE - 100.0), rng.range(100.0, SIZE - 100.0)];
            let column = model.column(p);
            for (s, z) in truth.surfaces(p).into_iter().enumerate() {
                let e: Option<Estimate> = match column.answer(s, f64::INFINITY) {
                    Answer::Known(e) | Answer::TooUncertain(e) => Some(e),
                    Answer::NoData => None,
                };
                if let Some(e) = e {
                    stats[s].add(e.value - z, e.sigma(), max_sigma);
                }
            }
        }
    }
    Ok(Bench { surfaces: names.into_iter().zip(stats).collect(), sites })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::logs::parse;

    #[test]
    fn field_has_the_stated_variance() {
        // Over many independent fields, the variance at one point is sd^2.
        let n = 4000;
        let vals: Vec<f64> = (0..n).map(|s| Field::new(&mut Rng::new(s), 0.0, 2.0, 300.0).at([123.0, 456.0])).collect();
        let var = vals.iter().map(|v| v * v).sum::<f64>() / n as f64;
        assert!((var - 4.0).abs() < 0.3, "variance {var}");
    }

    #[test]
    fn field_has_the_stated_correlation() {
        // Matérn 3/2: 0.05 at the range by construction; at a tenth of it, (1 + a) exp(-a) with a = 0.474.
        let n = 4000;
        let a = MATERN_RANGE_A / 10.0;
        for (h, expect) in [(30.0, (1.0 + a) * (-a).exp()), (300.0, 0.05)] {
            let mut sxy = 0.0;
            for s in 0..n {
                let f = Field::new(&mut Rng::new(s), 0.0, 1.0, 300.0);
                sxy += f.at([0.0, 0.0]) * f.at([h, 0.0]);
            }
            let rho = sxy / n as f64;
            assert!((rho - expect).abs() < 0.06, "h={h}: rho {rho}, expected {expect}");
        }
    }

    #[test]
    fn rough_field_has_exponential_correlation() {
        let n = 4000;
        for (h, expect) in [(30.0, (-0.3f64).exp()), (300.0, (-3.0f64).exp())] {
            let mut sxy = 0.0;
            for s in 0..n {
                let f = Field::with(&mut Rng::new(s), Smoothness::Rough, 0.0, 1.0, 300.0);
                sxy += f.at([0.0, 0.0]) * f.at([h, 0.0]);
            }
            let rho = sxy / n as f64;
            assert!((rho - expect).abs() < 0.06, "h={h}: rho {rho}, expected {expect}");
        }
    }

    #[test]
    fn surfaces_never_cross() {
        let t = Truth::new(5);
        let mut r = Rng::new(1);
        for _ in 0..2000 {
            let s = t.surfaces([r.range(0.0, SIZE), r.range(0.0, SIZE)]);
            assert!(s.windows(2).all(|w| w[1] <= w[0]));
        }
    }

    #[test]
    fn site_is_valid_input_and_reproducible() {
        let t = Truth::new(9);
        let site = t.site(9, 40);
        assert_eq!(site.holes.len(), 40);
        let text = site.to_csv("synthetic");
        assert_eq!(parse(&text).unwrap(), site);
        assert_eq!(Truth::new(9).site(9, 40), site);
        let reached = site.holes.iter().filter(|h| h.last_unit() == site.units.len() - 1).count();
        assert!(reached > 25 && reached < 40, "{reached} of 40 reached the basement");
    }
}
