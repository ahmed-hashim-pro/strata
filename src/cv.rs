//! Leave-one-hole-out cross-validation.
//!
//! Each hole in turn is removed, the whole model is refitted without it
//! (variograms included, so the hole cannot leak into its own prediction), and
//! every surface that hole determines is predicted at its location. If the
//! stated 95% intervals are honest, about 95% of those predictions contain the
//! logged value, and the mean squared standardised error is about 1.

use crate::logs::Site;
use crate::model::{Model, Surface, Z95};

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SurfaceStats {
    pub n: usize,
    pub sum_error: f64,
    pub sum_sq_error: f64,
    pub sum_sq_standardised: f64,
    pub inside95: usize,
    /// Predictions whose sigma exceeds the limit (the model would refuse).
    pub refused: usize,
    /// Of the predictions the model would give, how many were inside.
    pub answered_inside95: usize,
}

impl SurfaceStats {
    pub fn add(&mut self, error: f64, sigma: f64, max_sigma: f64) {
        self.n += 1;
        self.sum_error += error;
        self.sum_sq_error += error * error;
        let inside = error.abs() <= Z95 * sigma;
        self.sum_sq_standardised += if sigma > 0.0 {
            (error / sigma).powi(2)
        } else if error == 0.0 {
            0.0
        } else {
            f64::INFINITY
        };
        self.inside95 += usize::from(inside);
        if sigma > max_sigma {
            self.refused += 1;
        } else {
            self.answered_inside95 += usize::from(inside);
        }
    }

    pub fn merge(&mut self, o: &SurfaceStats) {
        self.n += o.n;
        self.sum_error += o.sum_error;
        self.sum_sq_error += o.sum_sq_error;
        self.sum_sq_standardised += o.sum_sq_standardised;
        self.inside95 += o.inside95;
        self.refused += o.refused;
        self.answered_inside95 += o.answered_inside95;
    }

    pub fn bias(&self) -> f64 {
        self.sum_error / self.n as f64
    }
    pub fn rmse(&self) -> f64 {
        (self.sum_sq_error / self.n as f64).sqrt()
    }
    /// Mean squared standardised error: 1 when the variances are right.
    pub fn msse(&self) -> f64 {
        self.sum_sq_standardised / self.n as f64
    }
    pub fn coverage(&self) -> f64 {
        self.inside95 as f64 / self.n as f64
    }
    pub fn answered_coverage(&self) -> Option<f64> {
        let answered = self.n - self.refused;
        (answered > 0).then(|| self.answered_inside95 as f64 / answered as f64)
    }
}

#[derive(Debug, Clone)]
pub struct Report {
    pub surfaces: Vec<(String, SurfaceStats)>,
    /// Holes whose removal left the model unbuildable.
    pub skipped: Vec<String>,
}

impl Report {
    pub fn total(&self) -> SurfaceStats {
        let mut t = SurfaceStats::default();
        for (_, s) in &self.surfaces {
            t.merge(s);
        }
        t
    }
}

pub fn leave_one_out(site: &Site, max_sigma: f64) -> Result<Report, String> {
    let full = Model::build(site)?;
    let mut surfaces: Vec<(String, SurfaceStats)> =
        (0..full.surface_count()).map(|s| (full.surface_name(s), SurfaceStats::default())).collect();
    let mut skipped = Vec::new();
    for (i, hole) in site.holes.iter().enumerate() {
        let Ok(model) = Model::build(&site.without_hole(i)) else {
            skipped.push(hole.id.clone());
            continue;
        };
        let column = model.column(hole.xy());
        let truth = std::iter::once(Some(hole.collar)).chain((0..site.units.len() - 1).map(|u| hole.base(u)));
        for (s, logged) in truth.enumerate() {
            if let (Some(logged), Surface::Known(e)) = (logged, column.surfaces[s]) {
                surfaces[s].1.add(e.value - logged, e.sigma(), max_sigma);
            }
        }
    }
    surfaces.retain(|(_, s)| s.n > 0);
    Ok(Report { surfaces, skipped })
}
