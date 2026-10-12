//! Ordinary kriging: the best linear unbiased estimate under a variogram, and
//! the variance of its error.
//!
//! The system is factored once; each query is one triangular solve.

use crate::linalg::Lu;
use crate::variogram::{dist, Variogram};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Estimate {
    pub value: f64,
    /// Kriging (error) variance, never negative.
    pub variance: f64,
}

impl Estimate {
    pub fn sigma(&self) -> f64 {
        self.variance.sqrt()
    }
}

#[derive(Debug, Clone)]
pub struct Kriging {
    points: Vec<[f64; 2]>,
    values: Vec<f64>,
    variogram: Variogram,
    lu: Lu,
}

impl Kriging {
    pub fn new(points: Vec<[f64; 2]>, values: Vec<f64>, variogram: Variogram) -> Result<Self, String> {
        let n = points.len();
        if n == 0 || n != values.len() {
            return Err("kriging needs at least one point and one value per point".into());
        }
        // [ Gamma 1 ] [ w  ]   [ gamma_0 ]
        // [ 1'    0 ] [ mu ] = [ 1       ]
        let m = n + 1;
        let mut a = vec![0.0; m * m];
        for i in 0..n {
            for j in 0..n {
                a[i * m + j] = variogram.gamma(dist(points[i], points[j]));
            }
            a[i * m + n] = 1.0;
            a[n * m + i] = 1.0;
        }
        let lu = Lu::factor(a, m).ok_or("kriging system is singular (two holes at the same place?)")?;
        Ok(Kriging { points, values, variogram, lu })
    }

    pub fn variogram(&self) -> &Variogram {
        &self.variogram
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Weights (summing to 1) and the Lagrange multiplier for a location.
    pub fn weights(&self, p: [f64; 2]) -> (Vec<f64>, f64) {
        let n = self.points.len();
        let mut b: Vec<f64> = self.points.iter().map(|&q| self.variogram.gamma(dist(p, q))).collect();
        b.push(1.0);
        let mut x = self.lu.solve(&b);
        let mu = x.pop().expect("system has n + 1 unknowns");
        debug_assert_eq!(x.len(), n);
        (x, mu)
    }

    pub fn predict(&self, p: [f64; 2]) -> Estimate {
        let (w, mu) = self.weights(p);
        let mut value = 0.0;
        let mut variance = mu;
        for (i, wi) in w.iter().enumerate() {
            value += wi * self.values[i];
            variance += wi * self.variogram.gamma(dist(p, self.points[i]));
        }
        Estimate { value, variance: variance.max(0.0) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;
    use crate::variogram::Kind;

    fn setup(seed: u64, nugget: f64) -> Kriging {
        let mut r = Rng::new(seed);
        let pts: Vec<[f64; 2]> = (0..25).map(|_| [r.range(0.0, 500.0), r.range(0.0, 500.0)]).collect();
        let vals: Vec<f64> = pts.iter().map(|p| 10.0 + 0.01 * p[0] + r.normal()).collect();
        let v = Variogram { kind: Kind::Exponential, nugget, sill: 2.0, range: 200.0 };
        Kriging::new(pts, vals, v).unwrap()
    }

    #[test]
    fn exact_at_the_data() {
        for nugget in [0.0, 0.3] {
            let k = setup(1, nugget);
            for i in 0..k.len() {
                let e = k.predict(k.points[i]);
                assert!((e.value - k.values[i]).abs() < 1e-8);
                assert!(e.variance < 1e-8);
            }
        }
    }

    #[test]
    fn weights_sum_to_one() {
        let k = setup(2, 0.1);
        let mut r = Rng::new(9);
        for _ in 0..50 {
            let (w, _) = k.weights([r.range(-200.0, 700.0), r.range(-200.0, 700.0)]);
            assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn variance_grows_away_from_the_data() {
        let k = setup(3, 0.0);
        let near = k.predict([k.points[0][0] + 1.0, k.points[0][1]]).variance;
        let far = k.predict([5000.0, 5000.0]).variance;
        assert!(near < far);
        // Far from everything, the error variance exceeds the sill (the mean
        // itself is estimated) and stays bounded.
        assert!(far > 2.0 && far < 2.0 * 2.0 + 1.0);
    }

    #[test]
    fn single_point_returns_its_value() {
        let v = Variogram { kind: Kind::Spherical, nugget: 0.0, sill: 1.0, range: 10.0 };
        let k = Kriging::new(vec![[0.0, 0.0]], vec![4.0], v).unwrap();
        let e = k.predict([3.0, 4.0]);
        assert!((e.value - 4.0).abs() < 1e-12);
        assert!((e.variance - 2.0 * v.gamma(5.0)).abs() < 1e-12);
    }

    #[test]
    fn coincident_points_are_refused() {
        let v = Variogram { kind: Kind::Exponential, nugget: 0.0, sill: 1.0, range: 10.0 };
        assert!(Kriging::new(vec![[1.0, 1.0], [1.0, 1.0]], vec![1.0, 2.0], v).is_err());
    }
}
