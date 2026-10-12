//! Dense LU factorisation with partial pivoting.
//!
//! A kriging system is factored once per variable and then solved for every
//! query point, so factor and solve are separate. The ordinary-kriging matrix
//! has a zero on its diagonal (the Lagrange row), which is why pivoting is
//! not optional here.

#[derive(Debug, Clone)]
pub struct Lu {
    n: usize,
    a: Vec<f64>,
    perm: Vec<usize>,
}

impl Lu {
    /// Factors the row-major `n`x`n` matrix `a`. Returns `None` if it is
    /// singular to working precision.
    pub fn factor(mut a: Vec<f64>, n: usize) -> Option<Lu> {
        assert_eq!(a.len(), n * n, "matrix is not {n}x{n}");
        let scale = a.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if scale == 0.0 || !scale.is_finite() {
            return None;
        }
        let mut perm: Vec<usize> = (0..n).collect();
        for k in 0..n {
            let mut p = k;
            for i in k + 1..n {
                if a[i * n + k].abs() > a[p * n + k].abs() {
                    p = i;
                }
            }
            if a[p * n + k].abs() <= 1e-13 * scale {
                return None;
            }
            if p != k {
                for j in 0..n {
                    a.swap(k * n + j, p * n + j);
                }
                perm.swap(k, p);
            }
            let pivot = a[k * n + k];
            for i in k + 1..n {
                let f = a[i * n + k] / pivot;
                a[i * n + k] = f;
                if f != 0.0 {
                    for j in k + 1..n {
                        a[i * n + j] -= f * a[k * n + j];
                    }
                }
            }
        }
        Some(Lu { n, a, perm })
    }

    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        assert_eq!(b.len(), n);
        let mut x: Vec<f64> = self.perm.iter().map(|&p| b[p]).collect();
        for i in 0..n {
            let mut s = x[i];
            for j in 0..i {
                s -= self.a[i * n + j] * x[j];
            }
            x[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = x[i];
            for j in i + 1..n {
                s -= self.a[i * n + j] * x[j];
            }
            x[i] = s / self.a[i * n + i];
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    fn matvec(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
        (0..n).map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum()).collect()
    }

    #[test]
    fn solves_random_systems() {
        let mut r = Rng::new(11);
        for n in [1, 2, 5, 17, 40] {
            let a: Vec<f64> = (0..n * n).map(|_| r.range(-1.0, 1.0)).collect();
            let x: Vec<f64> = (0..n).map(|_| r.range(-5.0, 5.0)).collect();
            let b = matvec(&a, &x, n);
            let got = Lu::factor(a, n).expect("random matrix is regular").solve(&b);
            for (g, e) in got.iter().zip(&x) {
                assert!((g - e).abs() < 1e-8, "n={n}: {g} vs {e}");
            }
        }
    }

    #[test]
    fn needs_pivoting_for_zero_diagonal() {
        // The shape of an ordinary-kriging system: zero top-left and bottom-right.
        let a = vec![0.0, 1.0, 1.0, 0.0];
        let x = Lu::factor(a, 2).unwrap().solve(&[3.0, 4.0]);
        assert_eq!(x, vec![4.0, 3.0]);
    }

    #[test]
    fn refuses_singular() {
        assert!(Lu::factor(vec![1.0, 2.0, 2.0, 4.0], 2).is_none());
        assert!(Lu::factor(vec![0.0; 9], 3).is_none());
    }
}
