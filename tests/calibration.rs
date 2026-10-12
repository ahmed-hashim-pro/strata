//! The claim: of the stated 95% intervals, about 95% contain the truth.
//! Checked against synthetic sites where the truth between the holes is
//! known, for smooth and for rough contacts. The model fits exponential and
//! spherical variograms; neither truth is drawn from those, so it never gets
//! the true shape for free.

use strata::cv::SurfaceStats;
use strata::synth::{bench, Smoothness};

fn total(smoothness: Smoothness) -> SurfaceStats {
    let b = bench(40, 40, 150, 1.5, smoothness).unwrap();
    let mut t = SurfaceStats::default();
    for (name, s) in &b.surfaces {
        // Every surface on its own, not just the pooled number.
        assert!((0.90..=0.985).contains(&s.coverage()), "{smoothness:?} {name}: {:.3}", s.coverage());
        assert!(s.bias().abs() < 0.15, "{smoothness:?} {name}: bias {:.3}", s.bias());
        t.merge(s);
    }
    t
}

#[test]
fn smooth_contacts_are_calibrated() {
    let t = total(Smoothness::Smooth);
    assert!((0.93..=0.97).contains(&t.coverage()), "coverage {:.3}", t.coverage());
    assert!((0.75..=1.25).contains(&t.msse()), "msse {:.3}", t.msse());
}

#[test]
fn rough_contacts_are_calibrated() {
    let t = total(Smoothness::Rough);
    assert!((0.92..=0.97).contains(&t.coverage()), "coverage {:.3}", t.coverage());
    assert!((0.75..=1.25).contains(&t.msse()), "msse {:.3}", t.msse());
}

#[test]
fn more_holes_means_smaller_errors() {
    let rmse = |holes| {
        let b = bench(15, holes, 150, 1.5, Smoothness::Smooth).unwrap();
        let mut t = SurfaceStats::default();
        b.surfaces.iter().for_each(|(_, s)| t.merge(s));
        (t.rmse(), t.refused as f64 / t.n as f64)
    };
    let (sparse, dense) = (rmse(15), rmse(80));
    assert!(dense.0 < 0.75 * sparse.0, "rmse {sparse:?} -> {dense:?}");
    assert!(dense.1 < sparse.1, "refused fraction {sparse:?} -> {dense:?}");
}
