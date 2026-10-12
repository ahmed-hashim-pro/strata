//! Properties that must hold for any site, checked over many random ones.

use strata::logs::Site;
use strata::model::{Answer, Model, Surface};
use strata::rng::Rng;
use strata::synth::{Smoothness, Truth, SIZE};

fn sites() -> impl Iterator<Item = (u64, Site)> {
    (1..=25u64).map(|seed| {
        let smooth = if seed % 2 == 0 { Smoothness::Smooth } else { Smoothness::Rough };
        let holes = 8 + (seed as usize * 7) % 50;
        (seed, Truth::with(seed, smooth).site(seed, holes))
    })
}

#[test]
fn layers_never_cross() {
    for (seed, site) in sites() {
        let model = Model::build(&site).unwrap();
        let mut r = Rng::new(seed);
        for _ in 0..400 {
            // Including well outside the holes, where kriging extrapolates.
            let p = [r.range(-500.0, SIZE + 500.0), r.range(-500.0, SIZE + 500.0)];
            let z: Vec<f64> = model
                .column(p)
                .surfaces
                .iter()
                .filter_map(|s| match s {
                    Surface::Known(e) => Some(e.value),
                    Surface::NoData => None,
                })
                .collect();
            assert!(z.windows(2).all(|w| w[1] <= w[0]), "seed {seed} at {p:?}: {z:?}");
        }
    }
}

#[test]
fn honours_every_hole_exactly() {
    for (seed, site) in sites() {
        let model = Model::build(&site).unwrap();
        for h in &site.holes {
            let c = model.column(h.xy());
            let Surface::Known(g) = c.surfaces[0] else { panic!() };
            assert!((g.value - h.collar).abs() < 1e-6, "seed {seed} {}", h.id);
            assert!(g.variance < 1e-9);
            for u in 0..site.units.len() - 1 {
                if let (Some(base), Surface::Known(e)) = (h.base(u), c.surfaces[u + 1]) {
                    assert!((e.value - base).abs() < 1e-6, "seed {seed} {} base of {u}: {} vs {base}", h.id, e.value);
                }
            }
        }
    }
}

#[test]
fn uncertainty_only_grows_downward() {
    for (seed, site) in sites() {
        let model = Model::build(&site).unwrap();
        let mut r = Rng::new(seed + 100);
        for _ in 0..200 {
            let c = model.column([r.range(0.0, SIZE), r.range(0.0, SIZE)]);
            let v: Vec<f64> = c
                .surfaces
                .iter()
                .filter_map(|s| if let Surface::Known(e) = s { Some(e.variance) } else { None })
                .collect();
            assert!(v.windows(2).all(|w| w[1] >= w[0]), "seed {seed}: {v:?}");
        }
    }
}

#[test]
fn a_looser_limit_never_answers_less() {
    for (seed, site) in sites().take(8) {
        let model = Model::build(&site).unwrap();
        let mut r = Rng::new(seed + 200);
        for _ in 0..200 {
            let c = model.column([r.range(-200.0, SIZE + 200.0), r.range(-200.0, SIZE + 200.0)]);
            let mut prev = 0;
            for limit in [0.1, 0.5, 1.0, 1.5, 2.0, 5.0, f64::INFINITY] {
                let a = c.answered(limit);
                assert!(a >= prev);
                prev = a;
                // Answered surfaces form a prefix: after the first refusal, nothing is known.
                for s in 0..c.surfaces.len() {
                    assert_eq!(matches!(c.answer(s, limit), Answer::Known(_)) && s < a, s < a);
                }
            }
        }
    }
}

#[test]
fn far_from_every_hole_the_model_refuses() {
    for (_, site) in sites() {
        let model = Model::build(&site).unwrap();
        let c = model.column([SIZE * 20.0, -SIZE * 20.0]);
        assert!(c.answered(1.5) < c.surfaces.len());
    }
}
