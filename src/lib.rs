//! strata: borehole logs to a layered subsurface model that says how sure it
//! is between the holes, and refuses to answer where it isn't.

// Index loops read like the maths they implement (LU, kriging systems, station
// runs); iterator chains would hide which row is which.
#![allow(clippy::needless_range_loop)]

pub mod cv;
pub mod kriging;
pub mod linalg;
pub mod logs;
pub mod model;
pub mod rng;
pub mod svg;
pub mod synth;
pub mod variogram;
