//! The CLI as a program: exit codes and output are its interface.

use std::process::{Command, Output};

fn strata(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_strata"))
        .args(args)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("runs")
}

fn stdout(o: &Output) -> String {
    String::from_utf8_lossy(&o.stdout).into_owned()
}

const SITE: &str = "examples/site.csv";

#[test]
fn validate_reports_the_fit() {
    let o = strata(&["validate", SITE]);
    assert_eq!(o.status.code(), Some(0));
    let s = stdout(&o);
    assert!(s.contains("40 holes, 5 units: topsoil > clay > sand > gravel > chalk"), "{s}");
    assert!(s.contains("clay thickness"));
    assert!(s.contains("leave-one-out"));
}

#[test]
fn query_answers_near_holes_and_exits_zero() {
    // BH16 is at (487.7, 435.7).
    let o = strata(&["query", SITE, "490", "440"]);
    assert_eq!(o.status.code(), Some(0), "{}", stdout(&o));
    assert!(!stdout(&o).contains("REFUSED"));
}

#[test]
fn query_refuses_far_away_and_exits_one() {
    let o = strata(&["query", SITE, "5000", "5000"]);
    assert_eq!(o.status.code(), Some(1));
    assert!(stdout(&o).contains("REFUSED"));
}

#[test]
fn query_limit_is_adjustable() {
    let strict = strata(&["query", SITE, "560", "480", "--max-sigma", "0.2"]);
    let loose = strata(&["query", SITE, "560", "480", "--max-sigma", "10"]);
    assert_eq!(strict.status.code(), Some(1));
    assert_eq!(loose.status.code(), Some(0));
}

#[test]
fn query_json() {
    let o = strata(&["query", SITE, "5000", "5000", "--json"]);
    let s = stdout(&o);
    assert!(s.starts_with("{\"x\":5000,\"y\":5000,\"surfaces\":["), "{s}");
    assert!(s.contains("\"status\":\"too_uncertain\""));
}

#[test]
fn bad_input_names_the_line_and_exits_two() {
    let o = strata(&["validate", "tests/data/typo.csv"]);
    assert_eq!(o.status.code(), Some(2));
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(e.contains("line 4: unit 'clya' is not in the units list"), "{e}");
}

#[test]
fn usage_errors_exit_two() {
    for args in [
        &[][..],
        &["frobnicate"],
        &["query", SITE, "1"],
        &["map", SITE, "--surface", "chalk"],
        &["query", SITE, "1", "2", "--max-sigma", "-1"],
        &["validate", SITE, "--bogus", "1"],
    ] {
        let o = strata(args);
        assert_eq!(o.status.code(), Some(2), "{args:?}");
    }
}

#[test]
fn pictures_are_svg() {
    let s = stdout(&strata(&["section", SITE, "--from", "40,470", "--to", "980,560", "--truth", "7"]));
    assert!(s.starts_with("<svg") && s.trim_end().ends_with("</svg>"));
    assert!(s.contains("true contacts"));
    let m = stdout(&strata(&["map", SITE, "--surface", "clay"]));
    assert!(m.contains("Where the model answers: base of clay"));
}

#[test]
fn cross_validation_check_passes_on_the_example() {
    let o = strata(&["cv", SITE, "--check"]);
    assert_eq!(o.status.code(), Some(0), "{}", stdout(&o));
}

#[test]
fn example_is_what_synth_generates() {
    let o = strata(&["synth", "--seed", "7", "--holes", "40"]);
    let committed = std::fs::read_to_string(SITE).unwrap();
    assert_eq!(stdout(&o), committed, "examples/site.csv is stale: run `make examples`");
}
