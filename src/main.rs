use std::fmt::Write as _;
use std::process::ExitCode;

use strata::cv::{self, SurfaceStats};
use strata::logs::{self, Site};
use strata::model::{Answer, Fitted, Model, DEFAULT_MAX_SIGMA, VARIANCE_SCALE_LIMITS, Z95};
use strata::svg::{self, MapOptions, SectionOptions};
use strata::synth::{self, Smoothness, Truth};

const USAGE: &str = "\
strata: borehole logs to a layered subsurface model with honest uncertainty

usage:
  strata validate <logs.csv>                     parse, fit, report the variograms and any holes the model contradicts
  strata query    <logs.csv> <x> <y> [--json]    elevations of every surface at a point, with 95% intervals
  strata section  <logs.csv> --from X,Y --to X,Y [--buffer M] [--truth SEED] [-o out.svg]
  strata map      <logs.csv> [--surface NAME] [--cell M] [--line X,Y:X,Y] [-o out.svg]
  strata cv       <logs.csv> [--check]           leave-one-hole-out cross-validation
  strata synth    [--seed N] [--holes N] [-o out.csv]
  strata bench    [--sites N] [--holes N] [--points N] [--rough] [--check]
  strata version

options:
  --max-sigma M   refuse any surface whose standard error exceeds M metres (default 1.5)

exit codes: 0 ok, 1 refused (query) or check failed (cv, bench), 2 bad input or usage";

struct Args {
    positional: Vec<String>,
    flags: Vec<(String, Option<String>)>,
}

impl Args {
    fn parse(raw: &[String]) -> Result<Args, String> {
        const SWITCHES: [&str; 3] = ["--json", "--check", "--rough"];
        let (mut positional, mut flags) = (Vec::new(), Vec::new());
        let mut it = raw.iter();
        while let Some(a) = it.next() {
            if a.starts_with('-') && a.len() > 1 && a.parse::<f64>().is_err() {
                if SWITCHES.contains(&a.as_str()) {
                    flags.push((a.clone(), None));
                } else {
                    let v = it.next().ok_or(format!("{a} needs a value"))?;
                    flags.push((a.clone(), Some(v.clone())));
                }
            } else {
                positional.push(a.clone());
            }
        }
        Ok(Args { positional, flags })
    }

    fn has(&self, name: &str) -> bool {
        self.flags.iter().any(|(n, _)| n == name)
    }

    fn get(&self, name: &str) -> Option<&str> {
        self.flags.iter().rev().find(|(n, _)| n == name).and_then(|(_, v)| v.as_deref())
    }

    fn num<T: std::str::FromStr>(&self, name: &str, default: T) -> Result<T, String> {
        match self.get(name) {
            None => Ok(default),
            Some(v) => v.parse().map_err(|_| format!("{name} '{v}' is not a valid number")),
        }
    }

    fn check_known(&self, allowed: &[&str]) -> Result<(), String> {
        for (n, _) in &self.flags {
            if !allowed.contains(&n.as_str()) {
                return Err(format!("unknown option {n}"));
            }
        }
        Ok(())
    }

    fn max_sigma(&self) -> Result<f64, String> {
        let m: f64 = self.num("--max-sigma", DEFAULT_MAX_SIGMA)?;
        if m.is_nan() || m <= 0.0 {
            return Err("--max-sigma must be positive".into());
        }
        Ok(m)
    }
}

fn point(s: &str) -> Result<[f64; 2], String> {
    let parts: Vec<&str> = s.split(',').collect();
    match parts.as_slice() {
        [x, y] => match (x.trim().parse::<f64>(), y.trim().parse::<f64>()) {
            (Ok(x), Ok(y)) if x.is_finite() && y.is_finite() => Ok([x, y]),
            _ => Err(format!("'{s}' is not X,Y")),
        },
        _ => Err(format!("'{s}' is not X,Y")),
    }
}

fn load(path: &str) -> Result<Site, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    logs::parse(&text).map_err(|e| format!("{path}: {e}"))
}

fn output(path: Option<&str>, content: &str) -> Result<(), String> {
    match path {
        None | Some("-") => {
            print!("{content}");
            Ok(())
        }
        Some(p) => std::fs::write(p, content).map_err(|e| format!("{p}: {e}")),
    }
}

fn one_file(args: &Args, cmd: &str) -> Result<String, String> {
    match args.positional.as_slice() {
        [f] => Ok(f.clone()),
        _ => Err(format!("{cmd} takes one logs file")),
    }
}

fn validate(args: &Args) -> Result<u8, String> {
    args.check_known(&[])?;
    let site = load(&one_file(args, "validate")?)?;
    let model = Model::build(&site)?;
    println!("{} holes, {} units: {}", site.holes.len(), site.units.len(), site.units.join(" > "));
    println!();
    println!("{:<20} {:>5}  {:>7}  variogram", "modelled as", "holes", "var. ×");
    let g = &model.ground;
    let scale = |f: &Fitted| {
        let at_limit = f.variance_scale == VARIANCE_SCALE_LIMITS.0 || f.variance_scale == VARIANCE_SCALE_LIMITS.1;
        format!("{:.2}{}", f.variance_scale, if at_limit { "*" } else { " " })
    };
    let mut clamped = false;
    clamped |= scale(g).ends_with('*');
    println!("{:<20} {:>5}  {:>7}  {}", "ground elevation", g.observations, scale(g), g.variogram());
    for (u, t) in model.thickness.iter().enumerate() {
        let name = format!("{} thickness", site.units[u]);
        match t {
            Ok(f) => {
                clamped |= scale(f).ends_with('*');
                println!("{:<20} {:>5}  {:>7}  {}", name, f.observations, scale(f), f.variogram())
            }
            Err(e) => println!("{:<20} {:>5}  {:>7}  not modelled: {e}", name, "-", "-"),
        }
    }
    println!("\nvar. × is the kriging variance correction measured by leave-one-out on each variable.");
    if clamped {
        println!(
            "* at its limit ({} or {}): the fitted variogram family fits this variable poorly; check `strata cv`.",
            VARIANCE_SCALE_LIMITS.0, VARIANCE_SCALE_LIMITS.1
        );
    }
    let stopped = site.holes.iter().filter(|h| h.last_unit() + 1 < site.units.len()).count();
    println!();
    println!("{stopped} holes stopped above the basement; their last unit's base is only known to be deeper.");
    let bad = model.contradictions(&site);
    if bad.is_empty() {
        println!("None of them contradicts the model.");
    } else {
        println!("{} contradict the model (it is confidently too deep there):", bad.len());
        for b in &bad {
            println!("  {b}");
        }
    }
    Ok(0)
}

fn query(args: &Args) -> Result<u8, String> {
    args.check_known(&["--json", "--max-sigma"])?;
    let max_sigma = args.max_sigma()?;
    let [file, x, y] = args.positional.as_slice() else {
        return Err("query takes <logs.csv> <x> <y>".into());
    };
    let p = point(&format!("{x},{y}"))?;
    let model = Model::build(&load(file)?)?;
    let column = model.column(p);
    let mut refused = false;
    let mut out = String::new();
    if args.has("--json") {
        let mut items = Vec::new();
        for s in 0..model.surface_count() {
            let name = model.surface_name(s);
            items.push(match column.answer(s, max_sigma) {
                Answer::Known(e) => format!(
                    "{{\"surface\":\"{name}\",\"status\":\"known\",\"elevation\":{:.3},\"sigma\":{:.3},\"low95\":{:.3},\"high95\":{:.3}}}",
                    e.value, e.sigma(), e.value - Z95 * e.sigma(), e.value + Z95 * e.sigma()
                ),
                Answer::TooUncertain(e) => {
                    refused = true;
                    format!("{{\"surface\":\"{name}\",\"status\":\"too_uncertain\",\"sigma\":{:.3},\"max_sigma\":{max_sigma}}}", e.sigma())
                }
                Answer::NoData => {
                    refused = true;
                    format!("{{\"surface\":\"{name}\",\"status\":\"no_data\"}}")
                }
            });
        }
        let _ = writeln!(out, "{{\"x\":{},\"y\":{},\"surfaces\":[{}]}}", p[0], p[1], items.join(","));
    } else {
        let _ = writeln!(out, "at ({}, {}), max sigma {max_sigma} m", p[0], p[1]);
        let _ = writeln!(out, "{:<20} {:>10}  {:>18}  {:>6}", "surface", "elevation", "95% interval", "sigma");
        for s in 0..model.surface_count() {
            let name = model.surface_name(s);
            let line = match column.answer(s, max_sigma) {
                Answer::Known(e) => format!(
                    "{name:<20} {:>10.2}  {:>8.2} .. {:<8.2}  {:>6.2}",
                    e.value,
                    e.value - Z95 * e.sigma(),
                    e.value + Z95 * e.sigma(),
                    e.sigma()
                ),
                Answer::TooUncertain(e) => {
                    refused = true;
                    format!("{name:<20} {:>10}  {:>18}  {:>6.2}", "REFUSED", "sigma over limit", e.sigma())
                }
                Answer::NoData => {
                    refused = true;
                    format!("{name:<20} {:>10}  {:>18}", "REFUSED", "too few holes")
                }
            };
            let _ = writeln!(out, "{line}");
        }
    }
    print!("{out}");
    Ok(u8::from(refused))
}

fn section(args: &Args) -> Result<u8, String> {
    args.check_known(&["--from", "--to", "--buffer", "--truth", "--max-sigma", "-o"])?;
    let site = load(&one_file(args, "section")?)?;
    let from = point(args.get("--from").ok_or("section needs --from X,Y")?)?;
    let to = point(args.get("--to").ok_or("section needs --to X,Y")?)?;
    if from == to {
        return Err("--from and --to are the same point".into());
    }
    let truth = match args.get("--truth") {
        None => None,
        Some(s) => Some(Truth::new(s.parse().map_err(|_| format!("--truth '{s}' is not a seed"))?)),
    };
    let model = Model::build(&site)?;
    let opts = SectionOptions {
        from,
        to,
        buffer: args.num("--buffer", 60.0)?,
        max_sigma: args.max_sigma()?,
        truth: truth.as_ref(),
    };
    output(args.get("-o"), &svg::section(&site, &model, &opts))?;
    Ok(0)
}

fn map(args: &Args) -> Result<u8, String> {
    args.check_known(&["--surface", "--cell", "--line", "--max-sigma", "-o"])?;
    let site = load(&one_file(args, "map")?)?;
    let model = Model::build(&site)?;
    let surface = match args.get("--surface") {
        None => model.surface_count() - 1,
        Some(n) => model.surface_by_name(n).ok_or_else(|| {
            format!(
                "--surface '{n}': expected 'ground' or a unit with a base ({})",
                model.units[..model.units.len() - 1].join(", ")
            )
        })?,
    };
    let line = match args.get("--line") {
        None => None,
        Some(l) => {
            let (a, b) = l.split_once(':').ok_or("--line is X,Y:X,Y")?;
            Some((point(a)?, point(b)?))
        }
    };
    let cell: f64 = args.num("--cell", 0.0)?;
    let opts = MapOptions { surface, cell: (cell > 0.0).then_some(cell), max_sigma: args.max_sigma()?, line };
    output(args.get("-o"), &svg::map(&site, &model, &opts))?;
    Ok(0)
}

fn stats_table(rows: &[(String, SurfaceStats)], max_sigma: f64) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "{:<18} {:>6} {:>7} {:>6} {:>6} {:>9}  {:>8} {:>13}",
        "surface", "n", "bias", "rmse", "msse", "inside95", "refused", "inside95 ans."
    );
    let mut total = SurfaceStats::default();
    let mut row = |name: &str, s: &SurfaceStats| {
        let ans = s.answered_coverage().map_or("-".into(), |c| format!("{:.1}%", 100.0 * c));
        let _ = writeln!(
            out,
            "{:<18} {:>6} {:>7.2} {:>6.2} {:>6.2} {:>8.1}%  {:>7.1}% {:>13}",
            name,
            s.n,
            s.bias(),
            s.rmse(),
            s.msse(),
            100.0 * s.coverage(),
            100.0 * s.refused as f64 / s.n as f64,
            ans
        );
    };
    for (name, s) in rows {
        row(name, s);
        total.merge(s);
    }
    row("all", &total);
    let _ =
        writeln!(out, "\nbias and rmse in metres; msse ~ 1 and inside95 ~ 95% when the stated intervals are honest;");
    let _ = write!(
        out,
        "refused = sigma over {max_sigma} m; the last column is coverage among the answers the model would give."
    );
    out
}

/// The interval check used by `cv --check` and `bench --check`.
fn calibrated(total: &SurfaceStats, lo: f64, hi: f64) -> bool {
    let c = total.coverage();
    (lo..=hi).contains(&c)
}

fn cross_validate(args: &Args) -> Result<u8, String> {
    args.check_known(&["--check", "--max-sigma"])?;
    let max_sigma = args.max_sigma()?;
    let site = load(&one_file(args, "cv")?)?;
    let report = cv::leave_one_out(&site, max_sigma)?;
    println!("leave-one-hole-out over {} holes (variograms refitted each time)\n", site.holes.len());
    println!("{}", stats_table(&report.surfaces, max_sigma));
    if !report.skipped.is_empty() {
        println!("\nskipped (model unbuildable without them): {}", report.skipped.join(", "));
    }
    let total = report.total();
    if args.has("--check") {
        // Leave-one-out has few predictions per site, so the band is wide.
        let ok = calibrated(&total, 0.85, 0.995);
        println!(
            "\ncheck: {:.1}% inside the 95% intervals: {}",
            100.0 * total.coverage(),
            if ok { "ok" } else { "FAIL (expected 85% - 99.5%)" }
        );
        return Ok(u8::from(!ok));
    }
    Ok(0)
}

fn synthesize(args: &Args) -> Result<u8, String> {
    args.check_known(&["--seed", "--holes", "-o"])?;
    let seed: u64 = args.num("--seed", 7)?;
    let holes: usize = args.num("--holes", 40)?;
    if holes == 0 || holes > 500 {
        return Err("--holes must be 1..500".into());
    }
    let site = Truth::new(seed).site(seed, holes);
    let comment = format!(
        "Synthetic site: strata synth --seed {seed} --holes {holes}\n{} m square; the true geology is reproducible from the seed.",
        synth::SIZE
    );
    output(args.get("-o"), &site.to_csv(&comment))?;
    Ok(0)
}

fn bench(args: &Args) -> Result<u8, String> {
    args.check_known(&["--sites", "--holes", "--points", "--check", "--rough", "--max-sigma"])?;
    let max_sigma = args.max_sigma()?;
    let sites: usize = args.num("--sites", 50)?;
    let holes: usize = args.num("--holes", 40)?;
    let points: usize = args.num("--points", 200)?;
    if sites == 0 || points == 0 || !(6..=500).contains(&holes) {
        return Err("--sites and --points must be positive, --holes 6..500".into());
    }
    let smoothness = if args.has("--rough") { Smoothness::Rough } else { Smoothness::Smooth };
    let b = synth::bench(sites, holes, points, max_sigma, smoothness)?;
    println!(
        "{} synthetic sites with {} contacts, {holes} holes and {points} test points each, scored against the truth\n",
        b.sites,
        if smoothness == Smoothness::Rough { "rough" } else { "smooth" }
    );
    println!("{}", stats_table(&b.surfaces, max_sigma));
    if args.has("--check") {
        let mut total = SurfaceStats::default();
        for (_, s) in &b.surfaces {
            total.merge(s);
        }
        let ok = calibrated(&total, 0.92, 0.98);
        println!(
            "\ncheck: {:.1}% inside the 95% intervals: {}",
            100.0 * total.coverage(),
            if ok { "ok" } else { "FAIL (expected 92% - 98%)" }
        );
        return Ok(u8::from(!ok));
    }
    Ok(0)
}

fn run(raw: &[String]) -> Result<u8, String> {
    let Some((cmd, rest)) = raw.split_first() else {
        return Err(USAGE.into());
    };
    let args = Args::parse(rest)?;
    match cmd.as_str() {
        "validate" => validate(&args),
        "query" => query(&args),
        "section" => section(&args),
        "map" => map(&args),
        "cv" => cross_validate(&args),
        "synth" => synthesize(&args),
        "bench" => bench(&args),
        "version" | "--version" => {
            println!("strata {}", env!("CARGO_PKG_VERSION"));
            Ok(0)
        }
        "help" | "--help" | "-h" => {
            println!("{USAGE}");
            Ok(0)
        }
        other => Err(format!("unknown command '{other}'\n\n{USAGE}")),
    }
}

/// Rust ignores SIGPIPE, so `strata validate f | head` would panic on the
/// first write after `head` exits. Restore the default: die quietly, as
/// every other filter in the pipe does.
#[cfg(unix)]
fn restore_sigpipe() {
    extern "C" {
        fn signal(signum: i32, handler: usize) -> usize;
    }
    const SIGPIPE: i32 = 13;
    const SIG_DFL: usize = 0;
    // SAFETY: setting a signal to its default disposition, before any threads exist.
    unsafe {
        signal(SIGPIPE, SIG_DFL);
    }
}

#[cfg(not(unix))]
fn restore_sigpipe() {}

fn main() -> ExitCode {
    restore_sigpipe();
    let raw: Vec<String> = std::env::args().skip(1).collect();
    match run(&raw) {
        Ok(code) => ExitCode::from(code),
        Err(e) => {
            eprintln!("strata: {e}");
            ExitCode::from(2)
        }
    }
}
