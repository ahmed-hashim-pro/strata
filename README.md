# strata

[![CI](https://github.com/ahmed-hashim-pro/strata/actions/workflows/ci.yml/badge.svg)](https://github.com/ahmed-hashim-pro/strata/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A command-line tool that turns borehole logs into a layered model of the
ground. For any point it gives the elevation of every contact with a 95%
interval. Where it can't give a useful interval, it refuses to answer.

Every subsurface model is a guess between the holes. A site investigation
might drill forty boreholes over a square kilometre, so most of the ground has
never been seen. Most tools still draw smooth, confident layers through it. The
section looks equally certain next to a borehole and 150 m from one, and
whoever designs the foundation or prices the excavation can't tell which parts
are measured and which are made up. `strata` makes the difference visible. It
also checks the intervals it states: of its 95% intervals, about 95% contain
the truth.

```
$ strata query examples/site.csv 490 440        # next to hole BH16
surface               elevation        95% interval   sigma
ground surface            63.63     63.28 .. 63.97       0.17
base of topsoil           62.64     62.28 .. 63.01       0.19
base of clay              59.34     58.65 .. 60.03       0.35
base of sand              56.79     55.99 .. 57.60       0.41
base of gravel            51.90     51.06 .. 52.74       0.43

$ strata query examples/site.csv 600 480        # 120 m from the nearest hole
surface               elevation        95% interval   sigma
ground surface            63.93     62.74 .. 65.13       0.61
base of topsoil           63.18     61.88 .. 64.48       0.66
base of clay              59.16     56.75 .. 61.57       1.23
base of sand              55.41     52.55 .. 58.26       1.45
base of gravel          REFUSED    sigma over limit    1.52
$ echo $?
1
```

![Cross-section A–B through six projected boreholes: coloured layers of topsoil, clay, sand, gravel and chalk between the holes, each contact with a grey 95% band; a hatched column between 520 m and 650 m where the deeper contacts are refused; red dashed lines for the true contacts of the synthetic site, which stay inside the grey bands](docs/section.svg)

*`strata section examples/site.csv --from 40,470 --to 980,560 --truth 7`.
The site is synthetic, so the true contacts are known and drawn in red. Each
hole shows its own log. Between the holes, the grey bands widen. In the
hatched column, 120 m from the nearest hole, the model gives the upper
contacts and refuses the deeper ones.*

![Plan map of the 1 km square site: blue cells shading from light (sure) near the boreholes to dark (unsure) between them, hatched where sigma for the base of gravel exceeds 1.5 m, which is 40% of the area; filled dots for holes that reached the base of gravel, hollow dots for the six that stopped above it; the section line A–B dashed across](docs/map.svg)

*`strata map examples/site.csv --line 40,470:980,560`. The answer to "where
do we need another borehole?" is the hatching.*

![Terminal recording: the two queries above, then strata bench scoring 50 synthetic sites against their true geology with 95.6% inside the 95% intervals, then a misspelt unit name rejected with its line number](docs/demo.gif)

## Quickstart

You need Rust 1.74 or newer. `strata` has no dependencies.

```bash
git clone https://github.com/ahmed-hashim-pro/strata && cd strata
cargo build --release
export PATH=$PWD/target/release:$PATH

strata validate examples/site.csv                # parse, fit, show the variograms
strata query    examples/site.csv 490 440        # a column at a point
strata section  examples/site.csv --from 40,470 --to 980,560 -o section.svg
strata map      examples/site.csv --surface clay -o clay.svg
strata cv       examples/site.csv                # leave-one-hole-out
strata bench                                     # 50 synthetic sites against the truth
```

`cargo test` runs the whole suite in a few seconds. `make check` adds
rustfmt, clippy, and the full 100-site calibration runs that CI does.

## The input

One CSV file of interval logs, preceded by the stratigraphic column:

```
# comments and blank lines are ignored
units: topsoil, clay, sand, gravel, chalk
hole,x,y,collar,unit,from,to
BH01,412.3,880.1,61.42,topsoil,0.0,0.7
BH01,412.3,880.1,61.42,clay,0.7,5.3
BH01,412.3,880.1,61.42,sand,5.3,10.1
...
```

- **`units`** lists the units from top to bottom. The last one is the
  basement: the model finds its top but never its base.
- **`collar`** is the ground elevation at the hole. `from` and `to` are depths
  below it, in the same units as `x` and `y` (metres in every example here).
- **A pinched-out unit** is a unit a hole skips between two units it did log.
  It counts as thickness 0 there, which is real information.
- **The unit a hole stops in** has an unknown base. The hole only shows that
  the base is deeper than the bottom of the hole. That is not used for
  fitting. `strata validate` checks it afterwards and lists any hole whose
  bottom lies more than 1.96σ below the predicted base, which would mean the
  model is confidently wrong there.

The parser fails closed. Each of these is an error that names the line:

- a misspelt unit
- a gap or an overlap between intervals
- a unit logged twice, or out of stratigraphic order
- a hole whose coordinates change between rows
- a hole's rows split up by another hole
- two holes at the same spot

```
$ strata validate tests/data/typo.csv
strata: tests/data/typo.csv: line 4: unit 'clya' is not in the units list (topsoil, clay, sand, chalk)
```

## How the model works

### Thicknesses, not elevations

The obvious approach is to interpolate each contact's elevation separately,
and it lets layers cross: base of sand above base of clay, between two holes
where neither does. `strata` interpolates the **ground surface** and each
unit's **thickness** (an isopach) instead. It clamps every thickness at zero
and stacks them:

```
base of unit k = ground - thickness(0) - ... - thickness(k)
```

Layers can't cross by construction, and a unit can thin out to nothing. The
property tests check this at thousands of random points on randomly generated
sites, including far outside the holes.

### Kriging, and why its intervals need checking

Each variable is interpolated by **ordinary kriging**. A variogram is fitted
by weighted least squares (exponential or spherical, nugget, sill, range) over
a fixed grid, so the same data always gives the same model. The kriging system
is factored once and solved for each query. Kriging gives the best linear
estimate under that variogram, along with its error variance. The variance of
a contact is the sum of the variances above it.

That variance assumes the fitted variogram is exactly right, and a variogram
fitted to forty holes isn't. The error goes in both directions, and on
smooth geology it is large:

| Truth (100 synthetic sites) | What goes wrong | Coverage before | After |
| --- | --- | --- | --- |
| Smooth contacts (most sediments) | Exponential and spherical models overstate short-range variation, so the intervals are too wide (msse 0.59) | 98.2% | 95.3% |
| Rough contacts (erosion surfaces) | The fitted sill comes out a little low, so the intervals are slightly narrow (msse 1.09) | 93.7% | 94.0% |

So each variable is **measured** before it's used. Each hole is kriged from
the others under the fitted variogram, and the variance is scaled by the mean
squared standardised error. The factor is clamped to [0.25, 4] so a handful
of holes can't make the model wildly sure or unsure. `strata validate` prints
every factor and marks any that hit the limit.

### Refusing

A surface is answered only if its standard error σ is at most `--max-sigma`
(default 1.5 m, so a 95% interval of about ±2.9 m). Errors only accumulate
downward, so once one contact is refused, every contact below it is refused
too. A unit that fewer than six holes got through isn't modelled at all, and
nothing below it is answered.

`strata query` exits 1 if anything was refused, so a script can't mistake a
partial answer for a full one.

## How the intervals are checked

Real boreholes can't show what lies between them, so the intervals are tested
two ways.

**Against the truth, on synthetic sites.** `strata synth` builds a site whose
geology is known everywhere. The ground surface and each thickness are
Gaussian random fields, drawn as a sum of random cosines. Boreholes are
sampled from it and logged to 0.1 m, the way a driller would, and about one
in five stops before the basement. `strata bench` builds the model from the
holes and scores it at 200 random points per site *between* the holes. Two
deliberate choices:

- **The truth isn't from the model's family.** Smooth sites use a Matérn 3/2
  covariance, which `strata` can't fit; rough sites use an exponential. Both
  must pass.
- **Every surface must pass on its own**, not just the pooled figure.

`make bench` runs 100 sites of each kind:

| 100 sites × 200 points | Smooth contacts | Rough contacts |
| --- | --- | --- |
| Inside the 95% interval, all surfaces | **95.3%** | **94.0%** |
| Worst single surface | 94.9% | 93.7% |
| Mean squared standardised error (1 = honest) | 0.92 | 1.07 |
| Bias | ≤ 0.01 m | ≤ 0.07 m |
| Refused at σ > 1.5 m, base of gravel | 32% | 95% |
| Inside the interval, among answers given | 95.6% | 93.7% |

The refusal rate is supposed to differ. Rough contacts really are less
predictable from 40 holes, and the model says so instead of drawing them just
as smoothly. CI runs both benches on every push and fails if the pooled
coverage leaves 92–98%. The test suite checks each surface on 40 sites. It
also checks that more holes give smaller errors and fewer refusals.

**On the holes themselves, by cross-validation.** `strata cv` removes each
hole in turn and refits everything without it, variograms and variance
corrections included, so the hole can't leak into its own prediction. Then it
predicts every contact that hole logged. This is the check you can run on a
real site:

```
$ strata cv examples/site.csv
surface                 n    bias   rmse   msse  inside95   refused inside95 ans.
ground surface         40   -0.02   0.58   0.95     97.5%      0.0%         97.5%
base of topsoil        40   -0.01   0.62   0.94     97.5%      0.0%         97.5%
base of clay           40   -0.02   1.10   0.80    100.0%      2.5%        100.0%
base of sand           35    0.04   1.41   0.96     94.3%     22.9%         96.3%
base of gravel         34    0.00   1.47   0.98     97.1%     29.4%         95.8%
all                   189   -0.00   1.08   0.92     97.4%     10.1%         97.6%
```

Forty holes give a noisy coverage figure: one miss moves it by 2.5 points.
`--check` fails outside 85–99.5%, which still catches a model that's badly
off.

## Commands

| Command | Does |
| --- | --- |
| `strata validate <logs>` | Parse, fit, print the variograms and variance corrections, list holes that contradict the model |
| `strata query <logs> <x> <y> [--json]` | Every surface at a point, with its 95% interval or the reason it was refused |
| `strata section <logs> --from X,Y --to X,Y [--buffer M] [--truth SEED]` | SVG cross-section; holes within `--buffer` (60 m) are projected onto it |
| `strata map <logs> [--surface NAME] [--cell M] [--line X,Y:X,Y]` | SVG plan of σ for one surface (default: the deepest), with the refused area hatched |
| `strata cv <logs> [--check]` | Leave-one-hole-out cross-validation |
| `strata synth [--seed N] [--holes N]` | A synthetic site; the truth is reproducible from the seed |
| `strata bench [--sites N] [--holes N] [--points N] [--rough] [--check]` | Score the model against synthetic truth |

Every command that refuses takes `--max-sigma M`. Pictures go to stdout
unless you pass `-o`. **Exit codes:** 0 ok; 1 refused (`query`) or a failed
`--check`; 2 bad input or usage.

## How it is tested

`cargo test` runs 46 tests:

- **Unit tests (27).**
  - LU solve against random systems, including the zero diagonal of a kriging
    matrix.
  - Kriging is exact at the data, its weights sum to 1, and its variance grows
    with distance.
  - Variogram shape and fitting.
  - The parser's 19 ways to reject a file.
  - The RNG's moments.
  - The synthetic fields' variance and correlation, measured against the
    covariance they claim.
- **Properties (5).** These run over 25 random sites, both smooth and rough,
  with 12 to 57 holes:
  - layers never cross, anywhere
  - every hole is honoured exactly
  - variance only grows downward
  - a looser limit never answers less, and answers are always a prefix of
    the column
  - far from all holes, the model refuses
- **Calibration (3).**
  - Smooth sites, 40 sites × 150 points: every surface inside 90–98.5%,
    pooled 93–97%, msse 0.75–1.25.
  - Rough sites: the same bounds.
  - More holes give lower error and fewer refusals.
- **CLI (11).** The real binary's exit codes and output, that it exits quietly
  when a pipe closes, and a check that
  `examples/site.csv` is exactly what `strata synth --seed 7` produces.

CI runs:

- the tests on Rust 1.74 and on stable
- rustfmt and clippy with warnings as errors
- both 100-site benches
- a check that `examples/` and `docs/` match what the current code generates

## What this does not do

- **Layer-cake stratigraphy only.** Each unit is a single layer, in one
  order. Faults, lenses, channels, repeated units and overturned beds can't be
  represented. The parser rejects a unit logged twice in one hole instead of
  guessing.
- **Thicknesses are kriged independently.** Adding their variances assumes
  their errors are uncorrelated. Real thicknesses often correlate (a thick
  clay over a thin sand). The leave-one-out correction and `strata cv` catch
  the effect on the intervals, but cokriging would model it.
- **Stationary variograms.** One variogram per variable for the whole site.
  On a site that's really two settings (a valley fill against a terrace),
  split the file.
- **No censored-data kriging.** A hole that stopped early contributes nothing
  to the unit it stopped in. It's checked against the model but not used to
  fit it.
- **The 95% is Gaussian.** Intervals are estimate ± 1.96σ. The benches show
  this holds for the fields tested here, not that it holds for every real
  site. Run `strata cv` on yours.
- **Depths are vertical.** Inclined holes need their intervals converted to
  true vertical depth first.

## Layout

```
src/
  logs.rs        input format and its strict parser
  variogram.rs   experimental variogram, model fitting
  kriging.rs     ordinary kriging, factored once per variable
  linalg.rs      LU with partial pivoting
  model.rs       ground + thicknesses, the leave-one-out correction, refusal
  cv.rs          leave-one-hole-out cross-validation
  synth.rs       synthetic sites with known truth, and the bench
  svg.rs         section and map
  rng.rs         SplitMix64, so every result is reproducible from a seed
  main.rs        the CLI
tests/           properties, calibration, CLI
examples/        site.csv (strata synth --seed 7)
docs/            the pictures above; `make demo` regenerates them
```

## License

MIT
