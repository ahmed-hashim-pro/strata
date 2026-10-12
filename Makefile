# Thin wrappers over cargo, plus the targets that regenerate what the README shows.
CARGO ?= cargo
STRATA := target/release/strata

.PHONY: build test lint bench check examples demo clean

build:
	$(CARGO) build --release

test:
	$(CARGO) test

lint:
	$(CARGO) fmt --check
	$(CARGO) clippy --all-targets -- -D warnings

# The calibration claim, at full size: 100 synthetic sites each, smooth and rough.
bench: build
	$(STRATA) bench --sites 100 --check
	$(STRATA) bench --sites 100 --rough --check

check: lint test bench

examples: build
	$(STRATA) synth --seed 7 --holes 40 -o examples/site.csv

# Re-render the README pictures and re-record docs/demo.gif from the real CLI.
# The GIF needs Docker (VHS); the binary is built on Debian bookworm glibc, as the VHS image is.
demo: build
	$(STRATA) section examples/site.csv --from 40,470 --to 980,560 --truth 7 -o docs/section.svg
	$(STRATA) map examples/site.csv --line 40,470:980,560 -o docs/map.svg
	docker run --rm -u $$(id -u):$$(id -g) -e HOME=/tmp -v $(CURDIR):/vhs \
	  -v $(CURDIR)/target/release:/opt/strata:ro ghcr.io/charmbracelet/vhs:v0.10.0 docs/demo.tape

clean:
	$(CARGO) clean
