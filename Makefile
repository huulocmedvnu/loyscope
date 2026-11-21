# loyscope: one-command reproduction
#   make test        build the Rust tools and run the estimator tests (seconds)
#   make figures     download the archived intermediate outputs, then recompute every statistic, table and figure (minutes)
#   make reproduce   full analysis from the public raw data: download, extract, fit, statistics, figures
#                    (Linux, about 200 GB of temporary disk, many hours; see workflow/config.sh)
SHELL := /bin/bash
RESULTS ?= $(CURDIR)/results
export LOY_RESULTS := $(RESULTS)

.PHONY: build test intermediate figures reproduce clean
build:
	cargo build --release
test: build
	cargo test --release
intermediate:
	bash workflow/fetch_intermediate.sh "$(RESULTS)"
figures: build intermediate
	bash workflow/06_tables.sh
	bash workflow/10_figures.sh
	@echo "figures in $(RESULTS)/figures, tables in $(RESULTS)/tables"
reproduce: build
	bash workflow/run_all.sh
clean:
	rm -rf $(RESULTS)/figures $(RESULTS)/tables
