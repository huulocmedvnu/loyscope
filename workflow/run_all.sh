#!/bin/bash
# full reproduction: build, resources, all datasets, robustness, tables
set -euo pipefail
cargo build --release
for s in 00_resources 01_cellxgene 02_onek1k_geo 03_kidney_geo 04_seaad_atac 05_robustness 06_tables; do echo "== $s"; bash workflow/$s.sh; done
