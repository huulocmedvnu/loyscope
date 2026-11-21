#!/bin/bash
# SEA-AD MTG multiome chromatin: per-nucleus chrY (MSY) / chrX / total ATAC counts for joining with RNA calls
set -euo pipefail; source workflow/config.sh
A=$LOY_DATA/seaad; mkdir -p $A
workflow/fastdl.sh https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad $A/atac.h5ad 8
$BIN/atacy $A/atac.h5ad $LOY_RESULTS/extract/seaad_atac && rm -f $A/atac.h5ad
