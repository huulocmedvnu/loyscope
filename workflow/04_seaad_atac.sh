#!/bin/bash
# SEA-AD MTG multiome chromatin: per-nucleus MSY, chrX and total ATAC counts (with matched autosomal peak sets), joined to RNA
set -euo pipefail; source workflow/config.sh
A=$LOY_DATA/seaad; mkdir -p $A
workflow/fastdl.sh https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad $A/atac.h5ad 8
$BIN/atacy $A/atac.h5ad $LOY_RESULTS/extract/seaad_atac && rm -f $A/atac.h5ad
# join with the RNA of the same nuclei (middle temporal gyrus RNA extracted in step 01 as seaad_mtg)
$BIN/multijoin $LOY_RESULTS/extract/seaad_mtg $LOY_RESULTS/extract/seaad_atac.atac.tsv.gz $LOY_RESULTS/extract/seaad_multi
