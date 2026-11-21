#!/bin/bash
# genome-wide expression of cells without versus with Y reads, matched on library size deciles (Immune Aging cohort),
# with the same contrast for an autosomal control split; requires the step 01 extraction of iaging
set -euo pipefail; source workflow/config.sh
u=$(awk -F'\t' '$1=="iaging"{print $4}' config/datasets.tsv); f=$LOY_DATA/iaging.h5ad
[ -s "$f" ] || workflow/fastdl.sh "$u" "$f" 8
E=$LOY_RESULTS/extract
$BIN/loydepth "$f" $E/iaging $E/iaging.fitamb.tsv $E/iaging
[ "${KEEP_H5AD:-0}" = 1 ] || rm -f "$f"
