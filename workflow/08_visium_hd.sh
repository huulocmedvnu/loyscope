#!/bin/bash
# Visium HD sections (config/visium_hd.tsv): stream the binned outputs, keep the 8 um filtered matrix, run spaty on tiles
# of 16 x 16 bins; the archive is not stored (about 12 GB per section)
set -euo pipefail; source workflow/config.sh
O=$LOY_RESULTS/spatial; mkdir -p "$O"; U=https://cf.10xgenomics.com/samples/spatial-exp
while IFS=$'\t' read -r v n; do d=$LOY_DATA/spatial/$n; mkdir -p "$d"
  [ -s "$d/m.h5" ] || { curl -sSL --retry 5 "$U/$v/$n/${n}_binned_outputs.tar.gz" | tar xzf - -C "$d" --wildcards "*square_008um/filtered_feature_bc_matrix.h5"
    mv "$(find "$d" -name filtered_feature_bc_matrix.h5 | head -1)" "$d/m.h5"; }
  "$BIN/spaty" "$d/m.h5" "$TSS" 16 "$O/$n.t16" > "$LOY_RESULTS/logs/spaty_$n.log" 2>&1
  [ "${KEEP_H5AD:-0}" = 1 ] || rm -rf "$d"
done < config/visium_hd.tsv
