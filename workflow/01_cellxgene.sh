#!/bin/bash
# per CELLxGENE dataset: download -> loyx (per-cell Y, pseudo-Y, X/Y panel, obs, pseudobulk) -> loyfit (with and without
# ambient term, all disease labels) -> loycell -> loybench -> loyde -> delete the h5ad. Datasets run JOBS at a time.
set -euo pipefail; source workflow/config.sh
one() { name=$1; url=$2; f=$LOY_DATA/$name.h5ad; o=$LOY_RESULTS/extract/$name; log=$LOY_RESULTS/logs/$name.log
  [ -s $o.fitamb.tsv ] && return 0
  workflow/fastdl.sh "$url" "$f" 8 || { echo "download failed: $name" >> $LOY_RESULTS/logs/status.log; return 0; }
  { $BIN/loyx "$f" "$TSS" "$o" && LOY_ALL_DISEASE=1 $BIN/loyfit "$o" && LOY_ALL_DISEASE=1 LOY_AMBIENT=1 $BIN/loyfit "$o" \
    && $BIN/loycell "$o" && $BIN/loybench "$o" && $BIN/loyde "$o"; } > $log 2>&1 && echo "ok $name" >> $LOY_RESULTS/logs/status.log || echo "failed $name" >> $LOY_RESULTS/logs/status.log
  [ "${KEEP_H5AD:-0}" = 1 ] || rm -f "$f"; }
export -f one
awk -F'\t' 'NR>1{print $1"\t"$4}' config/datasets.tsv | xargs -P $JOBS -L 1 -d '\n' bash -c 'IFS=$'"'"'\t'"'"' read n u <<< "$0"; one "$n" "$u"'
