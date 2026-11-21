#!/bin/bash
# OneK1K author matrices from GEO GSE196830 (75 pools) joined to CELLxGENE labels (barcode core + donor)
set -euo pipefail; source workflow/config.sh
G=$LOY_DATA/onek1k_geo; mkdir -p $G
while read f; do gsm=${f%%_*}; [ -s $G/$f ] || curl -sS --retry 5 -o $G/$f https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/$gsm/suppl/$f; done < config/onek1k_geo_files.txt
u=$(awk -F'\t' '$1=="onek1k"{print $4}' config/datasets.tsv); workflow/fastdl.sh "$u" $G/onek1k.h5ad 8
$BIN/obsdump $G/onek1k.h5ad $G/obs.tsv && rm -f $G/onek1k.h5ad
for f in $G/*RawCounts.csv.gz; do b=$(basename $f); gsm=${b%%_*}; pool=$(echo $b | sed -E 's/.*_Pool([0-9]+)_.*/\1/'); echo -e "$pool\t$f\t$(ls $G/${gsm}_*Individual_Barcodes.csv.gz)"; done > $G/manifest.tsv
o=$LOY_RESULTS/extract/onek1k_geo
$BIN/geox $G/manifest.tsv $G/obs.tsv "$TSS" $o && $BIN/loyfit $o && LOY_AMBIENT=1 $BIN/loyfit $o
