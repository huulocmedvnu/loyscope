#!/bin/bash
# published kidney LOY study: GEO RNA matrices (GSE151302, GSE195460, GSE220289); density-trough call vs ZINB
set -euo pipefail; source workflow/config.sh
G=$LOY_DATA/gb2024; mkdir -p $G
while IFS=$'\t' read g f; do gsm=${f%%_*}; [ -s $G/$f ] || curl -sS --retry 5 -o $G/$f https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/$gsm/suppl/$f; done < config/gb2024_rna.tsv
awk -F'\t' -v G=$G '{split($2,a,"_"); print $1"_"a[2]"\t"G"/"$2}' config/gb2024_rna.tsv > $G/manifest.tsv
o=$LOY_RESULTS/extract/gb2024
$BIN/tenx $G/manifest.tsv "$TSS" $o && LOY_ALL_DISEASE=1 $BIN/loyfit $o && LOY_ALL_DISEASE=1 LOY_AMBIENT=1 $BIN/loyfit $o
# chromatin of the same study: ATAC fragments of donor DN4 (male) and Control6 (female background), nuclei matched to the
# peak matrix of each library; fraction of nuclei without Y by the detection model on MSY fragment counts
for p in GSM5837798_DN4 GSM5837793_Control6; do gsm=${p%%_*}
  for s in fragments.tsv.gz filtered_peak_bc_matrix.h5; do [ -s $G/${p}_$s ] || curl -sS --retry 5 -o $G/${p}_$s https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/$gsm/suppl/${p}_$s; done
  $BIN/fragy $G/${p}_fragments.tsv.gz $G/${p}_filtered_peak_bc_matrix.h5 $LOY_RESULTS/extract/gb2024_${p#*_}.fragy.tsv; done
$BIN/atacmix $LOY_RESULTS/extract/gb2024_Control6.fragy.tsv $LOY_RESULTS/extract/gb2024_DN4.fragy.tsv $LOY_RESULTS/extract/gb2024_DN4.atacmix.json
