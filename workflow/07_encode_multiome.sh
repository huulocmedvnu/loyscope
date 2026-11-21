#!/bin/bash
# ENCODE 10x multiome libraries (config/encode_manifest.tsv): RNA LOY estimates and, with CHROMATIN=1 (default),
# chromatin counts of the same nuclei joined by barcode. Archives are streamed to LOY_DATA and deleted after extraction.
set -euo pipefail; source workflow/config.sh
O=$LOY_RESULTS/encode; W=$LOY_DATA/encode; mkdir -p "$O" "$W" "$LOY_RESULTS/logs"
CHROMATIN=${CHROMATIN:-1}
WL=resources/whitelist; mkdir -p $WL
for s in atac gex; do [ -s $WL/${s}_737K-arc-v1.txt.gz ] || curl -sSL -o $WL/${s}_737K-arc-v1.txt.gz \
  https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/${s}_737K-arc-v1.txt.gz; done
EMPTY=$W/empty_fragments.tsv.gz; [ -s "$EMPTY" ] || printf "" | gzip > "$EMPTY"
# female background: MSY fragments per fragment in female nuclei
if [ "$CHROMATIN" = 1 ]; then
  while IFS=$'\t' read -r acc url; do d=$W/fem_$acc; mkdir -p "$d"
    [ -s "$O/fem_$acc.tsv" ] || { curl -sSL --retry 5 "$url" | tar xzf - -C "$d"
      "$BIN/encmulti" --female "$(find "$d" -name fragments.tsv.gz | head -1)" "$O/fem_$acc.tsv"; rm -rf "$d"; }
  done < config/encode_female.tsv
  BG=$(cat "$O"/fem_*.tsv | awk -F'\t' '$1!="barcode"{t+=$2; y+=$3} END{printf "%.4e", y/t}')
fi
one() { acc=$1; stage=$4; tissue=$5; furl=$6; rurl=$7; d=$W/w_$acc; o=$O/$acc
  [ -s "$o.fitamb.tsv" ] && return 0
  mkdir -p "$d"; curl -sSL --retry 5 "$rurl" | tar xzf - -C "$d"
  frag=$EMPTY; if [ "$CHROMATIN" = 1 ]; then curl -sSL --retry 5 "$furl" | tar xzf - -C "$d"; frag=$(find "$d" -name fragments.tsv.gz | head -1); fi
  st=$(echo "$stage" | awk '{if($2=="year") print $1"-year-old stage"; else print "embryonic"}')
  ts=$(echo "$tissue" | sed 's/Homo sapiens //; s/ tissue.*//; s/ nuclear.*//' | tr ' ' '_')
  "$BIN/encmulti" "$frag" "$(dirname "$(find "$d" -path '*filtered/matrix.mtx' | head -1)")" $WL/atac_737K-arc-v1.txt.gz $WL/gex_737K-arc-v1.txt.gz "$TSS" "$acc" "$ts" "$st" "$o"
  LOY_AMBIENT=1 "$BIN/loyfit" "$o" 30
  [ "$CHROMATIN" = 1 ] && ATAC_BG=$BG RNA_BG=0 "$BIN/multijoin" "$o" "$o.atac.tsv.gz" "$o.multi" 30
  rm -rf "$d"; }
export -f one; export O W WL EMPTY CHROMATIN BIN TSS; export BG=${BG:-0}
while IFS=$'\t' read -r a b c d e f g; do one "$a" "$b" "$c" "$d" "$e" "$f" "$g" > "$LOY_RESULTS/logs/encode_$a.log" 2>&1 &
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done; done < config/encode_manifest.tsv
wait
# pooled heart RNA table and chromatin dosage
python3 analysis/encode_pool.py "$O" config/encode_manifest.tsv "$LOY_RESULTS/extract/encode_heart.fitamb.tsv"
[ "$CHROMATIN" = 1 ] && python3 analysis/encode_dosage.py "$LOY_RESULTS" || true
