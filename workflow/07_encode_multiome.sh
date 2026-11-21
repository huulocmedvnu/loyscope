#!/bin/bash
# ENCODE 10x multiome, male donors: RNA LOY (loyfit) vs chromatin LOY (multijoin) in the same nuclei. Stream, extract, delete.
cd "${LOY_ROOT:-.}"; mkdir -p encode logs
T="${TSS:-data/tss.tsv}"; WA=geo/wl/atac_737K-arc-v1.txt.gz; WG=geo/wl/gex_737K-arc-v1.txt.gz
echo "ENCODE START $(date)" >> logs/run.log
# 1. female background (chrY MSY fragments per fragment in female heart nuclei)
while IFS=$'\t' read acc url; do d=encode/fem_$acc; mkdir -p $d
  [ -s encode/fem_$acc.tsv ] || { curl -sSL --retry 5 -o $d/f.tar.gz "$url" && tar xzf $d/f.tar.gz -C $d && target/release/encmulti --female $(find $d -name fragments.tsv.gz | head -1) encode/fem_$acc.tsv; rm -rf $d; }
done < meta/encode_female.tsv
BG=$(cat encode/fem_*.tsv | awk -F"\t" "\$1!=\"barcode\"{t+=\$2; y+=\$3} END{printf \"%.4e\", y/t}"); echo "ATAC_BG $BG" >> logs/run.log
# 2. male multiome libraries, 3 at a time
one() { acc=$1; donor=$2; age=$3; stage=$4; tissue=$5; furl=$6; rurl=$7; d=encode/w_$acc; o=encode/$acc
  [ -s $o.multi.groups_cmp.tsv ] && return 0
  mkdir -p $d && curl -sSL --retry 5 -o $d/f.tar.gz "$furl" && curl -sSL --retry 5 -o $d/r.tar.gz "$rurl" && tar xzf $d/f.tar.gz -C $d && tar xzf $d/r.tar.gz -C $d || { echo "ENC DLFAIL $acc" >> logs/run.log; rm -rf $d; return 0; }
  st=$(echo "$stage" | awk "{if(\$2==\"year\") print \$1\"-year-old stage\"; else print \"embryonic\"}")
  ts=$(echo "$tissue" | sed "s/Homo sapiens //; s/ tissue.*//; s/ nuclear.*//" | tr " " "_")
  { target/release/encmulti $(find $d -name fragments.tsv.gz | head -1) $(dirname $(find $d -path "*filtered/matrix.mtx" | head -1)) $WA $WG $T $acc "$ts" "$st" $o \
    && LOY_AMBIENT=1 target/release/loyfit $o 30 \
    && ATAC_BG=$BG RNA_BG=0 target/release/multijoin $o $o.atac.tsv.gz $o.multi 30; } > logs/enc_$acc.log 2>&1 && echo "ENC ok $acc" >> logs/run.log || echo "ENC FAIL $acc" >> logs/run.log
  rm -rf $d; }
export -f one; export BG T WA WG
while IFS=$'\t' read -r a b c d e f g; do
  RAYON_NUM_THREADS=8 one "$a" "$b" "$c" "$d" "$e" "$f" "$g" &
  while [ $(jobs -rp | wc -l) -ge 3 ]; do sleep 5; done
done < meta/encode_manifest.tsv
wait
echo "ENCODE ALL DONE $(date)" >> logs/run.log
