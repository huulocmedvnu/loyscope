#!/bin/bash
# remaining ENCODE heart multiome libraries: RNA only (chromatin dosage not informative at this depth)
cd "${LOY_ROOT:-.}"
T="${TSS:-data/tss.tsv}"; WA=geo/wl/atac_737K-arc-v1.txt.gz; WG=geo/wl/gex_737K-arc-v1.txt.gz; E=encode/empty_fragments.tsv.gz
[ -s $E ] || printf "" | gzip > $E
rm -rf encode/w_*
echo "ENCODE RNA START $(date)" >> logs/run.log
one() { acc=$1; stage=$4; tissue=$5; rurl=$7; d=encode/w_$acc; o=encode/$acc
  [ -s $o.multi.groups_cmp.tsv ] && return 0; [ -s $o.fitamb.tsv ] && return 0
  mkdir -p $d && curl -sSL --retry 5 -o $d/r.tar.gz "$rurl" && tar xzf $d/r.tar.gz -C $d || { echo "ENC DLFAIL $acc" >> logs/run.log; rm -rf $d; return 0; }
  st=$(echo "$stage" | awk '{if($2=="year") print $1"-year-old stage"; else print "embryonic"}')
  ts=$(echo "$tissue" | sed "s/Homo sapiens //; s/ tissue.*//; s/ nuclear.*//" | tr " " "_")
  { target/release/encmulti $E $(dirname $(find $d -path "*filtered/matrix.mtx" | head -1)) $WA $WG $T $acc "$ts" "$st" $o && LOY_AMBIENT=1 target/release/loyfit $o 30; } > logs/encr_$acc.log 2>&1 && echo "ENCR ok $acc" >> logs/run.log || echo "ENCR FAIL $acc" >> logs/run.log
  rm -rf $d; }
export -f one; export T WA WG E
while IFS=$'\t' read -r a b c d e f g; do RAYON_NUM_THREADS=8 one "$a" "$b" "$c" "$d" "$e" "$f" "$g" & while [ $(jobs -rp | wc -l) -ge 3 ]; do sleep 5; done; done < meta/encode_manifest.tsv
wait; echo "ENCODE ALL DONE $(date)" >> logs/run.log
