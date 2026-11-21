#!/bin/bash
# parallel ranged download: fastdl.sh <url> <out> [parts=8]; falls back to single curl if size unknown
url=$1; out=$2; n=${3:-8}
size=$(curl -sIL "$url" | grep -i '^content-length' | tail -1 | tr -dc '0-9')
if [ -z "$size" ] || [ "$size" -lt 100000000 ]; then curl -sSL --retry 5 -o "$out" "$url"; exit $?; fi
part=$(( (size + n - 1) / n )); pids=()
for i in $(seq 0 $((n-1))); do s=$((i*part)); e=$((s+part-1)); [ $e -ge $size ] && e=$((size-1))
  ( for t in 1 2 3 4 5; do curl -sSL -r $s-$e -o "$out.part$i" "$url" && [ $(stat -c %s "$out.part$i") -eq $((e-s+1)) ] && exit 0; sleep 5; done; exit 1 ) & pids+=($!); done
ok=0; for p in "${pids[@]}"; do wait $p || ok=1; done
[ $ok -ne 0 ] && { echo "fastdl failed $url" >&2; exit 1; }
cat $(for i in $(seq 0 $((n-1))); do echo "$out.part$i"; done) > "$out" && rm -f "$out".part*
[ $(stat -c %s "$out") -eq $size ]
