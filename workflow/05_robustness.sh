#!/bin/bash
# offline re-analyses from saved extraction outputs: library-level grouping and ambient
# contamination fractions three and ten times the estimate from female cells
set -euo pipefail; source workflow/config.sh; E=$LOY_RESULTS/extract
for n in retina; do $BIN/loyregroup $E/$n $E/${n}_lib library_id - && LOY_ALL_DISEASE=1 LOY_AMBIENT=1 $BIN/loyfit $E/${n}_lib; done
for p in $E/*.cells.bin; do o=${p%.cells.bin}; case $(basename $o) in *_lib|onek1k) continue;; esac
  for s in 3 10; do LOY_ALL_DISEASE=1 LOY_AMBIENT=1 LOY_BG_SCALE=$s $BIN/loyfit $o 50 1; done; done
