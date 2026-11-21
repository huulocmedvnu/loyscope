#!/bin/bash
# offline re-analyses from saved extraction outputs: library-level grouping, RPS4Y1-only Y set, and ambient
# contamination fractions three and ten times the estimate from female cells
set -euo pipefail; source workflow/config.sh; E=$LOY_RESULTS/extract
for n in aida retina; do $BIN/loyregroup $E/$n $E/${n}_lib library_id - && $BIN/loyfit $E/${n}_lib; done
$BIN/loyregroup $E/retina $E/retina_rps4y1 - RPS4Y1 && $BIN/loyfit $E/retina_rps4y1
for p in $E/*.cells.bin; do o=${p%.cells.bin}; case $(basename $o) in *_lib|retina_rps4y1|onek1k) continue;; esac
  for s in 3 10; do LOY_ALL_DISEASE=1 LOY_AMBIENT=1 LOY_BG_SCALE=$s $BIN/loyfit $o 50 1; done; done
