#!/bin/bash
# offline re-analyses from saved extraction outputs: library-level grouping, RPS4Y1-only Y set
set -euo pipefail; source workflow/config.sh; E=$LOY_RESULTS/extract
for n in aida retina; do $BIN/loyregroup $E/$n $E/${n}_lib library_id - && $BIN/loyfit $E/${n}_lib; done
$BIN/loyregroup $E/retina $E/retina_rps4y1 - RPS4Y1 && $BIN/loyfit $E/retina_rps4y1
