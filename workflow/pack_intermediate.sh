#!/bin/bash
# pack the intermediate outputs needed by `make figures` (run after steps 1 to 9); writes loyscope_intermediate.tar.gz
set -euo pipefail; source workflow/config.sh
cd "$LOY_RESULTS"
find -L extract -type f \( -name '*.fitamb*.tsv' -o -name '*.fit.tsv' -o -name '*.fitrel.tsv' -o -name '*.de*.tsv' -o -name '*.cellde.tsv' \
  -o -name '*.bench.tsv' -o -name '*.rf_info.txt' -o -name '*.trough.tsv' -o -name '*.groups.tsv' -o -name '*.info.txt' \
  -o -name '*.groups_cmp.tsv' -o -name '*.fragy.tsv' -o -name '*.atacmix.json' -o -name '*.multi.*' \) > ../pack.list
find -L encode spatial -type f 2>/dev/null >> ../pack.list || true
tar czhf ../loyscope_intermediate.tar.gz -T ../pack.list && cd .. && sha256sum loyscope_intermediate.tar.gz 2>/dev/null > config_intermediate.sha256 || shasum -a 256 loyscope_intermediate.tar.gz > config_intermediate.sha256
ls -la loyscope_intermediate.tar.gz
