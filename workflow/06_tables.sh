#!/bin/bash
# summary tables (Python reads only numbers written by the Rust tools)
set -euo pipefail; source workflow/config.sh; E=$LOY_RESULTS/extract; T=$LOY_RESULTS/tables
python3 analysis/bench_summ.py $E/*.bench.tsv > $T/benchmark.txt
python3 analysis/cell_summ.py $E/*.cellde.tsv > $T/par_dosage.txt
python3 analysis/consolidate.py && python3 analysis/turnover.py > $T/turnover.txt && python3 analysis/ts_donor.py > $T/within_donor.txt && python3 analysis/bio345.py > $T/clonality_selection.txt
for d in "kidney:chronic kidney disease" "ckd_sn:chronic kidney disease" "ckd_sn:acute kidney injury" "ad_microglia:Alzheimer disease" "ad_mssm:Alzheimer disease" "ad_radc:Alzheimer disease" "dcm:dilated cardiomyopathy" "covid:COVID-19"; do
  n=${d%%:*}; dis=${d#*:}; python3 analysis/dis_summ.py $E/$n.fitamb.tsv "$dis" 5; done > $T/disease_contrasts.txt
# robustness of age slopes, meta-analysis, disease contrasts, threshold sensitivity, overdispersion
LOY_RESULTS=$LOY_RESULTS python3 analysis/robust_stats.py > $T/robust_stats.txt
python3 analysis/dn4_estimates.py > $T/dn4.txt
python3 analysis/depth_retention.py > $T/depth_retention.txt
[ -s $E/encode_heart.fitamb.tsv ] && LOY_RESULTS=$LOY_RESULTS python3 analysis/encode_dosage.py "$LOY_RESULTS" > $T/encode_dosage.txt || true
