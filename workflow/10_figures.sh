#!/bin/bash
# the numbers quoted in the manuscript (numbers.json, then text_numbers.json), all figures, the graphical abstract and the study flowchart
set -euo pipefail; source workflow/config.sh
cd analysis/figures
export LOY_RESULTS LOY_SPATIAL=$LOY_RESULTS/spatial
python3 manuscript_numbers.py > $LOY_RESULTS/tables/numbers.txt
python3 ../text_numbers.py > $LOY_RESULTS/tables/text_numbers.log
for s in fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig_measure figS3_spatial make_flowchart make_graphical_abstract; do python3 $s.py; done
