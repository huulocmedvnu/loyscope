#!/bin/bash
# all figures, the graphical abstract, the study flowchart and the numbers quoted in the manuscript
set -euo pipefail; source workflow/config.sh
cd analysis/figures
export LOY_RESULTS LOY_SPATIAL=$LOY_RESULTS/spatial
python3 manuscript_numbers.py > $LOY_RESULTS/tables/numbers.txt
for s in fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig_measure figS3_spatial make_flowchart make_graphical_abstract; do python3 $s.py; done
