#!/bin/bash
# gene -> chromosome, TSS (UCSC hg38 ncbiRefSeqSelect) and cell lifespans (Sender and Milo 2021, milo-lab/cellular_turnover)
set -euo pipefail; source workflow/config.sh
curl -sL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz | gunzip -c \
  | awk -F'\t' '{tss = ($4=="+") ? $5+1 : $6; print $13"\t"$3"\t"tss}' | sort -u -k1,1 > resources/tss.tsv
curl -sL -o resources/turnover_Summary.xlsx https://raw.githubusercontent.com/milo-lab/cellular_turnover/HEAD/Summary.xlsx
python3 analysis/lifespan_from_xlsx.py resources/turnover_Summary.xlsx resources/lifespan.tsv
wc -l resources/tss.tsv resources/lifespan.tsv
