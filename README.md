# loyscope

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24431454.svg)](https://doi.org/10.5281/zenodo.24431454)

Detection-aware estimation of mosaic loss of chromosome Y (LOY) from single-cell and single-nucleus RNA-seq.

A cell without reads from the male-specific region of chromosome Y is not necessarily a cell without chromosome Y:
at typical depths, intact cells frequently yield zero Y reads. loyscope estimates the fraction of cells that have lost Y
in each donor and cell type with a mixture model that accounts for sequencing depth and ambient RNA, and calibrates it
against expression-matched autosomal gene sets that cannot be lost (pseudo-Y) and against female cells.

## Model

For cell *i* with total UMI count *N_i* and Y-gene UMI count *y_i*, within a donor x cell type group:

    y_i ~ f * Poisson(rho * mu * N_i) + (1 - f) * NegBin(mean = mu * N_i, size = k)

*f* is the LOY fraction, *mu* the Y expression rate of intact cells, *k* the dispersion (selected by likelihood), and
*rho* the ambient Y level relative to intact cells, estimated from female cells of the same dataset (rho = 0 gives the
zero-inflated form). Parameters are fitted by EM; confidence intervals come from bootstrap resampling of cells.
Twenty pseudo-Y sets are fitted in the same way and must return f close to 0.

## Tools (Rust, `src/bin`)

| binary | input | purpose |
|---|---|---|
| `loyx` | CELLxGENE h5ad | streams the sparse matrix once; per-cell Y and pseudo-Y counts, X/Y gene panel, obs, pseudobulk |
| `geox` | GEO dense CSV pools | same records from author count matrices (OneK1K) |
| `tenx` | 10x `filtered_feature_bc_matrix.h5` | same records, sex inference, marker-based cell classes, density-trough comparator |
| `loyfit` | extraction outputs | ZINB / ambient-aware fits, bootstrap CI, pseudo-Y controls |
| `loybench` | extraction outputs | Random Forest caller (as in Nature 2025), spike-in and depth-thinning benchmarks |
| `loycell` | extraction outputs | per-cell posteriors; PAR dosage and X/Y gametolog contrasts with a pseudo-zero control |
| `loyde` | pseudobulk | genome-wide LOY vs Y+ contrasts in groups where zero-Y cells are mostly true LOY |
| `loyregroup` | extraction outputs | regroup by any obs column or recompute Y from a gene subset, without the h5ad |
| `atacy` | ATAC h5ad | per-nucleus chrY (outside PAR) / chrX / total chromatin counts |
| `obsdump` | h5ad | obs table export |
| `fragy` | ATAC fragment files | per-barcode total, X and MSY fragment counts, streamed |
| `encmulti` | ENCODE 10x multiome | RNA records for `loyfit` and chromatin counts (MSY, X, 20 matched autosomal control bin sets) of the same nuclei |
| `multijoin` | `encmulti` outputs | joins RNA and chromatin of each nucleus; per-group chromatin and RNA estimates and per-nucleus posteriors |
| `atacmix` | `fragy` outputs | fraction of nuclei without Y from MSY fragment counts, female background fixed |
| `loydepth` | h5ad + extraction | genome-wide contrast of cells without and with Y reads within library-size deciles, with an autosomal control split |
| `spaty` | Visium HD 8 um bins | tile-level test for a deficit of MSY reads against X gametologs, with matched autosomal controls |

Matrices are read in row blocks through HDF5 hyperslabs with a prefetching reader thread; memory stays near one block
plus per-group summaries.

## Reproducing the study

Requirements: Linux or macOS, Rust >= 1.75, HDF5 >= 1.10 (`apt install libhdf5-dev`, or `brew install hdf5` and
`export HDF5_DIR=$(brew --prefix hdf5)` on macOS), curl, librsvg (`rsvg-convert`, for the schematic figures), Python >= 3.10 with the packages in `requirements.txt`,
about 200 GB of free disk at peak (downloads are deleted after extraction).

```
cargo build --release
cargo test --release                # estimator recovery on simulated counts
bash workflow/run_all.sh          # or run the numbered steps in workflow/ one by one
```

Steps (each writes to `results/`, set in `workflow/config.sh`):

| step | script | what it does |
|---|---|---|
| 0 | `00_resources.sh` | gene coordinates (UCSC hg38 RefSeq Select) and cell lifespans (Sender and Milo 2021) |
| 1 | `01_cellxgene.sh` | every dataset in `config/datasets.tsv`: extraction, fits, PAR1 dosage, benchmarks, pseudobulk contrasts |
| 2 | `02_onek1k_geo.sh` | OneK1K author count matrices (GEO GSE196830) |
| 3 | `03_kidney_geo.sh` | kidney RNA (GSE151302, GSE195460, GSE220289), density-trough comparator, ATAC fragments of donor DN4 |
| 4 | `04_seaad_atac.sh` | middle temporal gyrus chromatin peaks joined to the RNA of the same nuclei |
| 5 | `05_robustness.sh` | grouping by library, RPS4Y1-only Y set, ambient contamination x3 and x10 |
| 7 | `07_encode_multiome.sh` | ENCODE heart multiome (`config/encode_manifest.tsv`), RNA and chromatin of the same nuclei |
| 8 | `08_visium_hd.sh` | Visium HD sections (`config/visium_hd.tsv`) |
| 9 | `09_depth_de.sh` | library-size-matched genome-wide contrast of cells without and with Y reads (Immune Aging) |
| 6 | `06_tables.sh` | summary tables and all statistics (`analysis/robust_stats.py` and the scripts in `analysis/`) |
| 10 | `10_figures.sh` | figures, graphical abstract, flowchart and `results/tables/numbers.json` |

`run_all.sh` runs them in this order. Steps 1 to 9 need the network and about 200 GB of temporary disk; downloads
are deleted after extraction unless `KEEP_H5AD=1`. Steps 5, 6 and 10 run from saved outputs in minutes.

### Where each result comes from

| manuscript item | produced by |
|---|---|
| Figure 1 (study design) | `analysis/figures/make_flowchart.py` |
| Figure 2 (measurement) | `analysis/figures/fig_measure.py` (pseudo-Y, thinning, calibration: `loyfit`, `loybench`) |
| Figure 3 (loss of the chromosome) | `fig3.py` (`loycell`; kidney chromatin `fragy` + `atacmix`, `analysis/dn4_estimates.py`) |
| Figure 4 (where loss occurs) | `fig5.py` (`analysis/robust_stats.py`: REML and Hartung-Knapp meta-analysis) |
| Figure 5 (systemic or local) | `fig6.py` (`analysis/ts_donor.py`, `turnover.py`, `bio345.py`) |
| Figure 6 (cells without Y) | `fig7.py` (`loyde`, `loydepth`) |
| Figure 7 (published estimates) | `fig4.py` (`resources/literature_search.tsv`, `robust_stats.py` disease contrasts) |
| Supplementary Figures 1 and 2 | `fig1.py`, `fig2.py` (`loybench`) |
| Supplementary Figure 3 | `figS3_spatial.py` (`spaty`) |
| Supplementary Table 2 | `resources/literature_search.tsv` (Europe PMC search, screened by hand) |
| Supplementary Table 3 | `analysis/depth_retention.py` |
| Supplementary Tables 4 and 5 | `analysis/robust_stats.py` (`results/tables/robust_*.tsv`) |
| Numbers in the text | `analysis/figures/manuscript_numbers.py` and `results/tables/*.txt` |

## Data

All inputs are public. `config/datasets.tsv` lists CELLxGENE dataset identifiers, download URLs and collection DOIs;
`config/geo.tsv` lists GEO and SEA-AD sources. Only cells with raw integer counts from 10x assays are analysed.

## Citation

See `CITATION.cff`.

## License

MIT (code). Input data remain under the licences of their providers.
