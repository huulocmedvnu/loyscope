# loyscope

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

Matrices are read in row blocks through HDF5 hyperslabs with a prefetching reader thread; memory stays near one block
plus per-group summaries.

## Reproducing the study

Requirements: Linux, Rust >= 1.75, HDF5 >= 1.10, curl, Python >= 3.10 with numpy and openpyxl (summary tables only),
about 200 GB of free disk at peak (downloads are deleted after extraction).

```
cargo build --release
bash workflow/run_all.sh          # or run the numbered steps in workflow/ one by one
```

Steps:

1. `00_resources.sh`  gene coordinates (UCSC hg38 RefSeq Select) and cell lifespans (Sender and Milo 2021)
2. `01_cellxgene.sh`  all datasets in `config/datasets.tsv`
3. `02_onek1k_geo.sh` OneK1K author matrices (GEO GSE196830)
4. `03_kidney_geo.sh` kidney matrices of a published LOY study (GSE151302, GSE195460, GSE220289)
5. `04_seaad_atac.sh` SEA-AD MTG chromatin accessibility
6. `05_robustness.sh` library-level grouping and single-gene Y set
7. `06_tables.sh`     summary tables in `results/tables`

Settings (threads, parallel datasets, paths) are in `workflow/config.sh` and can be overridden by environment variables.

## Data

All inputs are public. `config/datasets.tsv` lists CELLxGENE dataset identifiers, download URLs and collection DOIs;
`config/geo.tsv` lists GEO and SEA-AD sources. Only cells with raw integer counts from 10x assays are analysed.

## Citation

See `CITATION.cff`.

## License

MIT (code). Input data remain under the licences of their providers.
