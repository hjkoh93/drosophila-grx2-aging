# Age-Dependent Gene Expression Analysis in *Drosophila* Brain

Single-cell RNA-seq analysis of **Grx2** (Glutaredoxin 1, gene symbol: *Grx1*) expression changes with aging in *Drosophila melanogaster* brain.

## Datasets

| Dataset | Source | Cells | Ages (days) | Tissue |
|---------|--------|-------|-------------|--------|
| Davie et al. (2018) | [GSE107451](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107451) | 56,902 | 0, 1, 3, 6, 9, 15, 30, 50 | Whole brain |
| Lu et al. (2023) | [AFCA](https://www.flycellatlas.org/) | 289,981 | 5, 30, 50, 70 | Head |

## Repository Structure

```
├── scripts/
│   ├── 01_extract_and_normalize.py   # Grx1 extraction & CPM normalization
│   ├── 02_statistical_analysis.py     # Spearman, Kruskal-Wallis, Mann-Whitney U
│   └── 03_visualization.py            # Trajectory, heatmaps, dopaminergic analysis
├── requirements.txt
└── README.md
```

## Data Preparation

### Davie et al.
1. Download from GEO ([GSE107451](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107451)):
   - `GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv.tar`
   - `GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv`
2. Extract the tar file to obtain `genes.tsv`, `barcodes.tsv`, and `matrix.mtx`

### Lu et al.
1. Download `adata_head_S_v1.0.h5ad` from the [Fly Cell Atlas](https://www.flycellatlas.org/)

### Directory layout
Place downloaded files in a `data/` directory:
```
data/
├── davie/
│   ├── genes.tsv
│   ├── barcodes.tsv
│   ├── matrix.mtx
│   └── metadata.tsv
└── lu/
    └── adata_head_S_v1.0.h5ad
```

## Usage

```bash
# Install dependencies
pip install -r requirements.txt

# Step 1: Extract Grx1 expression and normalize to CPM
python scripts/01_extract_and_normalize.py

# Step 2: Run statistical analyses
python scripts/02_statistical_analysis.py

# Step 3: Generate figures
python scripts/03_visualization.py
```

All output figures are saved to the `output/` directory.

## Statistical Methods

- **Spearman's rank correlation**: age–expression association (robust to zero-inflated data)
- **Kruskal–Wallis H test**: multi-group comparison across ages
- **Mann–Whitney U test**: pairwise comparisons with Bonferroni correction
- **Effect sizes**: rank-biserial correlation (*r*) and Cohen's *d*
- **Two-part (hurdle) model**: (1) logistic regression on detection probability, (2) Spearman correlation among expressing cells only

## References

- Davie K, et al. A single-cell transcriptome atlas of the aging *Drosophila* brain. *Cell*. 2018;174(4):982-998.
- Lu T-C, et al. Aging Fly Cell Atlas identifies exhaustive aging features at cellular resolution. *Science*. 2023;380(6650):eadg0934.
