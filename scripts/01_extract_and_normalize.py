"""
01_extract_and_normalize.py
===========================
Extract Grx1 (Grx2) expression from two Drosophila brain scRNA-seq
datasets and normalize to CPM.

Input:
  - data/davie/genes.tsv, barcodes.tsv, matrix.mtx, metadata.tsv
  - data/lu/adata_head_S_v1.0.h5ad

Output:
  - output/davie_grx1.csv
  - output/lu_grx1.csv
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

# --------------- Paths (adjust as needed) ---------------
DAVIE_DIR = 'data/davie'
LU_H5AD   = 'data/lu/adata_head_S_v1.0.h5ad'
OUT_DIR    = 'output'
os.makedirs(OUT_DIR, exist_ok=True)

GENES_OF_INTEREST = ['Grx1']


# ========================================================
# Davie et al. — Matrix Market extraction
# ========================================================
def extract_davie():
    """Extract genes from Davie sparse matrix by scanning row indices."""
    print('=' * 60)
    print('Davie et al. — Extracting gene expression')
    print('=' * 60)

    genes = pd.read_csv(
        os.path.join(DAVIE_DIR, 'genes.tsv'),
        sep='\t', header=None, names=['gene_id', 'gene_name']
    )
    barcodes = pd.read_csv(
        os.path.join(DAVIE_DIR, 'barcodes.tsv'),
        sep='\t', header=None, names=['barcode']
    )
    metadata = pd.read_csv(
        os.path.join(DAVIE_DIR, 'metadata.tsv'), sep='\t'
    )

    n_cells = len(barcodes)

    # Find row indices (1-based) for each gene
    gene_rows = {}
    for g in GENES_OF_INTEREST:
        matches = genes[genes['gene_name'] == g]
        if len(matches) == 0:
            print(f'  WARNING: {g} not found in gene list')
            continue
        gene_rows[g] = matches.index[0] + 1  # 1-based
        print(f'  {g}: row {gene_rows[g]}')

    # Initialize expression arrays
    expr = {g: np.zeros(n_cells, dtype=np.float32) for g in gene_rows}

    # Scan matrix file — only store entries for target genes
    mtx_path = os.path.join(DAVIE_DIR, 'matrix.mtx')
    print(f'  Scanning {mtx_path} ...')
    target_rows = set(gene_rows.values())
    row_to_gene = {v: k for k, v in gene_rows.items()}

    with open(mtx_path, 'r') as f:
        for line in f:
            if line.startswith('%'):
                continue
            parts = line.strip().split()
            if len(parts) == 3:
                r, c, v = int(parts[0]), int(parts[1]), int(parts[2])
                if r in target_rows:
                    expr[row_to_gene[r]][c - 1] = v

    # Merge with metadata and compute CPM
    meta_sub = metadata[['new_barcode', 'Age', 'Genotype', 'nUMI', 'annotation']].rename(
        columns={'new_barcode': 'barcode'}
    )

    for gene_name, values in expr.items():
        df = pd.DataFrame({'barcode': barcodes['barcode'].values, 'raw': values})
        df = df.merge(meta_sub, on='barcode', how='inner')
        df['cpm'] = (df['raw'] / df['nUMI']) * 1e6
        out_path = os.path.join(OUT_DIR, f'davie_{gene_name.lower()}.csv')
        df.to_csv(out_path, index=False)
        n_expr = (df['raw'] > 0).sum()
        print(f'  {gene_name}: {len(df)} cells, {n_expr} expressing ({n_expr/len(df)*100:.1f}%)')
        print(f'    Saved: {out_path}')

    print()


# ========================================================
# Lu et al. — h5ad extraction
# ========================================================
def extract_lu():
    """Extract genes from Lu AnnData file in backed/chunked mode."""
    print('=' * 60)
    print('Lu et al. — Extracting gene expression')
    print('=' * 60)

    adata = sc.read_h5ad(LU_H5AD, backed='r')
    n = adata.shape[0]
    var_names = list(adata.var_names)

    # Find gene indices
    gene_indices = {}
    for g in GENES_OF_INTEREST:
        if g in var_names:
            gene_indices[g] = var_names.index(g)
            print(f'  {g}: var index {gene_indices[g]}')
        else:
            print(f'  WARNING: {g} not found')

    # Chunked extraction
    vals = {g: np.zeros(n, dtype=np.float32) for g in gene_indices}
    chunk_size = 50000

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        for g, idx in gene_indices.items():
            chunk = adata.X[start:end, idx]
            vals[g][start:end] = (
                chunk.toarray().flatten() if hasattr(chunk, 'toarray')
                else np.array(chunk).flatten()
            )
        print(f'  {end}/{n}')

    # Build DataFrames
    for gene_name, values in vals.items():
        df = pd.DataFrame({
            'age': adata.obs['age'].values.astype(int),
            'total_counts': adata.obs['total_counts'].values,
            'raw': values,
            'afca_annotation': adata.obs['afca_annotation'].values,
        })
        df['cpm'] = (df['raw'] / df['total_counts']) * 1e6
        out_path = os.path.join(OUT_DIR, f'lu_{gene_name.lower()}.csv')
        df.to_csv(out_path, index=False)
        n_expr = (df['raw'] > 0).sum()
        print(f'  {gene_name}: {len(df)} cells, {n_expr} expressing ({n_expr/len(df)*100:.1f}%)')
        print(f'    Saved: {out_path}')

    print()


# ========================================================
if __name__ == '__main__':
    extract_davie()
    extract_lu()
    print('Done. All CSV files saved to', OUT_DIR)
