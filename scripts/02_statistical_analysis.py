"""
02_statistical_analysis.py
==========================
Statistical analysis of age-dependent Grx1 (Grx2) expression:
  - Summary statistics per age group
  - Spearman correlation (whole-population)
  - Kruskal-Wallis H test
  - Pairwise Mann-Whitney U with Bonferroni correction

Input:  output/davie_grx1.csv, output/lu_grx1.csv
Output: printed results to stdout
"""

import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

OUT_DIR = 'output'


# ========================================================
# Summary statistics
# ========================================================
def summary_table(df, age_col, dataset, gene):
    print(f'\n{"="*65}')
    print(f'{dataset} — {gene} Summary Statistics')
    print(f'{"="*65}')
    summary = df.groupby(age_col).agg(
        n=('raw', 'size'),
        mean_cpm=('cpm', 'mean'),
        sem_cpm=('cpm', 'sem'),
        sd_cpm=('cpm', 'std'),
        median_cpm=('cpm', 'median'),
        pct_expr=('raw', lambda x: (x > 0).mean() * 100),
    ).reset_index()
    print(f'{"Age":>5} {"n":>7} {"Mean CPM":>10} {"SEM":>8} {"SD":>8} '
          f'{"Median":>8} {"%Expr":>7}')
    for _, r in summary.iterrows():
        print(f'{int(r[age_col]):>5}d {int(r["n"]):>7} {r["mean_cpm"]:>10.2f} '
              f'{r["sem_cpm"]:>8.2f} {r["sd_cpm"]:>8.2f} '
              f'{r["median_cpm"]:>8.2f} {r["pct_expr"]:>6.1f}%')
    return summary


# ========================================================
# Spearman correlation
# ========================================================
def spearman_test(df, age_col, dataset, gene):
    rho, p = stats.spearmanr(df[age_col], df['cpm'])
    print(f'\n{dataset} — {gene}')
    print(f'  Spearman (all cells): rho = {rho:.4f}, p = {p:.2e}, n = {len(df)}')
    return rho, p


# ========================================================
# Kruskal-Wallis
# ========================================================
def kruskal_test(df, age_col, dataset, gene):
    groups = [g['cpm'].values for _, g in df.groupby(age_col)]
    if len(groups) > 1:
        try:
            kw, kwp = stats.kruskal(*groups)
            print(f'  Kruskal-Wallis: H = {kw:.2f}, p = {kwp:.2e}')
            return kw, kwp
        except ValueError:
            print('  Kruskal-Wallis: insufficient variance')
    return None, None


# ========================================================
# Pairwise Mann-Whitney U
# ========================================================
def pairwise_mwu(df, age_col, pairs, dataset, gene):
    print(f'\n{"="*65}')
    print(f'{dataset} — {gene} Pairwise Mann-Whitney U')
    print(f'{"="*65}')
    n_comp = len(pairs)
    for young, old in pairs:
        g1 = df[df[age_col] == young]['cpm']
        g2 = df[df[age_col] == old]['cpm']
        n1, n2 = len(g1), len(g2)
        u_stat, mw_p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
        r_rb = 1 - (2 * u_stat) / (n1 * n2)
        pooled_std = np.sqrt(
            ((n1 - 1) * g1.std()**2 + (n2 - 1) * g2.std()**2) / (n1 + n2 - 2)
        )
        cohens_d = (g1.mean() - g2.mean()) / pooled_std if pooled_std > 0 else 0
        bonf_p = min(mw_p * n_comp, 1.0)
        fc = g1.mean() / g2.mean() if g2.mean() > 0 else float('inf')

        print(f'\n  {young}d (n={n1:,}) vs {old}d (n={n2:,})')
        print(f'    Mean CPM: {young}d = {g1.mean():.2f}, {old}d = {g2.mean():.2f}')
        print(f'    Median:   {young}d = {g1.median():.2f}, {old}d = {g2.median():.2f}')
        print(f'    FC ({young}d/{old}d) = {fc:.2f}x')
        print(f'    %Expr: {young}d = {(g1>0).mean()*100:.1f}%, '
              f'{old}d = {(g2>0).mean()*100:.1f}%')
        print(f'    U = {u_stat:,.0f}, p = {mw_p:.2e}, Bonf p = {bonf_p:.2e}')
        print(f'    r_rb = {r_rb:.4f}, Cohen\'s d = {cohens_d:.4f}')
    print(f'\n  (Bonferroni correction: {n_comp} comparisons)')


# ========================================================
# Main
# ========================================================
def run_analysis(csv_path, age_col, dataset, gene, pairs):
    df = pd.read_csv(csv_path)
    summary_table(df, age_col, dataset, gene)
    spearman_test(df, age_col, dataset, gene)
    kruskal_test(df, age_col, dataset, gene)
    pairwise_mwu(df, age_col, pairs, dataset, gene)


if __name__ == '__main__':
    davie_pairs = [(0, 15), (0, 30), (0, 50)]
    lu_pairs = [(5, 30), (5, 50), (5, 70)]

    run_analysis(f'{OUT_DIR}/davie_grx1.csv', 'Age',
                 'Davie et al.', 'Grx2', davie_pairs)
    run_analysis(f'{OUT_DIR}/lu_grx1.csv', 'age',
                 'Lu et al.', 'Grx2', lu_pairs)

    print('\n' + '=' * 65)
    print('Effect Size Interpretation')
    print('  Rank-biserial |r|: <0.1 negligible, 0.1-0.3 small, '
          '0.3-0.5 medium, >0.5 large')
    print('  Cohen\'s |d|:      <0.2 negligible, 0.2-0.5 small, '
          '0.5-0.8 medium, >0.8 large')
    print('=' * 65)
