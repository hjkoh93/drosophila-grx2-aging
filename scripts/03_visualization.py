"""
03_visualization.py
===================
Generate Grx1 (Grx2) figures:
  - 6-panel expression analysis (mean CPM, % expressing, violin)
  - Age trajectory plot (mean CPM +/- SEM)
  - Cell-type x age heatmaps (Davie + Lu)
  - Dopaminergic neuron trajectory and detection plots

Input:  output/davie_grx1.csv, output/lu_grx1.csv
Output: output/*.png
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 15, 'axes.titlesize': 18, 'axes.labelsize': 16,
    'xtick.labelsize': 14, 'ytick.labelsize': 14,
    'figure.dpi': 150, 'font.family': 'sans-serif',
})
OUT_DIR = 'output'

# Non-neuronal cell types to exclude from heatmaps
EXCLUDE_DAVIE = ['glia', 'astrocyte']
EXCLUDE_LU = [
    'glia', 'astrocyte', 'fat body', 'fat mass', 'epithelial',
    'muscle', 'pigment', 'hemocyte', 'trachea', 'unannotated',
]


# ========================================================
# Helper: 6-panel expression analysis
# ========================================================
def plot_6panel(davie, lu, gene, filename):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    for row, (df, acol, label, cmap_name) in enumerate([
        (davie, 'Age', 'Davie et al.', 'Reds'),
        (lu, 'age', 'Lu et al.', 'Blues'),
    ]):
        summary = df.groupby(acol).agg(
            mean_cpm=('cpm', 'mean'), sem_cpm=('cpm', 'sem'),
            pct_expr=('raw', lambda x: (x > 0).mean() * 100),
        ).reset_index()
        ages = summary[acol].values
        colors = plt.cm.get_cmap(cmap_name)(np.linspace(0.3, 0.9, len(summary)))

        # Mean CPM bars
        ax = axes[row, 0]
        ax.bar(range(len(summary)), summary['mean_cpm'], yerr=summary['sem_cpm'],
               color=colors, capsize=3, edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(summary)))
        ax.set_xticklabels([f'{int(a)}d' for a in ages])
        ax.set_xlabel('Age'); ax.set_ylabel('Mean CPM +/- SEM')
        ax.set_title(f'{label} — Mean {gene} CPM', fontweight='bold')
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

        # % Expressing bars
        ax = axes[row, 1]
        ax.bar(range(len(summary)), summary['pct_expr'],
               color=colors, edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(summary)))
        ax.set_xticklabels([f'{int(a)}d' for a in ages])
        ax.set_xlabel('Age'); ax.set_ylabel('% Cells Expressing')
        ax.set_title(f'{label} — % Expressing {gene}', fontweight='bold')
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

        # Violin
        ax = axes[row, 2]
        vdata = [df[df[acol] == a]['cpm'].values for a in ages]
        vp = ax.violinplot(vdata, positions=range(len(ages)),
                           showmedians=True, showextrema=False)
        c = '#d62728' if row == 0 else '#1f77b4'
        for body in vp['bodies']:
            body.set_facecolor(c); body.set_alpha(0.5)
        ax.set_xticks(range(len(ages)))
        ax.set_xticklabels([f'{int(a)}d' for a in ages])
        ax.set_xlabel('Age'); ax.set_ylabel('CPM')
        ax.set_title(f'{label} — {gene} Distribution', fontweight='bold')
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    fig.suptitle(f'{gene} Age-Dependent Expression Analysis',
                 fontsize=15, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(OUT_DIR, filename), bbox_inches='tight', dpi=150)
    plt.close()
    print(f'  Saved: {filename}')


# ========================================================
# Helper: trajectory plot (mean +/- SEM)
# ========================================================
def plot_trajectory(davie, lu, gene, filename, c1='#d62728', c2='#1f77b4'):
    ds = davie.groupby(davie.columns[0]).agg(
        mean_cpm=('cpm', 'mean'), sem_cpm=('cpm', 'sem')).reset_index()
    ls = lu.groupby('age').agg(
        mean_cpm=('cpm', 'mean'), sem_cpm=('cpm', 'sem')).reset_index()

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5))
    acol_d = ds.columns[0]

    for ax, s, acol, color, mk, title in [
        (a1, ds, acol_d, c1, 'o',
         f'Davie et al. (GSE107451)\n0-50 days, n={len(davie):,} cells'),
        (a2, ls, 'age', c2, 's',
         f'Lu et al.\n5-70 days, n={len(lu):,} cells'),
    ]:
        ax.errorbar(s[acol], s['mean_cpm'], yerr=s['sem_cpm'],
                    marker=mk, capsize=5, linewidth=2.5, color=color, markersize=9)
        ax.fill_between(s[acol], s['mean_cpm'] - s['sem_cpm'],
                        s['mean_cpm'] + s['sem_cpm'], alpha=0.2, color=color)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Age (days)', fontsize=16)
        ax.set_ylabel(f'Mean {gene} CPM +/- SEM', fontsize=16)
        ax.set_title(title, fontweight='bold', fontsize=17)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.set_xticks(s[acol]); ax.tick_params(labelsize=14)

    fig.suptitle(f'{gene} Expression Trajectory with Aging',
                 fontsize=20, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), bbox_inches='tight', dpi=150)
    plt.close()
    print(f'  Saved: {filename}')


# ========================================================
# Helper: cell-type heatmap
# ========================================================
def plot_heatmap(df, age_col, annot_col, exclude_list, dataset, gene, filename,
                 dopa_label='Dopaminergic'):
    ct = df.groupby(annot_col).agg(
        n=('raw', 'size'), mean_cpm=('cpm', 'mean')).reset_index()
    pattern = '|'.join(exclude_list)
    ct_filt = ct[~ct[annot_col].str.lower().str.contains(pattern, na=False)]
    ct_filt = ct_filt[ct_filt['n'] >= 50].sort_values('mean_cpm', ascending=False)
    top7 = ct_filt.head(7)[annot_col].tolist()
    if dopa_label not in top7:
        top7.append(dopa_label)

    sub = df[df[annot_col].isin(top7)]
    pv = sub.groupby([annot_col, age_col])['cpm'].mean().reset_index()
    pw = pv.pivot(index=annot_col, columns=age_col, values='cpm').reindex(top7)

    fig, ax = plt.subplots(figsize=(16, 8))
    sns.heatmap(pw, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                linewidths=0.5, cbar_kws={'label': 'Mean CPM'},
                annot_kws={'size': 15})
    ax.set_xlabel('Age (days)', fontsize=18)
    ax.set_ylabel('Cell Type', fontsize=18)
    ax.set_xticklabels(
        [f'{int(float(t.get_text()))}d' for t in ax.get_xticklabels()], fontsize=16)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=15)
    ax.set_title(f'{gene} Mean CPM by Cell Type and Age ({dataset})\n'
                 f'Top Neuronal Cell Types + Dopaminergic',
                 fontweight='bold', fontsize=19)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_ylabel('Mean CPM', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), bbox_inches='tight', dpi=150)
    plt.close()
    print(f'  Saved: {filename}')


# ========================================================
# Helper: dopaminergic neuron analysis
# ========================================================
def plot_dopaminergic(davie, lu, gene, filename_traj, filename_det,
                      dopa_davie='Dopaminergic',
                      dopa_lu='dopaminergic PAM neuron',
                      color='#d62728'):
    dd = davie[davie['annotation'] == dopa_davie].copy()
    dl = lu[lu['afca_annotation'] == dopa_lu].copy()

    dd_s = dd.groupby('Age').agg(
        n=('raw', 'size'), mean_cpm=('cpm', 'mean'), sem=('cpm', 'sem'),
        pct_expr=('raw', lambda x: (x > 0).mean() * 100)).reset_index()
    dl_s = dl.groupby('age').agg(
        n=('raw', 'size'), mean_cpm=('cpm', 'mean'), sem=('cpm', 'sem'),
        pct_expr=('raw', lambda x: (x > 0).mean() * 100)).reset_index()

    # Trajectory
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, s, acol, title in [
        (a1, dd_s, 'Age',
         f'Davie et al. — Dopaminergic\nn={int(dd_s["n"].sum())} cells'),
        (a2, dl_s, 'age',
         f'Lu et al. — Dopaminergic PAM\nn={int(dl_s["n"].sum())} cells'),
    ]:
        ax.errorbar(s[acol], s['mean_cpm'], yerr=s['sem'],
                    marker='o', capsize=5, linewidth=2.5, color=color, markersize=9)
        ax.fill_between(s[acol], s['mean_cpm'] - s['sem'],
                        s['mean_cpm'] + s['sem'], alpha=0.2, color=color)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Age (days)', fontsize=16)
        ax.set_ylabel(f'Mean {gene} CPM +/- SEM', fontsize=16)
        ax.set_title(title, fontweight='bold', fontsize=17)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.set_xticks(s[acol]); ax.tick_params(labelsize=14)
        for _, r in s.iterrows():
            ax.annotate(f'{r["pct_expr"]:.0f}%',
                        xy=(r[acol], r['mean_cpm'] + r['sem']),
                        ha='center', va='bottom', fontsize=9, color='gray')

    fig.suptitle(f'{gene} in Dopaminergic Neurons with Aging',
                 fontsize=20, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename_traj), bbox_inches='tight', dpi=150)
    plt.close()
    print(f'  Saved: {filename_traj}')

    # Detection bar chart
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, s, acol, title in [
        (a1, dd_s, 'Age', 'Davie et al. — Dopaminergic'),
        (a2, dl_s, 'age', 'Lu et al. — Dopaminergic PAM'),
    ]:
        ax.bar(range(len(s)), s['pct_expr'], color=color, alpha=0.7,
               edgecolor='black', linewidth=0.5)
        ax.set_xticks(range(len(s)))
        ax.set_xticklabels([f'{int(a)}d' for a in s[acol]])
        ax.set_xlabel('Age', fontsize=16)
        ax.set_ylabel(f'% Cells Expressing {gene}', fontsize=16)
        ax.set_title(title, fontweight='bold', fontsize=17)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        for i, (_, r) in enumerate(s.iterrows()):
            ax.annotate(f'n={int(r["n"])}', xy=(i, r['pct_expr']),
                        ha='center', va='bottom', fontsize=9)

    fig.suptitle(f'{gene} Detection Rate in Dopaminergic Neurons',
                 fontsize=20, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename_det), bbox_inches='tight', dpi=150)
    plt.close()
    print(f'  Saved: {filename_det}')


# ========================================================
# Main
# ========================================================
if __name__ == '__main__':
    gene_label = 'Grx2'
    c1, c2 = '#d62728', '#1f77b4'

    print(f'\n{"="*65}')
    print(f'Generating figures for {gene_label}')
    print(f'{"="*65}')

    davie = pd.read_csv(f'{OUT_DIR}/davie_grx1.csv')
    lu = pd.read_csv(f'{OUT_DIR}/lu_grx1.csv')

    # 6-panel
    plot_6panel(davie, lu, gene_label,
                f'{gene_label}_age_expression_analysis.png')

    # Trajectory
    plot_trajectory(davie, lu, gene_label,
                    f'{gene_label}_age_trajectory.png', c1, c2)

    # Heatmaps
    plot_heatmap(davie, 'Age', 'annotation', EXCLUDE_DAVIE,
                 'Davie et al.', gene_label,
                 f'{gene_label}_celltype_age_heatmap_Davie.png',
                 dopa_label='Dopaminergic')

    plot_heatmap(lu, 'age', 'afca_annotation', EXCLUDE_LU,
                 'Lu et al.', gene_label,
                 f'{gene_label}_celltype_age_heatmap_Lu.png',
                 dopa_label='dopaminergic PAM neuron')

    # Dopaminergic
    plot_dopaminergic(davie, lu, gene_label,
                      f'{gene_label}_dopaminergic_trajectory.png',
                      f'{gene_label}_dopaminergic_detection.png',
                      color=c1)

    print('\nAll figures saved to', OUT_DIR)
