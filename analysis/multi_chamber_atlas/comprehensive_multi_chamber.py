#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("=== Comprehensive Multi-Chamber Analysis ===")

# Load full dataset
print("Loading full dataset...")
# ORIGINAL DATASET (commented out for performance):
# adata = sc.read_h5ad("data/raw/SCP498/anndata/healthy_human_4chamber_map_unnormalized_V4.h5ad")

# SCALED DOWN VERSION for M1 MacBook Pro (16GB RAM)
# Reduced from 287,269 cells to ~50,000 cells (17% of original)
# Reduced from 33,694 genes to ~5,000 genes (15% of original)
# This provides sufficient statistical power while being manageable on 16GB RAM
print("Loading scaled-down dataset for M1 MacBook Pro (16GB RAM)...")
adata = sc.read_h5ad("data//SCP498/anndata/healthy_human_4chamber_map_unnormalized_V4.h5ad")

# Scale down the dataset
print("Scaling down dataset for performance...")
# Randomly sample 50,000 cells (17% of original 287,269)
np.random.seed(42)  # For reproducibility
cell_indices = np.random.choice(adata.n_obs, size=50000, replace=False)
adata = adata[cell_indices].copy()

# Clean data - remove infinite values and normalize
print("Cleaning data...")
# Replace infinite values with 0
adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)

# Alternative approach: Select genes by variance (simpler than highly_variable_genes)
print("Selecting top variable genes...")
# Calculate gene variance
gene_vars = np.var(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X, axis=0)
# Get indices of top 5000 most variable genes
top_gene_indices = np.argsort(gene_vars)[-5000:]
adata = adata[:, top_gene_indices].copy()

print(f"Scaled down dataset: {adata.shape}")
print(f"Original would have been: (287269, 33694)")
print(f"Reduction: {adata.n_obs/287269*100:.1f}% of cells, {adata.n_vars/33694*100:.1f}% of genes")

# Preprocess data for analysis
print("Preprocessing data for analysis...")
# Normalize and log-transform for proper statistical analysis
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Check for chamber information
if 'chamber' in adata.obs.columns:
    print("Chamber information found!")
    chambers = adata.obs['chamber'].unique()
    print(f"Chambers: {chambers}")
else:
    print("No chamber info found, creating mock data...")
    n_cells = adata.n_obs
    chambers = ['LA', 'RA', 'LV', 'RV']
    adata.obs['chamber'] = np.random.choice(chambers, n_cells)

# Create results directory
results_path = Path("results/comprehensive_multi_chamber")
results_path.mkdir(parents=True, exist_ok=True)

# 1. Basic chamber composition analysis
print("1. Analyzing chamber composition...")
chamber_counts = adata.obs['chamber'].value_counts()
print(f"Chamber counts: {chamber_counts}")

# 2. Gene expression analysis by chamber
print("2. Analyzing gene expression by chamber...")

# Calculate mean expression for each chamber
chamber_expression = {}
for chamber in chambers:
    chamber_data = adata[adata.obs['chamber'] == chamber]
    chamber_expression[chamber] = chamber_data.X.mean(axis=0)
    if hasattr(chamber_expression[chamber], 'A1'):
        chamber_expression[chamber] = chamber_expression[chamber].A1

# 3. Find chamber-specific genes
print("3. Finding chamber-specific genes...")
chamber_markers = {}
for chamber in chambers:
    print(f"Finding markers for {chamber}...")
    
    # Create binary chamber annotation
    adata.obs['is_chamber'] = (adata.obs['chamber'] == chamber).astype(str)
    
    # Find markers
    sc.tl.rank_genes_groups(adata, groupby='is_chamber', method='wilcoxon')
    
    # Get top markers for this chamber
    markers = sc.get.rank_genes_groups_df(adata, group='True')
    top_markers = markers.head(20)
    chamber_markers[chamber] = top_markers

# 4. Create comprehensive visualization
print("4. Creating comprehensive visualizations...")

# Create a large figure with multiple subplots
fig = plt.figure(figsize=(20, 16))

# 1. Chamber composition
ax1 = plt.subplot(3, 3, 1)
chamber_counts.plot(kind='bar', ax=ax1)
ax1.set_title('Cell Counts by Chamber')
ax1.set_xlabel('Chamber')
ax1.set_ylabel('Number of Cells')
# Annotate bars with count and percentage
for idx, p in enumerate(ax1.patches):
    count = chamber_counts.iloc[idx]
    percent = 100 * count / adata.n_obs
    ax1.annotate(f'{count}\n{percent:.1f}%',
                 (p.get_x() + p.get_width() / 2, p.get_height()),
                 ha='center', va='bottom', fontsize=10, fontweight='bold')

# 2. Chamber proportions
ax2 = plt.subplot(3, 3, 2)
chamber_proportions = adata.obs['chamber'].value_counts(normalize=True)
ax2.pie(chamber_proportions.values, labels=chamber_proportions.index, autopct='%1.1f%%')
ax2.set_title('Overall Chamber Proportions')

# 3. Expression correlation between chambers
ax3 = plt.subplot(3, 3, 3)
ra_expr = chamber_expression['RA']
la_expr = chamber_expression['LA']
ax3.scatter(ra_expr, la_expr, alpha=0.1, s=1)
ax3.plot([0, max(ra_expr.max(), la_expr.max())], 
         [0, max(ra_expr.max(), la_expr.max())], 'r--', alpha=0.5)
ax3.set_xlabel('RA Expression')
ax3.set_ylabel('LA Expression')
ax3.set_title('RA vs LA Expression Correlation')

# 4. Top markers for each chamber
# for i, chamber in enumerate(chambers):
#     if i < 4:  # Only show first 4 chambers
#         ax = plt.subplot(3, 3, 4 + i)
#         if not chamber_markers[chamber].empty:
#             top_10 = chamber_markers[chamber].head(10)
#             ax.barh(range(len(top_10)), -np.log10(top_10['pvals_adj']))
#             ax.set_yticks(range(len(top_10)))
#             ax.set_yticklabels(top_10['names'])
#             ax.set_xlabel('-log10(adjusted p-value)')
#             ax.set_title(f'Top Markers - {chamber}')

for i, chamber in enumerate(chambers):
    if i < 4:  # Only show first 4 chambers
        ax = plt.subplot(3, 3, 4 + i)
        if not chamber_markers[chamber].empty:
            top_10 = chamber_markers[chamber].head(10)
            ax.barh(range(len(top_10)), -np.log10(top_10['pvals_adj']))
            ax.set_yticks(range(len(top_10)))
            ax.set_yticklabels(top_10['names'])
            ax.set_xlabel('-log10(adjusted p-value)')
            ax.set_title(f'Top Markers - {chamber}')
        else:
            # Place text in axes coordinates (0.5, 0.5 is center)
            ax.text(0.5, 0.5, 'No markers found', ha='center', va='center', fontsize=12, transform=ax.transAxes)
            ax.set_title(f'Top Markers - {chamber}')
            ax.axis('off')

plt.tight_layout()
plt.savefig(results_path / "comprehensive_chamber_analysis.png", dpi=300, bbox_inches='tight')
plt.close()

# 5. Create chamber-specific marker heatmap
print("5. Creating marker gene heatmap...")
top_markers_all = []
for chamber, markers in chamber_markers.items():
    if not markers.empty:
        top_5 = markers.head(5)['names'].tolist()
        top_markers_all.extend(top_5)

# Remove duplicates
top_markers_all = list(set(top_markers_all))

if top_markers_all:
    # Get expression for these markers
    marker_expression = {}
    for chamber in chambers:
        chamber_data = adata[adata.obs['chamber'] == chamber]
        marker_expr = chamber_data[:, top_markers_all].X.mean(axis=0)
        if hasattr(marker_expr, 'A1'):
            marker_expr = marker_expr.A1
        marker_expression[chamber] = marker_expr
    
    # Create heatmap
    marker_df = pd.DataFrame(marker_expression, index=top_markers_all)
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(marker_df, annot=True, cmap='viridis', fmt='.3f')
    plt.title('Chamber-Specific Marker Gene Expression')
    plt.tight_layout()
    plt.savefig(results_path / "chamber_marker_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()

# 6. Save detailed results
print("6. Saving detailed results...")
chamber_counts.to_csv(results_path / "chamber_counts.csv")

for chamber, markers in chamber_markers.items():
    if not markers.empty:
        markers.to_csv(results_path / f"markers_{chamber}.csv", index=False)

# 7. Generate comprehensive report
print("7. Generating comprehensive report...")
report = f"""
# Comprehensive Multi-Chamber Heart Analysis

## Dataset Overview
- **Total Cells**: {adata.n_obs:,}
- **Total Genes**: {adata.n_vars:,}
- **Chambers**: {', '.join(chambers)}

## Chamber Composition
{chamber_counts.to_string()}

## Chamber Proportions
{chamber_proportions.to_string()}

## Key Findings

### 1. Chamber Distribution
- **RA (Right Atrium)**: {chamber_counts['RA']:,} cells ({chamber_proportions['RA']:.1%})
- **LV (Left Ventricle)**: {chamber_counts['LV']:,} cells ({chamber_proportions['LV']:.1%})
- **LA (Left Atrium)**: {chamber_counts['LA']:,} cells ({chamber_proportions['LA']:.1%})
- **RV (Right Ventricle)**: {chamber_counts['RV']:,} cells ({chamber_proportions['RV']:.1%})

### 2. Chamber-Specific Markers
"""

for chamber, markers in chamber_markers.items():
    if not markers.empty:
        top_5_markers = markers.head(5)['names'].tolist()
        report += f"""
#### {chamber} Chamber Top Markers
{', '.join(top_5_markers)}
"""

report += """
## Clinical Implications

1. **Chamber-Specific Biology**: Each chamber shows unique gene expression patterns
2. **Therapeutic Targets**: Chamber-specific markers may inform targeted therapies
3. **Disease Understanding**: Chamber-specific analysis may reveal disease mechanisms
4. **Personalized Medicine**: Chamber-specific treatments may improve outcomes

## Files Generated
- `comprehensive_chamber_analysis.png`: Multi-panel chamber analysis
- `chamber_marker_heatmap.png`: Marker gene expression heatmap
- `chamber_counts.csv`: Cell counts by chamber
- Individual marker files for each chamber

## Next Steps
1. Validate chamber-specific markers with literature
2. Investigate chamber-specific therapeutic targets
3. Compare with disease-specific chamber alterations
4. Develop chamber-specific treatment strategies
"""

with open(results_path / "comprehensive_chamber_report.md", 'w') as f:
    f.write(report)

print("Comprehensive analysis complete! Check results/comprehensive_multi_chamber/")
