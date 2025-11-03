import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def plot_qc_metrics(adata, save_path):
    """Plot quality control metrics"""
    print("Generating QC plots...")
    
    # Set up the plotting
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Violin plot: total counts per cell, number of genes per cell, and mitochondrial gene percentage
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, show=False)
    plt.savefig(save_path / "qc_violin_validation_gse145154.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Scatter plot: total counts vs mitochondrial gene percentage (detect contamination)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False)
    plt.savefig(save_path / "qc_scatter1_validation_gse145154.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Scatter plot: total counts vs number of genes per cell (detect doublets/high complexity)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False)
    plt.savefig(save_path / "qc_scatter2_validation_gse145154.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Scatter plot: number of genes per cell vs mitochondrial gene percentage (detect contamination)
    sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', show=False)
    plt.savefig(save_path / "qc_scatter3_validation_gse145154.png", dpi=300, bbox_inches='tight')
    plt.close()

# Filter cells and genes
def filter_cells_genes(adata, max_mt=20, min_genes=200, max_genes=5000):
    """Filter cells and genes based on QC metrics"""
    print("Filtering cells and genes...")
    
    # Store original counts
    n_cells_orig = adata.n_obs
    n_genes_orig = adata.n_vars
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    adata = adata[adata.obs.pct_counts_mt < max_mt, :]
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=3)
    
    print(f"Filtered {n_cells_orig - adata.n_obs} cells")
    print(f"Filtered {n_genes_orig - adata.n_vars} genes")
    print(f"Remaining: {adata.n_obs} cells, {adata.n_vars} genes")
    
    return adata

def normalize_data(adata):
    """Normalize the data"""
    print("Normalizing data...")
    
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log transform
    sc.pp.log1p(adata)
    
    return adata

def main():
    # Load preprocessed data
    project_root = Path(__file__).parent.parent.parent
    data_path = project_root / "data/processed/heart_data_preprocessed.h5ad"
    adata = sc.read_h5ad(data_path)
    
    # Create results directory
    results_path = Path("results/figures")
    results_path.mkdir(parents=True, exist_ok=True)
    
    # Plot QC metrics
    plot_qc_metrics(adata, results_path)
    
    # Filter cells and genes
    adata = filter_cells_genes(adata)
    
    # Normalize data
    adata = normalize_data(adata)
    
    # Save filtered and normalized data
    adata.write(project_root / "data/processed/heart_data_filtered_normalized_validation_gse145154.h5ad")
    
    print("Quality control complete!")

if __name__ == "__main__":
    main()