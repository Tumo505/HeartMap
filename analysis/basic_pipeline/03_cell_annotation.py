import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def find_variable_genes(adata, n_top_genes=2000):
    """Find highly variable genes"""
    print("Finding highly variable genes...")
    
    # This computes statistics (mean, variance, dispersion) for each gene and identifies the set of highly variable genes, which are most informative for downstream analysis.
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_top_genes)
    
    # Plot highly variable genes
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(Path("results/figures/highly_variable_genes_validation_gse145154.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return adata

def perform_pca(adata, n_comps=50):
    """Perform Principal Component Analysis(PCA)"""
    print("Performing PCA...")
    
    # Keep only highly variable genes for PCA
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
    
    # Plot PCA
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
    plt.savefig(Path("results/figures/pca_variance_validation_gse145154.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return adata

def compute_neighborhood_graph(adata, n_neighbors=10, n_pcs=40):
    """Compute neighborhood graph"""
    print("Computing neighborhood graph...")
    
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    return adata

def perform_clustering(adata, resolution=0.5):
    """Perform clustering (K-means as fallback if igraph not available)"""
    print("Performing clustering...")
    
    try:
        # Try Leiden clustering first
        sc.tl.leiden(adata, resolution=resolution)
        print("Used Leiden clustering")
    except ImportError:
        # Fallback to K-means clustering if igraph is not available
        print("igraph not available, using K-means clustering instead")
        # Use PCA components for K-means
        n_clusters = 10  # You can adjust this number
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(adata.obsm['X_pca'][:, :20])
        # Convert to categorical dtype for Scanpy compatibility
        adata.obs['leiden'] = pd.Categorical(cluster_labels.astype(str))
        print(f"Used K-means clustering with {n_clusters} clusters")
    
    return adata

def compute_umap(adata):
    """Compute UMAP embedding"""
    print("Computing UMAP...")
    
    sc.tl.umap(adata)
    
    return adata

def plot_clustering_results(adata, save_path):
    """Plot clustering results"""
    print("Plotting clustering results...")
    
    # Plot UMAP with clusters
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data',
               title='Leiden Clustering', frameon=False, save='_leiden_clustering_validation_gse145154.png')
    
    # If cell type annotations exist, plot them too
    if 'cell_type' in adata.obs.columns:
        sc.pl.umap(adata, color=['cell_type'], legend_loc='on data',
                   title='Cell Types', frameon=False, save='_cell_types_validation_gse145154.png')

def identify_marker_genes(adata, groupby='leiden', n_genes=25):
    """Identify marker genes for each cluster"""
    print("Identifying marker genes...")
    
    sc.tl.rank_genes_groups(adata, groupby, method='wilcoxon')
    
    # Plot marker genes
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False, show=False)
    plt.savefig(Path("results/figures/marker_genes_validation_gse145154.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save marker genes to file
    marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    marker_genes.to_csv(Path("results/tables/marker_genes_validation_gse145154.csv"))
    
    return adata

def main():
    # Load filtered and normalized data
    project_root = Path(__file__).parent.parent.parent
    data_path = project_root / "data/processed/heart_data_filtered_normalized.h5ad"
    adata = sc.read_h5ad(data_path)
    
    # Create results directories
    Path("results/figures").mkdir(parents=True, exist_ok=True)
    Path("results/tables").mkdir(parents=True, exist_ok=True)
    
    # Find highly variable genes using the full dataset
    # adata = find_variable_genes(adata) // This is the full dataset, use it if you have enough compute resources, otherwise use the below function

    # Subset to 10,000 cells, 1,000 genes (the full dataset is 287,269 cells and 33,694 genes(full dataset), and 262,003 cells and 28,550 genes after filtering)
    adata_subset = adata[:10000, :1000]  
    adata_subset = find_variable_genes(adata_subset)

    # Use the subset for downstream steps (temporary)
    adata = adata_subset
    
    # Perform PCA
    adata = perform_pca(adata)
    
    # Compute neighborhood graph
    adata = compute_neighborhood_graph(adata)
    
    # Perform clustering
    adata = perform_clustering(adata)
    
    # Compute UMAP
    adata = compute_umap(adata)
    
    # Plot results
    plot_clustering_results(adata, Path("results/figures"))
    
    # Identify marker genes
    adata = identify_marker_genes(adata)
    
    # Save annotated data
    adata.write(project_root / "data/processed/heart_data_annotated_validation_gse145154.h5ad")
    
    print("Cell annotation complete!")

if __name__ == "__main__":
    main()