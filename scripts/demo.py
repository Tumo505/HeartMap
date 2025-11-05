#!/usr/bin/env python3
"""
HeartMAP Quick Demo
This script demonstrates the basic functionality of the refactored HeartMAP system.
"""

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

# Import HeartMAP
from heartmap import HeartMapModel, Config
from heartmap.data import DataLoader

def create_demo_data():
    """Create mock single-cell heart data for demonstration"""
    print(" Creating mock heart data...")
    
    # Parameters
    n_cells = 500
    n_genes = 1000
    
    # Generate mock count matrix
    np.random.seed(42)
    X = np.random.negative_binomial(10, 0.3, size=(n_cells, n_genes))
    X = csr_matrix(X.astype(np.float32))
    
    # Create mock cell metadata
    chambers = np.random.choice(['RA', 'RV', 'LA', 'LV'], n_cells)
    cell_types = np.random.choice(['Cardiomyocyte', 'Fibroblast', 'Endothelial', 'Immune'], n_cells)
    
    obs = pd.DataFrame({
        'chamber': chambers,
        'cell_type': cell_types,
        'sample_id': np.random.choice(['Sample_001', 'Sample_002', 'Sample_003'], n_cells)
    })
    
    # Create gene metadata
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    var = pd.DataFrame(index=gene_names)
    
    # Create AnnData object
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    print(f" Created mock data: {n_cells} cells, {n_genes} genes")
    print(f"   Chambers: {list(obs['chamber'].unique())}")
    print(f"   Cell types: {list(obs['cell_type'].unique())}")
    
    return adata

def demo_configuration():
    """Demonstrate configuration management"""
    print("\n⚙️ Configuring HeartMAP...")
    
    # Create default configuration
    config = Config.default()
    
    # Customize for demo
    config.data.max_cells_subset = 200  # Small for demo
    config.data.max_genes_subset = 500
    config.data.test_mode = True  # Enable test mode
    config.model.model_type = "basic"
    
    print(f" Configuration created:")
    print(f"   Model type: {config.model.model_type}")
    print(f"   Max cells: {config.data.max_cells_subset}")
    print(f"   Max genes: {config.data.max_genes_subset}")
    print(f"   Test mode: {config.data.test_mode}")
    
    return config

def demo_data_loading(adata, config):
    """Demonstrate data loading and preprocessing"""
    print("\n Processing data...")
    
    # Initialize data loader
    data_loader = DataLoader(config)
    
    # Preprocess data
    processed_adata = data_loader.preprocess(adata)
    
    print(f" Data processing completed:")
    print(f"   Original: {adata.shape}")
    print(f"   Processed: {processed_adata.shape}")
    
    return processed_adata

def demo_analysis(processed_adata, config):
    """Demonstrate basic analysis"""
    print("\n Running HeartMAP analysis...")
    
    # Initialize HeartMAP model
    model = HeartMapModel(config)
    
    # For demo, we'll just show the model creation
    # Full analysis would be: results = model.analyze(processed_adata)
    print(" HeartMAP model initialized successfully")
    print(f"   Model type: {config.model.model_type}")
    print(f"   Configuration: {config.data.max_cells_subset} max cells")
    
    # Simulate some basic results
    results = {
        'processed_data': processed_adata,
        'n_cells_analyzed': processed_adata.n_obs,
        'n_genes_analyzed': processed_adata.n_vars,
        'chambers_detected': processed_adata.obs['chamber'].unique().tolist(),
        'cell_types_detected': processed_adata.obs['cell_type'].unique().tolist()
    }
    
    print(f" Analysis summary:")
    print(f"   Cells analyzed: {results['n_cells_analyzed']}")
    print(f"   Genes analyzed: {results['n_genes_analyzed']}")
    print(f"   Chambers: {results['chambers_detected']}")
    print(f"   Cell types: {results['cell_types_detected']}")
    
    return results

def demo_visualization(results):
    """Create simple visualization"""
    print("\n Creating visualizations...")
    
    adata = results['processed_data']
    
    # Chamber composition
    chamber_counts = adata.obs['chamber'].value_counts()
    cell_type_counts = adata.obs['cell_type'].value_counts()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Chamber distribution
    chamber_counts.plot(kind='bar', ax=ax1, color='lightblue')
    ax1.set_title('Cells per Chamber')
    ax1.set_xlabel('Chamber')
    ax1.set_ylabel('Number of Cells')
    ax1.tick_params(axis='x', rotation=45)
    
    # Cell type distribution  
    cell_type_counts.plot(kind='bar', ax=ax2, color='lightcoral')
    ax2.set_title('Cells per Type')
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('Number of Cells')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('figures/demo_results.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print(" Visualization saved to: figures/demo_results.png")

def main():
    """Main demo function"""
    print(" HeartMAP Quick Demo")
    print("=" * 50)
    
    try:
        # Create directories
        import os
        os.makedirs('figures', exist_ok=True)
        
        # Step 1: Create demo data
        adata = create_demo_data()
        
        # Step 2: Configure analysis
        config = demo_configuration()
        
        # Step 3: Process data
        processed_adata = demo_data_loading(adata, config)
        
        # Step 4: Run analysis
        results = demo_analysis(processed_adata, config)
        
        # Step 5: Visualize results
        demo_visualization(results)
        
        print("\n Demo completed successfully!")
        print("\nNext steps:")
        print("1. Try with real data: place .h5ad files in data/raw/")
        print("2. Use command line: heartmap data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad")
        print("3. Start web interface: python app.py")
        print("4. Explore notebooks: notebooks/quickstart_example.ipynb")
        print("5. Read documentation: DEPLOYMENT_GUIDE.md")
        
    except Exception as e:
        print(f" Demo failed: {e}")
        print("Please check the installation and try again.")
        
if __name__ == "__main__":
    main()
