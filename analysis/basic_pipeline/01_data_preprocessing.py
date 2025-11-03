import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import subprocess
# Set up scanpy
sc.settings.verbosity = 3
sc.set_figure_params(dpi=80, facecolor='white')

def load_and_preprocess_data(data_path):
    """Load and preprocess the heart dataset"""
    print("Loading data...")
    adata = sc.read_h5ad(data_path)
    
    # Basic info
    print(f"Dataset shape: {adata.shape}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    
    # Make gene names unique
    adata.var_names_make_unique()
    
    # Store raw data
    adata.raw = adata
    
    return adata

def basic_filtering(adata, min_genes=200, min_cells=3):
    """Basic filtering of cells and genes"""
    print("Performing basic filtering...")
    
    # Filter cells with too few genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Filter genes present in too few cells
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.shape}")
    return adata

def calculate_qc_metrics(adata):
    """Calculate quality control metrics"""
    print("Calculating QC metrics...")
    
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    
    # Hemoglobin genes
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], 
                              percent_top=None, log1p=False, inplace=True)
    
    return adata

#ensure data integrity before processing
def verify_data_integrity(data_dir, checksum_file):
    """Verify data integrity using SHA-256 checksums."""
    print("Verifying data integrity...")
    script_path = Path(__file__).parent.parent.parent / "utils/sha256_checksum.py"
    result = subprocess.run(
        ["python", script_path, "verify", data_dir, checksum_file],
        capture_output=True,
        text=True
    )
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError("Data integrity verification failed. Please check the checksums.")

def main():
    # Paths - use project root as reference
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "data/raw"
    checksum_file = project_root / "data/raw/checksums.txt"
    output_path = project_root / "data/processed/"
    output_path.mkdir(parents=True, exist_ok=True)

    # Verify data integrity
    verify_data_integrity(data_dir, checksum_file)
    
    # Load data
    data_path = project_root / "data/validation/gse145154_merged.h5ad"
    adata = load_and_preprocess_data(data_path)
    
    # Basic filtering
    adata = basic_filtering(adata)
    
    # Calculate QC metrics
    adata = calculate_qc_metrics(adata)
    
    # Save preprocessed data
    adata.write(output_path / "heart_data_preprocessed.h5ad")
    
    print("Preprocessing complete!")
    print(f"Final dataset shape: {adata.shape}")

if __name__ == "__main__":
    main()