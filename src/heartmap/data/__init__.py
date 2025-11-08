"""
Data processing utilities for HeartMAP
"""

import os
import hashlib
from pathlib import Path
from typing import Tuple, List, Union
import warnings

import scanpy as sc
import numpy as np
import anndata as ad
from scipy.sparse import issparse

from ..config import Config


class DataValidator:
    """Validate data integrity and format"""

    @staticmethod
    def verify_checksum(file_path: str, expected_checksum: str) -> bool:
        """Verify file checksum"""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest() == expected_checksum

    @staticmethod
    def validate_anndata(adata: ad.AnnData, check_qc_metrics: bool = True) -> Tuple[bool, List[str]]:
        """Validate AnnData object structure"""
        issues = []

        if adata.n_obs == 0:
            issues.append("No cells in dataset")
        if adata.n_vars == 0:
            issues.append("No genes in dataset")

        # Check for QC metrics only if requested (after they should be calculated)
        if check_qc_metrics:
            # scanpy creates these standard QC metric columns
            required_obs = ['n_genes_by_counts', 'total_counts']
            for col in required_obs:
                if col not in adata.obs.columns:
                    issues.append(f"Missing required obs column: {col}")

        # Check for NaN/inf values
        if issparse(adata.X):
            if not np.isfinite(adata.X.data).all():
                issues.append("Non-finite values in X matrix")
        else:
            if not np.isfinite(adata.X).all():
                issues.append("Non-finite values in X matrix")

        return len(issues) == 0, issues


class DataLoader:
    """Load and preprocess data"""

    def __init__(self, config: Config):
        self.config = config

    def load_raw_data(
        self, file_path: Union[str, Path], verify_integrity: bool = True
    ) -> ad.AnnData:
        """Load raw single-cell data"""
        file_path = Path(file_path)

        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")

        # Load data based on file format
        if file_path.suffix == '.h5ad':
            adata = sc.read_h5ad(file_path)
        elif file_path.suffix == '.h5':
            adata = sc.read_10x_h5(file_path, genome=None, gex_only=True)
        elif file_path.suffix == '.csv':
            adata = sc.read_csv(file_path).T  # Transpose to have genes as variables
        else:
            raise ValueError(f"Unsupported file format: {file_path.suffix}")

        # Validate data (skip QC metrics check for raw data)
        is_valid, issues = DataValidator.validate_anndata(adata, check_qc_metrics=False)
        if not is_valid:
            warnings.warn(f"Data validation issues: {'; '.join(issues)}")

        return adata

    def preprocess_basic(self, adata: ad.AnnData) -> ad.AnnData:
        """Basic preprocessing pipeline"""
        adata = adata.copy()

        # Make gene names unique
        adata.var_names_make_unique()

        # Store raw data
        adata.raw = adata

        # Basic filtering
        sc.pp.filter_cells(adata, min_genes=self.config.data.min_genes)
        sc.pp.filter_genes(adata, min_cells=self.config.data.min_cells)

        return adata

    def calculate_qc_metrics(self, adata: ad.AnnData) -> ad.AnnData:
        """Calculate quality control metrics"""
        adata = adata.copy()

        # Mitochondrial genes
        adata.var['mt'] = adata.var_names.str.startswith('MT-')

        # Ribosomal genes
        adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

        # Hemoglobin genes
        adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(
            adata,
            percent_top=None,
            log1p=False,
            inplace=True
        )
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=['mt', 'ribo', 'hb'],
            percent_top=None,
            log1p=False,
            inplace=True
        )

        return adata

    def scale_for_memory(self, adata: ad.AnnData) -> ad.AnnData:
        """Scale dataset for memory constraints"""
        if self.config.data.max_cells_subset and adata.n_obs > self.config.data.max_cells_subset:
            np.random.seed(self.config.data.random_seed)
            cell_indices = np.random.choice(
                adata.n_obs,
                size=self.config.data.max_cells_subset,
                replace=False
            )
            adata = adata[cell_indices].copy()

        if self.config.data.max_genes_subset and adata.n_vars > self.config.data.max_genes_subset:
            # Select most variable genes
            if issparse(adata.X):
                # For sparse matrices, convert to dense temporarily for variance calculation
                dense_subset = adata.X[:min(1000, adata.n_obs), :].toarray()
                gene_vars = np.var(dense_subset, axis=0)
            else:
                gene_vars = np.var(adata.X, axis=0)

            top_gene_indices = np.argsort(gene_vars)[-self.config.data.max_genes_subset:]
            adata = adata[:, top_gene_indices].copy()

        return adata

    def normalize_and_scale(self, adata: ad.AnnData) -> ad.AnnData:
        """Normalize and scale data"""
        adata = adata.copy()

        # Clean data - remove infinite values
        if issparse(adata.X):
            adata.X.data = np.nan_to_num(adata.X.data, nan=0, posinf=0, neginf=0)
        else:
            adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)

        # Normalize to target sum
        sc.pp.normalize_total(adata, target_sum=self.config.data.target_sum)

        # Log transform
        sc.pp.log1p(adata)

        # Sanitize after log1p (can create NaNs/Inf from edge cases)
        if issparse(adata.X):
            adata.X.data = np.nan_to_num(adata.X.data, nan=0, posinf=0, neginf=0)
        else:
            adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)

        return adata

    def preprocess(self, adata: ad.AnnData) -> ad.AnnData:
        """Complete preprocessing pipeline (convenience method)"""
        adata = self.preprocess_basic(adata)
        adata = self.scale_for_memory(adata)
        adata = self.normalize_and_scale(adata)
        return adata


class DataProcessor:
    """Main data processing class"""

    def __init__(self, config: Config):
        self.config = config
        self.loader = DataLoader(config)

    @staticmethod
    def _sanitize_before_pca(adata: ad.AnnData) -> ad.AnnData:
        """Ensure finite values and remove empty genes/cells before PCA."""
        adata = adata.copy()
        # Replace NaN/Inf with zeros
        if issparse(adata.X):
            import numpy as _np
            data = adata.X.data
            if data.size:
                adata.X.data = _np.nan_to_num(data, nan=0, posinf=0, neginf=0)
        else:
            adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)

        # Drop all-zero genes/cells to avoid zero-variance issues
        try:
            sc.pp.filter_genes(adata, min_counts=1)
            sc.pp.filter_cells(adata, min_counts=1)
        except Exception:
            pass
        return adata

    def process_from_raw(self, file_path: str, save_intermediate: bool = True) -> ad.AnnData:
        """Complete processing pipeline from raw data"""

        # Ensure processed data directory exists
        if save_intermediate:
            os.makedirs(self.config.paths.processed_data_dir, exist_ok=True)

        # Load raw data
        adata = self.loader.load_raw_data(file_path)

        # Basic preprocessing
        adata = self.loader.preprocess_basic(adata)

        if save_intermediate:
            adata.write(os.path.join(
                self.config.paths.processed_data_dir,
                "preprocessed.h5ad"
            ))

        # Calculate QC metrics
        adata = self.loader.calculate_qc_metrics(adata)

        # Validate data with QC metrics
        is_valid, issues = DataValidator.validate_anndata(adata, check_qc_metrics=True)
        if not is_valid:
            warnings.warn(f"Data validation issues after QC calculation: {'; '.join(issues)}")

        if save_intermediate:
            adata.write(os.path.join(
                self.config.paths.processed_data_dir,
                "qc_calculated.h5ad"
            ))

        # Scale for memory if needed
        if (self.config.data.max_cells_subset or
                self.config.data.max_genes_subset):
            adata = self.loader.scale_for_memory(adata)

            if save_intermediate:
                adata.write(os.path.join(
                    self.config.paths.processed_data_dir,
                    "scaled.h5ad"
                ))

        # Normalize and scale
        adata = self.loader.normalize_and_scale(adata)

        if save_intermediate:
            adata.write(os.path.join(
                self.config.paths.processed_data_dir,
                "normalized.h5ad"
            ))

        # Final sanitization before PCA (handles web deployment NaNs)
        adata = self._sanitize_before_pca(adata)

        # Compute PCA for dimensionality reduction
        sc.tl.pca(adata, svd_solver='arpack')

        # Compute neighborhood graph (required for clustering)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

        if save_intermediate:
            adata.write(os.path.join(
                self.config.paths.processed_data_dir,
                "processed_with_neighbors.h5ad"
            ))

        return adata

    def create_test_dataset(self, adata: ad.AnnData, n_cells: int = 1000) -> ad.AnnData:
        """Create small test dataset"""
        np.random.seed(self.config.data.random_seed)

        n_cells = min(n_cells, adata.n_obs)
        cell_indices = np.random.choice(adata.n_obs, size=n_cells, replace=False)

        return adata[cell_indices].copy()


# Import ligand-receptor database module
try:
    from .lr_database import get_ligand_receptor_pairs, LigandReceptorDatabase
    LR_DATABASE_AVAILABLE = True
except ImportError:
    LR_DATABASE_AVAILABLE = False
    warnings.warn("Ligand-receptor database module not available. Install liana for full functionality.")

# Export data processing classes
__all__ = [
    'DataValidator',
    'DataLoader',
    'DataProcessor',
    'get_ligand_receptor_pairs',
    'LigandReceptorDatabase',
    'LR_DATABASE_AVAILABLE'
]
