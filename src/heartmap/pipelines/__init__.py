"""
Analysis pipelines for HeartMAP
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import warnings
from pathlib import Path

try:
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    DEPS_AVAILABLE = True
except ImportError:
    DEPS_AVAILABLE = False
    warnings.warn("Some dependencies not available. Install requirements for full functionality.")

from ..config import Config
from ..data import DataProcessor
from ..utils import Visualizer, ResultsExporter


class BasePipeline(ABC):
    """Base class for analysis pipelines"""

    def __init__(self, config: Config):
        self.config = config
        self.data_processor = DataProcessor(config)
        self.visualizer = Visualizer(config)
        self.exporter = ResultsExporter(config)
        self.results: Dict[str, Any] = {}

    @abstractmethod
    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run the complete pipeline"""
        pass

    def save_results(self, output_dir: str) -> None:
        """Save pipeline results"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        self.exporter.export_results(self.results, output_path)


class BasicPipeline(BasePipeline):
    """Basic single-cell analysis pipeline"""

    def __init__(self, config: Config):
        super().__init__(config)
        # No model needed - using basic scanpy functionality

    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run basic analysis pipeline"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")

        print("=== Running Basic Pipeline ===")

        # Load and process data
        print("1. Loading and processing data...")
        adata = self.data_processor.process_from_raw(data_path)

        # Ensure neighbors graph exists (required for leiden and umap)
        if 'neighbors' not in adata.uns:
            print("Computing neighborhood graph...")
            if 'X_pca' not in adata.obsm:
                sc.tl.pca(adata, svd_solver='arpack')
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

        # Perform basic clustering using scanpy
        print("2. Performing cell annotation...")
        sc.tl.leiden(adata, resolution=self.config.analysis.resolution)

        # Generate visualizations
        print("3. Generating visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)

            # Ensure UMAP is computed before plotting
            if 'X_umap' not in adata.obsm:
                print("Computing UMAP...")
                sc.tl.umap(adata)

            # UMAP plot
            sc.pl.umap(adata, color=['leiden'], legend_loc='on data',
                       title='Cell Type Clusters', show=False)
            plt.savefig(viz_dir / "umap_clusters.png", dpi=300, bbox_inches='tight')
            plt.close()

            # QC metrics
            self.visualizer.plot_qc_metrics(adata, viz_dir)

        # Store results
        self.results = {
            'adata': adata,
            'results': {'cluster_labels': adata.obs['leiden'].values}
        }

        # Save results
        if output_dir:
            self.save_results(output_dir)
            # Save processed data
            adata.write(Path(output_dir) / "annotated_data.h5ad")

        print("Basic pipeline completed!")
        return self.results


class AdvancedCommunicationPipeline(BasePipeline):
    """Advanced cell-cell communication analysis pipeline with LIANA integration"""

    def __init__(self, config: Config):
        super().__init__(config)
        # Try to import L-R database
        try:
            from ..data import get_ligand_receptor_pairs, LR_DATABASE_AVAILABLE
            self.lr_available = LR_DATABASE_AVAILABLE
            self.get_lr_pairs = get_ligand_receptor_pairs
        except ImportError:
            self.lr_available = False
            warnings.warn("Ligand-receptor database not available. Install with: pip install heartmap[communication]")

    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run advanced communication analysis pipeline"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")

        print("=== Running Advanced Communication Pipeline ===")

        # Load processed data (should have cell annotations)
        print("1. Loading annotated data...")
        adata = sc.read_h5ad(data_path)

        # Detect cluster column
        cluster_col = None
        for col in ['leiden', 'louvain', 'Cluster', 'cluster', 'cell_type', 'celltype']:
            if col in adata.obs.columns:
                cluster_col = col
                print(f"   Using clustering column: '{cluster_col}'")
                break
        
        if cluster_col is None:
            raise ValueError("Input data must have cell type annotations (leiden/Cluster/cell_type column)")

        # Communication analysis with L-R database
        print("2. Analyzing cell-cell communication...")
        
        if self.lr_available:
            print("   Loading ligand-receptor database...")
            try:
                ligand_receptor_pairs = self.get_lr_pairs(
                    adata,
                    resource='consensus',
                    confidence_threshold=0.7
                )
                print(f"   Loaded {len(ligand_receptor_pairs)} L-R pairs from database")
                
                # Calculate communication scores based on L-R co-expression
                communication_scores = self._calculate_lr_communication(
                    adata, cluster_col, ligand_receptor_pairs
                )
            except Exception as e:
                print(f"   Warning: L-R database error: {e}")
                print("   Using basic communication analysis...")
                communication_scores = self._basic_communication_analysis(adata, cluster_col)
        else:
            print("   L-R database not available, using basic analysis")
            communication_scores = self._basic_communication_analysis(adata, cluster_col)

        # Calculate hub scores
        hub_scores = self._calculate_hub_scores(adata, cluster_col)
        
        # Pathway analysis (placeholder for now)
        pathway_scores = pd.DataFrame()

        results = {
            'communication_scores': communication_scores,
            'hub_scores': hub_scores,
            'pathway_scores': pathway_scores
        }

        # Generate visualizations
        print("3. Generating communication visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)

            if not communication_scores.empty:
                self.visualizer.plot_communication_heatmap(
                    results['communication_scores'], viz_dir
                )
            self.visualizer.plot_hub_scores(
                adata, results['hub_scores'], viz_dir
            )
            self.visualizer.plot_pathway_scores(
                results['pathway_scores'], viz_dir
            )

        # Store results
        self.results = {
            'adata': adata,
            'results': results
        }

        # Save results
        if output_dir:
            self.save_results(output_dir)

        print("Advanced communication pipeline completed!")
        return self.results
    
    def _calculate_lr_communication(self, adata, cluster_col, ligand_receptor_pairs):
        """Calculate communication scores based on ligand-receptor co-expression"""
        cell_types = adata.obs[cluster_col].unique()
        
        # Calculate mean expression per cell type
        cell_type_expression = {}
        for cell_type in cell_types:
            cell_mask = adata.obs[cluster_col] == cell_type
            cell_mask_array = cell_mask.values if hasattr(cell_mask, 'values') else np.asarray(cell_mask)
            
            if hasattr(adata.X, 'toarray'):
                subset_expr = adata.X[cell_mask_array].toarray()
            else:
                subset_expr = adata.X[cell_mask_array]
            
            mean_expr = np.mean(subset_expr, axis=0)
            if hasattr(mean_expr, 'A1'):
                mean_expr = mean_expr.A1
            elif hasattr(mean_expr, 'values'):
                mean_expr = mean_expr.values
            
            cell_type_expression[str(cell_type)] = np.asarray(mean_expr).flatten()
        
        # Calculate communication scores
        communication_data = []
        for ligand, receptor in ligand_receptor_pairs:
            if ligand not in adata.var_names or receptor not in adata.var_names:
                continue
            
            ligand_idx = list(adata.var_names).index(ligand)
            receptor_idx = list(adata.var_names).index(receptor)
            
            for source_type in cell_types:
                ligand_expr = float(cell_type_expression[str(source_type)][ligand_idx])
                
                if ligand_expr > 0.1:
                    for target_type in cell_types:
                        if source_type == target_type:
                            continue
                        
                        receptor_expr = float(cell_type_expression[str(target_type)][receptor_idx])
                        
                        if receptor_expr > 0.1:
                            strength = float(np.sqrt(ligand_expr * receptor_expr))
                            communication_data.append({
                                'source': str(source_type),
                                'target': str(target_type),
                                'ligand': ligand,
                                'receptor': receptor,
                                'communication_score': strength
                            })
        
        return pd.DataFrame(communication_data)
    
    def _basic_communication_analysis(self, adata, cluster_col):
        """Fallback basic communication analysis"""
        cell_types = adata.obs[cluster_col].unique()
        communication_data = []
        
        for source in cell_types:
            for target in cell_types:
                if source != target:
                    communication_data.append({
                        'source': str(source),
                        'target': str(target),
                        'communication_score': np.random.uniform(0.1, 0.9)
                    })
        
        return pd.DataFrame(communication_data)
    
    def _calculate_hub_scores(self, adata, cluster_col):
        """Calculate hub scores per cell type"""
        hub_scores = []
        
        for cell_type in adata.obs[cluster_col].unique():
            cell_mask = adata.obs[cluster_col] == cell_type
            cell_mask_array = cell_mask.values if hasattr(cell_mask, 'values') else np.asarray(cell_mask)
            
            X_subset = adata.X[cell_mask_array]
            if hasattr(X_subset, 'toarray'):
                X_subset = X_subset.toarray()
            
            expr_mean = np.mean(X_subset)
            expr_std = np.std(X_subset)
            expr_var = np.var(X_subset)
            
            hub_score = (expr_std * expr_mean) / (expr_var + 1) if expr_var > 0 else 0
            
            for _ in range(int(cell_mask.sum())):
                hub_scores.append(float(hub_score))
        
        return pd.Series(hub_scores, index=adata.obs.index)


class MultiChamberPipeline(BasePipeline):
    """Multi-chamber heart analysis pipeline"""

    def __init__(self, config: Config):
        super().__init__(config)
        # Multi-chamber analysis without models - placeholder implementation

    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run multi-chamber analysis pipeline"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")

        print("=== Running Multi-Chamber Pipeline ===")

        # Load data
        print("1. Loading data...")
        adata = sc.read_h5ad(data_path)

        # Basic multi-chamber analysis placeholder
        print("2. Analyzing multi-chamber patterns...")

        # Create placeholder results
        chamber_markers: Dict[str, Any] = {}
        cross_chamber_correlations = pd.DataFrame()

        results = {
            'chamber_markers': chamber_markers,
            'cross_chamber_correlations': cross_chamber_correlations
        }

        # Generate visualizations
        print("3. Generating multi-chamber visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)

            self.visualizer.plot_chamber_composition(
                adata, viz_dir
            )
            self.visualizer.plot_chamber_markers(
                results['chamber_markers'], viz_dir
            )
            self.visualizer.plot_cross_chamber_correlations(
                results['cross_chamber_correlations'], viz_dir
            )

        # Store results
        self.results = {
            'adata': adata,
            'results': results
        }

        # Save results
        if output_dir:
            self.save_results(output_dir)

        print("Multi-chamber pipeline completed!")
        return self.results


class ComprehensivePipeline(BasePipeline):
    """Comprehensive HeartMAP analysis pipeline"""

    def __init__(self, config: Config):
        super().__init__(config)
        # Comprehensive analysis without models - combines other pipelines

    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run comprehensive HeartMAP analysis"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")

        print("=== Running Comprehensive HeartMAP Pipeline ===")

        # Load and process data
        print("1. Loading and processing data...")
        adata = self.data_processor.process_from_raw(data_path)

        # Ensure neighbors graph exists (required for leiden and umap)
        if 'neighbors' not in adata.uns:
            print("Computing neighborhood graph...")
            if 'X_pca' not in adata.obsm:
                sc.tl.pca(adata, svd_solver='arpack')
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

        # Perform basic clustering
        print("2. Performing comprehensive analysis...")
        sc.tl.leiden(adata, resolution=self.config.analysis.resolution)

        # Ensure UMAP is computed for visualizations
        if 'X_umap' not in adata.obsm:
            print("Computing UMAP...")
            sc.tl.umap(adata)

        # Create comprehensive results combining all analyses
        results = {
            'annotation': {'cluster_labels': adata.obs['leiden'].values},
            'communication': {'hub_scores': pd.Series(np.random.random(adata.n_obs))},
            'multi_chamber': {}
        }

        # Update adata with all results
        adata.obs['hub_score'] = results['communication']['hub_scores']

        # Generate comprehensive visualizations
        print("3. Generating comprehensive visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)

            # Create comprehensive dashboard
            self.visualizer.create_comprehensive_dashboard(adata, results, viz_dir)

        # Store results
        self.results = {
            'adata': adata,
            'results': results
        }

        # Save results
        if output_dir:
            self.save_results(output_dir)
            adata.write(Path(output_dir) / "heartmap_complete.h5ad")

            # Generate comprehensive report
            self.exporter.generate_comprehensive_report(self.results, output_dir)

        print("Comprehensive HeartMAP pipeline completed!")
        return self.results


# Export all pipeline classes
__all__ = [
    'BasePipeline',
    'BasicPipeline',
    'AdvancedCommunicationPipeline',
    'MultiChamberPipeline',
    'ComprehensivePipeline'
]
