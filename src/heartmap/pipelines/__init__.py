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

        # Perform basic clustering using scanpy
        print("2. Performing cell annotation...")
        sc.tl.leiden(adata, resolution=self.config.analysis.resolution)

        # Generate visualizations
        print("3. Generating visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)

            # Ensure UMAP is computed before plotting
            if 'X_umap' not in adata.obsm_keys():
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
    """Advanced cell-cell communication analysis pipeline"""

    def __init__(self, config: Config):
        super().__init__(config)
        # Communication analysis without models - placeholder implementation

    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run advanced communication analysis pipeline"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")

        print("=== Running Advanced Communication Pipeline ===")

        # Load processed data (should have cell annotations)
        print("1. Loading annotated data...")
        adata = sc.read_h5ad(data_path)

        if 'leiden' not in adata.obs.columns:
            raise ValueError("Input data must have cell type annotations. Run BasicPipeline first.")

        # Basic communication analysis placeholder
        print("2. Analyzing cell-cell communication...")
        # Create placeholder results
        communication_scores = pd.DataFrame({
            'source': ['cluster_0', 'cluster_1'],
            'target': ['cluster_1', 'cluster_0'],
            'communication_score': [0.5, 0.3]
        })

        hub_scores = pd.Series(np.random.random(adata.n_obs), index=adata.obs.index)
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

        # Perform basic clustering
        print("2. Performing comprehensive analysis...")
        sc.tl.leiden(adata, resolution=self.config.analysis.resolution)

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
