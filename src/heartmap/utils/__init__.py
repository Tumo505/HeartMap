"""
Utility functions and classes for HeartMAP
"""

import hashlib
from pathlib import Path
from typing import Dict, Union
import warnings

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import scanpy as sc
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    warnings.warn("Plotting dependencies not available")

from ..config import Config


class Visualizer:
    """Visualization utilities for HeartMAP"""

    def __init__(self, config: Config):
        self.config = config

    def plot_qc_metrics(self, adata, save_dir: Path) -> None:
        """Plot quality control metrics"""
        if not PLOTTING_AVAILABLE:
            return

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Number of genes
        axes[0, 0].hist(adata.obs['n_genes'], bins=50, alpha=0.7)
        axes[0, 0].set_xlabel('Number of genes')
        axes[0, 0].set_ylabel('Number of cells')
        axes[0, 0].set_title('Genes per cell')

        # Total counts
        axes[0, 1].hist(adata.obs['total_counts'], bins=50, alpha=0.7)
        axes[0, 1].set_xlabel('Total counts')
        axes[0, 1].set_ylabel('Number of cells')
        axes[0, 1].set_title('UMI per cell')

        # Mitochondrial percentage
        if 'pct_counts_mt' in adata.obs.columns:
            axes[1, 0].hist(adata.obs['pct_counts_mt'], bins=50, alpha=0.7)
            axes[1, 0].set_xlabel('Mitochondrial %')
            axes[1, 0].set_ylabel('Number of cells')
            axes[1, 0].set_title('Mitochondrial gene %')

        # Scatter plot
        axes[1, 1].scatter(adata.obs['total_counts'], adata.obs['n_genes'], alpha=0.6, s=1)
        axes[1, 1].set_xlabel('Total counts')
        axes[1, 1].set_ylabel('Number of genes')
        axes[1, 1].set_title('Genes vs UMI')

        plt.tight_layout()
        plt.savefig(save_dir / "qc_metrics.png", dpi=300, bbox_inches='tight')
        plt.close()

    def plot_communication_heatmap(self, comm_scores: pd.DataFrame, save_dir: Path) -> None:
        """Plot cell-cell communication heatmap"""
        if not PLOTTING_AVAILABLE:
            return

        # Create pivot table
        pivot_df = comm_scores.pivot(
            index='source',
            columns='target',
            values='communication_score'
        )

        plt.figure(figsize=(10, 8))
        sns.heatmap(pivot_df, annot=True, cmap='viridis', fmt='.3f')
        plt.title('Cell-Cell Communication Specificity')
        plt.tight_layout()
        plt.savefig(save_dir / "communication_heatmap.png", dpi=300, bbox_inches='tight')
        plt.close()

    def plot_hub_scores(self, adata, hub_scores: pd.Series, save_dir: Path) -> None:
        """Plot communication hub scores"""
        if not PLOTTING_AVAILABLE:
            return

        # Add hub scores to adata for plotting
        adata.obs['hub_score'] = hub_scores

        # Ensure UMAP exists
        if 'X_umap' not in adata.obsm:
            import scanpy as sc
            sc.tl.umap(adata)

        sc.pl.umap(adata, color='hub_score', title='Communication Hub Score', show=False)
        plt.savefig(save_dir / "hub_scores.png", dpi=300, bbox_inches='tight')
        plt.close()

    def plot_pathway_scores(self, pathway_scores: pd.DataFrame, save_dir: Path) -> None:
        """Plot pathway activity scores"""
        if not PLOTTING_AVAILABLE or pathway_scores.empty:
            return

        plt.figure(figsize=(12, 8))
        sns.heatmap(pathway_scores, annot=True, cmap='Blues', fmt='.3f')
        plt.title('Pathway Activity by Cell Type')
        plt.tight_layout()
        plt.savefig(save_dir / "pathway_scores.png", dpi=300, bbox_inches='tight')
        plt.close()

    def plot_chamber_composition(self, adata, save_dir: Path) -> None:
        """Plot chamber composition"""
        if not PLOTTING_AVAILABLE:
            return

        if 'chamber' not in adata.obs.columns:
            return

        chamber_counts = adata.obs['chamber'].value_counts()

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Bar plot
        chamber_counts.plot(kind='bar', ax=ax1)
        ax1.set_title('Cell Counts by Chamber')
        ax1.set_xlabel('Chamber')
        ax1.set_ylabel('Number of Cells')

        # Pie chart
        ax2.pie(chamber_counts.values, labels=chamber_counts.index, autopct='%1.1f%%')
        ax2.set_title('Chamber Proportions')

        plt.tight_layout()
        plt.savefig(save_dir / "chamber_composition.png", dpi=300, bbox_inches='tight')
        plt.close()

    def plot_chamber_markers(self, chamber_markers: Dict, save_dir: Path) -> None:
        """Plot chamber-specific markers"""
        if not PLOTTING_AVAILABLE:
            return

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.flatten()

        for i, (chamber, markers) in enumerate(chamber_markers.items()):
            if i < 4 and not markers.empty:
                top_10 = markers.head(10)
                axes[i].barh(range(len(top_10)), -np.log10(top_10['pvals_adj']))
                axes[i].set_yticks(range(len(top_10)))
                axes[i].set_yticklabels(top_10['names'])
                axes[i].set_xlabel('-log10(adjusted p-value)')
                axes[i].set_title(f'Top Markers - {chamber}')

        plt.tight_layout()
        plt.savefig(save_dir / "chamber_markers.png", dpi=300, bbox_inches='tight')
        plt.close()

    def plot_cross_chamber_correlations(self, correlations: pd.DataFrame, save_dir: Path) -> None:
        """Plot cross-chamber correlations"""
        if not PLOTTING_AVAILABLE:
            return

        if correlations is None or getattr(correlations, 'empty', True):
            return

        plt.figure(figsize=(8, 6))
        sns.heatmap(correlations, annot=True, cmap='coolwarm', center=0, fmt='.3f')
        plt.title('Cross-Chamber Expression Correlations')
        plt.tight_layout()
        plt.savefig(save_dir / "cross_chamber_correlations.png", dpi=300, bbox_inches='tight')
        plt.close()

    def create_comprehensive_dashboard(self, adata, results: Dict, save_dir: Path) -> None:
        """Create comprehensive analysis dashboard"""
        if not PLOTTING_AVAILABLE:
            return

        # Create a large multi-panel figure
        plt.figure(figsize=(20, 16))

        # Panel 1: UMAP with clusters
        ax1 = plt.subplot(3, 3, 1)
        sc.pl.umap(adata, color='leiden', ax=ax1, show=False, frameon=False)
        ax1.set_title('Cell Type Clusters')

        # Panel 2: UMAP with hub scores
        ax2 = plt.subplot(3, 3, 2)
        sc.pl.umap(adata, color='hub_score', ax=ax2, show=False, frameon=False)
        ax2.set_title('Communication Hubs')

        # Panel 3: Chamber composition (if available)
        if 'chamber' in adata.obs.columns:
            ax3 = plt.subplot(3, 3, 3)
            chamber_counts = adata.obs['chamber'].value_counts()
            ax3.pie(chamber_counts.values, labels=chamber_counts.index, autopct='%1.1f%%')
            ax3.set_title('Chamber Distribution')

        # Additional panels for other analyses...

        plt.tight_layout()
        plt.savefig(save_dir / "comprehensive_dashboard.png", dpi=300, bbox_inches='tight')
        plt.close()


class ResultsExporter:
    """Export analysis results in various formats"""

    def __init__(self, config: Config):
        self.config = config

    def export_results(self, results: Dict, output_dir: Path) -> None:
        """Export all results to files"""
        output_dir.mkdir(parents=True, exist_ok=True)

        # Export tables
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(exist_ok=True)

        # Export specific result types
        if 'results' in results and 'marker_genes' in results['results']:
            marker_genes = results['results']['marker_genes']
            if marker_genes is not None:
                marker_genes.to_csv(tables_dir / "marker_genes.csv")

        if 'results' in results and 'communication_scores' in results['results']:
            comm_scores = results['results']['communication_scores']
            if comm_scores is not None:
                comm_scores.to_csv(tables_dir / "communication_scores.csv", index=False)

    def generate_comprehensive_report(self, results: Dict, output_dir: str) -> None:
        """Generate comprehensive analysis report"""
        output_path = Path(output_dir)

        # Get basic statistics
        adata = results.get('adata')
        if adata is None:
            return

        n_cells = adata.n_obs
        n_genes = adata.n_vars

        # Generate report
        report = f"""# HeartMAP Analysis Report

## Dataset Overview
- **Total Cells**: {n_cells:,}
- **Total Genes**: {n_genes:,}

## Analysis Components Completed
- ✅ Cell type annotation
- ✅ Cell-cell communication analysis
- ✅ Multi-chamber analysis

## Key Findings

### Cell Type Annotation
"""

        if 'leiden' in adata.obs.columns:
            cluster_counts = adata.obs['leiden'].value_counts()
            report += f"- **Number of cell types identified**: {len(cluster_counts)}\n"
            report += "- **Cell type distribution**:\n"
            for cluster, count in cluster_counts.head(5).items():
                pct = 100 * count / n_cells
                report += f"  - Cluster {cluster}: {count:,} cells ({pct:.1f}%)\n"

        if 'chamber' in adata.obs.columns:
            chamber_counts = adata.obs['chamber'].value_counts()
            report += "\n### Chamber Distribution\n"
            for chamber, count in chamber_counts.items():
                pct = 100 * count / n_cells
                report += f"- **{chamber}**: {count:,} cells ({pct:.1f}%)\n"

        report += """
### Communication Analysis
- Cell-cell communication patterns identified
- Communication hub cells detected
- Pathway activity scores calculated

## Files Generated
- `heartmap_complete.h5ad`: Complete processed dataset
- `heartmap_model.pkl`: Trained HeartMAP model
- `figures/`: All visualization outputs
- `tables/`: Exported data tables

## Next Steps
1. Validate findings with literature
2. Investigate specific cell type interactions
3. Apply model to new datasets
"""

        # Save report
        with open(output_path / "analysis_report.md", 'w') as f:
            f.write(report)


class ChecksumValidator:
    """Validate data integrity using checksums"""

    @staticmethod
    def calculate_sha256(file_path: str) -> str:
        """Calculate SHA-256 checksum of file"""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()

    @staticmethod
    def verify_file(file_path: str, expected_checksum: str) -> bool:
        """Verify file against expected checksum"""
        actual_checksum = ChecksumValidator.calculate_sha256(file_path)
        return actual_checksum == expected_checksum

    @staticmethod
    def generate_checksums(directory: Union[str, Path], output_file: str) -> None:
        """Generate checksums for all files in directory"""
        directory = Path(directory)
        checksums = {}

        for file_path in directory.rglob('*'):
            if file_path.is_file():
                rel_path = file_path.relative_to(directory)
                checksum = ChecksumValidator.calculate_sha256(str(file_path))
                checksums[str(rel_path)] = checksum

        # Save checksums
        with open(output_file, 'w') as f:
            for rel_path_str, checksum in checksums.items():
                f.write(f"{checksum}  {rel_path_str}\n")


def setup_logging(level: str = "INFO") -> None:
    """Setup logging configuration"""
    import logging

    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('heartmap.log')
        ]
    )


def create_logger(name: str):
    """Create a logger with the given name"""
    import logging
    return logging.getLogger(name)


# Export utility classes and functions
__all__ = [
    'Visualizer',
    'ResultsExporter',
    'ChecksumValidator',
    'setup_logging',
    'create_logger'
]
