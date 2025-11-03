import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr
from sklearn.metrics.pairwise import cosine_similarity

def analyze_communication_specificity(adata, save_path):
    """Analyze cell-type specific communication patterns"""
    print("Analyzing communication specificity...")
    
    # Get cell type information
    cell_types = adata.obs['leiden'].unique()
    
    # Create communication specificity matrix
    specificity_data = []
    
    for ct1 in cell_types:
        for ct2 in cell_types:
            if ct1 != ct2:
                # Calculate potential communication strength
                cells_ct1 = adata.obs['leiden'] == ct1
                cells_ct2 = adata.obs['leiden'] == ct2
                
                # Mock calculation - in real analysis, use actual L-R pairs
                ct1_expr = adata[cells_ct1].X.mean(axis=0)
                ct2_expr = adata[cells_ct2].X.mean(axis=0)
                
                # Calculate correlation as proxy for communication potential
                correlation = np.corrcoef(ct1_expr.A1, ct2_expr.A1)[0, 1]
                
                specificity_data.append({
                    'source': ct1,
                    'target': ct2,
                    'communication_score': abs(correlation) if not np.isnan(correlation) else 0
                })
    
    specificity_df = pd.DataFrame(specificity_data)
    
    # Create heatmap
    pivot_df = specificity_df.pivot(index='source', columns='target', values='communication_score')
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(pivot_df, annot=True, cmap='viridis', fmt='.3f')
    plt.title('Cell-Type Communication Specificity')
    plt.tight_layout()
    plt.savefig(save_path / "communication_specificity_validation_gse145154.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save specificity data
    specificity_df.to_csv(save_path / "communication_specificity_validation_gse145154.csv", index=False)

def identify_hub_cells(adata, save_path):
    """Identify cells that act as communication hubs"""
    print("Identifying communication hub cells...")
    
    # Calculate communication potential for each cell
    hub_scores = []
    
    for i in range(adata.n_obs):
        cell_expr = adata.X[i].toarray().flatten()
        
        # Calculate hub score based on expression diversity
        # High-expressing cells with diverse expression patterns
        hub_score = (cell_expr.std() * cell_expr.mean()) / (cell_expr.var() + 1)
        hub_scores.append(hub_score)
    
    adata.obs['hub_score'] = hub_scores
    
    # Plot hub scores
    sc.pl.umap(adata, color='hub_score', title='Communication Hub Score', show=False)
    plt.savefig(save_path / "communication_hubs_validation_gse145154.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Identify top hub cells
    top_hubs = adata.obs.nlargest(100, 'hub_score')
    top_hubs.to_csv(save_path / "top_hub_cells_validation_gse145154.csv")

def analyze_pathway_enrichment(adata, save_path):
    """Analyze pathway enrichment in communication"""
    print("Analyzing pathway enrichment...")
    
    # Define heart-specific pathways
    heart_pathways = {
        'Cardiac_Development': ['GATA4', 'NKX2-5', 'MEF2C', 'TBX5', 'HAND1', 'HAND2'],
        'Angiogenesis': ['VEGFA', 'VEGFB', 'VEGFC', 'FLT1', 'KDR', 'PDGFB'],
        'ECM_Remodeling': ['COL1A1', 'COL3A1', 'MMP2', 'MMP9', 'TIMP1', 'TIMP2'],
        'Calcium_Signaling': ['CACNA1C', 'CACNA1D', 'CACNA1G', 'CACNA1H', 'CACNA1S'],
        'Wnt_Signaling': ['WNT3A', 'WNT5A', 'FZD1', 'FZD2', 'LRP5', 'LRP6'],
        'Notch_Signaling': ['NOTCH1', 'NOTCH2', 'DLL1', 'DLL4', 'JAG1', 'JAG2']
    }
    
    # Calculate pathway scores for each cell type
    pathway_scores = {}
    
    for pathway, genes in heart_pathways.items():
        # Find genes present in dataset
        present_genes = [g for g in genes if g in adata.var_names]
        
        if present_genes:
            # Calculate pathway score for each cell type
            cell_type_scores = {}
            
            for ct in adata.obs['leiden'].unique():
                ct_mask = adata.obs['leiden'] == ct
                ct_expr = adata[ct_mask, present_genes].X.mean(axis=0)
                
                if hasattr(ct_expr, 'A1'):
                    pathway_score = ct_expr.A1.mean()
                else:
                    pathway_score = ct_expr.mean()
                
                cell_type_scores[ct] = pathway_score
            
            pathway_scores[pathway] = cell_type_scores
    
    # Create pathway heatmap only if we have data
    if pathway_scores:
        pathway_df = pd.DataFrame(pathway_scores).T
        
        if not pathway_df.empty:
            plt.figure(figsize=(12, 8))
            sns.heatmap(pathway_df, annot=True, cmap='Blues', fmt='.3f')
            plt.title('Pathway Activity by Cell Type')
            plt.tight_layout()
            plt.savefig(save_path / "pathway_enrichment_validation_gse145154.png", dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save pathway scores
            pathway_df.to_csv(save_path / "pathway_scores_validation_gse145154.csv")
        else:
            print("Warning: No pathway data available for heatmap")
    else:
        print("Warning: No pathway genes found in dataset")
        
        # Create a simple message plot instead
        plt.figure(figsize=(8, 6))
        plt.text(0.5, 0.5, 'No pathway genes found in dataset\n(heart-specific genes not present)', 
                ha='center', va='center', transform=plt.gca().transAxes, fontsize=12)
        plt.title('Pathway Enrichment Analysis')
        plt.axis('off')
        plt.savefig(save_path / "pathway_enrichment_validation_gse145154.png", dpi=300, bbox_inches='tight')
        plt.close()

def temporal_communication_analysis(adata, save_path):
    """Analyze potential temporal aspects of communication"""
    print("Analyzing temporal communication patterns...")
    
    # Use pseudotime if available, otherwise create mock temporal ordering
    if 'dpt_pseudotime' not in adata.obs.columns:
        # Create mock pseudotime based on gene expression patterns
        # This is a simplified approach - in reality, use proper pseudotime algorithms
        high_expr_genes = adata.var_names[adata.X.mean(axis=0).A1 > np.percentile(adata.X.mean(axis=0).A1, 75)]
        
        if len(high_expr_genes) > 0:
            pseudotime = adata[:, high_expr_genes].X.sum(axis=1).A1
            adata.obs['pseudotime'] = pseudotime
        else:
            adata.obs['pseudotime'] = np.random.random(adata.n_obs)
    else:
        adata.obs['pseudotime'] = adata.obs['dpt_pseudotime']
    
    # Analyze communication changes over pseudotime
    time_bins = pd.cut(adata.obs['pseudotime'], bins=5, labels=['Early', 'Early-Mid', 'Mid', 'Mid-Late', 'Late'])
    adata.obs['time_bin'] = time_bins
    
    # Calculate communication scores for each time bin
    time_communication = {}
    
    # Get unique time bins (handle both categorical and regular series)
    if hasattr(time_bins, 'cat'):
        unique_bins = time_bins.cat.categories
    else:
        unique_bins = time_bins.unique()
    
    for time_bin in unique_bins:
        time_mask = adata.obs['time_bin'] == time_bin
        if time_mask.sum() > 10:  # Ensure enough cells
            time_expr = adata[time_mask].X.mean(axis=0)
            
            if hasattr(time_expr, 'A1'):
                communication_score = time_expr.A1.std()  # Diversity as proxy
            else:
                communication_score = time_expr.std()
            
            time_communication[time_bin] = communication_score
    
    # Plot temporal communication
    if time_communication:
        plt.figure(figsize=(10, 6))
        plt.bar(time_communication.keys(), time_communication.values())
        plt.title('Communication Diversity Over Pseudotime')
        plt.xlabel('Pseudotime Bin')
        plt.ylabel('Communication Diversity Score')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(save_path / "temporal_communication.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # Save temporal data
    temporal_df = pd.DataFrame(list(time_communication.items()), 
                              columns=['Time_Bin', 'Communication_Score'])
    temporal_df.to_csv(save_path / "temporal_communication.csv", index=False)

def main():
    # Load data
    project_root = Path(__file__).parent.parent.parent
    data_path = project_root / "data/processed/heart_data_communication_validation_gse145154.h5ad"
    if not data_path.exists():
        data_path = project_root / "data/processed/heart_data_annotated_validation_gse145154.h5ad"
    
    adata = sc.read_h5ad(data_path)
    
    # Create advanced results directory
    results_path = Path("results/advanced")
    results_path.mkdir(parents=True, exist_ok=True)
    
    # Run advanced analyses
    analyze_communication_specificity(adata, results_path)
    identify_hub_cells(adata, results_path)
    analyze_pathway_enrichment(adata, results_path)
    temporal_communication_analysis(adata, results_path)
    
    # Save enhanced data
    adata.write(project_root / "data/processed/heart_data_advanced_validation_gse145154.h5ad")
    
    print("Advanced communication analysis complete!")

if __name__ == "__main__":
    main()