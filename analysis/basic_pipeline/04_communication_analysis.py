import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
from scipy.stats import norm

import sys
from pathlib import Path

# Add the project root to the Python path
project_root = Path(__file__).parent.parent.parent
sys.path.append(str(project_root))

from utils.statistical_tests import van_elteren_test
warnings.filterwarnings('ignore')

# import cell communication tools

# Check Python environment
import sys
print(f"Python executable: {sys.executable}")
print(f"Python path: {sys.path[:3]}...")

# Check if LIANA is available
LIANA_AVAILABLE = False
try:
    import liana
    import liana as li
    LIANA_AVAILABLE = True
    print(f"LIANA version {liana.__version__} successfully imported")
except ImportError as e:
    print(f"LIANA import failed: {e}")
    print("Install with: pip install liana")
except Exception as e:
    print(f"Unexpected error importing LIANA: {e}")
    print("Install with: pip install liana")

def _load_counts_adata(processed_adata, project_root: Path):
    """Try to load true counts from merged raw h5ad; fall back to .raw or X.

    Returns an AnnData with obs aligned to processed_adata and raw set.
    """
    counts_adata = None
    raw_path = project_root / 'data' / 'validation' / 'gse145154_merged.h5ad'
    try:
        if raw_path.exists():
            base = sc.read_h5ad(raw_path)
            # Align to processed obs
            common = processed_adata.obs_names.intersection(base.obs_names)
            if len(common) > 0:
                base = base[common].copy()
                processed_subset = processed_adata[common].copy()
            else:
                # If no overlap, just use processed
                base = processed_adata.copy()
                processed_subset = processed_adata
            base.obs = processed_subset.obs.copy()
            counts_adata = base
    except Exception:
        counts_adata = None

    if counts_adata is None:
        # Fallbacks
        if processed_adata.raw is not None:
            counts_adata = processed_adata.raw.to_adata()
            counts_adata.obs = processed_adata.obs.copy()
        else:
            counts_adata = processed_adata.copy()

    if counts_adata.raw is None:
        counts_adata.raw = counts_adata
    return counts_adata

def prepare_data_for_communication(adata, cell_type_col='leiden', sample_key='sample'):
    """Prepare data for cell communication analysis"""
    print("Preparing data for communication analysis...")
    
    # Ensure we have cell type annotations
    if cell_type_col not in adata.obs.columns:
        print(f"Warning: {cell_type_col} not found in observations")
        return None
    
    # Ensure we have a sample key for LIANA by-sample APIs
    if sample_key not in adata.obs.columns:
        print(f"'{sample_key}' not found in observations; creating a single-sample column.")
        adata.obs[sample_key] = 'sample_0'

    # Cast grouping keys to categorical
    try:
        adata.obs[cell_type_col] = adata.obs[cell_type_col].astype('category')
        adata.obs[sample_key] = adata.obs[sample_key].astype('category')
    except Exception:
        pass

    # Prefer true counts from merged raw if present
    project_root = Path(__file__).parent.parent.parent
    counts_adata = _load_counts_adata(adata, project_root)
    
    return counts_adata

def run_liana_analysis(adata, cell_type_col='leiden', sample_key='sample'):
    """Run LIANA cell communication analysis"""
    if not LIANA_AVAILABLE:
        print("LIANA not available. Skipping LIANA analysis.")
        return adata
    
    print("Running LIANA analysis...")
    
    try:
        # Use by-sample only when there are >= 2 samples; otherwise, use non-by-sample
        n_samples = adata.obs[sample_key].nunique() if sample_key in adata.obs.columns else 0
        if n_samples and n_samples >= 2:
            res = li.mt.rank_aggregate.by_sample(
                adata,
                groupby=cell_type_col,
                sample_key=sample_key,
                resource_name='consensus',
                n_perms=100,
                seed=42,
                verbose=True
            )
        else:
            res = li.mt.rank_aggregate(
                adata,
                groupby=cell_type_col,
                resource_name='consensus',
                n_perms=100,
                seed=42,
                verbose=True
            )
        # LIANA+ may return various structures; normalize to DataFrame
        stored = False
        try:
            import pandas as _pd  # type: ignore
            if res is not None:
                if isinstance(res, _pd.DataFrame):
                    adata.uns['liana_res'] = res
                    stored = True
                elif isinstance(res, dict):
                    # Common patterns: {'liana_res': df} or {'sampleA': df, ...}
                    if 'liana_res' in res and isinstance(res['liana_res'], _pd.DataFrame):
                        adata.uns['liana_res'] = res['liana_res']
                        stored = True
                    else:
                        dfs = [v for v in res.values() if isinstance(v, _pd.DataFrame) and not v.empty]
                        if dfs:
                            adata.uns['liana_res'] = _pd.concat(dfs, axis=0, ignore_index=True)
                            stored = True
                elif isinstance(res, list):
                    dfs = [x for x in res if isinstance(x, _pd.DataFrame) and not x.empty]
                    if dfs:
                        adata.uns['liana_res'] = _pd.concat(dfs, axis=0, ignore_index=True)
                        stored = True
        except Exception:
            pass
        
        # Check if results were stored
        if 'liana_res' in adata.uns:
            print(f"LIANA analysis completed successfully. Found {len(adata.uns['liana_res'])} interactions.")
        else:
            print("LIANA analysis completed but no results found in adata.uns['liana_res']")
            
    except Exception as e:
        print(f"Error running LIANA analysis: {e}")
        print("Continuing with mock analysis...")
    
    return adata

def analyze_communication_patterns(adata, save_path):
    """Analyze and visualize communication patterns"""
    if not LIANA_AVAILABLE or 'liana_res' not in adata.uns:
        print("No LIANA results available. Creating mock analysis...")
        create_mock_communication_analysis(adata, save_path)
        return
    
    print("Analyzing communication patterns...")
    
    # Get LIANA results
    liana_res = adata.uns['liana_res']

    # Normalize possible dict/list to DataFrame
    try:
        import pandas as _pd  # type: ignore
        if isinstance(liana_res, dict):
            if 'liana_res' in liana_res and isinstance(liana_res['liana_res'], _pd.DataFrame):
                liana_res = liana_res['liana_res']
            else:
                dfs = [v for v in liana_res.values() if isinstance(v, _pd.DataFrame) and not v.empty]
                liana_res = _pd.concat(dfs, axis=0, ignore_index=True) if dfs else _pd.DataFrame()
        elif isinstance(liana_res, list):
            dfs = [x for x in liana_res if isinstance(x, _pd.DataFrame) and not x.empty]
            liana_res = _pd.concat(dfs, axis=0, ignore_index=True) if dfs else _pd.DataFrame()
    except Exception:
        pass

    if liana_res is None or getattr(liana_res, 'empty', False):
        print("LIANA results empty. Falling back to mock analysis...")
        create_mock_communication_analysis(adata, save_path)
        return
    
    # Determine an appropriate ranking/score column
    rank_col = None
    if 'magnitude_rank' in liana_res.columns:
        rank_col = 'magnitude_rank'
    else:
        # Try any column ending with '_rank'
        rank_like = [c for c in liana_res.columns if c.endswith('_rank')]
        if rank_like:
            rank_col = rank_like[0]

    # Filter significant interactions using rank if available; otherwise fallback to magnitude/score
    if rank_col is not None:
        try:
            significant_interactions = liana_res[liana_res[rank_col] < 0.01]
        except Exception:
            significant_interactions = liana_res
    else:
        # Fallbacks: prefer higher magnitude or score
        if 'magnitude' in liana_res.columns:
            thr = liana_res['magnitude'].quantile(0.95)
            significant_interactions = liana_res[liana_res['magnitude'] >= thr]
        elif 'score' in liana_res.columns:
            thr = liana_res['score'].quantile(0.95)
            significant_interactions = liana_res[liana_res['score'] >= thr]
        else:
            significant_interactions = liana_res.head(100)

    if significant_interactions.empty:
        print("No significant interactions found. Falling back to mock analysis...")
        create_mock_communication_analysis(adata, save_path)
        return
    
    # Plot communication network
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create a heatmap of communication strength
    value_col = None
    for cand in ['magnitude_rank', 'magnitude', 'score']:
        if cand in significant_interactions.columns:
            value_col = cand
            break
    if value_col is None:
        value_col = significant_interactions.select_dtypes(include=[np.number]).columns.tolist()
        value_col = value_col[0] if value_col else None

    if value_col is None:
        print("No numeric value column to plot. Skipping heatmap.")
        return

    comm_matrix = significant_interactions.pivot_table(
        index='source', 
        columns='target', 
        values=value_col,
        aggfunc='mean'
    )
    
    sns.heatmap(comm_matrix, annot=True, cmap='viridis', ax=ax)
    plt.title('Cell-Cell Communication Strength')
    plt.tight_layout()
    plt.savefig(save_path / "communication_heatmap_validation_gse145154.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save significant interactions
    significant_interactions.to_csv(save_path / "significant_interactions_validation_gse145154.csv", index=False)

def create_mock_communication_analysis(adata, save_path):
    """Create mock communication analysis when LIANA is not available"""
    print("Creating mock communication analysis...")
    
    # Get cell type counts
    cell_types = adata.obs['leiden'].value_counts()
    
    # Create mock interaction data
    np.random.seed(42)
    n_interactions = 50
    
    mock_interactions = pd.DataFrame({
        'source': np.random.choice(cell_types.index, n_interactions),
        'target': np.random.choice(cell_types.index, n_interactions),
        'ligand': [f'LIGAND_{i}' for i in range(n_interactions)],
        'receptor': [f'RECEPTOR_{i}' for i in range(n_interactions)],
        'score': np.random.uniform(0, 1, n_interactions)
    })
    
    # Create communication heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    
    comm_matrix = mock_interactions.pivot_table(
        index='source',
        columns='target', 
        values='score',
        aggfunc='mean'
    ).fillna(0)
    
    sns.heatmap(comm_matrix, annot=True, cmap='viridis', ax=ax)
    plt.title('Mock Cell-Cell Communication Strength')
    plt.tight_layout()
    plt.savefig(save_path / "mock_communication_heatmap_validation_gse145154.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save mock interactions
    mock_interactions.to_csv(save_path / "mock_interactions_validation_gse145154.csv")

def analyze_ligand_receptor_pairs(adata, save_path):
    """Analyze ligand-receptor pairs"""
    print("Analyzing ligand-receptor pairs...")
    
    # Define some known heart-relevant ligand-receptor pairs
    heart_lr_pairs = {
        'VEGFA': ['FLT1', 'KDR'],
        'PDGFB': ['PDGFRB'],
        'TGFB1': ['TGFBR1', 'TGFBR2'],
        'IGF1': ['IGF1R'],
        'BMP2': ['BMPR1A', 'BMPR2'],
        'WNT3A': ['FZD1', 'FZD2'],
        'NOTCH1': ['DLL1', 'JAG1'],
        'CXCL12': ['CXCR4']
    }
    
    # Check expression of these pairs
    lr_expression = pd.DataFrame()
    
    for ligand, receptors in heart_lr_pairs.items():
        if ligand in adata.var_names:
            ligand_expr = adata[:, ligand].X.toarray().flatten()
            lr_expression[ligand] = ligand_expr
            
        for receptor in receptors:
            if receptor in adata.var_names:
                receptor_expr = adata[:, receptor].X.toarray().flatten()
                lr_expression[receptor] = receptor_expr
    
    # Plot expression patterns
    if not lr_expression.empty:
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create expression heatmap by cell type
        cell_types = adata.obs['leiden'].unique()
        expr_by_celltype = pd.DataFrame()
        
        for ct in cell_types:
            ct_mask = adata.obs['leiden'] == ct
            ct_expr = lr_expression[ct_mask].mean()
            expr_by_celltype[ct] = ct_expr
        
        sns.heatmap(expr_by_celltype, annot=True, cmap='Blues', ax=ax)
        plt.title('Ligand-Receptor Expression by Cell Type')
        plt.tight_layout()
        plt.savefig(save_path / "ligand_receptor_expression_validation_gse145154.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save expression data
        expr_by_celltype.to_csv(save_path / "ligand_receptor_expression_validation_gse145154.csv")

# implement statistical tests for communication analysis
def perform_statistical_test(data, group, strata=None, method="wilcoxon"):
    """
    Perform statistical tests (Wilcoxon or Van Elteren).
    
    Parameters:
        data (array-like): The data values.
        group (array-like): Group labels (e.g., 0 or 1 for two groups).
        strata (array-like, optional): Stratification labels for Van Elteren test.
        method (str): Statistical test method ("wilcoxon" or "van_elteren").
    
    Returns:
        float: Test statistic.
        float: p-value.
    """
    if method == "wilcoxon":
        from scipy.stats import ranksums
        return ranksums(data[group == 1], data[group == 0])
    elif method == "van_elteren":
        if strata is None:
            raise ValueError("Strata must be provided for Van Elteren test.")
        return van_elteren_test(data, group, strata)
    else:
        raise ValueError(f"Unknown method: {method}")

# Example data
data = np.array([1.2, 2.3, 1.8, 2.5, 3.1, 2.9, 1.7, 2.8])
group = np.array([0, 0, 1, 1, 0, 0, 1, 1])  # Two groups: 0 and 1
strata = np.array([1, 1, 1, 1, 2, 2, 2, 2])  # Two strata: 1 and 2

# Perform Wilcoxon rank-sum test
stat, p_value = perform_statistical_test(data, group, method="wilcoxon")
print(f"Wilcoxon Test: Statistic={stat}, p-value={p_value}")

# Perform Van Elteren test
stat, p_value = perform_statistical_test(data, group, strata=strata, method="van_elteren")
print(f"Van Elteren Test: Statistic={stat}, p-value={p_value}")

def main():
    # Load annotated data
    project_root = Path(__file__).parent.parent.parent
    data_path = project_root / "data/processed/heart_data_annotated_validation_gse145154.h5ad"
    adata = sc.read_h5ad(data_path)
    
    # Create results directories
    results_path = Path("results/communication")
    results_path.mkdir(parents=True, exist_ok=True)
    
    # Prepare data for communication analysis
    comm_adata = prepare_data_for_communication(adata, cell_type_col='leiden', sample_key='sample')
    
    if comm_adata is None:
        print("Could not prepare data for communication analysis")
        return
    
    # Run LIANA analysis
    print(f"LIANA_AVAILABLE: {LIANA_AVAILABLE}")
    if LIANA_AVAILABLE:
        comm_adata = run_liana_analysis(comm_adata, cell_type_col='leiden', sample_key='sample')
    else:
        print("Skipping LIANA analysis due to import failure")
    
    # Analyze communication patterns
    analyze_communication_patterns(comm_adata, results_path)
    
    # Analyze ligand-receptor pairs
    analyze_ligand_receptor_pairs(comm_adata, results_path)
    
    # Save communication analysis results
    comm_adata.write(project_root / "data/processed/heart_data_communication_validation_gse145154.h5ad")
    
    print("Cell-cell communication analysis complete!")

if __name__ == "__main__":
    main()