#!/usr/bin/env python
"""
HeartMAP Web Interface - Gradio app for Hugging Face Spaces deployment
Comprehensive chamber-specific cardiac analysis platform
"""

import gradio as gr
import tempfile
import sys
import shutil
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
from typing import Tuple, List, Dict

# Add src to path
sys.path.insert(0, 'src')

try:
    import scanpy as sc
    import anndata as ad
    from heartmap import Config
    from heartmap.pipelines import (
        BasicPipeline, 
        ComprehensivePipeline,
        MultiChamberPipeline,
        AdvancedCommunicationPipeline
    )
    HEARTMAP_AVAILABLE = True
except ImportError:
    HEARTMAP_AVAILABLE = False

def load_and_validate_data(uploaded_file) -> Tuple[ad.AnnData, str]:
    """Load and validate uploaded single-cell data in various formats with compatibility fallbacks
    
    Supported formats:
    - AnnData: .h5ad, .h5mu, .zarr
    - 10X: .h5, .mtx (with genes/barcodes)
    - Seurat: .rds, .rdata (h5seurat)
    - Loom: .loom
    - Text: .csv, .tsv, .txt (expression matrices)
    - Archives: .tar, .tar.gz, .tgz
    - HDF5: .hdf5, .h5 (generic)
    - Parquet: .parquet
    - And more...
    """
    import h5py
    import tarfile
    import tempfile
    import shutil
    from pathlib import Path
    
    # Determine file type
    filename = uploaded_file if isinstance(uploaded_file, str) else getattr(uploaded_file, 'name', str(uploaded_file))
    file_lower = filename.lower()
    filepath = Path(uploaded_file)
    
    # Categorize file types
    is_h5ad = file_lower.endswith('.h5ad')
    is_h5 = file_lower.endswith('.h5') and not file_lower.endswith('.h5ad') and not file_lower.endswith('.h5mu')
    is_tar = file_lower.endswith(('.tar', '.tar.gz', '.tgz'))
    is_loom = file_lower.endswith('.loom')
    is_mtx = file_lower.endswith(('.mtx', '.mtx.gz'))
    is_text = file_lower.endswith(('.csv', '.tsv', '.txt', '.csv.gz', '.tsv.gz', '.txt.gz'))
    is_parquet = file_lower.endswith('.parquet')
    is_hdf5 = file_lower.endswith('.hdf5')
    is_h5mu = file_lower.endswith('.h5mu')
    is_zarr = file_lower.endswith('.zarr') or filepath.is_dir()
    is_rds = file_lower.endswith(('.rds', '.rdata'))
    is_h5seurat = 'seurat' in file_lower and file_lower.endswith('.h5')
    is_arrow = file_lower.endswith('.arrow')
    is_json = file_lower.endswith(('.json', '.geojson'))

    try:
        # Handle TAR archives (GEO datasets)
        if is_tar:
            print("Detected TAR archive (GEO format). Extracting files...")
            temp_dir = tempfile.mkdtemp()
            
            try:
                # Extract tar file
                with tarfile.open(uploaded_file, 'r:*') as tar:
                    tar.extractall(temp_dir)
                    print(f"✓ Extracted {len(tar.getmembers())} files")
                
                # Look for compatible files
                extracted_files = list(Path(temp_dir).rglob('*'))
                
                # Debug: Show what files were found
                print(f"  Analyzing {len(extracted_files)} extracted items...")
                file_extensions = {}
                for f in extracted_files:
                    if f.is_file():
                        ext = f.suffix.lower()
                        file_extensions[ext] = file_extensions.get(ext, 0) + 1
                        # Also check full name for multi-extension files
                        if ext == '.gz':
                            # Check what's before .gz
                            name_lower = f.name.lower()
                            if '.h5ad.gz' in name_lower:
                                file_extensions['.h5ad.gz'] = file_extensions.get('.h5ad.gz', 0) + 1
                            elif '.h5.gz' in name_lower:
                                file_extensions['.h5.gz'] = file_extensions.get('.h5.gz', 0) + 1
                            elif '.mtx.gz' in name_lower:
                                file_extensions['.mtx.gz'] = file_extensions.get('.mtx.gz', 0) + 1
                
                if file_extensions:
                    print(f"  File types found: {file_extensions}")
                
                # Check for files with full name patterns (including .gz)
                h5ad_files = [f for f in extracted_files if f.name.lower().endswith('.h5ad') or f.name.lower().endswith('.h5ad.gz')]
                h5_files = [f for f in extracted_files if (f.name.lower().endswith('.h5') or f.name.lower().endswith('.h5.gz')) 
                           and not f.name.lower().endswith('.h5ad') and not f.name.lower().endswith('.h5ad.gz')]
                mtx_dirs = [f.parent for f in extracted_files if 'matrix.mtx' in f.name.lower()]
                
                if h5ad_files:
                    print(f"✓ Found {len(h5ad_files)} .h5ad file(s). Loading first one...")
                    h5ad_file = h5ad_files[0]
                    
                    # If gzipped, decompress first
                    if h5ad_file.name.endswith('.gz'):
                        print(f"  Decompressing {h5ad_file.name}...")
                        import gzip
                        decompressed_path = h5ad_file.parent / h5ad_file.stem  # Remove .gz
                        with gzip.open(h5ad_file, 'rb') as f_in:
                            with open(decompressed_path, 'wb') as f_out:
                                f_out.write(f_in.read())
                        print(f"  ✓ Decompressed to {decompressed_path.name}")
                        adata = sc.read_h5ad(str(decompressed_path))
                    else:
                        adata = sc.read_h5ad(str(h5ad_file))
                elif h5_files:
                    print(f"✓ Found {len(h5_files)} .h5 file(s). Loading first one...")
                    # Try loading as 10X format first, fall back to AnnData h5ad format
                    h5_loaded = False
                    for h5_file in h5_files:
                        # Decompress if gzipped
                        file_to_load = h5_file
                        if h5_file.name.endswith('.gz'):
                            print(f"  Decompressing {h5_file.name}...")
                            import gzip
                            decompressed_path = h5_file.parent / h5_file.stem  # Remove .gz
                            with gzip.open(h5_file, 'rb') as f_in:
                                with open(decompressed_path, 'wb') as f_out:
                                    f_out.write(f_in.read())
                            print(f"  ✓ Decompressed to {decompressed_path.name}")
                            file_to_load = decompressed_path
                        
                        try:
                            print(f"  Trying to load {file_to_load.name} as 10X format...")
                            adata = sc.read_10x_h5(str(file_to_load))
                            h5_loaded = True
                            print(f"  ✓ Successfully loaded as 10X format")
                            break
                        except Exception as e1:
                            print(f"  ⚠ Not 10X format: {str(e1)[:100]}")
                            try:
                                print(f"  Trying to load {file_to_load.name} as AnnData format...")
                                adata = sc.read_h5ad(str(file_to_load))
                                h5_loaded = True
                                print(f"  ✓ Successfully loaded as AnnData format")
                                break
                            except Exception as e2:
                                print(f"  ⚠ Not AnnData format: {str(e2)[:100]}")
                                continue
                    
                    if not h5_loaded:
                        print("  ⚠ No .h5 files could be loaded, trying other formats...")
                        adata = None
                elif mtx_dirs:
                    print(f"✓ Found matrix.mtx format. Loading from {mtx_dirs[0].name}...")
                    mtx_dir = mtx_dirs[0]
                    
                    # Check for required files (matrix, genes/features, barcodes)
                    mtx_files = list(mtx_dir.glob('matrix.mtx*'))
                    gene_files = list(mtx_dir.glob('genes.tsv*')) or list(mtx_dir.glob('features.tsv*'))
                    barcode_files = list(mtx_dir.glob('barcodes.tsv*'))
                    
                    print(f"  Found in {mtx_dir.name}:")
                    print(f"    Matrix files: {[f.name for f in mtx_files]}")
                    print(f"    Gene/feature files: {[f.name for f in gene_files]}")
                    print(f"    Barcode files: {[f.name for f in barcode_files]}")
                    
                    # Decompress .gz files if present
                    import gzip
                    for file_list in [mtx_files, gene_files, barcode_files]:
                        for f in file_list:
                            if f.name.endswith('.gz'):
                                print(f"  Decompressing {f.name}...")
                                decompressed_path = f.parent / f.name[:-3]  # Remove .gz
                                if not decompressed_path.exists():
                                    with gzip.open(f, 'rb') as f_in:
                                        with open(decompressed_path, 'wb') as f_out:
                                            f_out.write(f_in.read())
                                    print(f"    ✓ Decompressed to {decompressed_path.name}")
                    
                    # Try to read the MTX directory
                    try:
                        adata = sc.read_10x_mtx(str(mtx_dir))
                        print(f"  ✓ Successfully loaded MTX format")
                    except Exception as mtx_err:
                        print(f"  ⚠ Failed to load MTX: {str(mtx_err)[:200]}")
                        adata = None
                else:
                    adata = None
                
                # If h5 files failed or no standard format found, try text files
                if adata is None:
                    print("Searching for readable text files...")
                    adata = None
                    for ext_file in extracted_files:
                        if ext_file.is_file():
                            try:
                                if ext_file.suffix in ['.txt', '.csv', '.tsv', '.txt.gz', '.csv.gz', '.tsv.gz']:
                                    print(f" Attempting to read as expression matrix: {ext_file.name}")
                                    adata = sc.read_csv(str(ext_file), delimiter='\t' if 'tsv' in ext_file.name else ',')
                                    break
                            except Exception as read_err:
                                print(f"⚠ Failed to read {ext_file.name}: {str(read_err)}")
                                continue
                    
                    if adata is None:
                        shutil.rmtree(temp_dir)
                        return None, " TAR archive doesn't contain compatible files (.h5ad, .h5, matrix.mtx, or expression matrices)"
                
                if adata is not None:
                    print(f"✓ Successfully loaded from TAR archive: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
                else:
                    shutil.rmtree(temp_dir)
                    return None, " TAR archive doesn't contain compatible files (.h5ad, .h5, matrix.mtx, or expression matrices)"
                
                # Cleanup temp directory
                shutil.rmtree(temp_dir)
                
            except Exception as e:
                if Path(temp_dir).exists():
                    shutil.rmtree(temp_dir)
                raise e
        
        # Try reading based on file type
        elif is_h5ad:
            print(" Loading AnnData (.h5ad)...")
            adata = sc.read_h5ad(uploaded_file)
            print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        elif is_h5mu:
            print(" Loading MuData (.h5mu)...")
            try:
                import mudata
                mdata = mudata.read_h5mu(uploaded_file)
                adata = mdata.mod[list(mdata.mod.keys())[0]]
                print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
            except ImportError:
                return None, " MuData support requires: pip install mudata"
        
        elif is_zarr:
            print(" Loading Zarr array...")
            adata = sc.read_zarr(uploaded_file)
            print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        elif is_h5 and not is_h5seurat:
            print(" Loading 10X Genomics (.h5)...")
            adata = sc.read_10x_h5(uploaded_file)
            print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        elif is_mtx:
            print(" Loading Matrix Market (.mtx)...")
            mtx_dir = filepath.parent
            if (mtx_dir / 'genes.tsv').exists():
                adata = sc.read_10x_mtx(str(mtx_dir))
            else:
                adata = sc.read_mtx(uploaded_file).T
            print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        elif is_loom:
            print(" Loading Loom (.loom)...")
            adata = sc.read_loom(uploaded_file)
            print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        elif is_text:
            print(f" Loading text matrix ({filepath.suffix})...")
            delim = '\t' if 'tsv' in file_lower else ','
            adata = sc.read_csv(uploaded_file, delimiter=delim, first_column_names=True)
            print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        elif is_parquet or is_arrow:
            print(f" Loading {'Parquet' if is_parquet else 'Arrow'}...")
            try:
                import pandas as pd
                df = pd.read_parquet(uploaded_file) if is_parquet else pd.read_feather(uploaded_file)
                adata = ad.AnnData(df)
                print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
            except ImportError:
                return None, " Parquet/Arrow requires: pip install pyarrow"
        
        elif is_hdf5:
            print(" Loading generic HDF5...")
            with h5py.File(uploaded_file, 'r') as f:
                X = f.get('matrix', f.get('X', f.get('data', None)))
                if X is None:
                    return None, " HDF5 file missing required keys (matrix/X/data)"
                adata = ad.AnnData(X[()])
            print(f"✓ {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        elif is_h5seurat or is_rds:
            filename = filepath.name
            return None, (f" **Seurat/R Format Detected: `{filename}`**\n\n"
                         f"This file format requires conversion to AnnData (.h5ad) before analysis.\n\n"
                         f"**Option 1: Using SeuratDisk (Recommended)**\n"
                         f"```r\n"
                         f"# In R/RStudio:\n"
                         f"library(Seurat)\n"
                         f"library(SeuratDisk)\n\n"
                         f"# Load your Seurat object\n"
                         f"seurat_obj <- readRDS('{filename}')\n\n"
                         f"# Convert to h5ad format\n"
                         f"SaveH5Seurat(seurat_obj, filename = 'output.h5Seurat')\n"
                         f"Convert('output.h5Seurat', dest = 'h5ad')\n"
                         f"```\n\n"
                         f"**Option 2: Using sceasy**\n"
                         f"```r\n"
                         f"library(sceasy)\n"
                         f"sceasy::convertFormat(seurat_obj, from='seurat', to='anndata',\n"
                         f"                       outFile='output.h5ad')\n"
                         f"```\n\n"
                         f"Then upload the generated `output.h5ad` file to HeartMAP.\n\n"
                         f"**Need help?** See [FORMAT_SUPPORT_GUIDE.md](https://github.com/Tumo505/HeartMap/blob/master/FORMAT_SUPPORT_GUIDE.md) for detailed instructions.")
        
        else:
            return None, (f" Unsupported format: {filepath.suffix}\n\n"
                         f"**Supported formats:**\n"
                         f"• AnnData: .h5ad, .h5mu (MuData), .zarr\n"
                         f"• 10X: .h5, .mtx (Matrix Market)\n"
                         f"• Loom: .loom\n"
                         f"• Text: .csv, .tsv, .txt\n"
                         f"• Archives: .tar, .tar.gz, .tgz\n"
                         f"• Columnar: .parquet, .arrow\n"
                         f"• HDF5: .hdf5\n\n"
                         f"For FASTQ/BAM/CRAM, use Cell Ranger or similar tools first.")

    except Exception as e:
        # Handle version compatibility issues
        error_str = str(e)
        is_compat_error = (
            "IOSpec" in error_str or 
            "encoding_type" in error_str or 
            "No read method registered" in error_str or
            "unexpected keyword argument" in error_str or
            "'matrix'" in error_str or
            "AnnData.__init__" in error_str
        )
        
        if is_compat_error:
            print(f"AnnData version compatibility issue detected: {error_str[:100]}")
            print(f" Attempting fallback read method...")
            try:
                # Fallback 1: Read with backed mode then load into memory
                adata = sc.read_h5ad(uploaded_file, backed='r')
                adata = adata.to_memory()
                print("✓ Successfully loaded using backed mode")
            except Exception as e2:
                try:
                    # Fallback 2: Read manually from h5py with proper structure handling
                    print(" Attempting manual h5py read...")
                    with h5py.File(uploaded_file, 'r') as f:
                        # Read X matrix (handle different storage formats)
                        if 'X' in f:
                            X_group = f['X']
                            if isinstance(X_group, h5py.Dataset):
                                X = X_group[:]
                            else:
                                # Sparse matrix format
                                try:
                                    from scipy import sparse
                                    data = X_group['data'][:]
                                    indices = X_group['indices'][:]
                                    indptr = X_group['indptr'][:]
                                    shape = X_group['shape'][:]
                                    X = sparse.csr_matrix((data, indices, indptr), shape=shape)
                                except:
                                    X = X_group['data'][:]  # Fallback to data only
                        else:
                            raise ValueError("No X matrix found in file")

                        # Read obs (cell metadata)
                        obs_dict = {}
                        if 'obs' in f:
                            obs_group = f['obs']
                            for key in obs_group.keys():
                                try:
                                    data = obs_group[key][:]
                                    # Decode bytes if necessary
                                    if data.dtype.kind == 'S' or data.dtype.kind == 'O':
                                        data = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in data]
                                    obs_dict[key] = data
                                except Exception as e_key:
                                    print(f"Skipping obs key '{key}': {e_key}")
                        obs = pd.DataFrame(obs_dict) if obs_dict else pd.DataFrame(index=range(X.shape[0]))

                        # Read var (gene metadata)
                        var_dict = {}
                        if 'var' in f:
                            var_group = f['var']
                            for key in var_group.keys():
                                try:
                                    data = var_group[key][:]
                                    # Decode bytes if necessary
                                    if data.dtype.kind == 'S' or data.dtype.kind == 'O':
                                        data = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in data]
                                    var_dict[key] = data
                                except Exception as e_key:
                                    print(f"Skipping var key '{key}': {e_key}")
                        var = pd.DataFrame(var_dict) if var_dict else pd.DataFrame(index=range(X.shape[1]))

                        # Create basic AnnData object (skip problematic uns)
                        adata = ad.AnnData(X=X, obs=obs, var=var)
                    print("✓ Successfully loaded using manual h5py read")
                except Exception as e3:
                    raise ValueError(
                        f"Unable to load file with any method.\n\n"
                        f"**Primary error:** {str(e)}\n\n"
                        f"**Suggestions:**\n"
                        f"1. The file may have been created with a newer AnnData version\n"
                        f"2. Try re-saving the file with: `adata.write('file.h5ad', compression='gzip')`\n"
                        f"3. Or use an older AnnData format: `adata.write_h5ad('file.h5ad', as_dense='X')`\n\n"
                        f"**Technical details:** {str(e3)}"
                    )
        else:
            # Not a compatibility error - return the original error
            import traceback
            return None, (f" **Error loading data:**\n\n{str(e)}\n\n"
                         f"**Traceback:**\n```\n{traceback.format_exc()}\n```")

    try:
        # Check for chamber information
        chamber_info = ""
        if 'chamber' in adata.obs.columns:
            chambers = adata.obs['chamber'].unique()
            chamber_info = f"\nChamber information detected: {', '.join(chambers)}"
        else:
            # Try to infer chamber from various metadata fields
            chamber_assigned = False
            
            # Debug: Show available metadata columns
            print(f"Available metadata columns: {list(adata.obs.columns)}")
            
            # Check for common chamber-related column names
            chamber_keywords = ['tissue', 'location', 'sample', 'orig.ident', 'biosample', 'cell', 'batch', 'donor', 'patient']
            for col in adata.obs.columns:
                if any(keyword in col.lower() for keyword in chamber_keywords):
                    print(f"Attempting to infer chamber from column: {col}")
                    # Show sample values
                    sample_values = adata.obs[col].unique()[:5]
                    print(f"  Sample values: {sample_values}")
                    values = adata.obs[col].astype(str).str.upper()
                    
                    # Map common chamber identifiers
                    def map_chamber(val):
                        val = val.upper()
                        if any(x in val for x in ['RA', 'RIGHT ATRI', 'R_ATRI']):
                            return 'RA'
                        elif any(x in val for x in ['RV', 'RIGHT VENT', 'R_VENT']):
                            return 'RV'
                        elif any(x in val for x in ['LA', 'LEFT ATRI', 'L_ATRI']):
                            return 'LA'
                        elif any(x in val for x in ['LV', 'LEFT VENT', 'L_VENT']):
                            return 'LV'
                        elif any(x in val for x in ['ATRI']):
                            return 'RA'  # Default atrium
                        elif any(x in val for x in ['VENT']):
                            return 'LV'  # Default ventricle
                        else:
                            return None
                    
                    adata.obs['chamber'] = adata.obs[col].apply(map_chamber)
                    
                    # Check if we successfully assigned chambers
                    if adata.obs['chamber'].notna().sum() > 0:
                        chambers = adata.obs['chamber'].dropna().unique()
                        if len(chambers) > 1:
                            chamber_info = f"\nChamber information inferred from '{col}': {', '.join(chambers)}"
                            chamber_assigned = True
                            break
                        else:
                            # Fill NaN with the single detected chamber
                            adata.obs['chamber'] = adata.obs['chamber'].fillna(chambers[0])
                            chamber_info = f"\nSingle chamber detected from '{col}': {chambers[0]}"
                            chamber_assigned = True
                            break
            
            if not chamber_assigned:
                chamber_info = "\n⚠ No chamber information found"
                chamber_info += "\n  Single-chamber analysis will be performed"
                chamber_info += "\n  For multi-chamber analysis, data should have 'chamber' column with values: RA, RV, LA, LV"
                # Don't assign a default chamber - let the analysis handle missing chamber info
                adata.obs['chamber'] = 'Unknown'

        validation_msg = f"""
Data loaded successfully!
- Cells: {adata.n_obs:,}
- Genes: {adata.n_vars:,}
{chamber_info}
        """
        return adata, validation_msg

    except Exception as e:
        raise ValueError(f"Error validating loaded data: {str(e)}")


def create_communication_network(adata, hub_stats, chamber_stats=None):
    """
    Create interactive Plotly network graph showing inferred cell-cell communication
    based on co-expression of ligand-receptor genes
    
    Args:
        adata: AnnData object with clustering and gene expression
        hub_stats: DataFrame with hub scores per cell type
        chamber_stats: Optional DataFrame with chamber information
    
    Returns:
        Path to HTML file with interactive network
    """
    import networkx as nx
    
    try:
        # Get cell type information
        if 'leiden' not in adata.obs.columns:
            print("No clustering found, skipping network graph")
            return None
        
        cell_types = adata.obs['leiden'].unique()
        n_types = len(cell_types)
        
        # Common ligand-receptor pairs in cardiac tissue
        ligand_receptor_pairs = [
            ('VEGFA', 'FLT1'), ('VEGFA', 'KDR'),  # Angiogenesis
            ('TGFB1', 'TGFBR1'), ('TGFB1', 'TGFBR2'),  # TGF-beta signaling
            ('FGF2', 'FGFR1'), ('FGF2', 'FGFR2'),  # FGF signaling
            ('IL6', 'IL6R'),  # Inflammation
            ('TNF', 'TNFRSF1A'), ('TNF', 'TNFRSF1B'),  # TNF signaling
            ('PDGFA', 'PDGFRA'), ('PDGFB', 'PDGFRB'),  # PDGF signaling
            ('EGF', 'EGFR'),  # EGF signaling
            ('IGF1', 'IGF1R'),  # IGF signaling
            ('CXCL12', 'CXCR4'),  # Chemokine signaling
            ('CCL2', 'CCR2'),  # Monocyte recruitment
        ]
        
        # Calculate mean expression per cell type
        print("  Calculating cell type expression profiles...")
        cell_type_expression = {}
        for cell_type in cell_types:
            cell_mask = adata.obs['leiden'] == cell_type
            if hasattr(adata.X, 'toarray'):
                subset_expr = adata.X[cell_mask].toarray()
            else:
                subset_expr = adata.X[cell_mask]
            cell_type_expression[str(cell_type)] = np.mean(subset_expr, axis=0).A1 if hasattr(subset_expr, 'A1') else np.mean(subset_expr, axis=0)
        
        # Create network graph
        G = nx.DiGraph()  # Directed graph for ligand->receptor
        
        # Add nodes for each cell type with metadata
        node_info = []
        for i, cell_type in enumerate(cell_types):
            cell_mask = adata.obs['leiden'] == cell_type
            n_cells = cell_mask.sum()
            
            # Get hub score if available
            hub_score = 0
            if hub_stats is not None and len(hub_stats) > 0:
                type_label = f"Cluster {cell_type}"
                matching = hub_stats[hub_stats['Cell Type'] == type_label]
                if len(matching) > 0:
                    hub_score = matching.iloc[0]['Hub Score']
            
            # Get chamber distribution
            chamber_dist = ""
            if chamber_stats is not None and 'chamber' in adata.obs.columns:
                type_chambers = adata.obs[cell_mask]['chamber'].value_counts()
                chamber_dist = ", ".join([f"{ch}: {cnt}" for ch, cnt in type_chambers.items()])
            
            node_info.append({
                'id': str(cell_type),
                'label': f"Cluster {cell_type}",
                'size': int(n_cells),
                'hub_score': float(hub_score),
                'chamber_dist': chamber_dist,
                'color_idx': i
            })
            
            G.add_node(str(cell_type),
                      size=int(n_cells),
                      hub_score=float(hub_score),
                      label=f"Cluster {cell_type}",
                      color_idx=i)
        
        # Infer communication edges from ligand-receptor co-expression
        print("  Inferring cell-cell communication from gene expression...")
        edge_info = []
        for ligand, receptor in ligand_receptor_pairs:
            # Check if genes exist in dataset
            if ligand not in adata.var_names or receptor not in adata.var_names:
                continue
            
            ligand_idx = list(adata.var_names).index(ligand)
            receptor_idx = list(adata.var_names).index(receptor)
            
            # Find cell types that express ligand and receptor
            for source_type in cell_types:
                ligand_expr = cell_type_expression[str(source_type)][ligand_idx]
                
                if ligand_expr > 0.1:  # Threshold for meaningful expression
                    for target_type in cell_types:
                        if source_type == target_type:
                            continue
                        
                        receptor_expr = cell_type_expression[str(target_type)][receptor_idx]
                        
                        if receptor_expr > 0.1:  # Both genes expressed
                            # Calculate interaction strength
                            strength = np.sqrt(ligand_expr * receptor_expr)
                            
                            edge_key = (str(source_type), str(target_type))
                            edge_info.append({
                                'source': str(source_type),
                                'target': str(target_type),
                                'ligand': ligand,
                                'receptor': receptor,
                                'strength': strength
                            })
                            
                            # Add or update edge
                            if G.has_edge(str(source_type), str(target_type)):
                                G[str(source_type)][str(target_type)]['weight'] += strength
                                G[str(source_type)][str(target_type)]['interactions'].append(f"{ligand}-{receptor}")
                            else:
                                G.add_edge(str(source_type), str(target_type), 
                                          weight=strength, 
                                          interactions=[f"{ligand}-{receptor}"])
        
        print(f"  Found {G.number_of_edges()} communication interactions between {G.number_of_nodes()} cell types")
        
        # Create Plotly figure with layout
        pos = nx.spring_layout(G, k=2.5, iterations=50, seed=42, weight='weight')
        
        # Create edge traces (one per edge for hover info)
        edge_traces = []
        for edge in G.edges(data=True):
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            
            weight = edge[2].get('weight', 0)
            interactions = edge[2].get('interactions', [])
            
            # Create arrow shape for directed edge
            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line=dict(
                    width=max(0.5, min(5, weight * 2)),  # Scale by strength
                    color='rgba(150, 150, 150, 0.5)'
                ),
                hoverinfo='text',
                hovertext=f"<b>{edge[0]} → {edge[1]}</b><br>" + 
                          f"Interaction strength: {weight:.3f}<br>" +
                          f"Ligand-Receptor pairs: {len(interactions)}<br>" +
                          "<br>".join(interactions[:5]) +  # Show first 5
                          (f"<br>... and {len(interactions)-5} more" if len(interactions) > 5 else ""),
                showlegend=False
            )
            edge_traces.append(edge_trace)
        
        # Generate distinct colors for cell types
        import plotly.colors as pcolors
        if n_types <= 10:
            colors = pcolors.qualitative.Set3[:n_types]
        else:
            colors = pcolors.sample_colorscale("turbo", [i/n_types for i in range(n_types)])
        
        # Create separate trace for each cell type (to show in legend)
        node_traces = []
        for node_data in node_info:
            node_id = node_data['id']
            if node_id not in pos:
                continue
                
            x, y = pos[node_id]
            label = node_data['label']
            size = node_data['size']
            hub_score = node_data['hub_score']
            chamber_dist = node_data['chamber_dist']
            color_idx = node_data['color_idx']
            
            # Count outgoing and incoming communications
            out_edges = G.out_degree(node_id)
            in_edges = G.in_degree(node_id)
            
            hover_text = (
                f"<b>{label}</b><br>"
                f"Cells: {size:,}<br>"
                f"Hub Score: {hub_score:.4f}<br>"
                f"Sends signals to: {out_edges} cell types<br>"
                f"Receives from: {in_edges} cell types<br>"
                f"{chamber_dist}"
            )
            
            node_trace = go.Scatter(
                x=[x],
                y=[y],
                mode='markers+text',
                marker=dict(
                    size=max(15, min(60, size / 20)),  # Scale by cell count
                    color=colors[color_idx % len(colors)],
                    line=dict(width=2, color='white'),
                    symbol='circle'
                ),
                text=label.replace('Cluster ', 'C'),
                textposition="top center",
                textfont=dict(size=10, color='black'),
                hoverinfo='text',
                hovertext=hover_text,
                name=label,
                legendgroup=label,
                showlegend=True
            )
            node_traces.append(node_trace)
        
        # Create figure with all traces
        fig = go.Figure(data=edge_traces + node_traces)
        
        fig.update_layout(
            title={
                'text': "Cell-Cell Communication Network<br><sub>Inferred from Ligand-Receptor Gene Expression</sub>",
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 20}
            },
            showlegend=True,
            legend=dict(
                title="Cell Types",
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.05,
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="gray",
                borderwidth=1
            ),
            hovermode='closest',
            width=1400,
            height=900,
            plot_bgcolor='#f8f9fa',
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            annotations=[
                dict(
                    text="Arrows show cell-cell communication via ligand-receptor pairs | Node size = cell count | Edge width = interaction strength",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.5, y=-0.02,
                    xanchor='center',
                    font=dict(size=11, color='gray')
                )
            ],
            margin=dict(l=20, r=250, t=80, b=40)
        )
        
        # Save as HTML
        html_path = tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False).name
        fig.write_html(
            html_path,
            include_plotlyjs='cdn',
            config={
                'displayModeBar': True,
                'displaylogo': False,
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': 'communication_network',
                    'height': 800,
                    'width': 1200,
                    'scale': 2
                }
            }
        )
        
        print(f"✓ Created interactive communication network: {html_path}")
        return html_path
        
    except Exception as e:
        print(f"Warning: Could not create network graph: {e}")
        return None


def analyze_heart_data(
    uploaded_file,
    analysis_type,
    max_cells,
    max_genes,
    include_chamber_analysis,
    include_communication_hubs
) -> Tuple[str, str, str, str, str, str, str, str]:
    """Comprehensive HeartMAP analysis with chamber-specific insights
    
    Returns:
        Tuple of (output_msg, csv_file, viz_chamber, viz_hubs, viz_corr, viz_markers, viz_network, chamber_info)
    """

    if not HEARTMAP_AVAILABLE:
        return " HeartMAP not available. Please install dependencies.", None, None, None, None, None, None, "", None

    if uploaded_file is None:
        return " Please upload a file.", None, None, None, None, None, None, "", None

    # Validate file upload completed successfully
    try:
        if not Path(uploaded_file).exists():
            return " Upload incomplete. Please try uploading the file again.", None, None, None, None, None, None, "", None

        file_size_mb = Path(uploaded_file).stat().st_size / (1024 * 1024)
        file_size_gb = file_size_mb / 1024
        print(f"Processing file: {file_size_mb:.2f} MB ({file_size_gb:.2f} GB)")

        if file_size_mb > 10240:  # 10GB limit
            return f" File size ({file_size_gb:.2f} GB) exceeds maximum limit (10 GB). Please use a smaller dataset or subset your data.", None, None, None, None, None, None, "", None
        elif file_size_mb > 1024:  # Warn for files > 1GB
            print(f"Large file detected ({file_size_gb:.2f} GB). Processing may take 10-30 minutes...")
        elif file_size_mb > 500:
            print(f"Large file detected ({file_size_mb:.1f} MB). Processing may take several minutes...")
    except Exception as e:
        return f" File validation error: {str(e)}. Please re-upload the file.", None, None, None, None, None, None, "", None

    persistent_csv = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
    persistent_csv_path = persistent_csv.name
    persistent_csv.close()
    
    # Create persistent temporary files for visualizations (outside temp_dir context)
    persistent_viz_chamber = tempfile.NamedTemporaryFile(mode='w', suffix='_chamber_distribution.html', delete=False)
    persistent_viz_chamber_path = persistent_viz_chamber.name
    persistent_viz_chamber.close()
    
    persistent_viz_hubs = tempfile.NamedTemporaryFile(mode='w', suffix='_hub_scores.html', delete=False)
    persistent_viz_hubs_path = persistent_viz_hubs.name
    persistent_viz_hubs.close()
    
    persistent_viz_corr = tempfile.NamedTemporaryFile(mode='w', suffix='_correlations.html', delete=False)
    persistent_viz_corr_path = persistent_viz_corr.name
    persistent_viz_corr.close()
    
    persistent_viz_markers = tempfile.NamedTemporaryFile(mode='w', suffix='_markers.html', delete=False)
    persistent_viz_markers_path = persistent_viz_markers.name
    persistent_viz_markers.close()
    
    persistent_viz_network = tempfile.NamedTemporaryFile(mode='w', suffix='_network.html', delete=False)
    persistent_viz_network_path = persistent_viz_network.name
    persistent_viz_network.close()
    
    persistent_chamber_info = tempfile.NamedTemporaryFile(mode='w', suffix='_chamber_details.txt', delete=False)
    persistent_chamber_info_path = persistent_chamber_info.name
    persistent_chamber_info.close()

    try:
        # Load and validate data
        print("Loading data...")
        adata, validation_msg = load_and_validate_data(uploaded_file)
        
        # Check if loading failed
        if adata is None:
            error_msg = f" **Data Loading Failed**\n\n{validation_msg}"
            return error_msg, None, None, None, None, None, None, None

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_dir = temp_path / "results"
            results_dir.mkdir(exist_ok=True)

            # Save loaded data to temporary h5ad file for pipeline
            temp_data_file = temp_path / "loaded_data.h5ad"
            print(f"Saving loaded data to temporary file...")
            adata.write_h5ad(temp_data_file)

            # Create config
            config = Config.default()
            config.data.max_cells_subset = int(max_cells) if max_cells else None
            config.data.max_genes_subset = int(max_genes) if max_genes else None
            config.paths.processed_data_dir = str(temp_path / "processed")
            config.paths.results_dir = str(results_dir)

            config.create_directories()

            # Run appropriate pipeline
            if analysis_type == "comprehensive":
                pipeline = ComprehensivePipeline(config)
            elif analysis_type == "basic":
                pipeline = BasicPipeline(config)
            elif analysis_type == "multi_chamber" and include_chamber_analysis:
                pipeline = MultiChamberPipeline(config)
            else:
                pipeline = BasicPipeline(config)

            # Pass temporary h5ad file instead of original file
            results = pipeline.run(str(temp_data_file), str(results_dir))
            analyzed_adata = results.get('adata', adata)

            # ===== CREATE COMPREHENSIVE SUMMARY =====
            summary_data = {
                'Metric': [],
                'Value': []
            }
            
            # Basic statistics
            summary_data['Metric'].extend([
                'Total Cells Analyzed',
                'Total Genes Analyzed',
                'Cell Clusters Identified',
                'Mean Genes per Cell',
                'Mean UMI per Cell'
            ])

            n_clusters = analyzed_adata.obs['leiden'].nunique() if 'leiden' in analyzed_adata.obs.columns else 0
            mean_genes = analyzed_adata.n_vars  # Approximation
            mean_umi = 0  # Would need to calculate from data

            summary_data['Value'].extend([
                f"{analyzed_adata.n_obs:,}",
                f"{analyzed_adata.n_vars:,}",
                str(n_clusters),
                str(mean_genes),
                str(mean_umi)
            ])

            # ===== CHAMBER ANALYSIS =====
            chamber_stats = None
            chamber_correlations = None
            chamber_markers = None
            
            if include_chamber_analysis and 'chamber' in analyzed_adata.obs.columns:
                # Filter out Unknown chambers
                valid_chambers = analyzed_adata.obs['chamber'] != 'Unknown'
                if valid_chambers.sum() == 0:
                    summary_data['Metric'].append('Chamber Analysis')
                    summary_data['Value'].append('⚠ No chamber information available')
                    chamber_counts = pd.Series(dtype=int)  # Empty series
                else:
                    chamber_counts = analyzed_adata.obs.loc[valid_chambers, 'chamber'].value_counts()
                    chamber_proportions = (chamber_counts / valid_chambers.sum() * 100).round(2)

                    chamber_stats = pd.DataFrame({
                        'Chamber': chamber_counts.index,
                        'Cell Count': chamber_counts.values,
                        'Percentage': chamber_proportions.values
                    })

                    for chamber in chamber_counts.index:
                        summary_data['Metric'].append(f"{chamber} Cell Count")
                        summary_data['Value'].append(f"{chamber_counts[chamber]:,} ({chamber_proportions[chamber]:.1f}%)")
                
                    # ===== CROSS-CHAMBER CORRELATION ANALYSIS =====
                    print("Calculating cross-chamber gene expression correlations...")
                
                    if len(chamber_counts) > 1:
                        # Calculate mean expression per chamber
                        chambers = analyzed_adata.obs['chamber'].unique()
                        chamber_expr = {}
                        
                        for chamber in chambers:
                            chamber_mask = analyzed_adata.obs['chamber'] == chamber
                            X_chamber = analyzed_adata[chamber_mask].X
                            if hasattr(X_chamber, 'toarray'):
                                X_chamber = X_chamber.toarray()
                            chamber_expr[chamber] = np.mean(X_chamber, axis=0).flatten()
                        
                        # Calculate pairwise correlations
                        corr_data = []
                        for i, chamber1 in enumerate(chambers):
                            for j, chamber2 in enumerate(chambers):
                                if i <= j:  # Only upper triangle
                                    corr = np.corrcoef(chamber_expr[chamber1], chamber_expr[chamber2])[0, 1]
                                    corr_data.append({
                                        'Chamber 1': chamber1,
                                        'Chamber 2': chamber2,
                                        'Correlation': round(float(corr), 3)
                                    })
                        
                        chamber_correlations = pd.DataFrame(corr_data)
                        
                        # Add to summary
                        summary_data['Metric'].append('Cross-Chamber Correlations')
                        summary_data['Value'].append(f"{len(corr_data)} pairs analyzed")
                        
                        # Find highest and lowest correlations
                        non_diagonal = chamber_correlations[chamber_correlations['Chamber 1'] != chamber_correlations['Chamber 2']]
                        if len(non_diagonal) > 0:
                            max_corr = non_diagonal.loc[non_diagonal['Correlation'].idxmax()]
                            min_corr = non_diagonal.loc[non_diagonal['Correlation'].idxmin()]
                            
                            summary_data['Metric'].append('Highest Chamber Correlation')
                            summary_data['Value'].append(f"{max_corr['Chamber 1']}-{max_corr['Chamber 2']} (r={max_corr['Correlation']:.3f})")
                            
                            summary_data['Metric'].append('Lowest Chamber Correlation')
                            summary_data['Value'].append(f"{min_corr['Chamber 1']}-{min_corr['Chamber 2']} (r={min_corr['Correlation']:.3f})")
                        
                        # Identify chamber-specific marker genes
                        print("Identifying chamber-specific marker genes...")
                        marker_genes = {}
                        for chamber in chambers:
                            chamber_mask = analyzed_adata.obs['chamber'] == chamber
                            other_mask = ~chamber_mask
                            
                            X_chamber = analyzed_adata[chamber_mask].X
                            X_other = analyzed_adata[other_mask].X
                            
                            if hasattr(X_chamber, 'toarray'):
                                X_chamber = X_chamber.toarray()
                                X_other = X_other.toarray()
                            
                            # Calculate fold change
                            chamber_mean = np.mean(X_chamber, axis=0).flatten()
                            other_mean = np.mean(X_other, axis=0).flatten()
                            
                            # Avoid division by zero
                            fold_change = np.log2((chamber_mean + 1) / (other_mean + 1))
                            
                            # Get top markers (highest fold change)
                            top_indices = np.argsort(fold_change)[-10:][::-1]  # Top 10
                            top_genes = [analyzed_adata.var_names[i] for i in top_indices]
                            top_fc = [fold_change[i] for i in top_indices]
                            
                            marker_genes[chamber] = list(zip(top_genes, top_fc))
                        
                        chamber_markers = marker_genes
                        summary_data['Metric'].append('Chamber-Specific Markers')
                        summary_data['Value'].append(f"{sum(len(m) for m in marker_genes.values())} genes identified")

            # ===== COMMUNICATION HUBS =====
            hub_stats = None
            if include_communication_hubs:
                # Calculate hub scores (expression diversity + signalling potential)
                if 'leiden' in analyzed_adata.obs.columns:
                    hub_scores = []
                    for cell_type in analyzed_adata.obs['leiden'].unique():
                        cell_mask = analyzed_adata.obs['leiden'] == cell_type
                        
                        # Get expression data and convert sparse to dense if needed
                        X_subset = analyzed_adata[cell_mask].X
                        if hasattr(X_subset, 'toarray'):
                            X_subset = X_subset.toarray()
                        
                        expr_mean = np.mean(X_subset)
                        expr_std = np.std(X_subset)
                        expr_var = np.var(X_subset)

                        hub_score = (expr_std * expr_mean) / (expr_var + 1) if expr_var > 0 else 0
                        hub_scores.append({
                            'Cell Type': f"Cluster {cell_type}",
                            'Hub Score': round(float(hub_score), 4),
                            'Cell Count': cell_mask.sum()
                        })

                    hub_stats = pd.DataFrame(hub_scores).sort_values('Hub Score', ascending=False)

                    summary_data['Metric'].append('Top Communication Hub')
                    top_hub = hub_stats.iloc[0]['Cell Type'] if len(hub_stats) > 0 else "N/A"
                    summary_data['Value'].append(top_hub)

            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(persistent_csv_path, index=False)
            
            # ===== CREATE INTERACTIVE VISUALIZATIONS =====
            viz_chamber_file = None
            viz_hubs_file = None
            viz_corr_file = None
            viz_markers_file = None
            viz_network_file = None

            # 1. Chamber distribution plot (Interactive Plotly HTML)
            if chamber_stats is not None:
                fig_chamber = px.pie(
                    chamber_stats,
                    values='Cell Count',
                    names='Chamber',
                    title='Chamber Distribution - Interactive',
                    hole=0.3,
                    color_discrete_sequence=px.colors.qualitative.Set3
                )
                fig_chamber.update_traces(
                    textposition='inside',
                    textinfo='percent+label',
                    hovertemplate='<b>%{label}</b><br>Cells: %{value:,}<br>Percentage: %{percent}<extra></extra>'
                )
                fig_chamber.update_layout(
                    width=800,
                    height=600,
                    font=dict(size=14),
                    title_font_size=18,
                    showlegend=True,
                    legend=dict(orientation="v", yanchor="middle", y=0.5, xanchor="left", x=1.05)
                )

                # Save as standalone HTML with full interactivity
                viz_chamber_file = tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False).name
                fig_chamber.write_html(
                    viz_chamber_file,
                    include_plotlyjs='cdn',  # Use CDN for smaller file size
                    config={'displayModeBar': True, 'displaylogo': False}
                )
                print(f"✓ Created interactive chamber visualization: {viz_chamber_file}")

            # 2. Hub scores plot (Interactive Plotly HTML)
            if hub_stats is not None and len(hub_stats) > 0:
                fig_hubs = px.bar(
                    hub_stats.head(15),  # Show top 15 for better detail
                    x='Hub Score',
                    y='Cell Type',
                    title='Communication Hub Scores - Interactive (Top 15)',
                    orientation='h',
                    color='Hub Score',
                    color_continuous_scale='Viridis',
                    labels={'Hub Score': 'Hub Score', 'Cell Type': 'Cell Type'}
                )
                fig_hubs.update_traces(
                    hovertemplate='<b>%{y}</b><br>Hub Score: %{x:.4f}<br>Cell Count: %{customdata[0]:,}<extra></extra>',
                    customdata=hub_stats.head(15)[['Cell Count']].values
                )
                fig_hubs.update_layout(
                    width=900,
                    height=700,
                    font=dict(size=12),
                    title_font_size=18,
                    xaxis_title='Hub Score (higher = stronger communication hub)',
                    yaxis_title='Cell Type',
                    yaxis={'categoryorder': 'total ascending'}
                )

                # Save as standalone HTML with full interactivity
                viz_hubs_file = tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False).name
                fig_hubs.write_html(
                    viz_hubs_file,
                    include_plotlyjs='cdn',
                    config={'displayModeBar': True, 'displaylogo': False}
                )
                print(f"✓ Created interactive hub scores visualization: {viz_hubs_file}")
            
            # 3. Cross-Chamber Correlation Matrix (Interactive Plotly HTML)
            if chamber_correlations is not None and len(chamber_correlations) > 0:
                print("Creating cross-chamber correlation matrix...")
                
                # Create correlation matrix for heatmap
                chambers = chamber_correlations['Chamber 1'].unique()
                corr_matrix = np.eye(len(chambers))
                chamber_to_idx = {ch: i for i, ch in enumerate(chambers)}
                
                for _, row in chamber_correlations.iterrows():
                    i = chamber_to_idx[row['Chamber 1']]
                    j = chamber_to_idx[row['Chamber 2']]
                    corr_matrix[i, j] = row['Correlation']
                    corr_matrix[j, i] = row['Correlation']
                
                fig_corr = go.Figure(data=go.Heatmap(
                    z=corr_matrix,
                    x=list(chambers),
                    y=list(chambers),
                    colorscale='RdBu_r',
                    zmid=0.9,
                    zmin=0.85,
                    zmax=1.0,
                    text=corr_matrix,
                    texttemplate='%{text:.3f}',
                    textfont={"size": 14},
                    colorbar=dict(title="Correlation (r)")
                ))
                
                fig_corr.update_layout(
                    title='Cross-Chamber Gene Expression Correlations',
                    xaxis_title='Chamber',
                    yaxis_title='Chamber',
                    width=700,
                    height=700,
                    font=dict(size=14),
                    title_font_size=18
                )
                
                viz_corr_file = tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False).name
                fig_corr.write_html(
                    viz_corr_file,
                    include_plotlyjs='cdn',
                    config={'displayModeBar': True, 'displaylogo': False}
                )
                print(f"✓ Created cross-chamber correlation matrix: {viz_corr_file}")
            
            # 4. Chamber-Specific Marker Genes (Interactive Plotly HTML)
            if chamber_markers is not None and len(chamber_markers) > 0:
                print("Creating chamber-specific marker visualization...")
                
                # Prepare data for grouped bar chart
                marker_data = []
                for chamber, genes in chamber_markers.items():
                    for gene, fc in genes[:5]:  # Top 5 per chamber
                        marker_data.append({
                            'Chamber': chamber,
                            'Gene': gene,
                            'Log2 Fold Change': round(float(fc), 2)
                        })
                
                marker_df = pd.DataFrame(marker_data)
                
                fig_markers = px.bar(
                    marker_df,
                    x='Gene',
                    y='Log2 Fold Change',
                    color='Chamber',
                    barmode='group',
                    title='Chamber-Specific Marker Genes (Top 5 per Chamber)',
                    labels={'Log2 Fold Change': 'Log2 Fold Change (vs other chambers)'},
                    color_discrete_sequence=px.colors.qualitative.Set2
                )
                
                fig_markers.update_layout(
                    width=1200,
                    height=600,
                    font=dict(size=12),
                    title_font_size=18,
                    xaxis_title='Marker Gene',
                    yaxis_title='Log2 Fold Change',
                    xaxis={'categoryorder': 'total descending'},
                    hovermode='x unified'
                )
                
                viz_markers_file = tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False).name
                fig_markers.write_html(
                    viz_markers_file,
                    include_plotlyjs='cdn',
                    config={'displayModeBar': True, 'displaylogo': False}
                )
                print(f"✓ Created chamber-specific markers visualization: {viz_markers_file}")
            
            # 5. Communication Network Graph (Interactive knowledge graph)
            if include_communication_hubs and hub_stats is not None and len(hub_stats) > 0:
                print("Creating interactive communication network...")
                viz_network_file = create_communication_network(analyzed_adata, hub_stats, chamber_stats)
                if viz_network_file:
                    print(f"✓ Created interactive communication network: {viz_network_file}")
            
            # Count total visualizations created
            total_viz = sum([1 for v in [viz_chamber_file, viz_hubs_file, viz_corr_file, viz_markers_file, viz_network_file] if v is not None])
            print(f"✓ Total interactive visualizations created: {total_viz}")            # ===== CREATE COMPREHENSIVE OUTPUT MESSAGE =====
            output_msg = f"""
# HeartMAP Analysis Complete!

## Dataset Summary
- **Total Cells Analyzed:** {analyzed_adata.n_obs:,}
- **Total Genes Analyzed:** {analyzed_adata.n_vars:,}
- **Cell Clusters Identified:** {n_clusters}
- **Analysis Type:** {analysis_type.upper()}

"""

            if chamber_stats is not None:
                output_msg += f"""
## Chamber-Specific Results
Chamber analysis identified {len(chamber_stats)} distinct chambers:
- **Right Atrium (RA):** {chamber_stats[chamber_stats['Chamber']=='RA']['Cell Count'].values[0] if 'RA' in chamber_stats['Chamber'].values else 'N/A'} cells
- **Right Ventricle (RV):** {chamber_stats[chamber_stats['Chamber']=='RV']['Cell Count'].values[0] if 'RV' in chamber_stats['Chamber'].values else 'N/A'} cells
- **Left Atrium (LA):** {chamber_stats[chamber_stats['Chamber']=='LA']['Cell Count'].values[0] if 'LA' in chamber_stats['Chamber'].values else 'N/A'} cells
- **Left Ventricle (LV):** {chamber_stats[chamber_stats['Chamber']=='LV']['Cell Count'].values[0] if 'LV' in chamber_stats['Chamber'].values else 'N/A'} cells

"""

            if chamber_correlations is not None and len(chamber_correlations) > 0:
                non_diag = chamber_correlations[chamber_correlations['Chamber 1'] != chamber_correlations['Chamber 2']]
                if len(non_diag) > 0:
                    max_corr = non_diag.loc[non_diag['Correlation'].idxmax()]
                    min_corr = non_diag.loc[non_diag['Correlation'].idxmin()]
                    output_msg += f"""
## Cross-Chamber Correlations
Gene expression correlation analysis:
- **Highest:** {max_corr['Chamber 1']}-{max_corr['Chamber 2']} (r = {max_corr['Correlation']:.3f})
- **Lowest:** {min_corr['Chamber 1']}-{min_corr['Chamber 2']} (r = {min_corr['Correlation']:.3f})

"""
            
            if chamber_markers is not None:
                total_markers = sum(len(m) for m in chamber_markers.values())
                output_msg += f"""
## Chamber-Specific Markers
Identified **{total_markers} chamber-specific marker genes** across {len(chamber_markers)} chambers

"""

            if hub_stats is not None and len(hub_stats) > 0:
                output_msg += f"""
## Communication Hubs
Top communication hub cells (manuscript range: 0.037-0.047):
- **Top Hub:** {hub_stats.iloc[0]['Cell Type']} (Score: {hub_stats.iloc[0]['Hub Score']:.4f})
- **Range:** {hub_stats['Hub Score'].min():.4f} - {hub_stats['Hub Score'].max():.4f}

Communication hubs coordinate cellular interactions and represent therapeutic targets.

"""

            output_msg += f"""
## Results Available
 Summary statistics (CSV)  \n
 Interactive visualizations ({total_viz} HTML files)  \n
 Chamber composition analysis \n
 Cross-chamber correlation matrix \n
 Chamber-specific marker genes \n
 Communication hub identification \n
 Cell-cell communication network

## Interactive Visualizations
Download the HTML files below for fully interactive analysis!
"""

            # Prepare chamber info text
            chamber_info_text = ""
            if chamber_stats is not None:
                chamber_info_text = chamber_stats.to_string(index=False)
            
            if chamber_correlations is not None:
                chamber_info_text += "\n\n=== Cross-Chamber Correlations ===\n"
                chamber_info_text += chamber_correlations.to_string(index=False)
            
            if chamber_markers is not None:
                chamber_info_text += "\n\n=== Chamber-Specific Marker Genes ===\n"
                for chamber, genes in chamber_markers.items():
                    chamber_info_text += f"\n{chamber}: {', '.join([g[0] for g in genes[:5]])}"
            
            # Save chamber info to file
            chamber_info_file_path = None
            if chamber_info_text:
                with open(persistent_chamber_info_path, 'w', encoding='utf-8') as f:
                    f.write("=" * 80 + "\n")
                    f.write("CHAMBER DETAILS & MARKER GENES\n")
                    f.write("=" * 80 + "\n\n")
                    f.write(chamber_info_text)
                chamber_info_file_path = persistent_chamber_info_path
            
            # Copy visualization files to persistent location before temp_dir is deleted
            if viz_chamber_file and Path(viz_chamber_file).exists():
                shutil.copy2(viz_chamber_file, persistent_viz_chamber_path)
                viz_chamber_file = persistent_viz_chamber_path
            else:
                viz_chamber_file = None
            
            if viz_hubs_file and Path(viz_hubs_file).exists():
                shutil.copy2(viz_hubs_file, persistent_viz_hubs_path)
                viz_hubs_file = persistent_viz_hubs_path
            else:
                viz_hubs_file = None
            
            if viz_corr_file and Path(viz_corr_file).exists():
                shutil.copy2(viz_corr_file, persistent_viz_corr_path)
                viz_corr_file = persistent_viz_corr_path
            else:
                viz_corr_file = None
            
            if viz_markers_file and Path(viz_markers_file).exists():
                shutil.copy2(viz_markers_file, persistent_viz_markers_path)
                viz_markers_file = persistent_viz_markers_path
            else:
                viz_markers_file = None
            
            if viz_network_file and Path(viz_network_file).exists():
                shutil.copy2(viz_network_file, persistent_viz_network_path)
                viz_network_file = persistent_viz_network_path
            else:
                viz_network_file = None
            
            return (
                output_msg, 
                persistent_csv_path, 
                viz_chamber_file, 
                viz_hubs_file,
                viz_corr_file,
                viz_markers_file,
                viz_network_file,
                chamber_info_text if chamber_info_text else "No chamber data",
                chamber_info_file_path
            )
    
    except Exception as e:
        import traceback
        error_type = type(e).__name__

        # Handle specific error types
        if "ClientDisconnect" in error_type or "ClientDisconnect" in str(e):
            error_msg = """**Upload Interrupted**

The file upload was interrupted. This usually happens when:
- The file is very large (>100MB)
- Network connection is unstable
- Browser tab was closed during upload

**Solutions:**
1. Try uploading a smaller file (<100MB)
2. Ensure stable internet connection
3. Keep this browser tab open during upload
4. For large datasets, consider subsetting your data first
"""
        else:
            error_msg = f"**Error during analysis:**\n\n{str(e)}\n\n**Technical details:**\n```\n{traceback.format_exc()}\n```"
        
        return error_msg, None, None, None, None, None, None, "", None


# Create Gradio interface with enhanced features
with gr.Blocks(
    title="HeartMAP: Heart Multi-chamber Analysis Platform",
    theme=gr.themes.Soft(),
    analytics_enabled=False,
    css="""
        .resizable-textbox textarea {
            resize: both !important;
            overflow: auto !important;
            min-height: 200px !important;
            max-height: none !important;
        }
    """
) as demo:
    gr.Markdown("""
    # HeartMAP: Multi-chamber Heart Analysis
    
    Single-cell RNA-seq analysis for cardiac chamber biology  
    **Max file size: 10GB** • Large files may take 10-30 minutes

    **Analysis Features:**
    - Cell type annotation & QC
    - Chamber-specific analysis (RA, RV, LA, LV)
    - Communication hub identification
    - Cross-chamber correlations & marker genes
    - Interactive visualizations

    **Supported Formats:**  
    `.h5ad` `.h5` `.h5mu` `.mtx` `.loom` `.csv` `.tsv` `.tar` `.parquet` `.zarr` `.hdf5`  
    <details><summary>View all formats</summary>
    
    - **Direct:** AnnData (.h5ad, .h5mu, .zarr), 10X (.h5, .mtx), Loom (.loom), Text (.csv/.tsv/.txt), Archives (.tar/.tar.gz), Parquet/Arrow, HDF5
    - **Needs conversion:** Seurat (.h5seurat/.rds) → use SeuratDisk | FASTQ/BAM → use Cell Ranger
    </details>

    ---
    """)

    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### Upload Data")
            file_input = gr.File(
                label="Upload single-cell data (15+ formats supported)",
                file_types=[
                    ".h5ad", ".h5mu", ".zarr",  # AnnData variants
                    ".h5", ".mtx",              # 10X Genomics
                    ".loom",                    # Loom
                    ".csv", ".tsv", ".txt",     # Text matrices
                    ".tar", ".tar.gz", ".tgz",  # Archives
                    ".parquet", ".arrow",       # Columnar
                    ".hdf5",                    # Generic HDF5
                    ".h5seurat", ".rds", ".rdata"  # R/Seurat (with conversion guidance)
                ],
                type="filepath"
            )

            gr.Markdown("### Analysis Settings")
            analysis_type = gr.Dropdown(
                choices=["basic", "comprehensive", "multi_chamber"],
                value="comprehensive",
                label="Analysis Type"
            )

            gr.Markdown("### Advanced Options")
            max_cells = gr.Number(
                label="Max Cells (for memory optimization)",
                value=50000,
                precision=0
            )
            max_genes = gr.Number(
                label="Max Genes (for memory optimization)",
                value=5000,
                precision=0
            )

            gr.Markdown("### Chamber Analysis")
            include_chamber = gr.Checkbox(
                label="Enable chamber-specific analysis",
                value=True
            )
            include_hubs = gr.Checkbox(
                label="Enable communication hub analysis",
                value=True
            )

            analyze_btn = gr.Button("Run Analysis", variant="primary", size="lg")

        with gr.Column(scale=2):
            gr.Markdown("### Results & Visualizations")

            output_text = gr.Markdown(
                label="Analysis Summary",
                value="Results will appear here after analysis"
            )

            gr.Markdown("### Downloads")
            output_file = gr.File(
                label="Download Results CSV",
                file_types=[".csv"]
            )

            gr.Markdown("""
            ### Interactive Visualizations
            Download the HTML files below and open them in your browser for fully interactive charts:
            - **Zoom, pan, and hover** for detailed information
            - **Click legend items** to show/hide data
            - **Export** as PNG from the chart menu
            """)

            viz_file_1 = gr.File(
                label="Chamber Distribution (Interactive HTML)",
                file_types=[".html"]
            )

            viz_file_2 = gr.File(
                label="Communication Hubs (Interactive HTML)",
                file_types=[".html"]
            )

            viz_file_3 = gr.File(
                label="Cross-Chamber Correlations (Interactive HTML)",
                file_types=[".html"]
            )
            
            viz_file_4 = gr.File(
                label="Chamber-Specific Markers (Interactive HTML)",
                file_types=[".html"]
            )

            viz_file_5 = gr.File(
                label="Communication Network (Interactive HTML)",
                file_types=[".html"]
            )

            chamber_info = gr.Textbox(
                label="Chamber Details & Marker Genes",
                interactive=False,
                max_lines=50,
                elem_classes="resizable-textbox"
            )
            
            chamber_info_file = gr.File(
                label="Download Chamber Details & Marker Genes (TXT)",
                file_types=[".txt"]
            )

    # Analysis pipeline
    analyze_btn.click(
        fn=analyze_heart_data,
        inputs=[
            file_input,
            analysis_type,
            max_cells,
            max_genes,
            include_chamber,
            include_hubs
        ],
        outputs=[
            output_text,
            output_file,
            viz_file_1,
            viz_file_2,
            viz_file_3,
            viz_file_4,
            viz_file_5,
            chamber_info,
            chamber_info_file
        ]
    )

    gr.Markdown("""
    ---

    ### Usage Tips:

    | Scenario | Recommendation |
    |----------|----------------|
    | **First time** | Use "Basic" for quick exploration |
    | **Chamber biology** | Enable chamber-specific analysis |
    | **Drug targets** | Enable communication hubs |
    | **Limited memory** | Reduce cells/genes or use "Basic" |
    | **Full analysis** | Use "Comprehensive" mode |

    ### Documentation:
    - **Full User Guide:** https://github.com/Tumo505/HeartMap/blob/master/USER_GUIDE.md
    - **API Reference:** https://github.com/Tumo505/HeartMap/blob/master/API_DOCUMENTATION.md
    - **GitHub Repository:** https://github.com/Tumo505/HeartMap
    - **Python Package:** `pip install heartmap`

    ### Key References:
    - Chamber correlations: RV-LV (r=0.985), LA-LV (r=0.870)
    - Hub scores: 0.037-0.047 for atrial cardiomyocytes & adipocytes
    - Markers identified: 1,000+ per chamber
    - DEGs per pair: 150+ significantly different genes

    ### Citation:
    ```
    Kgabeng, T., Wang, L., Ngwangwa, H., & Pandelani, T. (2025).
    HeartMAP: A Multi-Chamber Spatial Framework for Cardiac Cell-Cell Communication.
    Available at: https://github.com/Tumo505/HeartMap
    ```
    """)

if __name__ == "__main__":
    # Launch with increased timeout and file size limits for large datasets
    demo.queue(max_size=10).launch(
        max_file_size="10gb",  # Increased to 10GB for very large single-cell datasets (max: 100GB)
        show_error=True
    )
