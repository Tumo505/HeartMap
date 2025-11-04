#!/usr/bin/env python
"""
Gradio app for Hugging Face Spaces deployment
"""

import gradio as gr
import tempfile
import sys
import shutil
import pandas as pd
from pathlib import Path

# Add src to path
sys.path.insert(0, 'src')

try:
    from heartmap import Config
    from heartmap.pipelines import BasicPipeline, ComprehensivePipeline
    HEARTMAP_AVAILABLE = True
except ImportError:
    HEARTMAP_AVAILABLE = False

def analyze_heart_data(uploaded_file, analysis_type, max_cells, max_genes):
    """Analyze uploaded heart data"""
    
    if not HEARTMAP_AVAILABLE:
        return "HeartMAP not available. Please install dependencies.", None, []
    
    if uploaded_file is None:
        return "Please upload a file.", None, []
    
    # Create persistent CSV file outside any context
    persistent_csv = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
    persistent_csv_path = persistent_csv.name
    persistent_csv.close()
    
    try:
        # Create temporary directory for results
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_dir = temp_path / "results"
            results_dir.mkdir(exist_ok=True)
            
            # Create config
            config = Config.default()
            config.data.max_cells_subset = int(max_cells) if max_cells else None
            config.data.max_genes_subset = int(max_genes) if max_genes else None
            config.paths.processed_data_dir = str(temp_path / "processed")
            config.paths.results_dir = str(results_dir)
            
            # Ensure all directories exist
            config.create_directories()
            
            # Create pipeline based on type
            if analysis_type == "comprehensive":
                pipeline = ComprehensivePipeline(config)
            elif analysis_type == "basic":
                pipeline = BasicPipeline(config)
            else:
                return f"Analysis type '{analysis_type}' not implemented yet. Use 'basic' or 'comprehensive'.", None, []
            
            # Run pipeline
            results = pipeline.run(str(uploaded_file), str(results_dir))
            
            # Create summary CSV
            summary_data = {
                'metric': [],
                'value': []
            }
            
            if 'adata' in results:
                adata = results['adata']
                summary_data['metric'].extend([
                    'Total Cells',
                    'Total Genes',
                    'Number of Clusters'
                ])
                summary_data['value'].extend([
                    adata.n_obs,
                    adata.n_vars,
                    adata.obs['leiden'].nunique() if 'leiden' in adata.obs.columns else 0
                ])
                
                # Add cluster distribution
                if 'leiden' in adata.obs.columns:
                    cluster_counts = adata.obs['leiden'].value_counts()
                    for cluster, count in cluster_counts.items():
                        summary_data['metric'].append(f'Cluster {cluster} Count')
                        summary_data['value'].append(count)
            
            # Write directly to persistent CSV file (outside temp context)
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(persistent_csv_path, index=False)
            
            # Collect visualization files and copy to persistent locations
            viz_dir = results_dir / "figures"
            visualization_files = []
            
            if viz_dir.exists():
                # Key visualization files to display
                viz_files = [
                    "umap_clusters.png",
                    "qc_metrics.png"
                ]
                
                # For comprehensive pipeline, also include dashboard
                if analysis_type == "comprehensive":
                    viz_files.append("comprehensive_dashboard.png")
                
                for viz_file in viz_files:
                    source_file = viz_dir / viz_file
                    if source_file.exists():
                        # Create persistent copy (binary mode for images)
                        persistent_viz = tempfile.NamedTemporaryFile(
                            mode='wb', suffix='.png', delete=False
                        )
                        persistent_viz_path = persistent_viz.name
                        persistent_viz.close()
                        
                        # Copy visualization to persistent location
                        shutil.copy2(source_file, persistent_viz_path)
                        visualization_files.append(persistent_viz_path)
            
            # Create summary message
            summary_msg = f"""Analysis Complete!

**Dataset Summary:**
- Total Cells: {summary_data['value'][0] if summary_data['value'] else 'N/A'}
- Total Genes: {summary_data['value'][1] if len(summary_data['value']) > 1 else 'N/A'}
- Clusters Identified: {summary_data['value'][2] if len(summary_data['value']) > 2 else 'N/A'}

**Results saved to:**
- Results CSV: heartmap_results.csv
- Visualizations: {len(visualization_files)} figure(s) displayed below

**Analysis Type:** {analysis_type}
"""
            
            # Return text, CSV, and visualization images
            # Return empty list if no visualizations (Gradio Gallery expects list)
            return summary_msg, persistent_csv_path, visualization_files if visualization_files else []
        
    except Exception as e:
        import traceback
        error_msg = f"❌ Error: {str(e)}\n\nTraceback:\n{traceback.format_exc()}"
        return error_msg, None, []

# Create Gradio interface
with gr.Blocks(title="HeartMAP: Heart Multi-chamber Analysis Platform") as demo:
    gr.Markdown("""
    # HeartMAP: Heart Multi-chamber Analysis Platform
    
    Upload single-cell RNA-seq data from heart tissue for comprehensive analysis.
    
    **Features:**
    - Cell type annotation and clustering
    - Quality control metrics
    - UMAP visualization
    - Cell-cell communication analysis (comprehensive mode)
    - Multi-chamber analysis
    
    **Input format:** AnnData (.h5ad) files
    """)
    
    with gr.Row():
        with gr.Column():
            file_input = gr.File(
                label="Upload single-cell data (.h5ad)",
                file_types=[".h5ad"],
                type="filepath"
            )
            analysis_type = gr.Dropdown(
                choices=["basic", "comprehensive"],
                value="basic",
                label="Analysis Type"
            )
            max_cells = gr.Number(
                label="Max Cells (for memory optimization)",
                value=10000,
                precision=0
            )
            max_genes = gr.Number(
                label="Max Genes (for memory optimization)",
                value=2000,
                precision=0
            )
            analyze_btn = gr.Button("Run Analysis", variant="primary")
        
        with gr.Column():
            output_text = gr.Textbox(
                label="Analysis Results",
                lines=10,
                interactive=False
            )
            output_file = gr.File(
                label="Download Results CSV",
                file_types=[".csv"]
            )
            output_gallery = gr.Gallery(
                label="Visualizations",
                show_label=True,
                elem_id="gallery",
                columns=2,
                rows=2,
                height="auto"
            )
    
    analyze_btn.click(
        fn=analyze_heart_data,
        inputs=[file_input, analysis_type, max_cells, max_genes],
        outputs=[output_text, output_file, output_gallery]
    )
    
    gr.Markdown("""
    ### Tips:
    - Start with **basic** analysis for faster results
    - Use **comprehensive** for full analysis including communication
    - Reduce max cells/genes if you encounter memory issues
    - Results CSV contains summary statistics
    """)

if __name__ == "__main__":
    # On Hugging Face Spaces, the app is already publicly accessible
    # No need for share=True (it's not supported anyway)
    # Public link: https://huggingface.co/spaces/Tumo505/HeartMAP
    demo.launch()
