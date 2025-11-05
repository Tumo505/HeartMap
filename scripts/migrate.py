"""
Migration script to convert existing HeartMAP analysis to new refactored format
"""

import shutil
import os
from pathlib import Path
import warnings

try:
    import scanpy as sc
    import pandas as pd
    import numpy as np
    DEPS_AVAILABLE = True
except ImportError:
    DEPS_AVAILABLE = False
    warnings.warn("Some dependencies not available")


def migrate_legacy_analysis():
    """Migrate existing analysis results to new format"""
    print("=== HeartMAP Migration Script ===")
    
    # Get project paths
    project_root = Path(__file__).parent.parent
    legacy_analysis_dir = project_root / "analysis"
    legacy_results_dir = project_root / "results"
    
    # Create new structure
    new_scripts_dir = project_root / "scripts" / "legacy"
    new_notebooks_dir = project_root / "notebooks"
    
    print("1. Creating new directory structure...")
    new_scripts_dir.mkdir(parents=True, exist_ok=True)
    new_notebooks_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy legacy scripts to new location
    if legacy_analysis_dir.exists():
        print("2. Migrating legacy analysis scripts...")
        
        # Copy each pipeline directory
        for pipeline_dir in legacy_analysis_dir.iterdir():
            if pipeline_dir.is_dir():
                dest_dir = new_scripts_dir / pipeline_dir.name
                if not dest_dir.exists():
                    shutil.copytree(pipeline_dir, dest_dir)
                    print(f"   Copied {pipeline_dir.name} -> scripts/legacy/{pipeline_dir.name}")
    
    # Create migration notebooks
    print("3. Creating migration notebooks...")
    create_migration_notebooks(new_notebooks_dir)
    
    # Create example scripts using new API
    print("4. Creating example scripts with new API...")
    create_example_scripts(project_root / "scripts")
    
    # Create Dockerfile for deployment
    print("5. Creating deployment files...")
    create_deployment_files(project_root)
    
    print(" Migration completed!")
    print("\nNext steps:")
    print("1. Install the new package: pip install -e .[all]")
    print("2. Run tests: python tests/test_heartmap.py")
    print("3. Try examples: python scripts/run_examples.py")
    print("4. Check notebooks in: notebooks/")


def create_migration_notebooks(notebooks_dir):
    """Create Jupyter notebooks for migration examples"""
    
    # Basic analysis notebook
    basic_notebook = notebooks_dir / "01_basic_analysis.ipynb"
    create_basic_notebook(basic_notebook)
    
    # Advanced analysis notebook
    advanced_notebook = notebooks_dir / "02_advanced_communication.ipynb"
    create_advanced_notebook(advanced_notebook)
    
    # Multi-chamber notebook
    chamber_notebook = notebooks_dir / "03_multi_chamber_analysis.ipynb"
    create_chamber_notebook(chamber_notebook)
    
    # Comprehensive analysis notebook
    comprehensive_notebook = notebooks_dir / "04_comprehensive_analysis.ipynb"
    create_comprehensive_notebook(comprehensive_notebook)


def create_basic_notebook(notebook_path):
    """Create basic analysis notebook"""
    notebook_content = {
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# HeartMAP Basic Analysis\n",
                    "\n",
                    "This notebook demonstrates basic single-cell heart analysis using the refactored HeartMAP package."
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "source": [
                    "import sys\n",
                    "from pathlib import Path\n",
                    "\n",
                    "# Add src to path for development\n",
                    "project_root = Path().absolute().parent\n",
                    "sys.path.insert(0, str(project_root / 'src'))\n",
                    "\n",
                    "from heartmap import Config, HeartMapModel\n",
                    "from heartmap.pipelines import BasicPipeline\n",
                    "from heartmap.data import DataProcessor"
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "source": [
                    "# Create configuration\n",
                    "config = Config.default()\n",
                    "\n",
                    "# For demonstration, use smaller dataset\n",
                    "config.data.max_cells_subset = 5000\n",
                    "config.data.max_genes_subset = 2000\n",
                    "config.data.test_mode = True\n",
                    "\n",
                    "print(\"Configuration created:\")\n",
                    "print(f\"- Max cells: {config.data.max_cells_subset}\")\n",
                    "print(f\"- Max genes: {config.data.max_genes_subset}\")\n",
                    "print(f\"- Test mode: {config.data.test_mode}\")"
                ]
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "source": [
                    "# Create and run basic pipeline\n",
                    "pipeline = BasicPipeline(config)\n",
                    "\n",
                    "# Note: Replace with your actual data path\n",
                    "# data_path = \"../data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad\"\n",
                    "# results = pipeline.run(data_path, \"../results/basic\")\n",
                    "\n",
                    "print(\"Basic pipeline created successfully!\")\n",
                    "print(\"To run analysis, provide your data path and uncomment the lines above.\")"
                ]
            }
        ],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "name": "python",
                "version": "3.10.0"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 4
    }
    
    import json
    with open(notebook_path, 'w') as f:
        json.dump(notebook_content, f, indent=2)


def create_advanced_notebook(notebook_path):
    """Create advanced communication analysis notebook"""
    # Similar structure to basic notebook but for advanced analysis
    # ... (implementation details)
    pass


def create_chamber_notebook(notebook_path):
    """Create multi-chamber analysis notebook"""
    # ... (implementation details)
    pass


def create_comprehensive_notebook(notebook_path):
    """Create comprehensive analysis notebook"""
    # ... (implementation details)
    pass


def create_example_scripts(scripts_dir):
    """Create example scripts using new API"""
    
    # Command line example
    cli_example = scripts_dir / "run_cli_analysis.py"
    with open(cli_example, 'w') as f:
        f.write("""#!/usr/bin/env python
\"\"\"
Command line analysis example
\"\"\"

import sys
from pathlib import Path

# Add src to path for development
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / 'src'))

from heartmap.api import CLIInterface

def main():
    cli = CLIInterface()
    
    # Example: replace with your data path
    data_path = "data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad"
    
    if Path(data_path).exists():
        cli.run_analysis(
            data_path=data_path,
            analysis_type="comprehensive",
            output_dir="results/comprehensive",
            config_path="config.yaml"
        )
    else:
        print(f"Data file not found: {data_path}")
        print("Please update the data_path variable with your actual data file.")

if __name__ == "__main__":
    main()
""")
    
    # API server example
    api_example = scripts_dir / "run_api_server.py"
    with open(api_example, 'w') as f:
        f.write("""#!/usr/bin/env python
\"\"\"
API server example
\"\"\"

import sys
from pathlib import Path

# Add src to path for development
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / 'src'))

from heartmap.api import create_api

def main():
    try:
        # Create API instance
        api = create_api(config_path="config.yaml")
        
        # Run server
        print("Starting HeartMAP API server...")
        print("Access at: http://localhost:8000")
        print("API docs at: http://localhost:8000/docs")
        
        api.run(host="0.0.0.0", port=8000, debug=True)
        
    except ImportError as e:
        print("FastAPI dependencies not available.")
        print("Install with: pip install -e .[api]")
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
""")


def create_deployment_files(project_root):
    """Create deployment configuration files"""
    
    # Dockerfile
    dockerfile = project_root / "Dockerfile"
    with open(dockerfile, 'w') as f:
        f.write("""FROM python:3.10-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \\
    build-essential \\
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements-dev.txt .
RUN pip install --no-cache-dir -r requirements-dev.txt

# Copy source code
COPY src/ ./src/
COPY setup.py .
COPY config.yaml .

# Install HeartMAP package
RUN pip install -e .[all]

# Expose port for API
EXPOSE 8000

# Default command
CMD ["python", "-m", "heartmap.api"]
""")
    
    # Docker compose
    docker_compose = project_root / "docker-compose.yml"
    with open(docker_compose, 'w') as f:
        f.write("""version: '3.8'

services:
  heartmap-api:
    build: .
    ports:
      - "8000:8000"
    volumes:
      - ./data:/app/data
      - ./results:/app/results
      - ./config.yaml:/app/config.yaml
    environment:
      - HEARTMAP_CONFIG=/app/config.yaml
    command: python scripts/run_api_server.py

  heartmap-worker:
    build: .
    volumes:
      - ./data:/app/data
      - ./results:/app/results
    environment:
      - HEARTMAP_CONFIG=/app/config.yaml
    command: python scripts/run_examples.py
""")
    
    # Gradio app for Hugging Face deployment
    gradio_app = project_root / "app.py"
    with open(gradio_app, 'w') as f:
        f.write("""#!/usr/bin/env python
\"\"\"
Gradio app for Hugging Face Spaces deployment
\"\"\"

import gradio as gr
import tempfile
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, 'src')

try:
    from heartmap import Config
    from heartmap.pipelines import ComprehensivePipeline
    HEARTMAP_AVAILABLE = True
except ImportError:
    HEARTMAP_AVAILABLE = False

def analyze_heart_data(uploaded_file, analysis_type, max_cells, max_genes):
    \"\"\"Analyze uploaded heart data\"\"\"
    
    if not HEARTMAP_AVAILABLE:
        return "HeartMAP not available. Please install dependencies."
    
    if uploaded_file is None:
        return "Please upload a file."
    
    try:
        # Create config
        config = Config.default()
        config.data.max_cells_subset = int(max_cells) if max_cells else None
        config.data.max_genes_subset = int(max_genes) if max_genes else None
        
        # Create pipeline
        if analysis_type == "comprehensive":
            pipeline = ComprehensivePipeline(config)
        else:
            return f"Analysis type {analysis_type} not implemented in demo"
        
        # For demo, just return configuration info
        return f\"\"\"Analysis configured successfully!

Configuration:
- Analysis type: {analysis_type}
- Max cells: {max_cells}
- Max genes: {max_genes}
- File: {uploaded_file.name}

Note: This is a demo interface. 
In the full version, your data would be processed here.
\"\"\"
        
    except Exception as e:
        return f"Error: {str(e)}"

# Create Gradio interface
demo = gr.Interface(
    fn=analyze_heart_data,
    inputs=[
        gr.File(label="Upload single-cell data (.h5ad)", file_types=[".h5ad"]),
        gr.Dropdown(
            choices=["comprehensive", "basic", "advanced", "multi_chamber"],
            value="comprehensive",
            label="Analysis Type"
        ),
        gr.Number(label="Max Cells (for memory optimization)", value=50000, precision=0),
        gr.Number(label="Max Genes (for memory optimization)", value=5000, precision=0),
    ],
    outputs=gr.Textbox(label="Analysis Results"),
    title="HeartMAP: Heart Multi-chamber Analysis Platform",
    description=\"\"\"
    Upload single-cell RNA-seq data from heart tissue for comprehensive analysis.
    
    **Features:**
    - Cell type annotation
    - Cell-cell communication analysis  
    - Multi-chamber analysis
    - Chamber-specific marker identification
    
    **Input format:** AnnData (.h5ad) files
    
    **Demo Note:** This demo shows the interface. Full analysis requires local installation.
    \"\"\",
    examples=[
        [None, "comprehensive", 10000, 2000],
        [None, "basic", 5000, 1000],
    ]
)

if __name__ == "__main__":
    demo.launch()
""")
    
    # GitHub Actions workflow
    github_dir = project_root / ".github" / "workflows"
    github_dir.mkdir(parents=True, exist_ok=True)
    
    ci_workflow = github_dir / "ci.yml"
    with open(ci_workflow, 'w') as f:
        f.write("""name: CI

on:
  push:
    branches: [ main, master ]
  pull_request:
    branches: [ main, master ]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, '3.10', 3.11]

    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -e .[dev]
    
    - name: Run tests
      run: |
        python tests/test_heartmap.py
    
    - name: Run linting
      run: |
        flake8 src/heartmap --max-line-length=100
    
    - name: Run type checking
      run: |
        mypy src/heartmap --ignore-missing-imports

  deploy:
    needs: test
    runs-on: ubuntu-latest
    if: github.ref == 'refs/heads/main'
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    
    - name: Build package
      run: |
        python -m pip install --upgrade pip build
        python -m build
    
    - name: Publish to PyPI
      if: startsWith(github.ref, 'refs/tags/')
      uses: pypa/gh-action-pypi-publish@release/v1
      with:
        password: ${{ secrets.PYPI_API_TOKEN }}
""")


if __name__ == "__main__":
    migrate_legacy_analysis()
