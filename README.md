# HeartMAP: Heart Multi-chamber Analysis Platform

[![PyPI version](https://badge.fury.io/py/heartmap.svg)](https://badge.fury.io/py/heartmap)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![CI Status](https://github.com/Tumo505/HeartMap/workflows/CI/badge.svg)](https://github.com/Tumo505/HeartMap/actions)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17512107.svg)](https://doi.org/10.5281/zenodo.17512107)

> **A production-ready Python package for comprehensive single-cell heart analysis with chamber-specific insights**

## Package Overview

HeartMAP is a specialized bioinformatics package that decodes cellular communication across all four chambers of the human heart. Unlike general single-cell tools, HeartMAP is purpose-built for cardiac biology, offering chamber-specific insights crucial for understanding heart function, disease, and therapeutic opportunities.

**Key Features:**
- **Production Ready**: Fully tested, documented, and deployed on PyPI
- **Multiple Interfaces**: CLI, Python API, REST API, and Web interface  
- **Easy Installation**: `pip install heartmap`
- **Configurable**: Works on 8GB+ RAM with memory optimization
- **Validated**: Tested on real human heart datasets
- **Comprehensive**: From basic QC to advanced communication analysis

## Quick Installation

```bash
# Install from PyPI
pip install heartmap

# Install with all features
pip install heartmap[all]

# Verify installation
python -c "import heartmap; print('✅ HeartMAP ready!')"
```

## Quick Start

### 30-Second Analysis
```bash
# Analyze your heart data with one command
heartmap your_heart_data.h5ad
```

### 2-Minute Python Analysis
```python
from heartmap import Config
from heartmap.pipelines import ComprehensivePipeline

# Quick analysis
config = Config.default()
pipeline = ComprehensivePipeline(config)
results = pipeline.run('your_data.h5ad', 'results/')

print("✅ Analysis complete! Check 'results/' directory.")
```

## Documentation

| Document | Description | Use When |
|----------|-------------|----------|
| **[User Guide](USER_GUIDE.md)** | Complete step-by-step tutorials | Learning HeartMAP |
| **[API Documentation](API_DOCUMENTATION.md)** | Full API reference | Programming with HeartMAP |
| **[Package README](PACKAGE_README.md)** | Package-specific documentation | Installing/using the package |
| **[Original README](README_ORIGINAL.md)** | Development documentation | Contributing to HeartMAP |

## What HeartMAP Can Do

### Analysis Pipelines

| Pipeline | Purpose | Output | Runtime |
|----------|---------|---------|---------|
| **Basic** | Quality control, cell typing | Cell annotations, QC metrics | 5-10 min |
| **Communication** | Cell-cell interactions | Communication networks, hubs | 10-15 min |
| **Multi-Chamber** | Chamber-specific analysis | Chamber markers, comparisons | 15-20 min |
| **Comprehensive** | Complete analysis | All of the above + reports | 20-30 min |

### Real-World Applications

```python
# Clinical Research: Chamber-specific targets
from heartmap.pipelines import MultiChamberPipeline
pipeline = MultiChamberPipeline(config)
results = pipeline.run('patient_data.h5ad')
lv_targets = results['chamber_markers']['LV']

# Drug Discovery: Communication pathways
from heartmap.pipelines import AdvancedCommunicationPipeline  
pipeline = AdvancedCommunicationPipeline(config)
results = pipeline.run('disease_data.h5ad')
drug_targets = results['communication_hubs']

# Education: Comparative analysis
results1 = pipeline.run('healthy_heart.h5ad')
results2 = pipeline.run('diseased_heart.h5ad')
```

## Performance

| Hardware | Dataset Size | Memory | Runtime | Status |
|----------|-------------|--------|---------|---------|
| 8GB RAM | 30K cells | ~6GB | 15 min | ✅ Recommended |
| 16GB RAM | 50K cells | ~12GB | 25 min | ✅ Optimal |
| 32GB RAM | 100K cells | ~24GB | 45 min | ✅ Production |

## Development

### For Contributors

```bash
# Development setup
git clone https://github.com/Tumo505/HeartMap.git
cd HeartMap
pip install -e .[dev]

# Run tests
python -m pytest tests/
python -m flake8 src/heartmap/
python -m mypy src/heartmap/
```

### For Package Users

The package is production-ready and maintained. See the [User Guide](USER_GUIDE.md) for complete usage instructions.

## Scientific Impact

- **Clinical**: Chamber-specific therapeutic strategies
- **Research**: First comprehensive multi-chamber communication atlas  
- **Education**: Accessible cardiac biology analysis platform
- **Industry**: Production-ready bioinformatics tool

## Use Cases

- **Pharmaceutical Research**: Drug target discovery, safety assessment
- **Clinical Cardiology**: Precision medicine, disease mechanisms
- **Basic Research**: Cardiac development, evolutionary biology
- **Computational Biology**: Method benchmarking, data integration

## Requirements

- **Python**: 3.8+ (tested on 3.8-3.11)
- **Memory**: 8GB+ RAM (configurable)
- **Storage**: 2GB+ for package and results
- **OS**: Linux, macOS, Windows

## Contributing

We welcome contributions! See our [development documentation](README_ORIGINAL.md) for details.

## Citation

```bibtex
@software{heartmap2025,
  title={HeartMAP: Heart Multi-chamber Analysis Platform},
  author={Kgabeng, Tumo and Wang, Lulu and Ngwangwa, Harry and Pandelani, Thanyani},
  year={2025},
  url={https://github.com/Tumo505/HeartMap},
  version={1.0.0}
}
```

## License

Apache 2.0 License - see [LICENSE](LICENSE) for details.

## Support

- **Documentation**: [User Guide](USER_GUIDE.md) | [API Docs](API_DOCUMENTATION.md)
- **Community**: [GitHub Discussions](https://github.com/Tumo505/HeartMap/discussions)
- **Issues**: [GitHub Issues](https://github.com/Tumo505/HeartMap/issues)
- **Contact**: 28346416@mylife.unisa.ac.za

---

**HeartMAP: Production-ready cardiac single-cell analysis for researchers worldwide**

*Install today: `pip install heartmap`*

### Supporting Files

```
├── scripts/                 # Setup and utility scripts
│   ├── setup.sh            # Automated environment setup
│   ├── validate.py         # Installation validation
│   ├── migrate.py          # Legacy code migration
│   └── deploy_huggingface.sh # HuggingFace deployment
├── tests/                   # Comprehensive test suite
├── notebooks/               # Jupyter notebook examples
├── config.yaml             # Default configuration
├── setup.py                # Package installation
├── Dockerfile              # Container deployment
├── docker-compose.yml      # Multi-service orchestration
└── app.py                  # Gradio web interface
```

## Quick Start

### Option 1: Automated Setup (Recommended)

```bash
# Clone the repository
git clone https://github.com/Tumo505/HeartMap.git
cd HeartMap

# Run automated setup script
./scripts/setup.sh

# Activate the environment
source heartmap_env/bin/activate  # Linux/Mac
# OR: heartmap_env\Scripts\activate  # Windows

# Validate installation
python scripts/validate.py

# Start analyzing!
heartmap data/raw/your_data.h5ad --analysis-type comprehensive
```

### Option 2: Manual Installation

```bash
# Create virtual environment
python3 -m venv heartmap_env
source heartmap_env/bin/activate

# Install dependencies
pip install -r requirements-dev.txt

# Install HeartMAP in development mode
pip install -e .[all]

# Validate installation
python scripts/validate.py
```

### Option 3: Package Installation

```bash
# Install from PyPI (when available)
pip install heartmap[all]

# Or install specific features
pip install heartmap[communication]  # Communication analysis only
pip install heartmap[api]            # API features only
```

## Usage Examples

### 1. Command Line Interface

```bash
# Basic analysis
heartmap data/raw/heart_data.h5ad

# Comprehensive analysis with custom output
heartmap data/raw/heart_data.h5ad \
    --analysis-type comprehensive \
    --output-dir results/comprehensive \
    --config my_config.yaml

# Specific analysis types
heartmap data/raw/heart_data.h5ad --analysis-type annotation
heartmap data/raw/heart_data.h5ad --analysis-type communication  
heartmap data/raw/heart_data.h5ad --analysis-type multi-chamber

# Memory-optimized for large datasets
heartmap data/raw/large_dataset.h5ad \
    --analysis-type comprehensive \
    --config config_large.yaml
```

### 2. Python API

```python
from heartmap import Config, HeartMapModel
from heartmap.pipelines import ComprehensivePipeline

# Load and customize configuration
config = Config.from_yaml('config.yaml')
config.data.max_cells_subset = 50000  # Optimize for your memory
config.data.max_genes_subset = 5000

# Option A: Use full HeartMAP model
model = HeartMapModel(config)
results = model.analyze('data/raw/heart_data.h5ad')

# Option B: Use specific pipeline
pipeline = ComprehensivePipeline(config)
results = pipeline.run('data/raw/heart_data.h5ad', 'results/')

# Save model for reuse
model.save('models/my_heartmap_model')

# Load and reuse saved model
loaded_model = HeartMapModel.load('models/my_heartmap_model')
new_results = loaded_model.predict(new_data)
```

### 3. REST API

```bash
# Start API server
python scripts/run_api_server.py
# Server available at http://localhost:8000
# API docs at http://localhost:8000/docs

# Use the API
curl -X POST "http://localhost:8000/analyze" \
     -H "Content-Type: multipart/form-data" \
     -F "file=@data/raw/heart_data.h5ad" \
     -F "analysis_type=comprehensive"

# Check available models
curl http://localhost:8000/models

# Update configuration
curl -X POST "http://localhost:8000/config" \
     -H "Content-Type: application/json" \
     -d '{"data": {"max_cells_subset": 30000}}'
```

### 4. Web Interface (Gradio)

```bash
# Start web interface
python app.py
# Access at http://localhost:7860

# Features:
# - Upload .h5ad files
# - Select analysis type
# - Configure memory settings
# - Download results
```

### 5. Jupyter Notebooks

```bash
# Install Jupyter
pip install jupyter

# Start notebook server
jupyter lab

# Open example notebooks:
# - notebooks/01_basic_analysis.ipynb
# - notebooks/02_advanced_communication.ipynb
# - notebooks/03_multi_chamber_analysis.ipynb
# - notebooks/04_comprehensive_analysis.ipynb
```

## Configuration

HeartMAP uses YAML configuration files for easy customization:

```yaml
# config.yaml or my_config.yaml
data:
  min_genes: 200
  min_cells: 3
  max_cells_subset: 50000        # Adjust based on your RAM
  max_genes_subset: 5000         # Reduce for faster analysis
  target_sum: 10000.0
  n_top_genes: 2000
  random_seed: 42
  test_mode: false               # Set true for quick testing

analysis:
  n_components_pca: 50
  n_neighbors: 10
  n_pcs: 40
  resolution: 0.5
  n_marker_genes: 25
  use_leiden: true
  use_liana: true                # Cell-cell communication

model:
  model_type: "comprehensive"
  save_intermediate: true
  use_gpu: false                 # Set true if GPU available
  batch_size: null
  max_memory_gb: null            # Auto-detect memory

paths:
  data_dir: "data"
  raw_data_dir: "data/raw"
  processed_data_dir: "data/processed"
  results_dir: "results"
  figures_dir: "figures"
  models_dir: "models"
```

### Memory Optimization Guidelines

| System RAM | max_cells_subset | max_genes_subset | Use Case |
|------------|------------------|------------------|----------|
| 8GB        | 10,000          | 2,000           | Laptop/Desktop |
| 16GB       | 30,000          | 4,000           | Workstation |
| 32GB       | 50,000          | 5,000           | Server |
| 64GB+      | 100,000+        | 10,000+         | HPC/Cloud |

## 🔬 Analysis Components

### 1. Basic Pipeline
- **Data Preprocessing**: Quality control, normalization, scaling
- **Cell Type Annotation**: Clustering and cell type identification  
- **Basic Visualization**: UMAP, t-SNE, cluster plots
- **Quality Metrics**: Cell and gene filtering statistics

### 2. Advanced Communication Analysis
- **Cell-Cell Communication**: Ligand-receptor interaction analysis
- **Communication Hubs**: Identification of key signaling cells
- **Pathway Enrichment**: Cardiac development and disease pathways
- **Network Analysis**: Communication network topology

### 3. Multi-Chamber Atlas
- **Chamber-Specific Analysis**: RA, RV, LA, LV specific patterns
- **Marker Identification**: Chamber-specific biomarkers
- **Cross-Chamber Correlations**: Inter-chamber relationship analysis
- **Comparative Analysis**: Chamber-to-chamber differences

### 4. Comprehensive Pipeline
- **Integrated Analysis**: All components combined
- **Advanced Visualizations**: Multi-panel figures and dashboards
- **Comprehensive Reports**: Automated result summaries
- **Model Persistence**: Save complete analysis state

## Deployment Guide

### Local Development Setup

#### Prerequisites
- Python 3.8+ (recommended: Python 3.10)
- Git
- Docker (optional, for containerized deployment)
- 8GB+ RAM (16GB+ recommended for larger datasets)

#### Quick Setup

```bash
# Clone and setup environment
git clone https://github.com/Tumo505/HeartMap.git
cd HeartMap

# Create virtual environment
python -m venv heartmap_env
source heartmap_env/bin/activate  # On Windows: heartmap_env\Scripts\activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements-dev.txt
pip install -e .[all]

# Configure the platform
cp config.yaml my_config.yaml
# Edit my_config.yaml to match your system resources

# Test installation
python scripts/validate.py
python scripts/run_examples.py
```

### Docker Deployment

```bash
# Single service
docker build -t heartmap .
docker run -p 8000:8000 -v $(pwd)/data:/app/data heartmap

# Multi-service with docker-compose
docker-compose up
# Services:
# - API server: http://localhost:8000  
# - Gradio interface: http://localhost:7860
# - Worker processes for batch analysis
```

### Hugging Face Spaces Deployment

```bash
# Prepare deployment files
./scripts/deploy_huggingface.sh

# Upload to your Hugging Face Space:
# 1. Create new Space at https://huggingface.co/new-space
# 2. Choose Gradio SDK
# 3. Upload generated files:
#    - app.py (Gradio interface)
#    - requirements.txt (Dependencies)
#    - src/ (Source code)
#    - config.yaml (Configuration)
```

### Cloud Platforms

#### AWS Deployment
```bash
# ECS deployment
docker build -t heartmap .
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin $ECR_URI
docker tag heartmap:latest $ECR_URI/heartmap:latest
docker push $ECR_URI/heartmap:latest

# Lambda deployment for serverless
sam build
sam deploy --guided
```

#### Google Cloud Platform
```bash
# Cloud Run deployment
gcloud builds submit --tag gcr.io/$PROJECT_ID/heartmap
gcloud run deploy --image gcr.io/$PROJECT_ID/heartmap --platform managed
```

#### Azure
```bash
# Container Instances
az container create --resource-group myResourceGroup \
    --name heartmap --image myregistry.azurecr.io/heartmap:latest
```

## Scientific Results

### Chamber Distribution
- **RA (Right Atrium):** 28.4% of cells
- **LV (Left Ventricle):** 27.0% of cells  
- **LA (Left Atrium):** 26.4% of cells
- **RV (Right Ventricle):** 18.2% of cells

### Chamber-Specific Markers
- **RA:** NPPA, MIR100HG, MYL7, MYL4, PDE4D
- **RV:** NEAT1, MYH7, FHL2, C15orf41, PCDH7
- **LA:** NPPA, ELN, MYL7, EBF2, RORA
- **LV:** CD36, LINC00486, FHL2, RP11-532N4.2, MYH7

### Cross-Chamber Correlations
- **RV vs LV:** r = 0.985 (highest correlation)
- **RA vs LA:** r = 0.960
- **LA vs LV:** r = 0.870 (lowest correlation)

## Testing & Validation

### Run Tests

```bash
# Full test suite
python tests/test_heartmap.py

# Validation suite
python scripts/validate.py

# Example analysis with mock data
python scripts/demo.py
```

### Performance Benchmarks

| Dataset Size | Memory Usage | Processing Time | Output |
|-------------|--------------|-----------------|--------|
| 10K cells   | 2GB RAM     | 5 minutes      | Complete analysis |
| 50K cells   | 8GB RAM     | 15 minutes     | Complete analysis |
| 100K cells  | 16GB RAM    | 30 minutes     | Complete analysis |

## Use Cases

### Research Applications
- **Interactive Analysis**: Jupyter notebooks for exploration
- **Batch Processing**: Command-line analysis of multiple datasets
- **Pipeline Integration**: Python API for custom workflows
- **Collaborative Research**: Shared configurations and models

### Production Applications
- **Web Services**: REST API for applications
- **Public Access**: Hugging Face Spaces for community use
- **Microservices**: Containerized deployment in cloud
- **High-Throughput**: Scalable analysis for large cohorts

### Educational Use
- **Teaching Platform**: Web interface for students
- **Reproducible Science**: Containerized environments
- **Method Comparison**: Multiple analysis approaches
- **Best Practices**: Clean, documented codebase

## Data Integrity & Reproducibility

### SHA-256 Checksums
- **Purpose**: Ensure data file integrity during storage/transfer
- **Implementation**: Automatic verification before analysis
- **Usage**: `python utils/sha256_checksum.py verify data/raw data/raw/checksums.txt`

### Fixed Random Seeds
- **Purpose**: Ensure reproducible results across runs
- **Implementation**: Fixed seeds in all stochastic processes
- **Scope**: Random sampling, clustering, mock data generation

### Examples of Reproducible Components

1. **Random Sampling**:
   ```python
   np.random.seed(42)
   cell_indices = np.random.choice(adata.n_obs, size=50000, replace=False)
   ```

2. **Clustering**:
   ```python
   kmeans = KMeans(n_clusters=n_clusters, random_state=42)
   ```

3. **LIANA Analysis**:
   ```python
   li.mt.rank_aggregate.by_sample(
       adata, groupby=cell_type_col, resource_name='consensus',
       n_perms=100, seed=42, verbose=True
   )
   ```

## 🔧 Development & Contributing

### Development Setup

```bash
# Development setup
git clone https://github.com/Tumo505/HeartMap.git
cd HeartMap
./scripts/setup.sh
source heartmap_env/bin/activate

# Install development dependencies
pip install -e .[dev]

# Run tests before committing
python tests/test_heartmap.py
python scripts/validate.py

# Code quality checks
black src/ tests/ scripts/  # Code formatting
flake8 src/ --max-line-length=100  # Linting
mypy src/heartmap --ignore-missing-imports  # Type checking
```

### Adding New Features

1. **New Models**: Inherit from `BaseModel` in `src/heartmap/models/`
2. **New Pipelines**: Inherit from `BasePipeline` in `src/heartmap/pipelines/`
3. **API Endpoints**: Add to `src/heartmap/api/rest.py`
4. **Configuration**: Extend dataclasses in `src/heartmap/config/`

### Testing Guidelines

```python
# Test structure
tests/
├── test_config.py       # Configuration management
├── test_data.py         # Data processing
├── test_models.py       # Analysis models
├── test_pipelines.py    # Analysis pipelines
└── test_api.py          # API interfaces
```

## Troubleshooting

### Common Issues

#### Memory Errors
```bash
# Reduce dataset size in config
data:
  max_cells_subset: 10000
  max_genes_subset: 2000
```

#### Import Errors
```bash
# Reinstall with all dependencies
pip install -e .[all]

# Check Python path
python -c "import heartmap; print(heartmap.__file__)"
```

#### Data Loading Issues
```bash
# Verify data format
python -c "import scanpy as sc; print(sc.read_h5ad('data/raw/your_file.h5ad'))"

# Check file permissions
ls -la data/raw/
```

#### Performance Issues
```bash
# Enable test mode for quick validation
test_mode: true  # in config.yaml

# Use GPU acceleration (if available)
model:
  use_gpu: true
```

### Getting Help

1. **Validation**: `python scripts/validate.py`
2. **Logs**: Check logs in `results/` directory
3. **Test Mode**: Set `test_mode: true` in configuration
4. **Mock Data**: `python scripts/demo.py` for testing
5. **Documentation**: Comprehensive docstrings in source code

## Requirements

### System Requirements
- **Python**: 3.8+ (recommended: 3.10)
- **Memory**: 8GB+ recommended (configurable)
- **Storage**: 5GB+ for data and results
- **OS**: Linux, macOS, Windows

### Dependencies
- **Core**: scanpy, pandas, numpy, scipy, scikit-learn
- **Visualization**: matplotlib, seaborn, plotly
- **Communication**: liana, cellphonedb (optional)
- **API**: fastapi, uvicorn (optional)
- **Web**: gradio (optional)

## Clinical Applications

- **Personalized Medicine**: Chamber-specific treatment strategies
- **Drug Development**: Chamber-specific therapeutic targets
- **Disease Understanding**: Chamber-specific disease mechanisms
- **Biomarker Discovery**: Chamber and communication-specific markers

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Install development dependencies: `pip install -e .[dev]`
4. Run tests: `python tests/test_heartmap.py`
5. Submit a pull request

## 📄 License

Ssee the [LICENSE](LICENSE) file for details.

## Citation

```bibtex
@software{heartmap2025,
  title={HeartMAP: A Multi-Chamber Spatial Framework for Cardiac Cell-Cell Communication},
  author={Kgabeng, Tumo and Wang, Lulu and Ngwangwa, Harry and Pandelani, Thanyani},
  year={2025},
  url={https://github.com/Tumo505/HeartMap},
  version={1.0.0},
  doi={10.5281/zenodo.16745118}
}
```

## Support

- **Issues**: [GitHub Issues](https://github.com/Tumo505/HeartMap/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Tumo505/HeartMap/discussions)
- **Email**: 28346416@mylife.unisa.ac.za

## Acknowledgments

- Department of Mechanical, Bioresources and Biomedical Engineering, University of South Africa
- Department of Engineering, Reykjavik University
- Single Cell Portal (SCP498) for providing the heart dataset
- The open-source scientific Python community

---

**HeartMAP is now production-ready and available for research, deployment, and collaboration!**

Whether you're a researcher exploring cardiac biology, a developer building applications, or an educator teaching single-cell analysis, HeartMAP provides the tools and flexibility you need.
