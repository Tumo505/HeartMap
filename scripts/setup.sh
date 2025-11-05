#!/bin/bash

# HeartMAP Quick Start Setup Script
# This script automates the local development setup

set -e  # Exit on any error

echo " HeartMAP Quick Start Setup"
echo "=============================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

# Check if Python 3.8+ is available
check_python() {
    print_step "Checking Python version..."
    
    if command -v python3 &> /dev/null; then
        PYTHON_VERSION=$(python3 --version | cut -d' ' -f2)
        PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d'.' -f1)
        PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d'.' -f2)
        
        if [ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -ge 8 ]; then
            print_status "Python $PYTHON_VERSION found ✓"
            PYTHON_CMD="python3"
        else
            print_error "Python 3.8+ required, found $PYTHON_VERSION"
            exit 1
        fi
    elif command -v python &> /dev/null; then
        PYTHON_VERSION=$(python --version | cut -d' ' -f2)
        PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d'.' -f1)
        PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d'.' -f2)
        
        if [ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -ge 8 ]; then
            print_status "Python $PYTHON_VERSION found ✓"
            PYTHON_CMD="python"
        else
            print_error "Python 3.8+ required, found $PYTHON_VERSION"
            exit 1
        fi
    else
        print_error "Python not found. Please install Python 3.8+ first."
        exit 1
    fi
}

# Create virtual environment
setup_venv() {
    print_step "Setting up virtual environment..."
    
    if [ ! -d "heartmap_env" ]; then
        $PYTHON_CMD -m venv heartmap_env
        print_status "Virtual environment created ✓"
    else
        print_warning "Virtual environment already exists"
    fi
    
    # Activate virtual environment
    if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
        # Windows
        source heartmap_env/Scripts/activate
    else
        # Unix-like (Linux, macOS)
        source heartmap_env/bin/activate
    fi
    
    print_status "Virtual environment activated ✓"
    
    # Upgrade pip
    pip install --upgrade pip
    print_status "Pip upgraded ✓"
}

# Install dependencies
install_dependencies() {
    print_step "Installing dependencies..."
    
    # Install requirements
    if [ -f "requirements-dev.txt" ]; then
        pip install -r requirements-dev.txt
        print_status "Development requirements installed ✓"
    else
        print_warning "requirements-dev.txt not found, installing basic requirements"
        pip install scanpy pandas numpy scipy scikit-learn matplotlib seaborn anndata plotly networkx tqdm statsmodels pyyaml
    fi
    
    # Install package in development mode
    if [ -f "setup.py" ]; then
        pip install -e .[all]
        print_status "HeartMAP package installed in development mode ✓"
    else
        print_warning "setup.py not found, skipping package installation"
    fi
}

# Run migration
run_migration() {
    print_step "Running migration script..."
    
    if [ -f "scripts/migrate.py" ]; then
        $PYTHON_CMD scripts/migrate.py
        print_status "Migration completed ✓"
    else
        print_warning "Migration script not found, skipping"
    fi
}

# Create sample configuration
create_config() {
    print_step "Creating sample configuration..."
    
    if [ ! -f "my_config.yaml" ]; then
        cp config.yaml my_config.yaml 2>/dev/null || cat > my_config.yaml << EOF
# HeartMAP Configuration
data:
  min_genes: 200
  min_cells: 3
  max_cells_subset: 50000    # Adjust based on your RAM
  max_genes_subset: 5000     # Reduce for faster analysis
  target_sum: 10000.0
  n_top_genes: 2000
  random_seed: 42
  test_mode: false

analysis:
  n_components_pca: 50
  n_neighbors: 10
  n_pcs: 40
  resolution: 0.5
  n_marker_genes: 25
  use_leiden: true
  use_liana: true

model:
  model_type: "comprehensive"
  save_intermediate: true
  use_gpu: false
  batch_size: null
  max_memory_gb: null

paths:
  data_dir: "data"
  raw_data_dir: "data/raw"
  processed_data_dir: "data/processed"
  results_dir: "results"
  figures_dir: "figures"
  models_dir: "models"
EOF
        print_status "Sample configuration created: my_config.yaml ✓"
    else
        print_warning "Configuration file already exists: my_config.yaml"
    fi
}

# Run tests
run_tests() {
    print_step "Running tests..."
    
    if [ -f "tests/test_heartmap.py" ]; then
        $PYTHON_CMD tests/test_heartmap.py
        print_status "Tests completed ✓"
    else
        print_warning "Test file not found, skipping tests"
    fi
}

# Create data directories
create_directories() {
    print_step "Creating directory structure..."
    
    mkdir -p data/raw
    mkdir -p data/processed
    mkdir -p results
    mkdir -p figures
    mkdir -p models
    mkdir -p notebooks
    
    print_status "Directory structure created ✓"
}

# Download example data (optional)
download_example_data() {
    print_step "Setting up example data..."
    
    if [ ! -f "data/raw/example_heart_data.h5ad" ]; then
        print_warning "No example data available."
        print_warning "To use HeartMAP, place your .h5ad files in data/raw/"
        print_warning "Example: data/raw/your_heart_data.h5ad"
    fi
}

# Main setup function
main() {
    echo ""
    print_step "Starting HeartMAP setup..."
    echo ""
    
    check_python
    setup_venv
    install_dependencies
    create_directories
    run_migration
    create_config
    
    # Optional steps
    echo ""
    print_step "Running optional steps..."
    run_tests
    download_example_data
    
    echo ""
    echo " Setup completed successfully!"
    echo ""
    echo "Next steps:"
    echo "1. Activate environment: source heartmap_env/bin/activate  (Linux/Mac) or heartmap_env\\Scripts\\activate (Windows)"
    echo "2. Place your .h5ad data files in: data/raw/"
    echo "3. Run analysis: heartmap data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad --analysis-type comprehensive"
    echo "4. Or start API server: python scripts/run_api_server.py"
    echo "5. Or start Gradio interface: python app.py"
    echo ""
    echo "Configuration file: my_config.yaml (edit to customize)"
    echo "Documentation: DEPLOYMENT_GUIDE.md"
    echo ""
}

# Run setup
main "$@"
