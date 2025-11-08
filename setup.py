"""
Setup script for HeartMAP package
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README with UTF-8 encoding
readme_path = Path(__file__).parent / "README.md"
long_description = ""
if readme_path.exists():
    with open(readme_path, encoding='utf-8') as f:
        long_description = f.read()

# Read version with UTF-8 encoding
version_path = Path(__file__).parent / "src" / "heartmap" / "__init__.py"
version = "1.1.0"
if version_path.exists():
    with open(version_path, encoding='utf-8') as f:
        for line in f:
            if line.startswith("__version__"):
                version = line.split("=")[1].strip().strip('"').strip("'")
                break

setup(
    name="heartmap",
    version=version,
    author="Tumo Kgabeng, Lulu Wang, Harry Ngwangwa, Thanyani Pandelani",
    author_email="28346416@mylife.unisa.ac.za",
    description="Heart Multi-chamber Analysis Platform for single-cell RNA-seq",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Tumo505/HeartMap",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "scanpy>=1.9.0",
        "pandas>=1.5.0",
        "numpy>=1.21.0",
        "scipy>=1.9.0",
        "scikit-learn>=1.1.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "anndata>=0.8.0",
        "plotly>=5.0.0",
        "networkx>=2.8.0",
        "python-igraph>=0.10.0",
        "leidenalg>=0.8.0",
        "tqdm",
        "statsmodels",
        "pyyaml>=6.0",
    ],
    extras_require={
        "communication": [
            "liana>=0.1.0",
            "cellphonedb>=3.0.0",
            "omnipath>=1.0.0",
        ],
        "api": [
            "fastapi>=0.100.0",
            "uvicorn>=0.23.0",
            "pydantic>=2.0.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.0.0",
            "types-PyYAML>=6.0.0",
            "jupyter>=1.0.0",
            "notebook>=6.0.0",
        ],
        "all": [
            "liana>=0.1.0",
            "cellphonedb>=3.0.0",
            "omnipath>=1.0.0",
            "fastapi>=0.100.0",
            "uvicorn>=0.23.0",
            "pydantic>=2.0.0",
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.0.0",
            "types-PyYAML>=6.0.0",
            "jupyter>=1.0.0",
            "notebook>=6.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "heartmap=heartmap.api:run_cli",
            "heartmap-api=heartmap.api:create_api",
        ],
    },
    include_package_data=True,
    package_data={
        "heartmap": ["data/*.json", "config/*.yaml"],
    },
    keywords=[
        "single-cell",
        "RNA-seq",
        "heart",
        "cell-communication",
        "bioinformatics",
        "spatial-transcriptomics",
    ],
    project_urls={
        "Documentation": "https://github.com/Tumo505/HeartMap/wiki",
        "Source": "https://github.com/Tumo505/HeartMap",
        "Tracker": "https://github.com/Tumo505/HeartMap/issues",
    },
)
