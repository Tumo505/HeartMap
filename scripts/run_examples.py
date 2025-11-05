"""
Example script demonstrating HeartMAP usage
"""

import os
from pathlib import Path
import warnings

# Import HeartMAP components
from heartmap import Config, HeartMapModel
from heartmap.pipelines import ComprehensivePipeline, BasicPipeline
from heartmap.api import CLIInterface


def example_basic_analysis():
    """Example of basic HeartMAP analysis"""
    print("=== Basic HeartMAP Analysis Example ===")
    
    # Create configuration
    config = Config.default()
    
    # For demonstration, use smaller dataset
    config.data.max_cells_subset = 5000
    config.data.max_genes_subset = 2000
    config.data.test_mode = True
    
    # Update paths
    base_dir = Path(__file__).parent.parent
    config.update_paths(str(base_dir))
    config.create_directories()
    
    # Create basic pipeline
    pipeline = BasicPipeline(config)
    
    # Run analysis
    # data_path = "data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad"
    # results = pipeline.run(data_path, "results/basic")
    
    print("Basic analysis pipeline created successfully!")
    return pipeline


def example_comprehensive_analysis():
    """Example of comprehensive HeartMAP analysis"""
    print("=== Comprehensive HeartMAP Analysis Example ===")
    
    # Load configuration from file (if it exists)
    config_path = Path(__file__).parent.parent / "config.yaml"
    
    if config_path.exists():
        config = Config.from_yaml(str(config_path))
    else:
        config = Config.default()
    
    # Create comprehensive pipeline
    pipeline = ComprehensivePipeline(config)
    
    print("Comprehensive pipeline created successfully!")
    return pipeline


def example_model_training():
    """Example of training HeartMAP model"""
    print("=== HeartMAP Model Training Example ===")
    
    # Create configuration
    config = Config.default()
    config.data.test_mode = True
    
    # Create model
    model = HeartMapModel(config)
    
    print("HeartMAP model created successfully!")
    print("To train the model, use: model.fit(your_adata)")
    print("To make predictions, use: model.predict(new_adata)")
    
    return model


def example_cli_usage():
    """Example of CLI usage"""
    print("=== CLI Usage Example ===")
    
    # Create CLI interface
    cli = CLIInterface()
    
    print("CLI interface created!")
    print("To run analysis from command line:")
    print("heartmap data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad --analysis-type comprehensive --output-dir results")
    
    return cli


def example_config_management():
    """Example of configuration management"""
    print("=== Configuration Management Example ===")
    
    # Create default config
    config = Config.default()
    
    # Modify settings
    config.data.min_genes = 500
    config.data.min_cells = 5
    config.analysis.resolution = 0.8
    config.model.use_gpu = False
    
    # Save to file
    output_dir = Path(__file__).parent.parent / "examples"
    output_dir.mkdir(exist_ok=True)
    
    config.save_yaml(str(output_dir / "example_config.yaml"))
    config.save_json(str(output_dir / "example_config.json"))
    
    print("Configuration saved to:")
    print(f"- {output_dir / 'example_config.yaml'}")
    print(f"- {output_dir / 'example_config.json'}")
    
    # Load from file
    loaded_config = Config.from_yaml(str(output_dir / "example_config.yaml"))
    print("Configuration loaded successfully!")
    
    return config, loaded_config


def main():
    """Run all examples"""
    print("HeartMAP Examples")
    print("=" * 50)
    
    try:
        # Basic analysis
        basic_pipeline = example_basic_analysis()
        print()
        
        # Comprehensive analysis
        comprehensive_pipeline = example_comprehensive_analysis()
        print()
        
        # Model training
        model = example_model_training()
        print()
        
        # CLI usage
        cli = example_cli_usage()
        print()
        
        # Configuration management
        config, loaded_config = example_config_management()
        print()
        
        print(" All examples completed successfully!")
        print("\nNext steps:")
        print("1. Prepare your single-cell RNA-seq data in h5ad format")
        print("2. Run: heartmap your_data.h5ad --analysis-type comprehensive")
        print("3. Check results in the output directory")
        
    except Exception as e:
        print(f" Error running examples: {e}")
        print("Make sure all dependencies are installed: pip install -e .[all]")


if __name__ == "__main__":
    main()
