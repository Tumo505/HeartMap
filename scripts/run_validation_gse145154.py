#!/usr/bin/env python3
from pathlib import Path
import sys

# Add src to path for development mode
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / 'src'))

from heartmap.api import CLIInterface

def main():
	data_path = project_root / 'data' / 'validation' / 'gse145154_merged.h5ad'
	output_dir = project_root / 'results' / 'validation_gse145154'
	config_path = project_root / 'config.yaml'

	cli = CLIInterface()
	cli.run_analysis(
		data_path=str(data_path),
		analysis_type='basic',
		output_dir=str(output_dir),
		config_path=None
	)

if __name__ == '__main__':
	main()
