#!/usr/bin/env python3

import os
import gzip
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import anndata as ad


def find_sample_prefixes(raw_dir: Path) -> List[str]:
	"""Identify unique sample prefixes based on *_barcodes.tsv.gz files."""
	prefixes = []
	for barcodes_file in raw_dir.glob("*_barcodes.tsv.gz"):
		prefix = barcodes_file.name.replace("_barcodes.tsv.gz", "")
		features_file = raw_dir / f"{prefix}_features.tsv.gz"
		matrix_file = raw_dir / f"{prefix}_matrix.mtx.gz"
		if features_file.exists() and matrix_file.exists():
			prefixes.append(prefix)
	return sorted(prefixes)


def read_10x_triplet(raw_dir: Path, prefix: str) -> ad.AnnData:
	"""Read a 10x-style triplet with a prefix and return AnnData.

	Files expected:
	- {prefix}_matrix.mtx.gz
	- {prefix}_features.tsv.gz
	- {prefix}_barcodes.tsv.gz
	"""
	matrix_path = raw_dir / f"{prefix}_matrix.mtx.gz"
	features_path = raw_dir / f"{prefix}_features.tsv.gz"
	barcodes_path = raw_dir / f"{prefix}_barcodes.tsv.gz"

	# Read matrix (sparse)
	with gzip.open(matrix_path, "rb") as fh:
		X = mmread(fh)
	# 10x MTX is genes x cells; transpose to cells x genes and ensure CSR
	if sp.issparse(X):
		X = X.T.tocsr()
	else:
		X = np.asarray(X).T

	# Read features
	with gzip.open(features_path, "rt") as fh:
		features = pd.read_csv(
			fh,
			sep="\t",
			header=None,
			comment=None,
			dtype=str
		)
	# 10x features.tsv columns: id, name, type (sometimes 2 or 3 cols)
	if features.shape[1] == 1:
		# Fallback: single column of gene names
		gene_ids = features.iloc[:, 0].values
		gene_names = features.iloc[:, 0].values
	else:
		gene_ids = features.iloc[:, 0].values
		gene_names = features.iloc[:, 1].values

	var = pd.DataFrame(index=pd.Index(gene_names, name="gene"))
	var["feature_id"] = gene_ids

	# Read barcodes
	with gzip.open(barcodes_path, "rt") as fh:
		barcodes = pd.read_csv(fh, sep="\t", header=None, comment=None, dtype=str).iloc[:, 0].values

	obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
	obs["sample"] = prefix

	adata = ad.AnnData(X=X, obs=obs, var=var)
	adata.obs_names_make_unique()
	adata.var_names_make_unique()
	return adata


def merge_samples(adatas: List[ad.AnnData]) -> ad.AnnData:
	"""Concatenate AnnData objects along observations (cells)."""
	if not adatas:
		raise ValueError("No AnnData objects to merge.")
	adata_merged = ad.concat(adatas, axis=0, join="outer", label="batch", keys=[a.obs["sample"].iloc[0] for a in adatas])
	return adata_merged


def main():
	project_root = Path(__file__).parent.parent
	raw_dir = project_root / "data" / "validation" / "GSE145154_RAW"
	output_dir = project_root / "data" / "validation"
	output_dir.mkdir(parents=True, exist_ok=True)
	output_file = output_dir / "gse145154_merged.h5ad"

	if not raw_dir.exists():
		raise FileNotFoundError(f"Raw directory not found: {raw_dir}")

	prefixes = find_sample_prefixes(raw_dir)
	if not prefixes:
		raise RuntimeError(f"No 10x triplets found in {raw_dir}")

	print(f"Found {len(prefixes)} samples. Reading...")
	adatas: List[ad.AnnData] = []
	for i, prefix in enumerate(prefixes, start=1):
		print(f"[{i}/{len(prefixes)}] Reading {prefix}...")
		adata = read_10x_triplet(raw_dir, prefix)
		adatas.append(adata)

	print("Merging samples...")
	adata_merged = merge_samples(adatas)
	print(f"Merged shape: {adata_merged.shape}")

	print(f"Writing {output_file}...")
	adata_merged.write_h5ad(output_file)
	print("Done.")


if __name__ == "__main__":
	main()
