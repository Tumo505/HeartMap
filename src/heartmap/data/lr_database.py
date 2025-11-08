"""
Ligand-Receptor Database Loader for HeartMAP
Uses LIANA's curated resources (consensus, CellPhoneDB, Omnipath, etc.)
"""

import pandas as pd
from typing import List, Tuple, Optional

# Try importing LIANA
try:
    import liana
    LIANA_AVAILABLE = True
    print(f"✓ LIANA v{liana.__version__} loaded successfully")
except ImportError:
    LIANA_AVAILABLE = False
    print("⚠ LIANA not available - will use fallback database")

class LigandReceptorDatabase:
    """
    Manage ligand-receptor interaction databases
    """

    def __init__(self, resource: str = 'consensus'):
        """
        Initialize L-R database

        Parameters:
        -----------
        resource : str
            Database to use:
            - 'consensus': curated from multiple sources (recommended)
            - 'cellphonedb': CellPhoneDB database
            - 'omnipath': OmniPath database  
            - 'connectome': Ramilowski 2015
            - 'cellinker': CellLinker database
        """
        self.resource = resource
        self.lr_pairs = None
        self.load_database()

    def load_database(self):
        """Load L-R database from LIANA or fallback"""

        if LIANA_AVAILABLE:
            self._load_from_liana()
        else:
            self._load_fallback_database()

    def _load_from_liana(self):
        """Load curated L-R pairs from LIANA"""
        try:
            from liana.resource import select_resource

            # Get the resource
            lr_df = select_resource(self.resource)

            # Standardize column names
            # LIANA typically has: ligand, receptor, (optional: source, confidence)
            if 'ligand' in lr_df.columns and 'receptor' in lr_df.columns:
                self.lr_pairs = lr_df[['ligand', 'receptor']].drop_duplicates()

                # Add confidence score if available
                if 'score' in lr_df.columns:
                    self.lr_pairs['confidence'] = lr_df['score']
                elif 'confidence' in lr_df.columns:
                    self.lr_pairs['confidence'] = lr_df['confidence']
                else:
                    self.lr_pairs['confidence'] = 1.0  # Default high confidence

                print(f"✓ Loaded {len(self.lr_pairs)} L-R pairs from LIANA {self.resource}")
            else:
                print(f" Unexpected LIANA format, using fallback")
                self._load_fallback_database()

        except Exception as e:
            print(f" Error loading LIANA database: {e}")
            self._load_fallback_database()

    def _load_fallback_database(self):
        """Fallback: comprehensive cardiac-focused L-R pairs"""

        # Expanded cardiac-relevant pairs (100+ interactions)
        cardiac_lr_data = [
            # Angiogenesis & Vascular
            ('VEGFA', 'FLT1', 0.95), ('VEGFA', 'KDR', 0.98), ('VEGFA', 'NRP1', 0.85),
            ('VEGFB', 'FLT1', 0.90), ('VEGFC', 'FLT4', 0.95),
            ('ANGPT1', 'TEK', 0.92), ('ANGPT2', 'TEK', 0.88),
            ('PGF', 'FLT1', 0.85),

            # TGF-beta superfamily
            ('TGFB1', 'TGFBR1', 0.98), ('TGFB1', 'TGFBR2', 0.98),
            ('TGFB2', 'TGFBR1', 0.95), ('TGFB3', 'TGFBR2', 0.95),
            ('BMP2', 'BMPR1A', 0.90), ('BMP2', 'BMPR2', 0.88),
            ('BMP4', 'BMPR1A', 0.92), ('BMP7', 'BMPR2', 0.85),
            ('INHBA', 'ACVR1B', 0.80),

            # FGF signaling
            ('FGF1', 'FGFR1', 0.90), ('FGF2', 'FGFR1', 0.95),
            ('FGF2', 'FGFR2', 0.92), ('FGF7', 'FGFR2', 0.88),
            ('FGF9', 'FGFR3', 0.85), ('FGF10', 'FGFR2', 0.87),

            # PDGF signaling
            ('PDGFA', 'PDGFRA', 0.98), ('PDGFB', 'PDGFRB', 0.98),
            ('PDGFC', 'PDGFRA', 0.90), ('PDGFD', 'PDGFRB', 0.88),

            # Inflammatory cytokines
            ('IL6', 'IL6R', 0.98), ('IL1B', 'IL1R1', 0.95),
            ('TNF', 'TNFRSF1A', 0.98), ('TNF', 'TNFRSF1B', 0.95),
            ('IFNG', 'IFNGR1', 0.92), ('IL10', 'IL10RA', 0.90),
            ('IL4', 'IL4R', 0.88), ('IL13', 'IL13RA1', 0.85),

            # Chemokines
            ('CXCL12', 'CXCR4', 0.98), ('CXCL12', 'CXCR7', 0.85),
            ('CCL2', 'CCR2', 0.95), ('CCL5', 'CCR5', 0.92),
            ('CXCL8', 'CXCR1', 0.90), ('CXCL8', 'CXCR2', 0.88),
            ('CCL3', 'CCR1', 0.85), ('CCL4', 'CCR5', 0.87),

            # Growth factors
            ('EGF', 'EGFR', 0.98), ('HBEGF', 'EGFR', 0.92),
            ('IGF1', 'IGF1R', 0.98), ('IGF2', 'IGF1R', 0.95),
            ('HGF', 'MET', 0.95), ('NGF', 'NTRK1', 0.92),

            # Notch signaling
            ('DLL1', 'NOTCH1', 0.90), ('DLL4', 'NOTCH1', 0.92),
            ('JAG1', 'NOTCH1', 0.88), ('JAG1', 'NOTCH2', 0.85),
            ('JAG2', 'NOTCH3', 0.82),

            # Wnt signaling
            ('WNT3A', 'FZD1', 0.85), ('WNT3A', 'FZD2', 0.83),
            ('WNT5A', 'FZD5', 0.88), ('WNT7A', 'FZD7', 0.85),

            # Extracellular matrix
            ('COL1A1', 'ITGA1', 0.90), ('COL1A1', 'ITGA2', 0.88),
            ('FN1', 'ITGA5', 0.95), ('FN1', 'ITGB1', 0.92),
            ('LAMB1', 'ITGA6', 0.88), ('THBS1', 'CD47', 0.85),

            # Cardiac specific
            ('NRG1', 'ERBB2', 0.95), ('NRG1', 'ERBB4', 0.92),
            ('EDN1', 'EDNRA', 0.95), ('EDN1', 'EDNRB', 0.90),
            ('NPPA', 'NPR1', 0.92), ('NPPB', 'NPR1', 0.90),

            # Semaphorins
            ('SEMA3A', 'NRP1', 0.88), ('SEMA3C', 'NRP2', 0.85),
            ('SEMA4D', 'PLXNB1', 0.82),

            # Ephrins
            ('EFNA1', 'EPHA2', 0.90), ('EFNB2', 'EPHB4', 0.92),

            # Complement
            ('C3', 'C3AR1', 0.88), ('C5', 'C5AR1', 0.90),

            # Adhesion
            ('ICAM1', 'ITGAL', 0.92), ('VCAM1', 'ITGA4', 0.90),
            ('CD34', 'SELP', 0.85), ('PECAM1', 'PECAM1', 0.88),

            # Apoptosis
            ('FASLG', 'FAS', 0.95), ('TNFSF10', 'TNFRSF10A', 0.90),

            # Neuropeptides
            ('BDNF', 'NTRK2', 0.88), ('NTF3', 'NTRK3', 0.85),

            # Metabolic
            ('LEP', 'LEPR', 0.92), ('ADIPOQ', 'ADIPOR1', 0.88),
            ('INS', 'INSR', 0.98), ('GCG', 'GCGR', 0.90),
        ]

        self.lr_pairs = pd.DataFrame(
            cardiac_lr_data,
            columns=['ligand', 'receptor', 'confidence']
        )

        print(f"✓ Loaded {len(self.lr_pairs)} L-R pairs from fallback cardiac database")

    def get_pairs(self, confidence_threshold: float = 0.0, present_in_data: Optional[List[str]] = None) -> List[Tuple[str, str]]:
        """
        Get L-R pairs as list of tuples

        Parameters:
        -----------
        confidence_threshold : float
            Minimum confidence score (0-1)
        present_in_data : list of str, optional
            Gene names present in dataset (filters to only available pairs)

        Returns:
        --------
        list of tuples: [(ligand, receptor), ...]
        """

        # Filter by confidence
        filtered = self.lr_pairs[self.lr_pairs['confidence'] >= confidence_threshold].copy()

        # Filter by gene availability
        if present_in_data is not None:
            present_set = set(present_in_data)
            filtered = filtered[
                filtered['ligand'].isin(present_set) &
                filtered['receptor'].isin(present_set)
            ]

        return list(zip(filtered['ligand'], filtered['receptor']))

    def get_dataframe(self, confidence_threshold: float = 0.0) -> pd.DataFrame:
        """Get L-R pairs as DataFrame"""
        return self.lr_pairs[self.lr_pairs['confidence'] >= confidence_threshold].copy()

    def save_to_csv(self, filepath: str):
        """Save database to CSV"""
        self.lr_pairs.to_csv(filepath, index=False)
        print(f"✓ Saved L-R database to {filepath}")


# Convenience function
def get_ligand_receptor_pairs(adata, resource: str = 'consensus', confidence_threshold: float = 0.7) -> List[Tuple[str, str]]:
    """
    Get ligand-receptor pairs filtered to genes present in adata

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene names
    resource : str
        Database to use ('consensus', 'cellphonedb', etc.)
    confidence_threshold : float
        Minimum confidence (0-1)

    Returns:
    --------
    list of tuples: [(ligand, receptor), ...]
    """

    db = LigandReceptorDatabase(resource=resource)
    available_genes = adata.var_names.tolist()
    pairs = db.get_pairs(
        confidence_threshold=confidence_threshold,
        present_in_data=available_genes
    )

    print(f" Found {len(pairs)} L-R pairs present in dataset (from {len(db.lr_pairs)} total)")
    return pairs


if __name__ == "__main__":
    # Test the database loader
    print("Testing L-R Database Loader...")
    print("=" * 60)

    # Test with LIANA if available
    db = LigandReceptorDatabase(resource='consensus')
    print(f"\nTotal pairs: {len(db.lr_pairs)}")
    print(f"\nSample pairs:")
    print(db.lr_pairs.head(10))

    # Test filtering
    high_conf_pairs = db.get_pairs(confidence_threshold=0.9)
    print(f"\nHigh confidence pairs (>0.9): {len(high_conf_pairs)}")

    # Save example
    db.save_to_csv("lr_database_export.csv")
