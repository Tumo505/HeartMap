"""
HeartMAP: Heart Multi-chamber Analysis Platform

A comprehensive framework for analyzing cell-cell communication across all four chambers
of the human heart using single-cell RNA sequencing data.
"""

__version__ = "2.0.0"
__author__ = "Tumo Kgabeng, Lulu Wang, Harry Ngwangwa, Thanyani Pandelani"
__email__ = "28346416@mylife.unisa.ac.za"

# from .models import HeartMapModel
from .pipelines import BasicPipeline, AdvancedCommunicationPipeline, MultiChamberPipeline
from .config import Config

__all__ = [
    "BasicPipeline",
    "AdvancedCommunicationPipeline",
    "MultiChamberPipeline",
    "Config"
]
