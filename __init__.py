# -*- coding: utf-8 -*-
"""ShapG2P: biomarker pathway enrichment with PPI network topology and SHAP.

Usage:
    from shapg2p import score_pathways

    scores = score_pathways(['TP53', 'ATM', 'APOE', ...])   # gene list
    scores = score_pathways('biomarkers.csv')               # or a file path
"""
from .core import score_pathways

__version__ = '0.1.1'
__all__ = ['score_pathways']
