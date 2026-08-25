# -*- coding: utf-8 -*-
"""ShapG2P: biomarker pathway enrichment with PPI network topology and SHAP.

用法:
    from shapg2p import score_pathways

    scores = score_pathways(['TP53', 'ATM', 'APOE', ...])   # 基因列表
    scores = score_pathways('biomarkers.csv')               # 或文件路径
"""
from .core import score_pathways

__version__ = '0.1.1'
__all__ = ['score_pathways']
