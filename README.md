# ShapG2P

ShapG2P: a strategy for biomarker pathway enrichment with PPI network topology and SHAP analysis.

Given a list of gene symbols (from a file or passed directly as a list), it returns a `{pathway: SHAP score}` dictionary sorted by score in descending order.

## Method

- Each gene is featurized by its PPI network distance to the genes of each pathway (`exp(-d/2)` similarity);
- An XGBoost classifier distinguishes biomarkers from background genes (global 1:1 undersampling, SEED=42);
- The mean absolute SHAP contribution (`pred_contribs`) of each pathway feature is used as the pathway SHAP score.

Bundled data: STRING human PPI network (17,613 genes), KEGG / Hallmark / WikiPathway pathways.

## Install

```bash
pip install shapg2p
```

Or for local development:

```bash
cd shapg2p_pkg
pip install -e .
```

## Usage

```python
from shapg2p import score_pathways

genes = ['TP53', 'ATM', 'APOE', 'SOD1', 'CDKN2A', 'BRCA1', 'BRCA2', 'PTEN',
         'RB1', 'MYC', 'EGFR', 'KRAS', 'PIK3CA', 'AKT1', 'MTOR', 'BCL2',
         'IL6', 'TNF', 'IL1B', 'IL10', 'CXCL8', 'TGFB1', 'VEGFA', 'MMP9',
         'COL1A1', 'SPP1', 'CD44', 'STAT3', 'NFKB1', 'HIF1A', 'SIRT1',
         'FOXO3', 'NFE2L2', 'KEAP1', 'TERT', 'LMNA', 'WRN', 'PINK1',
         'PARK7', 'SNCA']

# from a gene symbol list
scores = score_pathways(genes)

# or from a file (CSV/TSV/TXT; common gene columns such as "gene symbol"
# are auto-detected, STRING_ID is also supported)
# scores = score_pathways('my_biomarkers.csv')

print(scores)
# {'p53 signaling pathway': 0.42, 'Alzheimer disease': 0.31, ...}
```

Each call takes about 1–3 minutes (one XGBoost training + full-genome SHAP computation).

## Output

Returns `dict[str, float]`: pathway name → mean |SHAP| score, sorted descending. Only pathways with score > 0 are returned (irrelevant pathways are omitted). Duplicate pathway names across databases take the maximum score.

## Dependencies

numpy / pandas / scipy / scikit-learn / xgboost (installed automatically with `pip install shapg2p`).
