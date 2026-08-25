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

# 1. pass a list of gene symbols
scores = score_pathways(['TP53', 'ATM', 'APOE', 'SOD1', 'CDKN2A'])

# 2. pass a file path (CSV/TSV/TXT with a common gene column, e.g. "gene symbol")
scores = score_pathways('my_biomarkers.csv')

# 3. pass a comma/space/newline-delimited string
scores = score_pathways('TP53, ATM, APOE')

# Output: dict sorted by SHAP score descending
print(scores)
# {'p53 signaling pathway': 0.42, 'Alzheimer disease': 0.31, ...}
```

Each call takes about 1–3 minutes (one XGBoost training + full-genome SHAP computation).

## Output

Returns `dict[str, float]`: pathway name → mean |SHAP| score, sorted descending. Only pathways with score > 0 are returned (irrelevant pathways are omitted). Duplicate pathway names across databases take the maximum score.

## Publish to PyPI (for distribution)

### 1. Register an account and create an API token

- Register: https://pypi.org/account/register/
- Create a token: https://pypi.org/manage/account/token/ (Scope: "Entire account")
- The token looks like `pypi-AgEIcHlwaS5vcmcC...`, shown only once — save it.

### 2. Build

```bash
pip install --upgrade build twine
cd shapg2p_pkg
python -m build
```

This produces `dist/shapg2p-0.1.1.tar.gz` and `dist/shapg2p-0.1.1-py3-none-any.whl`.

### 3. Upload to TestPyPI first (optional but recommended)

```bash
python -m twine upload --repository testpypi dist/*
# Username: __token__    Password: your API token

# Verify installation
pip install --index-url https://test.pypi.org/simple/ shapg2p
```

### 4. Upload to PyPI

```bash
python -m twine upload dist/*
```

Enter `__token__` as the username and your API token as the password. After a successful upload:

```bash
pip install shapg2p
```

### 5. Release a new version

Bump `version` in `pyproject.toml` (e.g. `0.1.1`), then `python -m build` and `twine upload dist/*` again. PyPI does not allow re-uploading the same version number.

## Dependencies

numpy / pandas / scipy / scikit-learn / xgboost (installed automatically with `pip install shapg2p`).
