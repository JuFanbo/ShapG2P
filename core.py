# -*- coding: utf-8 -*-
"""ShapG2P core: given a list of genes (file path or gene symbol list),
return a {pathway: SHAP score} dictionary sorted in descending order.

Method: each gene is featurized by its PPI network distance to the genes
of each pathway (1/(d+1)-like similarity, exp(-d/2)); an XGBoost classifier
distinguishes biomarkers from background genes (global 1:1 undersampling,
SEED=42); SHAP (pred_contribs) mean absolute contribution of each pathway
feature is used as the pathway SHAP score.
"""
import csv
import json
import os

import numpy as np
import pandas as pd
import xgboost as xgb
from xgboost import DMatrix

__all__ = ['score_pathways']

SEED = 42
XGB_PARAMS = dict(
    n_estimators=300,
    learning_rate=0.1,
    max_depth=4,
    tree_method='hist',
    random_state=SEED,
    eval_metric='logloss',
)

# Final model = concatenation of the following single features (same as the paper workflow)
FINAL_FEATS = ['kegg_min', 'hallmark_min', 'hallmark_mean',
               'wikipathway_min', 'wikipathway_mean']
DB_NAME = {'kegg': 'KEGG', 'hallmark': 'Hallmark',
           'wikipathway': 'WikiPathway'}
DIST_SIGMA = 2.0
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

GENE_COLS = ['string name', 'gene symbol', 'gene', 'symbol', 'marker',
             'marker symbol', 'gene name', 'hgnc symbol', 'sym',
             'gene id symbol']


def _dist_to_sim(X):
    """Distance matrix -> similarity matrix (0-1); inf (disconnected) -> 0."""
    out = np.exp(-X / DIST_SIGMA)
    out[np.isinf(X)] = 0.0
    return out


def _load_json(name):
    with open(os.path.join(DATA_DIR, name), encoding='utf-8') as f:
        return json.load(f)


def _open_csv(path):
    """Try opening a CSV with utf-8-sig, then cp1252."""
    for enc in ('utf-8-sig', 'cp1252'):
        try:
            f = open(path, encoding=enc, newline='')
            f.readline()
            f.seek(0)
            return f, csv.DictReader(f)
        except UnicodeDecodeError:
            f.close()
    return open(path, encoding='cp1252', errors='replace', newline=''), None


def _load_pos_genes(path, gene_set):
    """Read the biomarker CSV (common gene columns / STRING_ID / first column fallback)."""
    ensp2sym = {}
    with open(os.path.join(DATA_DIR, 'info_gene.csv'),
              encoding='utf-8') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) > 2 and row[2] not in ensp2sym:
                ensp2sym[row[2].replace('9606.', '')] = row[1]

    f, reader = _open_csv(path)
    if reader is None:
        raise ValueError(f'Unable to decode file: {path}')
    norm = lambda s: s.strip().lower().replace(' ', '') \
        .replace('.', '').replace('_', '')
    fields = reader.fieldnames or []
    col = next((c for c in fields if norm(c) in GENE_COLS), None)
    id_col = 'STRING_ID' if 'STRING_ID' in fields else None
    if col is None and id_col is None:
        col = fields[0]

    pos = set()
    for row in reader:
        sym = (row.get(col) or '').strip()
        if sym in ('', 'No Data') and id_col:
            sym = ensp2sym.get((row.get(id_col) or '').strip()
                               .replace('9606.', ''), '')
        if sym in gene_set:
            pos.add(sym)
    f.close()
    return pos


def _parse_input(genes, gene_set):
    """Normalize input: file path / gene list / delimited string -> in-network gene set."""
    if isinstance(genes, str) and os.path.isfile(genes):
        return _load_pos_genes(genes, gene_set)
    if isinstance(genes, str):
        items = [g.strip() for g in
                 genes.replace(',', ' ').replace('\t', ' ').replace(';', ' ')
                 .split()]
    else:
        items = [str(g).strip() for g in genes if str(g).strip()]
    return set(items) & gene_set


def _load_features():
    """Return the final feature matrix (17613 x N) and pathway metadata [(db, stat, pathway)]."""
    parts, meta = [], []
    for m in FINAL_FEATS:
        db, stat = m.rsplit('_', 1)
        X = np.load(os.path.join(DATA_DIR, f'{m}_dist.npz'))['arr_0']
        parts.append(_dist_to_sim(X.astype(np.float64)).astype(np.float32))
        with open(os.path.join(DATA_DIR, f'{db}_pathways.json'),
                  encoding='utf-8') as f:
            meta.extend((db, stat, pw) for pw in json.load(f))
    return np.hstack(parts), meta


def score_pathways(genes, verbose=True):
    """Score pathways with SHAP for a given list of genes.

    Parameters
    ----------
    genes : list[str] or str
        Gene symbol list, or a path to a CSV/TSV/TXT file containing gene
        symbols (common gene columns such as "gene symbol"/"gene"/"symbol"
        are auto-detected), or a comma/space/newline-delimited string.
    verbose : bool
        Whether to print progress.

    Returns
    -------
    dict[str, float]
        {pathway: mean_abs_shap}, sorted in descending order;
        only pathways with score > 0 are kept.
    """
    gene_list = _load_json('genes.json')
    gene_set = set(gene_list)

    pos = _parse_input(genes, gene_set)
    if not pos:
        raise ValueError('No input gene matched the network; '
                         'please check your gene symbols (HGNC symbols recommended)')
    if verbose:
        print(f'Network genes: {len(gene_set)}, matched biomarkers: {len(pos)}')

    X_all, meta = _load_features()
    y_all = np.array([1 if g in pos else 0 for g in gene_list], np.int32)

    # Global 1:1 undersampling (keep all positives, sample negatives to match, SEED=42)
    rng = np.random.default_rng(SEED)
    pos_idx = np.where(y_all == 1)[0]
    neg_idx = np.where(y_all == 0)[0]
    n = min(len(pos_idx), len(neg_idx))
    if len(pos_idx) > len(neg_idx):
        keep = rng.choice(pos_idx, size=n, replace=False)
        sample_idx = np.sort(np.concatenate([keep, neg_idx]))
    else:
        keep = rng.choice(neg_idx, size=n, replace=False)
        sample_idx = np.sort(np.concatenate([pos_idx, keep]))
    y = y_all[sample_idx]
    if verbose:
        print(f'After undersampling: pos={int(y.sum())}, neg={int((y == 0).sum())}, '
              f'{X_all.shape[1]} features, training ...')

    clf = xgb.XGBClassifier(**XGB_PARAMS)
    clf.fit(X_all[sample_idx], y)

    # Compute SHAP on all genes
    contribs = clf.get_booster().predict(DMatrix(X_all), pred_contribs=True)
    sv = contribs[:, :-1]

    rank = pd.DataFrame(meta, columns=['dataset', 'stat', 'pathway'])
    rank['mean_abs_shap'] = np.abs(sv).mean(axis=0)
    rank['dataset'] = rank['dataset'].map(DB_NAME)

    # Pathway-level aggregation: max score for duplicated names, keep score > 0, sort descending
    result = (rank.groupby('pathway')['mean_abs_shap'].max()
                  .loc[lambda s: s > 0]
                  .sort_values(ascending=False).to_dict())
    if verbose:
        top = list(result.items())[:5]
        print('Top 5 pathways: ' + ', '.join(
            f'{pw}={v:.4f}' for pw, v in top))
    return result
