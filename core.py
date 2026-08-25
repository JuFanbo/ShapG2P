# -*- coding: utf-8 -*-
"""ShapG2P 核心: 输入基因列表 (文件或 gene symbol list),
输出 {通路名: SHAP 分数} 字典, 按分数从大到小。

方法: 每个基因用其到各通路基因的 PPI 网络距离 (1/(d+1)-类相似度,
exp(-d/2)) 作为特征, XGBoost 区分 biomarker 与背景基因 (全局 1:1
欠采样, SEED=42), 用 SHAP (pred_contribs) 得到每条通路特征对分类的
平均绝对贡献, 即为通路 SHAP 分数。
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

# final 模型 = 下列单特征拼接 (与论文工作流一致)
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
    """距离矩阵 → 相似度矩阵 (0-1), inf (不连通) 映射为 0。"""
    out = np.exp(-X / DIST_SIGMA)
    out[np.isinf(X)] = 0.0
    return out


def _load_json(name):
    with open(os.path.join(DATA_DIR, name), encoding='utf-8') as f:
        return json.load(f)


def _open_csv(path):
    """按 utf-8-sig -> cp1252 顺序尝试打开 CSV。"""
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
    """读标志物 CSV (兼容常见基因列名 / STRING_ID / 第一列兜底)。"""
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
        raise ValueError(f'无法解码文件: {path}')
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
    """输入归一化: 文件路径 / 基因列表 / 分隔符字符串 -> 网络内基因集合。"""
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
    """返回 final 特征矩阵 (17613 x N) 与通路元数据 [(db, stat, pathway)]。"""
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
    """对输入基因列表做 SHAP 通路评分。

    参数
    ----
    genes : list[str] 或 str
        基因符号列表, 或包含基因符号的 CSV/TSV/TXT 文件路径
        (文件需含常见基因列, 如 "gene symbol"/"gene"/"symbol" 等),
        或逗号/空格/换行分隔的字符串。
    verbose : bool
        是否打印进度。

    返回
    ----
    dict[str, float]
        {通路名: mean_abs_shap}, 按 SHAP 分数从大到小排序;
        仅保留分数 > 0 的通路。
    """
    gene_list = _load_json('genes.json')
    gene_set = set(gene_list)

    pos = _parse_input(genes, gene_set)
    if not pos:
        raise ValueError('输入中未匹配到任何网络内基因, '
                         '请检查基因符号 (建议使用 HGNC 标准符号)')
    if verbose:
        print(f'网络基因: {len(gene_set)}, 匹配 biomarker: {len(pos)}')

    X_all, meta = _load_features()
    y_all = np.array([1 if g in pos else 0 for g in gene_list], np.int32)

    # 全局 1:1 欠采样 (正样本全留, 负样本随机抽至等量, SEED=42)
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
        print(f'欠采样后: pos={int(y.sum())}, neg={int((y == 0).sum())}, '
              f'特征 {X_all.shape[1]} 维, 训练中 ...')

    clf = xgb.XGBClassifier(**XGB_PARAMS)
    clf.fit(X_all[sample_idx], y)

    # SHAP 在全部基因上计算
    contribs = clf.get_booster().predict(DMatrix(X_all), pred_contribs=True)
    sv = contribs[:, :-1]

    rank = pd.DataFrame(meta, columns=['dataset', 'stat', 'pathway'])
    rank['mean_abs_shap'] = np.abs(sv).mean(axis=0)
    rank['dataset'] = rank['dataset'].map(DB_NAME)

    # 通路级聚合: 同名通路取最大分数, 仅保留分数 > 0, 按分数降序
    result = (rank.groupby('pathway')['mean_abs_shap'].max()
                  .loc[lambda s: s > 0]
                  .sort_values(ascending=False).to_dict())
    if verbose:
        top = list(result.items())[:5]
        print('Top 5 通路: ' + ', '.join(
            f'{pw}={v:.4f}' for pw, v in top))
    return result
