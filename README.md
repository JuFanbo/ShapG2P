# ShapG2P

ShapG2P: a strategy for biomarker pathway enrichment with PPI network topology and SHAP analysis.

输入一个基因符号列表（文件或直接传 list），输出 `{通路名: SHAP 分数}` 字典，按分数从大到小排序。

## 方法

- 每个基因以其到各通路基因的 PPI 网络距离（`exp(-d/2)` 相似度）为特征；
- XGBoost 分类器区分 biomarker 与背景基因（全局 1:1 欠采样，SEED=42）；
- 用 SHAP（`pred_contribs`）计算每条通路特征的平均绝对贡献 = 通路 SHAP 分数。

内置数据：STRING 人类 PPI 网络（17613 基因）、KEGG / Hallmark / WikiPathway 通路。

## 安装

```bash
pip install shapg2p
```

或本地开发安装：

```bash
cd shapg2p_pkg
pip install -e .
```

## 用法

```python
from shapg2p import score_pathways

# 方式 1: 直接传基因符号列表
scores = score_pathways(['TP53', 'ATM', 'APOE', 'SOD1', 'CDKN2A'])

# 方式 2: 传文件路径 (CSV/TSV/TXT, 含常见基因列如 gene symbol / gene / symbol)
scores = score_pathways('my_biomarkers.csv')

# 方式 3: 传逗号/空格/换行分隔的字符串
scores = score_pathways('TP53, ATM, APOE')

# 输出: 字典, 按 SHAP 分数从大到小
print(scores)
# {'p53 signaling pathway': 0.42, 'Alzheimer disease': 0.31, ...}
```

运行约需 1–3 分钟（一次 XGBoost 训练 + 全基因 SHAP 计算）。

## 输出说明

返回 `dict[str, float]`：key 为通路名，value 为该通路的 mean |SHAP| 分数，按分数降序；
仅包含分数 > 0 的通路（无关通路不返回）。同名通路（出现在多个数据库中）取最大分数。

## 上传到 PyPI（供他人 pip install）

### 1. 注册账号并创建 API token

- 注册：https://pypi.org/account/register/
- 创建 token：https://pypi.org/manage/account/token/ （Scope 选整个账号即可）
- token 形如 `pypi-AgEIcHlwaS5vcmcC...`，只显示一次，保存好

### 2. 构建

```bash
pip install --upgrade build twine
cd shapg2p_pkg
python -m build
```

生成 `dist/shapg2p-0.1.0.tar.gz` 和 `dist/shapg2p-0.1.0-py3-none-any.whl`。

### 3. 先传 TestPyPI 验证（可选但推荐）

```bash
python -m twine upload --repository testpypi dist/*
# 用户名输入: __token__   密码输入: pypi-xxx (你的 API token)

# 验证安装
pip install --index-url https://test.pypi.org/simple/ shapg2p
```

### 4. 上传正式 PyPI

```bash
python -m twine upload dist/*
```

同样输入 `__token__` + API token。上传成功后即可：

```bash
pip install shapg2p
```

### 5. 更新版本

修改 `pyproject.toml` 的 `version`（如 `0.1.1`），重新 `python -m build` 并 `twine upload dist/*`。
PyPI 不允许重复上传相同版本号。

## 依赖

numpy / pandas / scipy / scikit-learn / xgboost（`pip install shapg2p` 时自动安装）。
