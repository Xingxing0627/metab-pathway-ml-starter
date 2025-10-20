# Metabolic Pathway ML — Starter Kit

**Goal**: 从0起步，搭建“酶促代谢产物预测 + 多步路径搜索 + HRMS 验证”的最小可运行（MVP）工程骨架。

## 快速开始

```bash
# 1) 创建conda环境（或mamba）
conda env create -f environment.yml
conda activate metab-ml

# 2) 可选：安装预提交钩子
pre-commit install

# 3) 运行环境自检
python scripts/00_check_env.py

# 4) 合同/配置检查
python scripts/10_schema_contract_check.py --schema configs/schema_contract.yaml --example data/examples/minimal_events.csv

# 5) 预处理 & 填充环境上下文（示例）
python scripts/04_context_env_fill.py --in data/examples/minimal_events.csv --out data/interim/events_filled.csv --env configs/env_fields.json

# 6) 训练Lite占位（示例，不真正训练）
python scripts/80_train_lite.py --config configs/project.yaml

# 7) 路径搜索占位（Beam skeleton）
python scripts/33_beam_expand.py --candidates data/examples/candidates_v0.csv --out data/processed/pathway_tree.json

# 8) HRMS 验证占位
python scripts/90_hrms_validation.py --in data/processed/pathway_tree.json --out reports/metrics/hrms_fdr.csv
```

**重要**：本仓库是“骨架”，核心算法留有 TODO。你可以逐步替换 `scripts/` 与 `src/` 中的占位实现。

---

## 目录结构
```
metab-pathway-ml-starter/
  configs/                 # 项目配置 & 合同（contracts）
  data/
    raw/                   # 原始数据（只读）
    external/              # 外部来源数据（数据库导出等）
    interim/               # 中间产物
    processed/             # 最终分析/建模产物
    examples/              # 示例最小数据
  notebooks/               # 可选：探索笔记本
  reports/
    figures/
    metrics/
  scripts/                 # CLI 脚本（管道入口）
  src/metabml/             # 代码库（可被脚本调用）
  tests/                   # 单元/集成测试（占位）
  environment.yml          # conda 环境
  requirements.txt         # 纯pip备选
  Makefile / justfile      # 常用任务封装
  MODEL_CARD.md            # 模型卡模板
  PROVENANCE.yml           # 数据谱系/复现记录模板
```

## 最小MVP里程碑（首周）
- ✅ 环境可创建、脚本可运行、合同可校验
- ✅ `minimal_events.csv` 可通过 `04_context_env_fill.py` 产出 `events_filled.csv`
- ✅ `33_beam_expand.py` 能读入 candidates 并生成一个**小**的 `pathway_tree.json`
- ✅ `90_hrms_validation.py` 能输出一条示范 q-value 曲线（随机/诱饵占位）

> 之后逐步替换 Lite/Plus 模型与真实 HRMS 流程。
