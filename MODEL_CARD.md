# Model Card (Template)

## Overview
- **Model**: Metabolic event & product prediction
- **Intended Use**: Research
- **Out-of-scope**: Clinical or regulatory decisions

## Data
- Sources: (填写Rhea/ChEBI/UniProt等)
- Splits: scaffold / enzyme / species / homology
- OOD: 未见酶/物种

## Training
- Objectives: occur(T1), center(T2), retrieval(T3)
- Calibration: temperature scaling per stratum

## Evaluation
- Metrics: AUROC, IoU, Recall@k, ECE, FDR
- Risk: 数据泄漏、候选库漂移
