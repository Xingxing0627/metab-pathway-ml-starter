#!/usr/bin/env bash
set -euo pipefail

IN_EVENTS="data/examples/minimal_events.csv"
FILLED="data/interim/events_filled.csv"
PACK="data/interim/pack_v1.npz"
MODEL="data/processed/lite_clf_rdkit.joblib"
CANDS="data/external/candidates_demo_30.csv"
TREE="data/processed/pathway_tree.json"
REPORT_MD="reports/pathway_best.md"
REPORT_CSV="reports/pathway_best.csv"
PEAKS="data/hrms/example_peaks.csv"

echo "[1/5] fill context"
python scripts/04_context_env_fill.py --in "$IN_EVENTS" --out "$FILLED" --env configs/env_fields.json

echo "[2/5] pack RDKit"
python scripts/21_schema_adapter.py --in_csv "$FILLED" --out_npz "$PACK" --T 1.5

echo "[3/5] train baseline"
python scripts/80_train_lite.py --config configs/project.yaml --pack "$PACK" --model_out "$MODEL"

echo "[4/5] beam search"
python scripts/33_beam_expand.py --candidates "$CANDS" --out "$TREE" --start_smiles "CC(=O)O" --k0 12 --k_decay 1 --max_depth 4 --no_loop

echo "[5/5] render + HRMS"
python scripts/34_render_path.py --tree "$TREE" --candidates "$CANDS" --md_out "$REPORT_MD" --csv_out "$REPORT_CSV"
python scripts/91_hrms_validate_ppm.py --tree "$TREE" --features_csv "$PEAKS" --ppm 10 --adducts "[M+H]+,[M+Na]+,[M+K]+,[M+NH4]+" --min_intensity 800 --out_hits reports/metrics/hrms_hits.csv --out_q reports/metrics/hrms_qvalues.csv

echo "[OK] done. See $REPORT_MD and reports/metrics/"
python scripts/93_plot_qcurve.py --q_csv reports/metrics/hrms_qvalues.csv --out_png reports/metrics/qcurve.png
