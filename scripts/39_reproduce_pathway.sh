set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PEAKS="${PEAKS:-data/hrms/synthetic_peaks.csv}"
ADDUCTS="${ADDUCTS:-[M+H]+,[M+Na]+,[M+K]+,[M+NH4]+}"
PPM="${PPM:-10}"
MINI="${MINI:-800}"

echo "[INFO] using PEAKS=$PEAKS  ADDUCTS=$ADDUCTS  PPM=$PPM  MIN_INT=$MINI"
python scripts/04_context_env_fill.py --in data/examples/minimal_events.csv --out data/interim/events_filled.csv --env configs/env_fields.json
python scripts/21_schema_adapter.py --in_csv data/interim/events_filled.csv --out_npz data/interim/pack_v1.npz --T 1.5
python scripts/80_train_lite.py --config configs/project.yaml --pack data/interim/pack_v1.npz --model_out data/processed/lite_clf_rdkit.joblib
python scripts/33_beam_expand.py --candidates data/external/candidates_demo_30.csv --out data/processed/pathway_tree.json --start_smiles "CC(=O)O" --k0 12 --k_decay 1 --max_depth 4 --no_loop
python scripts/34_render_path.py --tree data/processed/pathway_tree.json --candidates data/external/candidates_demo_30.csv --md_out reports/pathway_best.md --csv_out reports/pathway_best.csv
python scripts/92_export_theo_mz.py --tree data/processed/pathway_tree.json --adducts "[M+H]+,[M+Na]+,[M+K]+,[M+NH4]+,[M-H]-" --out reports/pathway_theo_mz.csv
python scripts/91_hrms_validate_ppm.py --tree data/processed/pathway_tree.json --features_csv "$PEAKS" --ppm "$PPM" --adducts "$ADDUCTS" --min_intensity "$MINI" --out_hits reports/metrics/hrms_hits.csv --out_q reports/metrics/hrms_qvalues.csv
python scripts/93_plot_qcurve.py --q_csv reports/metrics/hrms_qvalues.csv --out_png reports/metrics/qcurve.png || true
GITHASH="$(git rev-parse --short HEAD 2>/dev/null || echo "no-git")"
DATE="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
MD5_EVENTS="$(md5sum data/examples/minimal_events.csv | awk '{print $1}')"
MD5_CANDS="$(md5sum data/external/candidates_demo_30.csv | awk '{print $1}')"
MD5_PEAKS="$(md5sum "$PEAKS" | awk '{print $1}')"
cat > PROVENANCE.yml <<YAML
run_at_utc: $DATE
git_hash: $GITHASH
inputs:
  events_csv: data/examples/minimal_events.csv
  candidates_csv: data/external/candidates_demo_30.csv
  peaks_csv: $PEAKS
input_md5:
  events: $MD5_EVENTS
  candidates: $MD5_CANDS
  peaks: $MD5_PEAKS
params:
  ppm: $PPM
  adducts: "$ADDUCTS"
  min_intensity: $MINI
artifacts:
  pathway_tree: data/processed/pathway_tree.json
  report_md: reports/pathway_best.md
  report_csv: reports/pathway_best.csv
  theo_mz_csv: reports/pathway_theo_mz.csv
  hrms_hits_csv: reports/metrics/hrms_hits.csv
  hrms_qvalues_csv: reports/metrics/hrms_qvalues.csv
  qcurve_png: reports/metrics/qcurve.png
YAML

echo "[OK] reproduce complete. See PROVENANCE.yml and reports/"
