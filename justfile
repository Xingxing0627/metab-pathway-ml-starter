setup:
	conda env create -f environment.yml
	pre-commit install

check:
	python scripts/00_check_env.py

data:
	python scripts/04_context_env_fill.py --in data/examples/minimal_events.csv --out data/interim/events_filled.csv --env configs/env_fields.json

lite:
	python scripts/80_train_lite.py --config configs/project.yaml

beam:
	python scripts/33_beam_expand.py --candidates data/examples/candidates_v0.csv --out data/processed/pathway_tree.json

hrms:
	python scripts/90_hrms_validation.py --in data/processed/pathway_tree.json --out reports/metrics/hrms_fdr.csv
