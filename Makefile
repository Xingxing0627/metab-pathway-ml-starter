.PHONY: setup check data lite beam hrms pos neg bundle

setup:
	conda env create -f environment.yml || true
	pre-commit install || true
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
.PHONY: pos neg bundle
pos:
	PEAKS=data/hrms/synthetic_peaks.csv ADDUCTS='[M+H]+,[M+Na]+,[M+K]+,[M+NH4]+' PPM=10 MINI=800 ./scripts/39_reproduce_pathway.sh
neg:
	PEAKS=data/hrms/real_neg_std.csv ADDUCTS='[M-H]-' PPM=20 MINI=0 ./scripts/39_reproduce_pathway.sh
.PHONY: help
help:
	@echo make pos|neg|bundle|check|data|lite|beam|hrms
bundle: pos
	python scripts/11_release_packager.py --out release_bundle.tar.gz
