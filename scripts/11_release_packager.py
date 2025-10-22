# -*- coding: utf-8 -*-
import argparse, hashlib, json, tarfile
from pathlib import Path

KEEP = [
  "scripts", "configs", "docs", "reports/pathway_best.md",
  "reports/metrics/hrms_qvalues.csv", "reports/metrics/hrms_fdr.csv",
  "PROVENANCE.yml", "docs/run_recipe.yaml", "docs/schema_contract.yaml"
]

def md5(p: Path) -> str:
    h=hashlib.md5()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(8192), b""):
            h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out", default="release_bundle.tar.gz")
    args=ap.parse_args()

    files=[]
    for pat in KEEP:
        for p in Path(".").glob(pat):
            if p.is_file():
                files.append(p)

    manifest=[]
    for p in files:
        manifest.append({"path": str(p), "md5": md5(p)})

    Path("reports").mkdir(parents=True, exist_ok=True)
    Path("docs").mkdir(parents=True, exist_ok=True)
    Path("reports/manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    with tarfile.open(args.out, "w:gz") as tar:
        for p in files+[Path("reports/manifest.json")]:
            tar.add(str(p), arcname=str(p))
    print("[OK] wrote", args.out, "with", len(files), "files")

if __name__=="__main__":
    main()
