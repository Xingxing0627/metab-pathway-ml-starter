# -*- coding: utf-8 -*-
import argparse, json, pandas as pd
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="data/processed/pathway_tree.json")
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--md_out", default="reports/pathway_best.md")
    ap.add_argument("--csv_out", default="reports/pathway_best.csv")
    args = ap.parse_args()

    t = json.load(open(args.tree, "r", encoding="utf-8"))
    cands = pd.read_csv(args.candidates)
    idxs = t.get("best_path_indices", [])
    rows=[]
    for step,i in enumerate(idxs, start=1):
        smi = str(cands.loc[i, "product_smiles"]) if "product_smiles" in cands.columns else ""
        name = str(cands.loc[i, "name"]) if "name" in cands.columns else ""
        tag = str(cands.loc[i, "rule_tag"]) if "rule_tag" in cands.columns else ""
        rows.append({"step":step,"product_smiles":smi,"name":name,"rule_tag":tag})
    df = pd.DataFrame(rows)
    Path(args.md_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.csv_out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.csv_out, index=False)

    with open(args.md_out, "w", encoding="utf-8") as f:
        f.write(f"# Best Pathway\n\n")
        f.write(f"- start_smiles: `{t.get('start_smiles','')}`\n")
        f.write(f"- allow_escape_elems: {t.get('allow_escape_elems',[])}\n")
        f.write(f"- allow_escape_molecules: {t.get('allow_escape_molecules',[])}\n")
        f.write(f"- beams_per_depth: {t.get('beams_per_depth',[])}\n\n")
        for r in rows:
            f.write(f"{r['step']:02d}. {r['name'] or r['product_smiles']}  ({r['rule_tag']})\n")
    print("[OK] wrote", args.md_out, "and", args.csv_out)

if __name__ == "__main__":
    main()
