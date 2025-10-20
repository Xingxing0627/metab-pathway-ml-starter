# -*- coding: utf-8 -*-
import argparse, pandas as pd
from rdkit import Chem

def check_smiles(col):
    bad = []
    for i, smi in enumerate(col.astype(str)):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            bad.append((i, smi))
    return bad

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--col", default="substrate_smiles")
    ap.add_argument("--out_bad", default="reports/metrics/bad_smiles.csv")
    args = ap.parse_args()
    df = pd.read_csv(args.csv)
    bad = check_smiles(df[args.col])
    if bad:
        import csv, os
        os.makedirs("reports/metrics", exist_ok=True)
        with open(args.out_bad, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["row_idx","bad_smiles"])
            w.writerows(bad)
        print(f"[WARN] {len(bad)} invalid SMILES written to {args.out_bad}")
    else:
        print("[OK] all SMILES valid")

if __name__ == "__main__":
    main()
