# -*- coding: utf-8 -*-
import argparse, json, math
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

PROTON = 1.007276466812
NA = 22.989218
ADDUCTS = {
    "[M+H]+":  lambda m: (m + PROTON),
    "[M+Na]+": lambda m: (m + NA),
    "[M-H]-":  lambda m: (m - PROTON),
}

def exact_mass(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return float("nan")
    return float(Descriptors.ExactMolWt(m))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--adducts", default="[M+H]+,[M+Na]+,[M-H]-")
    ap.add_argument("--out", default="reports/pathway_theo_mz.csv")
    args = ap.parse_args()

    t = json.load(open(args.tree, "r", encoding="utf-8"))
    smiles_list = t.get("best_path_smiles", [])
    names = [a.strip() for a in args.adducts.split(",") if a.strip()]
    funcs = [(a, ADDUCTS[a]) for a in names if a in ADDUCTS]

    rows=[]
    for step, smi in enumerate(smiles_list, start=1):
        m = exact_mass(smi)
        if math.isnan(m): continue
        for a,fn in funcs:
            rows.append({"step":step, "smiles":smi, "adduct":a, "theo_mz":fn(m)})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print("[OK] wrote", args.out)

if __name__ == "__main__":
    main()
