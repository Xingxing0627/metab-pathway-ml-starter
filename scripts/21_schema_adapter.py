# -*- coding: utf-8 -*-
"""
Schema adapter:
- Read events CSV
- Build X_env (pH, temp_C, redox one-hot-like prob vector of len 4)
- Build X_mol (RDKit Morgan FP 1024 from product/main_product/substrate_smiles)
- y_occur from 'occur'
- Robust handling for redox_prob (accept: JSON str, dict, list/tuple; column name 'redox_prob' or legacy 'env_redox_prob')
"""

import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import ConvertToNumpyArray

ORDER_REDOX = ["aerobic","anoxic","anaerobic","unknown"]

def parse_prob(x):
    """Return length-4 vector aligned to ORDER_REDOX."""
    if x is None or (isinstance(x,str) and not x.strip()):
        d = {"unknown":1.0}
    elif isinstance(x,(list,tuple,np.ndarray)):
        try:
            vec = [float(v) for v in x]
            vec = (vec + [0.0,0.0,0.0,0.0])[:4]
            s = sum(vec)
            return [v/s for v in vec] if s>0 else [0.0,0.0,0.0,1.0]
        except Exception:
            d = {"unknown":1.0}
    elif isinstance(x,dict):
        d = x
    else:
        s = str(x).strip()
        d = None
        # try JSON
        try:
            obj = json.loads(s)
            if isinstance(obj,dict):
                d = obj
            elif isinstance(obj,(list,tuple)):
                vec = [float(v) for v in obj]
                vec = (vec + [0.0,0.0,0.0,0.0])[:4]
                ss = sum(vec)
                return [v/ss for v in vec] if ss>0 else [0.0,0.0,0.0,1.0]
        except Exception:
            pass
        if d is None:
            # try literal
            try:
                import ast
                obj = ast.literal_eval(s)
                if isinstance(obj,dict):
                    d = obj
                elif isinstance(obj,(list,tuple)):
                    vec = [float(v) for v in obj]
                    vec = (vec + [0.0,0.0,0.0,0.0])[:4]
                    ss = sum(vec)
                    return [v/ss for v in vec] if ss>0 else [0.0,0.0,0.0,1.0]
            except Exception:
                d = {"unknown":1.0}
    vec = [float(d.get(k,0.0)) for k in ORDER_REDOX]
    ss = sum(vec)
    return [v/ss for v in vec] if ss>0 else [0.0,0.0,0.0,1.0]

def temp_scale_prob(vec, T=1.5):
    """Temperature scaling in logit space, then renorm."""
    eps = 1e-12
    v = [max(eps, min(1.0-eps, float(x))) for x in vec]
    logits = [math.log(p/(1.0-p)) for p in v]
    scaled = [l/max(T,eps) for l in logits]
    probs = [1.0/(1.0+math.exp(-l)) for l in scaled]
    s = sum(probs)
    return [p/s for p in probs] if s>0 else [0.0,0.0,0.0,1.0]

def pick_smiles_row(row):
    for col in ("main_product","product_smiles","substrate_smiles"):
        if col in row and isinstance(row[col], str) and row[col].strip():
            return row[col].strip()
    return ""

def morgan_fp_1024(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        arr = np.zeros((1024,), dtype=np.int8)
        return arr
    fp = AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=1024)
    arr = np.zeros((1024,), dtype=np.int8)
    ConvertToNumpyArray(fp, arr)
    return arr

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True)
    ap.add_argument("--out_npz", required=True)
    ap.add_argument("--T", type=float, default=1.5)
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)

    # redox column name compatibility
    redox_col = "redox_prob" if "redox_prob" in df.columns else ("env_redox_prob" if "env_redox_prob" in df.columns else None)
    redox_vecs = []
    if redox_col:
        for val in df[redox_col].tolist():
            vec = parse_prob(val)
            vec = temp_scale_prob(vec, T=args.T)
            redox_vecs.append(vec)
    else:
        redox_vecs = [[0.0,0.0,0.0,1.0] for _ in range(len(df))]

    # numeric env features
    def to_float(x, default):
        try:
            return float(x)
        except Exception:
            return float(default)

    pH = [to_float(v, 7.0) for v in df["pH"].tolist()] if "pH" in df.columns else [7.0]*len(df)
    tempC = [to_float(v, 25.0) for v in df["temp_C"].tolist()] if "temp_C" in df.columns else [25.0]*len(df)

    # X_env = [pH, temp_C, redox(4)]
    X_env = np.array([[pH[i], tempC[i]] + redox_vecs[i] for i in range(len(df))], dtype=np.float32)

    # X_mol from product/main_product/substrate_smiles
    smiles = [pick_smiles_row(r) for _, r in df.iterrows()]
    X_mol = np.stack([morgan_fp_1024(s) for s in smiles], axis=0).astype(np.int8)

    # y
    y = df["occur"].fillna(0).astype(int).clip(0,1).to_numpy()

    Path(args.out_npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out_npz, X_env=X_env, X_mol=X_mol, y_occur=y)
    print(f"[OK] wrote {args.out_npz} with X_env={X_env.shape}, X_mol={X_mol.shape}, y={y.shape}")

if __name__ == "__main__":
    main()
