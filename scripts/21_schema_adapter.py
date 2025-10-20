# -*- coding: utf-8 -*-
"""
Schema adapter:
- Read single-step events CSV
- Normalize env distributions (temperature scaling) + entropy
- Encode pH bin & temperature
- Build RDKit Morgan fingerprints for substrates
- Pack to NPZ for training
"""
import argparse, json, ast
import numpy as np
import pandas as pd
from pathlib import Path

# RDKit imports for fingerprints
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import ConvertToNumpyArray

def load_prob(s):
    if isinstance(s, str):
        try: return ast.literal_eval(s)
        except Exception: return None
    return s

def temp_scale_prob(p, T=1.5):
    """Temperature-scale a probability vector via logit -> /T -> softmax-like normalize."""
    import math
    eps = 1e-12
    p = [max(eps, min(1-eps, float(v))) for v in p]
    logit = [math.log(v/(1.0-v)) for v in p]
    mx = max(logit)/max(T, eps)
    ex = [math.exp((lv/max(T,eps)) - mx) for lv in logit]
    s = sum(ex)
    return [float(v/s) for v in ex]

def entropy(p):
    eps = 1e-12
    arr = np.clip(np.asarray(p, dtype=float), eps, 1.0)
    return float(-np.sum(arr * np.log(arr)))

from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
_morgan_gen_cache = {}
def morgan_fp(smiles: str, nbits: int = 1024, radius: int = 2):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros((nbits,), dtype="float32")
    key = (nbits, radius)
    gen = _morgan_gen_cache.get(key)
    if gen is None:
        gen = GetMorganGenerator(radius=radius, fpSize=nbits)
        _morgan_gen_cache[key] = gen
    bv = gen.GetFingerprint(mol)
    arr = np.zeros((nbits,), dtype=np.int8)
    ConvertToNumpyArray(bv, arr)
    return arr.astype("float32")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True)
    ap.add_argument("--out_npz", required=True)
    ap.add_argument("--T", type=float, default=1.5, help="temperature for env_redox_prob scaling")
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)

    # --- Build X_env: [scaled redox probs..., entropy, ph_bin, temp_c]
    X_env_rows = []
    ph_col = "env_ph_bin"
    temp_col = "temp_c"
    p_col = "env_redox_prob"
    for _, row in df.iterrows():
        p = load_prob(row.get(p_col))
        if p is None:
            p = [1/3, 1/3, 1/3]  # default 3-level redox
        ps = temp_scale_prob(p, T=args.T)
        ent = entropy(ps)
        ph = float(row.get(ph_col, 0))
        tc = float(row.get(temp_col, 25.0))
        X_env_rows.append(list(ps) + [ent, ph, tc])
    X_env = np.asarray(X_env_rows, dtype="float32")

    # --- Build X_mol with RDKit Morgan fingerprints (substrate_smiles)
    smi_col = "substrate_smiles"
    fps = [morgan_fp(str(smi), nbits=1024, radius=2) for smi in df[smi_col].astype(str).tolist()]
    X_mol = np.vstack(fps).astype("float32")

    # --- Labels (occur)
    y = df.get("occur_label", pd.Series([0]*len(df))).astype(int).values

    # Save
    Path(args.out_npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out_npz, X_env=X_env, X_mol=X_mol, y_occur=y)
    print(f"[OK] wrote {args.out_npz} with X_env={X_env.shape}, X_mol={X_mol.shape}, y={y.shape}")

if __name__ == "__main__":
    main()
