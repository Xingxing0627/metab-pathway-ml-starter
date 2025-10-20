# -*- coding: utf-8 -*-
import argparse, json, math, random
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors

PROTON = 1.007276466812
NA = 22.989218
K = 38.963158
NH4 = 18.033823
ADDUCTS = {
    "[M+H]+":   lambda m: (m + PROTON),
    "[M+Na]+":  lambda m: (m + NA),
    "[M+K]+":   lambda m: (m + K),
    "[M+NH4]+": lambda m: (m + NH4),
    "[M-H]-":   lambda m: (m - PROTON),
}

def exact_mass(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return float("nan")
    return float(Descriptors.ExactMolWt(m))

def bh(q):
    n=len(q); 
    if n==0: return []
    idx=sorted(range(n), key=lambda i:q[i], reverse=True)
    out=[0.0]*n; m=1.0
    for r,i in enumerate(idx, start=1):
        v=min(m, q[i]*n/(n-r+1)); m=v; out[i]=v
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--features_csv", required=True)
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--adducts", default="[M+H]+,[M+Na]+,[M+K]+,[M+NH4]+,[M-H]-")
    ap.add_argument("--n_decoy", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--out_fdr", default="reports/metrics/hrms_fdr.csv")
    ap.add_argument("--out_png", default="reports/metrics/fdr_curve.png")
    args=ap.parse_args()

    t=json.load(open(args.tree,"r",encoding="utf-8"))
    smiles=t.get("best_path_smiles",[])
    if not smiles:
        Path(args.out_fdr).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["step","smiles","p_decoy","q_value"]).to_csv(args.out_fdr,index=False)
        print("[WARN] empty path"); return

    peaks=pd.read_csv(args.features_csv)
    if not {"mz","intensity"}<=set(peaks.columns):
        raise ValueError("features_csv must have columns mz,intensity")
    peaks=peaks[["mz","intensity"]].dropna().sort_values("mz").reset_index(drop=True)

    names=[a.strip() for a in args.adducts.split(",") if a.strip()]
    funcs=[(a,ADDUCTS[a]) for a in names if a in ADDUCTS]

    # theoretical list
    theo=[]
    for step,smi in enumerate(smiles, start=1):
        m=exact_mass(smi)
        if math.isnan(m): continue
        for a,fn in funcs:
            theo.append((step,smi,a,fn(m)))

    # best hit per (step,smi) using real peaks
    per_step=[]
    for step,smi,a,mz in theo:
        lo=mz*(1-args.ppm*1e-6); hi=mz*(1+args.ppm*1e-6)
        sub=peaks[(peaks.mz>=lo)&(peaks.mz<=hi)]
        if sub.empty: continue
        ppm=(sub.mz-mz).abs()/mz*1e6
        best_idx=ppm.idxmin()
        per_step.append((step,smi,float(ppm.loc[best_idx])))

    if not per_step:
        Path(args.out_fdr).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["step","smiles","p_decoy","q_value"]).to_csv(args.out_fdr,index=False)
        print("[OK] no hits under window; wrote empty FDR"); return

    # decoy simulation: jitter theo m/z by uniform in [-ppm, ppm] ppm
    random.seed(args.seed)
    decoy_ppm=[]
    for _ in range(args.n_decoy):
        # sample a random theoretical mz as anchor
        step,smi,a,mz=random.choice(theo)
        jitter=mz*(1+random.uniform(-args.ppm,args.ppm)*1e-6)
        lo=mz*(1-args.ppm*1e-6); hi=mz*(1+args.ppm*1e-6)
        sub=peaks[(peaks.mz>=lo)&(peaks.mz<=hi)]
        if sub.empty:
            continue
        ppm=(sub.mz-jitter).abs()/mz*1e6
        decoy_ppm.append(float(ppm.min()))
    if not decoy_ppm:
        decoy_ppm=[args.ppm]  # safe fallback

    # p-value per step as decoy CDF at observed ppm
    decoy_sorted=sorted(decoy_ppm)
    def ecdf(x):
        import bisect
        i=bisect.bisect_right(decoy_sorted, x)
        return i/len(decoy_sorted)

    obs_p=[min(1.0, ecdf(ppm)) for (_,_,ppm) in per_step]
    q=bh(obs_p)

    out=pd.DataFrame({
        "step":[s for (s,_,_) in per_step],
        "smiles":[sm for (_,sm,_) in per_step],
        "p_decoy":obs_p,
        "q_value":q
    }).sort_values("step")
    Path(args.out_fdr).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_fdr, index=False)
    print(f"[OK] wrote {args.out_fdr} rows={len(out)}")

    # plot FDR curve
    plt.figure()
    plt.plot(sorted(q), marker="o")
    plt.xlabel("rank"); plt.ylabel("q-value (decoy)"); plt.title("Decoy-based q-values")
    plt.tight_layout(); plt.savefig(args.out_png, dpi=150)
    print(f"[OK] wrote {args.out_png}")
if __name__=="__main__":
    main()
