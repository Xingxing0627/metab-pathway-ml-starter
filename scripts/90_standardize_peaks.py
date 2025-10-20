# -*- coding: utf-8 -*-
import argparse, pandas as pd
from pathlib import Path

MZ_CAND = ["mz","m/z","mass","m.z.","m z","masses"]
I_CAND  = ["intensity","abundance","height","ion_intensity","peak_intensity","signal"]

def pick_col(cols, cands):
    low = [c.strip().lower() for c in cols]
    for name in cands:
        if name in low:
            return cols[low.index(name)]
    # fallback: first numeric-like
    return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--in", dest="fin", required=True)
    ap.add_argument("--out", dest="fout", required=True)
    args=ap.parse_args()

    # auto-detect delimiter
    df = pd.read_csv(args.fin, sep=None, engine="python")
    mz_col = pick_col(df.columns, MZ_CAND)
    it_col = pick_col(df.columns, I_CAND)
    if mz_col is None or it_col is None:
        raise SystemExit(f"[ERR] cannot find mz/intensity columns in {list(df.columns)}")

    out = df[[mz_col, it_col]].rename(columns={mz_col:"mz", it_col:"intensity"})
    out = out.dropna().astype({"mz":"float64","intensity":"float64"})
    Path(args.fout).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.fout, index=False)
    print(f"[OK] wrote standardized peaks -> {args.fout} rows={len(out)}")
if __name__=="__main__":
    main()
