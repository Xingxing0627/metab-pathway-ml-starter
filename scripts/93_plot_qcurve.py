# -*- coding: utf-8 -*-
import argparse, pandas as pd
import matplotlib.pyplot as plt
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--q_csv", required=True)
    ap.add_argument("--out_png", default="reports/metrics/qcurve.png")
    args=ap.parse_args()
    df=pd.read_csv(args.q_csv)
    if df.empty or "q_value" not in df.columns:
        print("[WARN] empty q-values"); return
    y=df["q_value"].sort_values().values
    x=list(range(1, len(y)+1))
    plt.figure()
    plt.plot(x, y, marker="o")
    plt.xlabel("rank"); plt.ylabel("q-value"); plt.title("HRMS q-value curve")
    plt.tight_layout(); plt.savefig(args.out_png, dpi=150)
    print("[OK] wrote", args.out_png)
if __name__=="__main__": main()
