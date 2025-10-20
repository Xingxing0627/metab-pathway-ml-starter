import argparse, yaml, json, ast, sys, pandas as pd
from pathlib import Path

def approx_equal(a, b, eps=1e-6):
    return abs(a-b) < eps

def load_prob(a):
    if isinstance(a, str):
        return ast.literal_eval(a)
    return a

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schema", required=True)
    ap.add_argument("--example", required=True)
    args = ap.parse_args()

    schema = yaml.safe_load(open(args.schema, "r"))
    df = pd.read_csv(args.example)

    # Field presence and basic type checks
    required = [f["name"] for f in schema["fields"] if f.get("required")]
    missing = [f for f in required if f not in df.columns]
    if missing:
        print("[FAIL] Missing required columns:", missing); sys.exit(1)

    # Probability rule
    if any("env_redox_prob" in c for c in df.columns):
        for i,row in df.iterrows():
            p = load_prob(row.get("env_redox_prob"))
            if p is None: continue
            s = sum(p)
            if not approx_equal(s, 1.0, 1e-3):
                print(f"[WARN] Row {i} redox_prob sum!=1: {s}")

    print("[OK] Basic schema checks passed for", args.example)

if __name__ == "__main__":
    main()
