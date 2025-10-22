# -*- coding: utf-8 -*-
"""
Context/Env filling with robust redox probability handling.

- Accept either column name 'env_redox_prob' (legacy) or 'redox_prob' (contract).
- Accept values as: JSON string ({"aerobic":0.7,...}), python dict, list/tuple,
  or empty. They are normalized to an ordered vector [aerobic, anoxic, anaerobic, unknown].
- Apply temperature scaling in logit space (T>=1). Write a JSON string column
  '<col>_scaled' alongside original.
"""

import argparse
import json
import pandas as pd
from pathlib import Path

ORDER_REDOX = ["aerobic", "anoxic", "anaerobic", "unknown"]

def parse_prob(x):
    """Return a 4-d list aligned to ORDER_REDOX."""
    if x is None or (isinstance(x, str) and not x.strip()):
        d = {"unknown": 1.0}
    elif isinstance(x, (list, tuple)):
        try:
            vec = [float(v) for v in x]
            # If a plain vector is given, assume it is already aligned; pad/trim to 4
            vec = (vec + [0.0, 0.0, 0.0, 0.0])[:4]
            s = sum(vec)
            return [v / s for v in vec] if s > 0 else [0.0, 0.0, 0.0, 1.0]
        except Exception:
            d = {"unknown": 1.0}
    elif isinstance(x, dict):
        d = x
    else:
        # string -> try json then literal
        s = str(x).strip()
        d = None
        try:
            obj = json.loads(s)
            if isinstance(obj, dict):
                d = obj
            elif isinstance(obj, (list, tuple)):
                vec = [float(v) for v in obj]
                vec = (vec + [0.0, 0.0, 0.0, 0.0])[:4]
                ssum = sum(vec)
                return [v / ssum for v in vec] if ssum > 0 else [0.0, 0.0, 0.0, 1.0]
        except Exception:
            pass
        if d is None:
            try:
                import ast as _ast
                obj = _ast.literal_eval(s)
                if isinstance(obj, dict):
                    d = obj
                elif isinstance(obj, (list, tuple)):
                    vec = [float(v) for v in obj]
                    vec = (vec + [0.0, 0.0, 0.0, 0.0])[:4]
                    ssum = sum(vec)
                    return [v / ssum for v in vec] if ssum > 0 else [0.0, 0.0, 0.0, 1.0]
            except Exception:
                d = {"unknown": 1.0}
    vec = [float(d.get(k, 0.0)) for k in ORDER_REDOX]
    ssum = sum(vec)
    if ssum > 0:
        vec = [v / ssum for v in vec]
    else:
        vec = [0.0, 0.0, 0.0, 1.0]
    return vec

def temp_scale_prob(vec, T=1.5):
    """Temperature scale probabilities in logit space and renormalize."""
    import math
    eps = 1e-12
    v = [max(eps, min(1.0 - eps, float(x))) for x in vec]
    logits = [math.log(p / (1.0 - p)) for p in v]
    scaled = [l / max(T, eps) for l in logits]
    # back to prob via sigmoid, then renorm to sum 1 (defensive)
    probs = [1.0 / (1.0 + math.exp(-l)) for l in scaled]
    ssum = sum(probs)
    return [p / ssum for p in probs] if ssum > 0 else [0.0, 0.0, 0.0, 1.0]

def to_json_from_vec(vec):
    return json.dumps({k: float(vec[i]) for i, k in enumerate(ORDER_REDOX)})

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="fin", required=True)
    ap.add_argument("--out", dest="fout", required=True)
    ap.add_argument("--env", dest="env_json", default=None, help="unused placeholder for now")
    ap.add_argument("--T", type=float, default=1.5)
    args = ap.parse_args()

    df = pd.read_csv(args.fin)

    # Backward compatibility for column name
    source_col = None
    if "env_redox_prob" in df.columns:
        source_col = "env_redox_prob"
    elif "redox_prob" in df.columns:
        source_col = "redox_prob"
    else:
        # If no column at all, create unknown prior
        df["redox_prob"] = [to_json_from_vec([0.0, 0.0, 0.0, 1.0]) for _ in range(len(df))]
        source_col = "redox_prob"

    # Parse -> scale -> write side column
    vecs = df[source_col].apply(parse_prob).tolist()
    scaled = [temp_scale_prob(v, T=args.T) for v in vecs]
    out_col = f"{source_col}_scaled"
    df[out_col] = [to_json_from_vec(v) for v in scaled]

    Path(args.fout).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.fout, index=False)
    print(f"[OK] wrote {args.fout}, rows={len(df)}")

if __name__ == "__main__":
    main()
