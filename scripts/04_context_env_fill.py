# -*- coding: utf-8 -*-
"""
Context/Env filling with robust redox probability handling and rule-based conflicts.

Features:
- Accept either column name 'redox_prob' (contract) or legacy 'env_redox_prob'.
- If redox prior is missing/empty, fill from a context dictionary:
    docs/context_env_fields.json  (or path via --dict / legacy --env)
  using 'priors' keyed by microbiome_type, otherwise default to {"unknown":1.0}.
- Write temperature-scaled redox vector into 'redox_prob_scaled' (logit T-scaling).
- Apply simple conflict rules from 'conflicts' section to set 'env_conflict' (0/1).
- Keep outputs contract-friendly; if only legacy col exists, also write 'redox_prob'.

CLI:
  --in / --out: input and output CSV
  --dict: path to context dict JSON (optional)
  --env:  legacy alias for --dict (optional)
  --T:    temperature for scaling (default 1.5)
"""

import argparse
import json
from pathlib import Path
import math
import pandas as pd

ORDER_REDOX = ["aerobic", "anoxic", "anaerobic", "unknown"]

def parse_prob(x):
    """Return a 4-d list aligned to ORDER_REDOX."""
    if x is None or (isinstance(x, str) and not x.strip()):
        d = {"unknown": 1.0}
    elif isinstance(x, (list, tuple)):
        try:
            vec = [float(v) for v in x]
            vec = (vec + [0.0, 0.0, 0.0, 0.0])[:4]
            s = sum(vec)
            return [v / s for v in vec] if s > 0 else [0.0, 0.0, 0.0, 1.0]
        except Exception:
            d = {"unknown": 1.0}
    elif isinstance(x, dict):
        d = x
    else:
        s = str(x).strip()
        d = None
        # try JSON
        try:
            obj = json.loads(s)
            if isinstance(obj, dict):
                d = obj
            elif isinstance(obj, (list, tuple)):
                vec = [float(v) for v in obj]
                vec = (vec + [0.0, 0.0, 0.0, 0.0])[:4]
                ss = sum(vec)
                return [v / ss for v in vec] if ss > 0 else [0.0, 0.0, 0.0, 1.0]
        except Exception:
            pass
        if d is None:
            # try literal
            try:
                import ast as _ast
                obj = _ast.literal_eval(s)
                if isinstance(obj, dict):
                    d = obj
                elif isinstance(obj, (list, tuple)):
                    vec = [float(v) for v in obj]
                    vec = (vec + [0.0, 0.0, 0.0, 0.0])[:4]
                    ss = sum(vec)
                    return [v / ss for v in vec] if ss > 0 else [0.0, 0.0, 0.0, 1.0]
            except Exception:
                d = {"unknown": 1.0}
    vec = [float(d.get(k, 0.0)) for k in ORDER_REDOX]
    ssum = sum(vec)
    return [v / ssum for v in vec] if ssum > 0 else [0.0, 0.0, 0.0, 1.0]

def temp_scale_prob(vec, T=1.5):
    """Temperature scale probabilities in logit space and renormalize."""
    eps = 1e-12
    v = [max(eps, min(1.0 - eps, float(x))) for x in vec]
    logits = [math.log(p / (1.0 - p)) for p in v]
    scaled = [l / max(T, eps) for l in logits]
    probs = [1.0 / (1.0 + math.exp(-l)) for l in scaled]
    ssum = sum(probs)
    return [p / ssum for p in probs] if ssum > 0 else [0.0, 0.0, 0.0, 1.0]

def vec_to_json(vec):
    return json.dumps({k: float(vec[i]) for i, k in enumerate(ORDER_REDOX)})

def load_context_dict(path_hint):
    """Load context dictionary JSON if available; return {} if not."""
    cand = None
    if path_hint:
        cand = Path(path_hint)
    else:
        # fallback to docs/context_env_fields.json if exists
        fallback = Path("docs/context_env_fields.json")
        cand = fallback if fallback.exists() else None
    if cand and cand.exists():
        try:
            return json.loads(cand.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}

def get_prior_vec(ctx, microbiome_type_value):
    """Return a prior vector based on microbiome type and ctx['priors']."""
    try:
        mt = str(microbiome_type_value).strip()
    except Exception:
        mt = ""
    priors = ctx.get("priors", {})
    if mt and mt in priors and isinstance(priors[mt], dict):
        d = priors[mt].get("redox_prob", {})
        if isinstance(d, dict) and d:
            vec = [float(d.get(k, 0.0)) for k in ORDER_REDOX]
            s = sum(vec)
            if s > 0:
                return [v / s for v in vec]
    # default unknown
    return [0.0, 0.0, 0.0, 1.0]

def apply_conflicts(df, ctx):
    """Set env_conflict = 1 if any rule matches; else 0. Non-destructive."""
    rules = ctx.get("conflicts", [])
    if not isinstance(rules, list) or not rules:
        if "env_conflict" not in df.columns:
            df["env_conflict"] = 0
        return df
    # Start with zeros; if already exists, OR into it
    if "env_conflict" not in df.columns:
        df["env_conflict"] = 0
    for rule in rules:
        cond = rule.get("if", {})
        if not isinstance(cond, dict) or not cond:
            continue
        # Build a boolean mask that all specified columns equal the given values (as strings)
        mask = pd.Series([True] * len(df))
        for col, val in cond.items():
            sval = None if val is None else str(val)
            if col in df.columns:
                mask = mask & (df[col].astype(str) == ("" if sval is None else sval))
            else:
                mask = mask & False
        # Set conflicts
        df.loc[mask, "env_conflict"] = 1
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="fin", required=True)
    ap.add_argument("--out", dest="fout", required=True)
    ap.add_argument("--dict", dest="ctx_json", default=None, help="context dict JSON (default tries docs/context_env_fields.json)")
    ap.add_argument("--env", dest="legacy_env", default=None, help="legacy alias for --dict")
    ap.add_argument("--T", type=float, default=1.5)
    args = ap.parse_args()

    df = pd.read_csv(args.fin)

    # Load context dictionary if provided or available
    ctx_path = args.ctx_json or args.legacy_env
    ctx = load_context_dict(ctx_path)

    # Resolve redox source column
    redox_col = None
    if "redox_prob" in df.columns:
        redox_col = "redox_prob"
    elif "env_redox_prob" in df.columns:
        # Promote legacy column to contract column as well
        df["redox_prob"] = df["env_redox_prob"]
        redox_col = "redox_prob"
    else:
        redox_col = "redox_prob"

    # Fill missing/empty redox values using priors
    if redox_col not in df.columns:
        df[redox_col] = ""

    # Vectorize existing or fill from priors
    mt_col = "microbiome_type" if "microbiome_type" in df.columns else None
    vecs = []
    for i, row in df.iterrows():
        raw = row.get(redox_col, "")
        have = False
        if isinstance(raw, str):
            have = raw.strip() != ""
        else:
            have = raw is not None
        if have:
            vec = parse_prob(raw)
        else:
            mt = row.get(mt_col, "") if mt_col else ""
            vec = get_prior_vec(ctx, mt)
        vecs.append(vec)

    # Temperature scaling and write scaled column
    scaled = [temp_scale_prob(v, T=args.T) for v in vecs]
    df["redox_prob"] = [vec_to_json(v) for v in vecs]
    df["redox_prob_scaled"] = [vec_to_json(v) for v in scaled]

    # Apply conflicts if any
    df = apply_conflicts(df, ctx)

    Path(args.fout).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.fout, index=False)
    print(f"[OK] wrote {args.fout}, rows={len(df)}")

if __name__ == "__main__":
    main()
