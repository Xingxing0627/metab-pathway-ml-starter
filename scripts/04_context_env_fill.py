import argparse, json, ast, pandas as pd
from pathlib import Path

def load_prob(s):
    if isinstance(s, str):
        try: return ast.literal_eval(s)
        except Exception: return None
    return s

def softmax(x):
    import math
    ex = [math.exp(v) for v in x]
    s = sum(ex)
    return [v/s for v in ex]

def temp_scale_prob(p, T=1.5):
    # numerically stable via log->/T->softmax
    import math
    eps=1e-12
    p=[max(eps, min(1-eps, v)) for v in p]
    logit=[math.log(v/(1-v)) for v in p]
    scaled=[lv/T for lv in logit]
    # convert back to prob via softmax-like normalizer
    mx=max(scaled)
    ex=[math.exp(v-mx) for v in scaled]
    s=sum(ex)
    return [v/s for v in ex]

def entropy(p):
    import math
    eps=1e-12
    return -sum(max(eps,v)*math.log(max(eps,v)) for v in p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="fin", required=True)
    ap.add_argument("--out", dest="fout", required=True)
    ap.add_argument("--env", dest="fenv", required=True)
    ap.add_argument("--T", type=float, default=1.5)
    args = ap.parse_args()

    env = json.load(open(args.fenv, "r"))
    df = pd.read_csv(args.fin)

    # Normalize and enrich env distributions
    redox_levels = env["redox_levels"]
    probs_scaled = []
    ent = []
    for p in df["env_redox_prob"].apply(load_prob):
        if p is None:
            probs_scaled.append([1/len(redox_levels)]*len(redox_levels))
            ent.append(0.0)
        else:
            ps = temp_scale_prob(p, T=args.T)
            probs_scaled.append(ps)
            ent.append(entropy(ps))
    df[[f"redox_{lv}" for lv in redox_levels]] = pd.DataFrame(probs_scaled, index=df.index)
    df["redox_entropy"] = ent

    df.to_csv(args.fout, index=False)
    print(f"[OK] wrote {args.fout}, rows={len(df)}")

if __name__ == "__main__":
    main()
