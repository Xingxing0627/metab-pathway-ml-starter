# -*- coding: utf-8 -*-
import argparse, json, math, time
from pathlib import Path
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.DataStructs import TanimotoSimilarity

RDLogger.DisableLog("rdApp.info")
RDLogger.DisableLog("rdApp.warning")

# ----- Fingerprints / similarity -----
_GEN = {}
def _gen(sz=1024, r=2):
    k = (sz, r); g = _GEN.get(k)
    if g is None:
        g = GetMorganGenerator(radius=r, fpSize=sz)
        _GEN[k] = g
    return g

def fp(s):
    m = Chem.MolFromSmiles(s)
    return None if m is None else _gen().GetFingerprint(m)

def tanimoto(a, b):
    fa, fb = fp(a), fp(b)
    return 0.0 if (fa is None or fb is None) else float(TanimotoSimilarity(fa, fb))

def base_score(prev_smi, cand_smi, step):
    return max(1e-8, tanimoto(prev_smi, cand_smi) * (1.0 / (1.0 + 0.2 * step)))

# ----- Element book-keeping -----
BALANCE_ELEMS = ("C", "H", "O", "Cl")  # include Cl in count-level balance

def ecount(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return {}
    d = {}
    for a in m.GetAtoms():
        z = a.GetSymbol()
        d[z] = d.get(z, 0) + 1
    return d

def eset(smi): 
    return set(ecount(smi).keys())

def norm_smiles(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m else s

def allowed_mol_vecs(names):
    lib = {
        "H2O": {"H": 2, "O": 1},
        "CO2": {"C": 1, "O": 2},
        "[H]Cl": {"H": 1, "Cl": 1},
        "[H][H]": {"H": 2},
        "ClCl": {"Cl": 2},
    }
    out = []
    for n in names:
        n = n.strip()
        if not n: 
            continue
        if n in lib:
            out.append(lib[n])
            continue
        m = Chem.MolFromSmiles(n)
        if m is not None:
            out.append(ecount(n))
    return out

def nonneg_combo(diff, vecs):
    # exact representation on C/H/O/Cl by nonnegative integer combo of vecs
    if not vecs:
        return all(v == 0 for v in diff.values())
    keys = set(k for v in vecs for k in v.keys())
    consider = [k for k in BALANCE_ELEMS if (k in diff) or (k in keys)]
    for k in consider:
        if diff.get(k, 0) > 0 and all(v.get(k, 0) == 0 for v in vecs):
            return False
    lim = 32
    V = vecs[:3]
    L = len(V)
    if L == 1:
        v = V[0]
        for n0 in range(lim):
            if all(n0 * v.get(k, 0) == diff.get(k, 0) for k in consider):
                return True
        return False
    if L == 2:
        v0, v1 = V
        for n0 in range(lim):
            for n1 in range(lim):
                if all(n0 * v0.get(k, 0) + n1 * v1.get(k, 0) == diff.get(k, 0) for k in consider):
                    return True
        return False
    v0, v1, v2 = V[0], V[1], V[2]
    for n0 in range(lim // 2):
        for n1 in range(lim // 2):
            for n2 in range(lim // 2):
                if all(n0 * v0.get(k, 0) + n1 * v1.get(k, 0) + n2 * v2.get(k, 0) == diff.get(k, 0) for k in consider):
                    return True
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--start_smiles", default="CC(=O)O")
    ap.add_argument("--allow_escape_elems", default="H,O")  # Cl handled by molecule balance
    ap.add_argument("--allow_escape_molecules", default="H2O,CO2,[H]Cl,[H][H],ClCl")
    ap.add_argument("--k0", type=int, default=10)
    ap.add_argument("--k_decay", type=int, default=2)
    ap.add_argument("--max_depth", type=int, default=3)
    ap.add_argument("--max_nodes", type=int, default=2000)
    ap.add_argument("--max_time_s", type=int, default=30)
    ap.add_argument("--no_loop", action="store_true", default=True)
    args = ap.parse_args()

    cands = pd.read_csv(args.candidates)
    if "product_smiles" not in cands.columns:
        raise ValueError("candidates CSV must contain 'product_smiles'")

    # rule_tag weights (tune freely)
    rule_w = {
        "phaseI_oxidation": 1.30,
        "phaseI_hydroxylation": 1.25,
        "phaseII_esterification": 1.10,
        "phaseI_reduction": 1.00,
        "phaseI_halogenation": 0.60,
        "byproduct": 0.80,
    }
    cands["__rule_w"] = cands.get("rule_tag", pd.Series(["phaseI_reduction"] * len(cands))).map(rule_w).fillna(1.0)

    allow_elems = set(x.strip() for x in args.allow_escape_elems.split(",") if x.strip())
    esc_names = [x.strip() for x in args.allow_escape_molecules.split(",") if x.strip()]
    esc_vecs = allowed_mol_vecs(esc_names)

    cand_smi = cands["product_smiles"].astype(str).tolist()
    cand_set = [eset(s) for s in cand_smi]

    start_set = eset(args.start_smiles)
    start_norm = norm_smiles(args.start_smiles)

    t0 = time.time(); expanded = 0
    beams = [[{"path": [], "logp": 0.0, "visited": {start_norm}}]]
    kept = [1]

    for depth in range(args.max_depth):
        k = max(1, args.k0 - depth * args.k_decay)
        curr = beams[-1]; pool = []
        for n in curr:
            used = set(n["path"]) if args.no_loop else set()
            prev = cand_smi[n["path"][-1]] if n["path"] else args.start_smiles
            prev_c = ecount(prev)
            prev_norm = norm_smiles(prev)

            for j in range(len(cand_smi)):
                if j in used: 
                    continue
                cs = cand_smi[j]; cs_norm = norm_smiles(cs)
                if cs_norm in n["visited"]:
                    continue
                # element set constraint
                if cand_set[j] - start_set - allow_elems:
                    continue
                # monotonicity for non-whitelist elements
                prod_c = ecount(cs)
                ok = True
                for el, v in prev_c.items():
                    if el not in allow_elems and prod_c.get(el, 0) > v:
                        ok = False; break
                if not ok: 
                    continue
                # count-level conservation via escaping molecules on C/H/O/Cl
                diff = {}
                for el in set(list(prev_c.keys()) + list(prod_c.keys()) + list(BALANCE_ELEMS)):
                    d = prev_c.get(el, 0) - prod_c.get(el, 0)
                    if d < 0 and el not in allow_elems:
                        ok = False; break
                    diff[el] = max(0, d) if el in BALANCE_ELEMS else 0
                if not ok:
                    continue
                if not nonneg_combo(diff, esc_vecs):
                    continue

                sc = base_score(prev, cs, depth) * float(cands.loc[j, "__rule_w"])
                new_visited = set(n["visited"]); new_visited.add(cs_norm)
                pool.append({"path": n["path"] + [j], "logp": n["logp"] + math.log(sc), "visited": new_visited})

        pool.sort(key=lambda x: x["logp"], reverse=True)
        beams.append(pool[:k]); kept.append(len(beams[-1])); expanded += len(pool)
        if time.time() - t0 > args.max_time_s: print("[INFO] stop by time budget"); break
        if expanded > args.max_nodes: print("[INFO] stop by node budget"); break
        if not beams[-1]: print("[INFO] beam exhausted at depth", depth); break

    last_non_empty = next((b for b in reversed(beams) if b), beams[0])
    best = max(last_non_empty, key=lambda x: x["logp"])

    out = {
        "start_smiles": args.start_smiles,
        "allow_escape_elems": sorted(list(allow_elems)),
        "allow_escape_molecules": esc_names,
        "best_path_indices": best["path"],
        "best_path_smiles": [cand_smi[i] for i in best["path"]] if best["path"] else [],
        "best_path_score": best["logp"],
        "beams_per_depth": kept,
        "stats": {
            "expanded": expanded, "levels": len(beams) - 1,
            "k0": args.k0, "k_decay": args.k_decay,
            "max_depth": args.max_depth, "max_nodes": args.max_nodes,
            "max_time_s": args.max_time_s, "no_loop": bool(args.no_loop)
        }
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"[OK] wrote {args.out}  best_len={len(best['path'])}  kept={kept}")

if __name__ == "__main__":
    main()
