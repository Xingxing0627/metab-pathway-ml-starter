# -*- coding: utf-8 -*-
# Minimal trainable baseline: [X_mol | X_env] -> LogisticRegression (binary T1)
import argparse, yaml, numpy as np, joblib
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.model_selection import StratifiedKFold

def safe_auc(y_true, y_prob):
    y_true = np.asarray(y_true)
    if len(np.unique(y_true)) < 2:
        return float("nan")
    try:
        return float(roc_auc_score(y_true, y_prob))
    except Exception:
        return float("nan")

def train_eval_small(X, y):
    n = len(y)
    # Very small N: leave-one-out AUC; fallback to train acc if class collapses
    if n <= 8:
        probs = np.zeros(n, dtype=float)
        for i in range(n):
            idx_tr = [j for j in range(n) if j != i]
            idx_va = [i]
            if len(np.unique(y[idx_tr])) < 2:
                clf = LogisticRegression(max_iter=200, n_jobs=1)
                clf.fit(X, y)
                return clf, float("nan"), float(accuracy_score(y, clf.predict(X)))
            clf = LogisticRegression(max_iter=200, n_jobs=1)
            clf.fit(X[idx_tr], y[idx_tr])
            probs[i] = clf.predict_proba(X[idx_va])[:,1][0]
        auc = safe_auc(y, probs)
        clf = LogisticRegression(max_iter=200, n_jobs=1)
        clf.fit(X, y)
        acc = accuracy_score(y, clf.predict(X))
        return clf, auc, float(acc)
    else:
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        val_probs = np.zeros(n, dtype=float)
        for tr, va in skf.split(X, y):
            if len(np.unique(y[tr])) < 2:
                continue
            clf = LogisticRegression(max_iter=500, n_jobs=1)
            clf.fit(X[tr], y[tr])
            val_probs[va] = clf.predict_proba(X[va])[:,1]
        auc = safe_auc(y, val_probs)
        clf = LogisticRegression(max_iter=500, n_jobs=1)
        clf.fit(X, y)
        acc = accuracy_score(y, clf.predict(X))
        return clf, auc, float(acc)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--pack", default="data/interim/pack_v0.npz")
    ap.add_argument("--model_out", default="data/processed/lite_clf.joblib")
    ap.add_argument("--metrics_out", default="reports/metrics/lite_metrics.json")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config, "r"))
    print("[INFO] project:", cfg["project_name"])

    pack = np.load(args.pack)
    X = np.hstack([pack["X_mol"], pack["X_env"]])   # (N, d_mol + d_env)
    y = pack["y_occur"].astype(int)

    clf, auc, acc = train_eval_small(X, y)
    print(f"[OK] Lite baseline trained. CV-ish AUC={auc}, train Acc={acc}")

    Path(args.model_out).parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(clf, args.model_out)

    Path(args.metrics_out).parent.mkdir(parents=True, exist_ok=True)
    import json
    with open(args.metrics_out, "w", encoding="utf-8") as f:
        json.dump({"auc": auc, "train_acc": acc, "n": int(len(y))}, f, ensure_ascii=False, indent=2)

    print("[OK] saved model ->", args.model_out)
    print("[OK] saved metrics ->", args.metrics_out)

if __name__ == "__main__":
    main()
