# -*- coding: utf-8 -*-
# Build an extended candidate library (~60 rows) for demo runs.
import csv, random
from pathlib import Path

random.seed(42)

SEEDS = [
    # acids / aldehydes / alcohols / simple carbonyls
    ("CC(=O)O", "phaseII_esterification", "acetic acid", 0, -2.0),
    ("O=CC(=O)O", "phaseI_oxidation", "glyoxylic acid", 0, -5.0),
    ("OCC(=O)O", "phaseI_oxidation", "glycolic acid", 0, -4.0),
    ("C=O", "phaseI_oxidation", "formaldehyde", 0, -3.0),
    ("O=CO", "phaseI_oxidation", "formic acid", 0, -4.5),
    ("CO", "phaseI_reduction", "methanol", 0, -0.8),
    ("CCO", "phaseI_reduction", "ethanol", 0, -1.0),
    ("O=CCO", "phaseI_oxidation", "glycolaldehyde", 0, -2.5),
    ("O=CCC(=O)O", "phaseI_oxidation", "3-oxopropionic acid", 0, -1.8),
    ("CC(=O)Cl", "byproduct", "acetyl chloride", 0, 2.0),
    ("CC(=O)OC", "phaseII_esterification", "methyl acetate", 0, -1.2),
    ("CC(=O)OCC", "phaseII_esterification", "ethyl acetate", 0, -1.1),
    ("CCC(=O)O", "phaseII_esterification", "propionic acid", 0, -1.6),
    ("CCCC(=O)O", "phaseII_esterification", "butyric acid", 0, -1.6),
    ("O=CC(=O)OC", "phaseII_esterification", "methyl glyoxylate", 0, -1.5),
    ("O=CC(=O)OCC", "phaseII_esterification", "ethyl glyoxylate", 0, -1.4),
    ("O=COC", "phaseII_esterification", "methyl formate", 0, -1.3),
    ("O=COCC", "phaseII_esterification", "ethyl formate", 0, -1.2),
    ("CC(=O)CO", "phaseI_reduction", "hydroxyacetone", 0, -1.0),
    ("O=CCO", "phaseI_oxidation", "glycolaldehyde", 0, -2.5),
    # phenyl & simple aromatics
    ("c1ccccc1", "phaseI_oxidation", "benzene", 0, 1.0),
    ("c1ccccc1O", "phaseI_hydroxylation", "phenol", 0, -1.5),
    ("CCc1ccccc1", "phaseI_oxidation", "ethylbenzene", 0, -0.3),
    ("CC(=O)c1ccccc1", "phaseI_oxidation", "acetophenone", 0, -0.6),
    ("OCc1ccccc1", "phaseI_hydroxylation", "benzyl alcohol", 0, -0.8),
    ("O=CCc1ccccc1", "phaseI_oxidation", "phenylacetaldehyde", 0, -0.7),
    ("O=COc1ccccc1", "phaseII_esterification", "phenyl formate", 0, -0.5),
    ("CC(=O)Oc1ccccc1", "phaseII_esterification", "phenyl acetate", 0, -0.4),
    # small diols / acids
    ("OCCO", "phaseI_hydroxylation", "ethylene glycol", 0, -1.0),
    ("OCC(O)CO", "phaseI_hydroxylation", "glycerol", 0, -1.2),
    ("OC(=O)CO", "phaseI_oxidation", "glycolic acid isomer", 0, -3.8),
    ("OOC=O", "byproduct", "peroxyformic-like", 0, 1.5),
    # halogenated (limited)
    ("ClCC(=O)O", "byproduct", "chloroacetic acid", 0, 0.8),
    ("ClC=O", "byproduct", "formyl chloride", 0, 2.0),
    ("Clc1ccccc1", "phaseI_halogenation", "chlorobenzene", 0, 1.2),
    ("Clc1ccccc1O", "phaseI_halogenation", "p-chlorophenol (isomer-free)", 0, 0.8),
]

# Expand with simple transforms to reach ~60:
MORE = [
    ("O=CC(=O)O", "phaseI_oxidation", "glyoxylic acid (dup-weight)", 0, -5.0),
    ("O=CO", "phaseI_oxidation", "formic acid (dup-weight)", 0, -4.5),
    ("C=O", "phaseI_oxidation", "formaldehyde (dup-weight)", 0, -3.0),
    ("CO", "phaseI_reduction", "methanol (dup-weight)", 0, -0.8),
    ("CCO", "phaseI_reduction", "ethanol (dup-weight)", 0, -1.0),
    ("OCC(=O)O", "phaseI_oxidation", "glycolic acid (dup-weight)", 0, -4.0),
    ("O=COC(=O)O", "phaseII_esterification", "carbonic monoester", 0, -0.5),
    ("O=CC(=O)OC", "phaseII_esterification", "methyl glyoxylate (dup)", 0, -1.5),
    ("CC(=O)O", "phaseII_esterification", "acetic acid (dup)", 0, -2.0),
    ("O=COC", "phaseII_esterification", "methyl formate (dup)", 0, -1.3),
    ("O=COCC", "phaseII_esterification", "ethyl formate (dup)", 0, -1.2),
    ("CC(=O)OC", "phaseII_esterification", "methyl acetate (dup)", 0, -1.2),
    ("CC(=O)OCC", "phaseII_esterification", "ethyl acetate (dup)", 0, -1.1),
    ("CCC(=O)O", "phaseII_esterification", "propionic acid (dup)", 0, -1.6),
    ("CCCC(=O)O", "phaseII_esterification", "butyric acid (dup)", 0, -1.6),
    ("OCCO", "phaseI_hydroxylation", "ethylene glycol (dup)", 0, -1.0),
    ("OCC(O)CO", "phaseI_hydroxylation", "glycerol (dup)", 0, -1.2),
    ("OC(=O)CO", "phaseI_oxidation", "glycolic acid isomer (dup)", 0, -3.8),
    ("ClCC(=O)O", "byproduct", "chloroacetic acid (dup)", 0, 0.8),
    ("Clc1ccccc1O", "phaseI_halogenation", "chlorophenol (dup)", 0, 0.8),
    # a few small extras
    ("O=COC(=O)OC", "phaseII_esterification", "dimethyl carbonate", 0, -0.4),
    ("O=COC(=O)OCC", "phaseII_esterification", "methyl ethyl carbonate", 0, -0.3),
    ("CCOC(=O)OCC", "phaseII_esterification", "diethyl carbonate", 0, -0.2),
    ("CC(O)C=O", "phaseI_oxidation", "hydroxypropanal", 0, -1.0),
    ("O=CC(O)CO", "phaseI_oxidation", "glyceraldehyde", 0, -1.4),
    ("OC(=O)C(O)O", "phaseI_oxidation", "tartronic-like", 0, -2.0)
]

rows = SEEDS + MORE
# pad if <60 by duplicating low-energy species
while len(rows) < 60:
    rows.append(random.choice(SEEDS))

Path("data/external").mkdir(parents=True, exist_ok=True)
out = "data/external/candidates_demo_60.csv"
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["product_smiles","rule_tag","name","reversible","dg_proxy"])
    for smi, tag, name, rev, dg in rows[:60]:
        w.writerow([smi, tag, name, rev, dg])
print("[OK] wrote", out, "rows=", 60)
