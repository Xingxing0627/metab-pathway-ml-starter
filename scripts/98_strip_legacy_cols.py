# -*- coding: utf-8 -*-
import pandas as pd, sys
p=sys.argv[1]
df=pd.read_csv(p)
drop=[c for c in ("env_redox_prob",) if c in df.columns]
if drop: df.drop(columns=drop).to_csv(p,index=False)
