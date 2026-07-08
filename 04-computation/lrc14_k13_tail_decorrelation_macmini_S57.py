"""mac-mini-S57 (HYP-5297/THM-657): k=13 tail (diam>=76) safe via mu>=(7/6)E[W], E[W]->iid."""
import numpy as np
from math import gcd
from functools import reduce
import random, collections
random.seed(1013); TH=1/7
def EW_mu(E,x):
    Ea=np.array(sorted(E),float);ph=np.mod(np.outer(x,Ea),1.0);ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return float(np.maximum(g-TH,0).sum(axis=1).mean()), float((g.max(axis=1)>TH).mean())
GRID=800_000;x=(np.arange(GRID)+0.5)/GRID; mP=0.0565; iid=(6/7)**13
res=collections.defaultdict(list)
for _ in range(400):
    D=random.choice([76,90,120,200,400,1000]); E=sorted(random.sample(range(D+1),13)); E=[e-E[0] for e in E]
    g=reduce(gcd,[E[i+1]-E[i] for i in range(12)]);d=max(E)//max(g,1)
    if d<76: continue
    res[d if d<76 else D].append(EW_mu(E,x))
for D in sorted(res):
    ews=[r[0] for r in res[D]]
    print(f"diam~{D}: n={len(res[D])} minE[W]={min(ews):.4f} mean={np.mean(ews):.4f} (7/6)min/mP={7/6*min(ews)/mP:.1f}x")
print(f"iid E[W]=(6/7)^13={iid:.4f}; tail SAFE if min E[W]>= {mP*6/7:.4f}")
