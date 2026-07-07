import numpy as np
from itertools import combinations
from math import gcd
from functools import reduce
import sys

GRID=20000
xs=(np.arange(GRID)+0.5)/GRID
def mu_17(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    allg=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
    return float((allg.max(axis=1)>1/7).mean())
def norm_shape(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]; g=reduce(gcd,E[1:]) if len(E)>1 else 1
    return tuple(e//g for e in E) if g and g>0 else tuple(E)

print("=== CENSUS of mu_17(E) over cluster shapes (affine classes), opus-S131 ===\n")
for k,Smax in [(8,13),(9,12),(10,11)]:
    ap=tuple(range(k)); mu_ap=mu_17(list(ap))
    seen={}
    for S in range(k-1, Smax+1):
        for interior in combinations(range(1,S), k-2):
            sh=norm_shape((0,)+interior+(S,))
            if len(sh)==k and sh not in seen:
                seen[sh]=mu_17(list(sh))
    mn=min(seen.values())
    minimizers=[sh for sh,v in seen.items() if v<mn+2e-4]
    below=[(sh,v) for sh,v in seen.items() if v<mu_ap-2e-4]
    order=sorted(seen.items(), key=lambda kv:kv[1])
    print(f"k={k}: {len(seen)} affine shapes (spread<={Smax}). mu_17(AP {ap})={mu_ap:.4f}")
    print(f"   global min={mn:.4f}; AP is min? {abs(mn-mu_ap)<2e-4}; #minimizers={len(minimizers)} e.g.{minimizers[:3]}")
    print(f"   #strictly below AP: {len(below)}  {[(s,round(v,4)) for s,v in below[:3]]}")
    print(f"   5 smallest shapes: {[(s,round(v,4)) for s,v in order[:5]]}", flush=True)
    print()
