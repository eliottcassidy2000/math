import numpy as np
from itertools import combinations
from math import gcd
from functools import reduce
GRID=100000; xs=(np.arange(GRID)+0.5)/GRID
def Emaxgap(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
    return float(g.max(axis=1).mean())
def norm(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]; g=reduce(gcd,E[1:]) if len(E)>1 else 1
    return tuple(e//g for e in E) if g else tuple(E)
print("=== is E[maxgap] AP-minimized? census (opus-S133) ===")
for k,Smax in [(8,13),(9,12),(10,11)]:
    ap=Emaxgap(list(range(k))); seen={}
    for S in range(k-1,Smax+1):
        for it in combinations(range(1,S),k-2):
            sh=norm((0,)+it+(S,))
            if len(sh)==k and sh not in seen: seen[sh]=Emaxgap(list(sh))
    mn=min(seen.values()); below=[s for s,v in seen.items() if v<ap-3e-4]
    print(f"  k={k}: {len(seen)} shapes. E[maxgap(AP)]={ap:.4f}; global min={mn:.4f}; #below AP={len(below)}  {below[:2]}")
print("  => AP-minimality of E[maxgap] holds if #below AP = 0 (like mu_17). Then E[maxgap(E)]>=E[maxgap(AP)]>1/7.")
