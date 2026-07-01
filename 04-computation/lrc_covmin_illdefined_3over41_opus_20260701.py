"""What IS the covering-min at n=14? Is it 2/27 (keep-core) or is there a covering set with M in (1/14, 2/27)?
e.g. 3/41=[0;13,1,2]=0.0732. Search covering 13-sets for the smallest M>1/14."""
import numpy as np, random
from fractions import Fraction as F
from math import gcd
def M_exact(S,Q=130):
    Sa=np.array(S); bn,bd=0,1
    for q in range(2,Q+1):
        A=np.outer(Sa,np.arange(1,q))%q; kk=int(np.minimum(A,q-A).min(axis=0).max())
        if kk*bd>bn*q: bn,bd=kk,q
    g=gcd(bn,bd) or 1; return bn//g,bd//g
def covering(S): return all(any(s%q==0 for s in S) for q in range(2,14))
floor=F(1,14); kc=F(2,27)
best=(1,1); bestS=None; rng=random.Random(0); tested=0
# structured: drop d core, add multiples (targeted) + random covering sets
for trial in range(400000):
    d=rng.choice([0,1,1,2,2,3])
    core=list(range(1,13))
    drop=rng.sample(core,d); keep=[c for c in core if c not in drop]
    # add 1+d elements from a pool (multiples-rich) to cover q=13 and dropped q's
    pool=list(range(13,46))
    add=rng.sample(pool,1+d)
    S=sorted(set(keep+add))
    if len(S)!=13 or not covering(S): continue
    tested+=1
    bn,bd=M_exact(S); v=F(bn,bd)
    if v>floor and v<F(best[0],best[1]):
        best=(bn,bd); bestS=S
g=gcd(*best) or 1
print(f"n=14: searched {tested} covering sets; smallest M>1/14 found = {best[0]//g}/{best[1]//g} = {best[0]/best[1]:.6f}")
print(f"   set = {bestS}")
print(f"   keep-core 2/27={float(kc):.6f}; floor 1/14={float(floor):.6f}; gap candidate 3/41={3/41:.6f}")
print(f"   => {'2/(2n-1) is the min (no smaller found)' if F(*best)>=kc else 'SMALLER than 2/27 exists -- covering-min < 2/(2n-1)!'}")
