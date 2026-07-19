from math import gcd
from functools import reduce
from itertools import combinations
import random, time
# gap at p/q exceeds num/den  <=>  min(vp mod q, q-vp mod q)*den > num*q  for all v
NUM,DEN=3,41          # threshold 3/41
def crit(V):
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b); d=abs(a-b)
        if d>1: D.add(d)
    for v in V: D.add(2*v)
    return sorted(q for q in D if q>=2)
def beats(V):
    """True iff some critical p/q has gap > 3/41 (pure integer, early exit)."""
    for q in crit(V):
        thr=NUM*q
        for p in range(1,q):
            ok=True
            for v in V:
                r=(v*p)%q
                m=r if r<q-r else q-r
                if m*DEN<=thr: ok=False; break
            if ok: return True
    return False
print("FAST INTEGER SCAN: primitive families with M <= 3/41")
print("  (known: two extremals at 1/14, and {1..11,13,36} at 3/41)")
print()
# sanity: the three known families must NOT beat 3/41
for nm,V in [("{1..13}",list(range(1,14))),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("{1..11,13,36}",[1,2,3,4,5,6,7,8,9,10,11,13,36]),
             ("{1..12,14}",list(range(1,13))+[14])]:
    print(f"    sanity {nm:16s} beats 3/41? {beats(V)}   (first three should be False)")
print()
random.seed(3982)
for N in [20,30,45,70,110]:
    surv=[]; n=0; t0=time.time()
    while time.time()-t0 < 60:
        V=sorted(random.sample(range(1,N+1),13))
        if reduce(gcd,V)!=1: continue
        n+=1
        if not beats(V): surv.append(list(V))
    print(f"  cap N={N:3d}: {n:7d} primitive families, {len(surv)} with M <= 3/41  [{time.time()-t0:.0f}s]")
    for V in surv[:2]: print(f"      {V}")
