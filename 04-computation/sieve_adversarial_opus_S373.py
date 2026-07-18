# opus-2026-07-17-S373 -- IS THE DENOMINATOR BOUNDED ON *PRIMITIVE* FAMILIES?
#
# FIRST, THE OBSTRUCTION I MUST RULE OUT (dilation has bitten three times:
# MISTAKE-154, THM-1055, MISTAKE-156).  For k*{1,...,13} the lonely set of
# {1,...,13} is the single point 1/14, so lonely t satisfy k t = 1/14 mod 1,
# i.e. t = (1+14m)/(14k): minimal denominator 14k, UNBOUNDED in k.  So there is
# NO uniform Q over all families.  Dilation invariance (THM-1050) reduces to
# PRIMITIVE families, and the real question is whether Q is bounded there.
# This run hunts adversarially for a primitive family needing a large one.
from math import gcd
from functools import reduce
import random

def lonely_at(V,p,q):
    for v in V:
        r=(v*p)%q
        if min(r,q-r)*14 < q: return False
    return True
def min_den(V,Qmax=600):
    for q in range(1,Qmax+1):
        for p in range(1,q):
            if gcd(p,q)!=1: continue
            if lonely_at(V,p,q): return q
    return None
def prim(V): return reduce(gcd,V)==1

print("(4) CONFIRM THE DILATION OBSTRUCTION (why we must restrict to primitive)")
for k in [1,2,3,5,7]:
    V=[k*i for i in range(1,14)]
    print(f"    {k}*{{1..13}}  gcd={reduce(gcd,V):3d}   min denominator = {min_den(V,800)}")

print()
print("(5) ADVERSARIAL HUNT: maximise the minimal denominator over PRIMITIVE families")
random.seed(3731)
best=(0,None)
for trial in range(60):
    V=sorted(random.sample(range(1,300),13))
    if not prim(V): continue
    cur=min_den(V) or 999
    for step in range(120):
        W=list(V); i=random.randrange(13)
        W[i]=max(1,W[i]+random.choice([-3,-2,-1,1,2,3]))
        W=sorted(set(W))
        if len(W)!=13 or not prim(W): continue
        d=min_den(W) or 999
        if d>cur: V,cur=W,d
    if cur>best[0]: best=(cur,list(V))
print(f"    worst primitive family found: min denominator = {best[0]}")
print(f"      V = {best[1]}")

print()
print("(6) DISTRIBUTION OF THE MINIMAL DENOMINATOR (primitive, 600 families)")
from collections import Counter
C=Counter()
for _ in range(600):
    V=sorted(random.sample(range(1,500),13))
    if not prim(V): continue
    C[min_den(V) or 999]+=1
for q in sorted(C): print(f"      q = {q:3d}: {C[q]:4d}")
print(f"    MAXIMUM over the sample: {max(C)}")
