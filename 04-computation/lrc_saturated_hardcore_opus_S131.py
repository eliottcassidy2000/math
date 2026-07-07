from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
random.seed(31)

def distZ(num,den):
    r=num%den; return F(min(r,den-r),den)
def reach_M(V):
    V=[abs(v) for v in V if v!=0]; n=len(V); dens=set()
    for i in range(n):
        dens.add(2*V[i])
        for j in range(i+1,n):
            dens.add(V[i]+V[j])
            if V[i]!=V[j]: dens.add(abs(V[i]-V[j]))
    best=F(0)
    for d in dens:
        if d==0: continue
        for m in range(1,d):
            mn=min(distZ(v*m,d) for v in V)
            if mn>best: best=mn
    return best
def saturated(V):  # for each q in 2..14, some speed divisible by q
    return all(any(v%q==0 for v in V) for q in range(2,15))

thr=F(1,14)
print("=== SATURATED 13-family census: the genuine sieve-hard core of LRC(14) (opus-S131) ===")
print("    saturated = a multiple of EVERY q in {2..14} => no single sieve t=1/q works\n")
# construct a base saturated skeleton, fill randomly
base=[14,13,11,9,8,5,12,10]   # covers 2..14 (14->2,7,14; 9->3,9; 8->4,8; 5->5; 12->6,12; 10->10; 13;11)
assert saturated(base), "base not saturated"
below=[]; near=[]; tested=0; minM=F(1); Mdist={}
for _ in range(6000):
    extra=random.sample([x for x in range(1,40) if x not in base], 5)
    V=base+extra
    if len(set(V))<13: continue
    if reduce(gcd,V)!=1: continue
    if not saturated(V): continue
    M=reach_M(V); tested+=1
    if M<minM: minM=M
    b=round(float(M),3); Mdist[b]=Mdist.get(b,0)+1
    if M<thr: below.append((tuple(sorted(V)),M))
    elif M<thr+F(1,120): near.append((tuple(sorted(V)),M))
# also try WIDER saturated families (bigger multiples) -- closer to the true extremals
for _ in range(6000):
    b2=[14*random.randint(1,3),13,11,9*random.randint(1,2),8,5*random.randint(1,3),12,10]
    extra=random.sample(range(1,60),5); V=list(set(b2+extra))
    if len(V)!=13 or not saturated(V) or reduce(gcd,V)!=1: continue
    M=reach_M(V); tested+=1
    if M<minM: minM=M
    if M<thr: below.append((tuple(sorted(V)),M))
    elif M<thr+F(1,120): near.append((tuple(sorted(V)),M))

print(f"tested {tested} SATURATED families")
print(f"  min M found = {minM} = {float(minM):.6f}  (threshold 1/14={1/14:.6f})")
print(f"  M < 1/14 (would REFUTE LRC14): {len(below)}   {[ (v,str(m)) for v,m in below[:5]]}")
print(f"  M in [1/14, 1/14+0.008): {len(near)} near-tight saturated families")
for V,M in sorted(near,key=lambda z:z[1])[:8]: print(f"     M={float(M):.5f}={M}  {V}")
print(f"  M distribution (rounded): {dict(sorted(Mdist.items())[:12])}")
print(f"\n  => saturated families are the sieve-hard core; min M >= 1/14? {minM>=thr}")
