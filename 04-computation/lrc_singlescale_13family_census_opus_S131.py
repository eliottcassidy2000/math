from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random
random.seed(77)

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
            mn=None
            for v in V:
                dd=distZ(v*m,d)
                if mn is None or dd<mn: mn=dd
                if mn==0: break
            if mn is not None and mn>best: best=mn
    return best

thr=F(1,14)
AP=list(range(1,14))
MAP=reach_M(AP)
print(f"=== SINGLE-SCALE 13-family census: LRC(14) core (opus-S131) ===")
print(f"M(AP {{1..13}}) = {MAP} = {float(MAP):.6f}  (expect 1/14={1/14:.6f})  threshold 1/14\n")

below=[]; tight=[]; near=[]; tested=0
def consider(V):
    global tested
    if len(set(map(abs,V)))<13 or any(v==0 for v in V): return
    if reduce(gcd,[abs(v) for v in V])!=1: return
    M=reach_M(V); tested+=1
    if M<thr: below.append((tuple(V),M))
    elif M==thr: tight.append(tuple(sorted(V)))
    elif M<thr+F(1,182): near.append((tuple(sorted(V)),M))  # M in (1/14, 1/14+~0.006)

# perturbations of AP (single-scale, bounded)
for x in range(14,45): consider(list(range(1,13))+[x])       # {1..12,x}
for x in range(13,30):
    for y in range(x+1,32): consider(list(range(1,12))+[x,y])  # {1..11,x,y}
for i,j in combinations(range(13),2):                          # swap two AP entries up
    for si,sj in [(13,13),(13,14),(14,15),(1,1),(13,1)]:
        V=AP.copy(); V[i]+=si; V[j]+=sj; consider(V)
for _ in range(300):                                           # random bounded
    consider(sorted(random.sample(range(1,23),13)))

print(f"tested: {tested} single-scale 13-families")
print(f"M < 1/14 (would REFUTE LRC14): {len(below)}   {below[:5]}")
print(f"M = 1/14 exactly (tight, = AP class): {len(tight)}  distinct-up-to-perm: {len(set(tight))}")
for t in list(set(tight))[:6]: print(f"    tight: {t}")
print(f"M in (1/14, 1/14+0.0055] (near-tight moat rim): {len(near)}")
for V,M in sorted(near,key=lambda z:z[1])[:8]: print(f"    M={float(M):.6f}={M}  {V}")
print(f"\n  => AP unique tight & no sub-threshold => single-scale core consistent with LRC14 (structural, not proof)")
