# opus-2026-07-19-S394 -- HYP-7800: THE FREIMAN / DOUBLING ANGLE ON THE GAP.
#
# klein's THM-1004/1005 are titled "the INVERSE/rigidity theorem" -- and an
# inverse theorem for additive structure is FREIMAN's theorem: small doubling
# |A+A| <= K|A| forces A into a generalised AP.  The LRC extremals ARE APs
# ({1,...,13}) or near-APs ({1,...,11,13,24}), and an AP has the MINIMUM
# possible doubling |A+A| = 2|A|-1.
#
# THE FREIMAN-SHAPED CONJECTURE for the two-element extremal set:
#     LARGE DOUBLING  =>  M BOUNDED AWAY FROM THE FLOOR 1/14.
# Equivalently: M near 1/14 forces small doubling, i.e. near-AP structure.
# That is exactly the shape klein proved at n=13 by hand (radius <= 2); a
# doubling hypothesis would replace Hamming radius by an ADDITIVE measure and
# so cover ALL perturbations at once, not just bounded-radius ones.
from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def M(V):
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b)
        if a!=b: D.add(abs(a-b))
    for v in V: D.add(2*v)
    best=F(0)
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best=g
    return best
def doubling(V):
    S={a+b for a in V for b in V}
    return F(len(S), len(V))            # |V+V| / |V|
FLOOR=F(1,14); NEXT=F(1,13)

print("(1) DOUBLING OF THE EXTREMALS AND NEAR-EXTREMALS")
print("    minimum possible for a 13-set: |V+V| = 2*13-1 = 25, ratio 25/13 = 1.923")
print("    family                 |V+V|  ratio    M          at floor?")
for nm,V in [("{1,...,13}",list(range(1,14))),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("{1,...,12,14}",list(range(1,13))+[14]),
             ("{1..11,13,25}",[1,2,3,4,5,6,7,8,9,10,11,13,25]),
             ("{2,...,14}",list(range(2,15))),
             ("AP d=3",[1+3*i for i in range(13)]),
             ("odd {1,3,..,25}",[2*i+1 for i in range(13)])]:
    S={a+b for a in V for b in V}
    m=M(V)
    print(f"    {nm:20s} {len(S):5d}  {float(doubling(V)):.3f}   {str(m):9s}  "
          f"{'YES' if m==FLOOR else ''}")

print()
print("(2) DOES LARGE DOUBLING FORCE M AWAY FROM THE FLOOR?")
print("    binned over random primitive 13-families")
random.seed(394)
bins={}
for _ in range(90):
    V=sorted(random.sample(range(1,60),13))
    if reduce(gcd,V)!=1: continue
    d=float(doubling(V)); m=M(V)
    k=round(d*2)/2
    bins.setdefault(k,[]).append(float(m))
print("    doubling  n    min M      median M    all >= 1/13?")
for k in sorted(bins):
    b=sorted(bins[k])
    if len(b)<3: continue
    print(f"    {k:7.1f} {len(b):4d}  {b[0]:.6f}   {b[len(b)//2]:.6f}    "
          f"{'yes' if b[0]>=float(NEXT) else 'NO'}")

print()
print("(3) THE SHARP TEST: is there a doubling threshold K with")
print("    doubling(V) > K  =>  M(V) >= 1/13 (i.e. off the floor and past the gap)?")
allpts=[(float(doubling(V)),float(M(V))) for V in
        [sorted(random.sample(range(1,60),13)) for _ in range(60)]]
viol=[(d,m) for d,m in allpts if m<float(NEXT)]
print(f"    families with M < 1/13: {len(viol)} of {len(allpts)}")
if viol:
    print(f"      their doublings: {sorted(round(d,3) for d,_ in viol)[:8]}")
    print(f"      max doubling among them: {max(d for d,_ in viol):.3f}")
    print(f"      => any valid threshold K must exceed {max(d for d,_ in viol):.3f}")
else:
    print("      none -- every sampled family already clears 1/13")
print()
print("    extremal doublings for comparison: "
      f"{float(doubling(list(range(1,14)))):.3f} (AP), "
      f"{float(doubling([1,2,3,4,5,6,7,8,9,10,11,13,24])):.3f} (second tight)")
