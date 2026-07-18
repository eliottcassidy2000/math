from fractions import Fraction as F
from math import gcd
import random
random.seed(3)
def covers214(fam): return all(any(v%q==0 for v in fam) for q in range(2,15))
def M_exact(fam,Qcap=None):
    if Qcap is None: Qcap=2*max(fam)+2
    best=F(0)
    for q in range(2,Qcap+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
# min M over 12-sets W that COVER 2..14 -- is it bounded away from 1/13? (a gap => stability closes it)
print("min M over 12-sets covering 2..14 (gap above 1/13 => the 'W covers 14' case closes):")
best=(F(1),None)
found=0
for _ in range(60000):
    Vmax=random.choice([26,28,30,35,40])
    W=sorted(random.sample(range(2,Vmax+1),12))
    if not covers214(W): continue
    found+=1
    M=M_exact(W)
    if M<best[0]: best=(M,W)
M,W=best
print(f"  checked {found} covering 12-sets; MIN M = {M}={float(M):.4f}  (1/13={1/13:.4f}, 1/12={1/12:.4f}) at W={W}")
print(f"  gap above 1/13: {float(M-F(1,13)):.4f}  => {'GAP (delta>=this)' if M>F(1,13) else 'no gap'}")
# and: can adding vmax to such a W drag M(V) below 1/13? (the counterexample route)
print("\n  can a covering 12-set W + vmax give M(V)<1/13? (=> inverse theorem counterexample)")
def M_lt_13(fam):
    Q=2*max(fam)+2
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if 13*m>=q: return False
    return True
cnt=0; hit=0
for _ in range(40000):
    Vmax=random.choice([26,30,40])
    W=sorted(random.sample(range(2,Vmax+1),12))
    if not covers214(W): continue
    cnt+=1
    for vmax in random.sample(range(2,60),8):
        if vmax in W: continue
        V=sorted(W+[vmax])
        if len(V)==13 and M_lt_13(V):
            hit+=1; print(f"    FOUND: W(cov 2..14)={W} + vmax={vmax} => M<1/13! W dilated-AP? {W==[W[0]*i for i in range(1,13)]}")
            break
print(f"  covering-12-set bases tried {cnt}; W-covers-14 families with M(V)<1/13: {hit}")
