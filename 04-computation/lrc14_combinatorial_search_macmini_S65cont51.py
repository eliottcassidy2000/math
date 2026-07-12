#!/usr/bin/env python3
"""cont.51: search for combinatorial simplifications. (1) sector-count N generating function
G(z) = E[z^N] for consec -- does it factor / have a clean form? (2) the moment structure --
is E[N(7-N)] related to a pair-covering count with a closed form? (3) the p_j distribution
across k -- OEIS-searchable pattern?"""
from fractions import Fraction as F
def pdist(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return p
print("(1) sector-count distribution p_j for consec {0..k-1}, k=7..12 (looking for structure):")
for k in range(7,13):
    p=pdist(list(range(k)))
    nz=[(j,p[j]) for j in range(8) if p[j]!=0]
    print(f"  k={k}: " + "  ".join(f"p{j}={p[j]}" for j,pv in nz))
print()
print("(2) common denominators + numerator patterns (OEIS candidates):")
for k in range(7,13):
    p=pdist(list(range(k)))
    from math import gcd
    dens=[pv.denominator for pv in p if pv!=0]
    L=1
    for d in dens: L=L*d//gcd(L,d)
    nums=[int(pv*L) for pv in p]
    print(f"  k={k}: LCD={L}, numerators (N=0..7) = {nums}")
print()
print("(3) E[N] and E[N(7-N)] closed-form check vs iid + the '7-sector coupon' formula:")
from math import comb
for k in range(7,13):
    p=pdist(list(range(k)))
    EN=sum(j*p[j] for j in range(8))
    J=sum(p[j]*j*(7-j) for j in range(8))
    # iid coupon: E[N] = 6*(6/7)^(k-1) (6 nonzero sectors, k-1 nonzero phases)
    EN_iid=6*F(6,7)**(k-1)
    print(f"  k={k}: E[N]={float(EN):.4f} (iid {float(EN_iid):.4f}), J={float(J):.4f}, J/E[N]={float(J/EN):.3f}")
