#!/usr/bin/env python3
"""The bounded-core dual is degree <=4 (solvable) and IS the bimodality functional.
Worpitzky bridge: miss-PGF (deg 6, value/p_t basis) <-> dual g (deg<=4, binomial/S_r basis).
L_y = sum_t g(t) p_t; consec maximizes it = three-gap bimodality. kind-mendel-2026-06-27-S10."""
from fractions import Fraction as F
from math import gcd, floor, comb
from functools import reduce
from itertools import combinations
def gl(xs): return reduce(lambda a,b:a*b//gcd(a,b),[x for x in xs if x],1)
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def miss_p(E):
    pos=[e for e in E if e]; D=7*gl(pos); bset=set([0,D])
    for e in pos:
        step=D//(7*e); x=0
        while x<=D: bset.add(x); x+=step
    bps=sorted(bset); p=[F(0)]*7
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid2=a+b; hit={(7*((e*mid2)%(2*D)))//(2*D) for e in E}
        N=len([j for j in range(1,7) if j not in hit]); p[N]+=F(b-a,D)
    return p
# THM-534 duals g(t) as functions of t in {0..6}
duals={8: lambda t:F((t-1)*(t-2)*(t-4)*(t-5),40),
       9: lambda t:F(-(t-2)*(t-3)*(t-6),36), 10: lambda t:F(-(t-2)*(t-3)*(t-6),36),
       11:lambda t:F((t-3)*(t-4),12),12:lambda t:F((t-3)*(t-4),12),13:lambda t:F((t-3)*(t-4),12)}
deg={8:4,9:3,10:3,11:2,12:2,13:2}
caps={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
print("=== dual g(t) values on {0..6} (degree <=4 = SOLVABLE; note bimodal: mass at t=0,6) ===")
for k in [8,9,10,11]:
    g=duals[k]; vals=[g(t) for t in range(7)]
    print(f"  k={k} (deg {deg[k]}): g = {[str(v) for v in vals]}  -> L_y = sum g(t) p_t")
print()
print("=== L_y = sum_t g(t) p_t for consec_k; verify <= cap_k AND consec maximizes ===")
for k in [8,9,10]:
    g=duals[k]; cap=caps[k]
    p=miss_p(list(range(k))); Ly=sum(g(t)*p[t] for t in range(7))
    # search bounded-spread E for max L_y
    best=(F(-1),None)
    for sub in combinations(range(1,17),k-1):
        E=[0]+list(sub)
        if gall(E)!=1: continue
        pe=miss_p(E); l=sum(g(t)*pe[t] for t in range(7))
        if l>best[0]: best=(l,E)
    print(f"  k={k}: L_y(consec)={float(Ly):.5f}  cap={float(cap):.5f} (<=cap:{Ly<=cap})  "
          f"max L_y over bounded={float(best[0]):.5f} at {'consec' if best[1]==list(range(k)) else best[1]} (consec=max:{best[1]==list(range(k))})")
print()
print("=== k=8 dual IS the bimodality functional: L_y = p_0 + p_6 + (1/10) p_3 ===")
g=duals[8]; p=miss_p(list(range(8)))
print(f"  g(0)={g(0)}, g(6)={g(6)}, g(3)={g(3)}, others={[str(g(t)) for t in (1,2,4,5)]}")
print(f"  L_y(consec_8) = p0+p6+p3/10 = {p[0]}+{p[6]}+{p[3]}/10 = {p[0]+p[6]+p[3]/10} = {float(p[0]+p[6]+p[3]/10):.5f}")
print(f"  bimodality(p0+p6)={float(p[0]+p[6]):.5f}; degree-4 dual rewards EXTREMES => consec (three-gap, most bimodal) wins")
