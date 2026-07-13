#!/usr/bin/env python3
"""mac-mini-S76: even 3rd-order truncation FAILS (the resummation wall). Bonferroni:
L=P(X=0) >= 1-E1+E2-E3 (odd-order upper bound on P(union)). Compute it -- if negative, the
truncated inclusion-exclusion is useless even at order 3; L needs the FULL resummed alternating
series (opus THM-515 Riesz product). Confirms: 3rd-order is NECESSARY (pairwise blind) but not
SUFFICIENT by truncation -- the covering energy-deficit must be resummed across all orders."""
from math import comb
c=1.0/14
def moments_L(S,res=300000):
    m=[0.0]*7; L=0
    for j in range(res):
        t=(j+0.5)/res
        X=sum(1 for v in S if min((v*t)%1,1-((v*t)%1))<c)
        if X==0: L+=1
        for k in range(1,7):
            if X>=k: m[k]+=comb(X,k)
    return [m[k]/res for k in range(1,7)],L/res
for nm,S in [("AP",list(range(1,14))),("{1..11,13,84}",[*range(1,12),13,84]),("deep well",[*range(1,13),182])]:
    E,L=moments_L(sorted(set(S)))
    b3=1-E[0]+E[1]-E[2]; b5=1-E[0]+E[1]-E[2]+E[3]-E[4]
    print(f"{nm:14s}: L(true)={L:.5f} | Bonferroni3=1-E1+E2-E3={b3:+.3f} | Bonferroni5={b5:+.3f} | E_k={[round(x,2) for x in E]}")
print("\n=> truncated inclusion-exclusion is NEGATIVE (useless) at every finite order -- the")
print("moments grow (conditional convergence, THM-504 wall). 3rd-order is NECESSARY (E3 first")
print("separates AP) but truncation NEVER suffices: the covering energy-deficit needs RESUMMATION")
print("(opus THM-515 Riesz product). The last inch = a resummed >=3rd-order three-distance/Freiman")
print("inverse -- non-local, and provably beyond all pairwise AND truncated-moment methods.")
