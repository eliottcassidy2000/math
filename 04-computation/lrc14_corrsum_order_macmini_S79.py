#!/usr/bin/env python3
"""mac-mini-S79: decompose the residue corrsum = Sum_{j>=2} INT e_j(g) by ORDER j (e_j =
elementary symmetric poly of the 13 recentred safe-indicators g_i). corrsum(AP)=-1; which
orders drive it there? If high orders matter, the residue is irreducibly the full resummation
(=> THM-730's Schur/order-3 alone cannot close it). Honest structural diagnostic."""
c=1.0/14
def g_vals(S,t):
    # g_i = -1 on danger (||v t||<c), +1/6 on safe
    return [(-1.0 if min((v*t)%1,1-((v*t)%1))<c else 1.0/6) for v in S]
def esym(x):
    # elementary symmetric e_0..e_k via product expansion
    e=[1.0]+[0.0]*len(x)
    for xi in x:
        for j in range(len(x),0,-1):
            e[j]+=e[j-1]*xi
    return e  # e[0..k]
def corrsum_by_order(S,res=300000):
    k=len(S); acc=[0.0]*(k+1)
    for jj in range(res):
        t=(jj+0.5)/res
        e=esym(g_vals(S,t))
        for j in range(k+1): acc[j]+=e[j]
    return [a/res for a in acc]  # INT e_j, j=0..k
from math import comb
print("INT e_j(g) by order j (corrsum = Sum_{j>=2} INT e_j; L=(6/7)^13*(1+corrsum)):\n")
for nm,S in [("AP {1..13}",list(range(1,14))),("{1..11,13,84}",[*range(1,12),13,84]),("deep well",[*range(1,13),182])]:
    S=sorted(set(S)); E=corrsum_by_order(S)
    cs=sum(E[2:])  # corrsum
    print(f"{nm}: corrsum={cs:+.4f}  (L={ (6/7)**13*(1+cs):.5f})")
    print("   INT e_j: " + "  ".join(f"e{j}={E[j]:+.3f}" for j in range(2,9)))
    print("            " + "  ".join(f"e{j}={E[j]:+.3f}" for j in range(9,14)))
    # cumulative partial sums Sum_{2..J}
    part=[sum(E[2:J+1]) for J in range(2,14)]
    print("   partial corrsum(2..J): " + " ".join(f"{p:+.2f}" for p in part))
print("\nREAD: if the partial sums Sum_{2..J} do NOT settle until J large (all 13 orders matter),")
print("the residue corrsum>-1 is IRREDUCIBLY the full resummation. THM-730 (order-3 Schur) is one")
print("term; the AP's corrsum=-1 is a MULTI-ORDER conspiracy => no bounded-order proof. Honest:")
print("the residue is LRC(14); this diagnoses WHY (high-order cancellation), not a way through.")
