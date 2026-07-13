#!/usr/bin/env python3
"""mac-mini-S76: THE LAST INCH via ADDITIVE ENERGY at the SL(3) level. opus: L(S)=P(safe)=
Sum_k(-1)^k E_k, E_k = k-fold danger-overlap = |T|=k additive energy; AP maximizes E_3,E_4
(perfect cancellation L=0) but NOT E_2 (Sidon-max, cyclotomic floor). So the AP-blanket/last-inch
lives at |T|>=3, invisible to pairwise methods -- explaining why union bound/2nd moment fail.
Freiman: covering (non-AP) => additive-energy DEFICIT >= order|A|. TEST: for near-AP covering
families, does the TRIPLE-energy deficit E_3(AP)-E_3(S) track the loneliness L(S)>0, and is E_2
flat (Sidon-blind)?"""
from fractions import Fraction as F
c=F(1,14)
def nn_lt(x,thr):  # ||x||<thr ?
    x%=1; return min(x,1-x)<thr
def energies_and_L(S, res=200004):
    # E_k = E[ C(X,k) ], X=#dangerous runners; L=P(X=0). exact-ish via fine grid.
    from math import comb
    cnt=[0]*6; L=0
    cf=float(c)
    for j in range(res):
        t=(j+0.5)/res
        X=sum(1 for v in S if min((v*t)%1,1-((v*t)%1))<cf)
        if X==0: L+=1
        for k in range(1,6):
            if X>=k: cnt[k]+=comb(X,k)
    return [cnt[k]/res for k in range(1,6)], L/res
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
print("family                          | E2    E3    E4   | L=P(safe) | note")
print("-"*82)
AP=list(range(1,14))
fams={
 "AP {1..13} (non-cov)": AP,
 "{1..11,13,84} (cov,12->84)": [*range(1,12),13,84],
 "{1..11,13,168}": [*range(1,12),13,168],
 "deep well {1..12,182}": [*range(1,13),182],
 "{2..14} (cov)": list(range(2,15)),
 "{1..10,12,13,154}(11->154)": [*range(1,11),12,13,154],
 "lacunary {1,2,4,8,..}": [1,2,4,8,16,32,64,128,256,512,1024,2048,4096],
}
E_AP=None
for nm,S in fams.items():
    S=sorted(set(S))
    if len(S)!=13: 
        print(f"  {nm}: not 13"); continue
    E,L=energies_and_L(S)
    if nm.startswith("AP"): E_AP=E
    d3 = (E_AP[1]-E[1]) if E_AP else 0.0
    cov = "COV" if is_cov(S) else "non-cov"
    print(f"  {nm:30s} | {E[0]:.2f}  {E[1]:.2f}  {E[2]:.2f} | {L:.5f} | {cov} E3def={d3:+.2f}")
print("\nKEY CHECKS: (1) E2 ~flat across AP vs covering (Sidon-blind, pair=cyclotomic floor)?")
print("(2) E3 deficit (AP-cov) POSITIVE for covering AND tracks L>0? => last inch = SL(3) energy")
print("deficit, Freiman lever (non-AP => E3 deficit => L>0), pairwise methods provably blind.")
