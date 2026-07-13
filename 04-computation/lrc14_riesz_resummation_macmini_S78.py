#!/usr/bin/env python3
"""mac-mini-S78: the RIESZ RESUMMATION structure -- Schur deficit => L(S)>0 (= covering LRC14).
Rewrite L exactly WITHOUT the divergent E_k truncation. Safe indicator f_i=1-1_{D_i}, mean 6/7.
Set g_i = f_i/(6/7)-1 (mean 0, g_i=-1 on danger, +1/6 on safe, |g_i|<=1). Then
   L = (6/7)^k * INT prod_i (1+g_i) = (6/7)^k * [1 + Sum_{|T|>=2} corr_T],  corr_T=INT prod_T g_i.
This is a FINITE (2^k-term) BOUNDED sum -- no conditional-convergence issue. L>0 <=> Sum_{|T|>=2}
corr_T > -1. The AP drives it to -1 (perfect anti-correlation, L->0); covering (Schur deficit)
should keep it > -1. TEST the split: (a) dissociated/loose families: corr small, L~(6/7)^k>0
(the LOOSE case, correlation-bounded); (b) near-AP: corr->-1, L->0 (the residual = knife-edge)."""
c=1.0/14; m=6.0/7  # safe-mean
def L_and_corrsum(S,res=400000):
    k=len(S); L=0
    for j in range(res):
        t=(j+0.5)/res
        if all(min((v*t)%1,1-((v*t)%1))>=c for v in S): L+=1
    L/=res
    corrsum = L/(m**k)-1   # = Sum_{|T|>=2} corr_T
    return L, corrsum, (m**k)
def Schur(A):
    Aset=set(A); return sum(1 for a in A for b in A if a+b in Aset)
from math import comb
print(f"independent value (6/7)^13 = {m**13:.5f}. L>0 <=> corrsum := L/(6/7)^13 - 1 > -1.\n")
print(f"{'set':26s} | Schur | deficit | L        | corrsum(>-1?) | regime")
print("-"*88)
fams={
 "AP {1..13}": list(range(1,14)),
 "{1..11,13,84}": [*range(1,12),13,84],
 "deep well {1..12,182}": [*range(1,13),182],
 "{2..14}": list(range(2,15)),
 "spread(coprime-ish)": [3,17,29,41,58,73,91,112,140,171,199,233,281],
 "dissociated {2^i}": [1,2,4,8,16,32,64,128,256,512,1024,2048,4096],
}
for nm,S in fams.items():
    S=sorted(set(S)); 
    if len(S)!=13: continue
    L,cs,indep=L_and_corrsum(S); T=Schur(S); dfc=comb(13,2)-T
    reg="AP/knife-edge" if cs<-0.9 else ("near-AP" if cs<-0.5 else "LOOSE (decorr)")
    print(f"{nm:26s} | {T:5d} | {dfc:7d} | {L:.6f} | {cs:+.4f}      | {reg}")
print("\nSTRUCTURE: L=(6/7)^13*(1+corrsum), a FINITE bounded sum (no divergence). L>0 <=> corrsum>-1.")
print("LOOSE (dissociated, large Schur deficit): corrsum near 0 => L~(6/7)^13>0 -- correlation-bounded,")
print("the tractable half (= kps E_grid/kissing, opus dissociated). NEAR-AP (small deficit): corrsum")
print("-> -1 => L->0 -- the knife-edge residual = the open covering LRC(14). The resummation SPLITS")
print("exactly as the covering-min does; the near-AP corrsum>-1 is the irreducible open statement.")
