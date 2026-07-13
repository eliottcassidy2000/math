#!/usr/bin/env python3
"""mac-mini-S84: pursue the DIRECT METRIC RESUMMATION of the middle-order cancellation.
Hypothesis: the +-20 middle-order terms in L=Sum_k(-1)^k E_k are the TRIVIAL binomial cancellation
of the all-dangerous tail. E_k=Sum_x C(x,k) p_x, so L=Sum_k(-1)^k E_k = Sum_x p_x Sum_k(-1)^k C(x,k)
= Sum_x p_x (1-1)^x = p_0 (each x>=1 contributes 0 BINOMIALLY). So the 'resummation wall' (S79) is
an ARTIFACT of the moment expansion; the direct object is L=p_0 (safe measure). TEST this, and
confirm the direct-L Fourier = the resonance sum (corrsum) = loose/near-AP split (no new lever)."""
import numpy as np
from math import comb
c=1.0/14
def pdist(S,res=600000):
    S=sorted(set(S)); k=len(S)
    ts=(np.arange(res)+0.5)/res; X=np.zeros(res,dtype=int)
    for v in S:
        r=(v*ts)%1.0; d=np.minimum(r,1-r); X+=(d<c).astype(int)
    return np.array([(X==x).mean() for x in range(k+1)])
print("(1) the middle-order +-20 = binomial cancellation of the all-dangerous tail:\n")
for nm,S in [("AP {1..13}",list(range(1,14))),("{1..11,13,84}",[*range(1,12),13,84])]:
    p=pdist(S); k=len(p)-1
    # E_k = sum_x C(x,k) p_x
    E=[sum(comb(x,j)*p[x] for x in range(k+1)) for j in range(k+1)]
    L_moment=sum((-1)**j*E[j] for j in range(k+1))
    # per-x contribution to L: p_x * sum_j (-1)^j C(x,j) = p_x*(1-1)^x = 0 for x>=1, p_0 for x=0
    perx=[p[x]*sum((-1)**j*comb(x,j) for j in range(x+1)) for x in range(k+1)]
    print(f"{nm}: E_6={E[6]:.1f} E_7={E[7]:.1f} (the '+-20') | L(moment sum)={L_moment:.5f} = p_0={p[0]:.5f}")
    print(f"   per-danger-level x contribution to L (= p_x*(1-1)^x): {[round(v,4) for v in perx]}")
    print(f"   => only x=0 contributes (p_0); every x>=1 cancels BINOMIALLY. The +-20 at k=6,7 are")
    print(f"      C(13,6)*p_13={comb(13,6)*p[13]:.1f}, C(13,7)*p_13={comb(13,7)*p[13]:.1f} -- the tail, cancelling.\n")
print("(2) VERDICT: the middle-order cancellation is TRIVIAL (binomial), an artifact of the moment")
print("expansion. There is NO hidden analytic 'resummation' -- the direct metric object is L=p_0,")
print("the safe measure. Its Fourier is the resonance sum corrsum (HYP-6430): loose (few resonances,")
print("provable density floor) + near-AP (many resonances, the residue). So the 'direct metric")
print("resummation' returns the RAW conjecture L=p_0>0 for covering = LRC(14) -- no new structure.")
print("Clarification: the S79 'middle-order wall' was partly illusory; the real (and only) content")
print("is that the safe set is nonempty for covering. This door, too, is the conjecture itself.")
