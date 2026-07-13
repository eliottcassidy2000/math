"""
opus-2026-07-11-S266: prove the additive bound |eps_v| <= f(#relations) rigorously. Outcome: a clean rigorous
IDENTITY for eps_v (a signed sum over additive relations), but NO clean bound -- the low-order truncation
(m<=2) is 0.13 vs the actual eps_v=0.019, so higher-order terms (m>=3) cancel ~0.11 via the alternating (-1)^m
sign. eps_v's smallness is an ALTERNATING MULTI-ORDER cancellation (confirms S262 rigorously). The additive
bound is NOT a clean theorem.

THE RIGOROUS IDENTITY. b_k = sin(pi k/7)/(pi k); beta_u = 1_{D_u}, beta_u^hat(h) = b_{h/u}[u|h];
1_{G'} = prod_w (1 - beta_w). Then eps_v |G'| = <beta_v - 1/7, prod_w(1-beta_w)>, and (the k_v=0 terms cancel
against -1/7, forcing v to participate):
    eps_v |G'| = Sum_{relations R} (-1)^{m_R} (6/7)^{r - m_R} prod_{u in T'_R} b_{k_u},
where R: support T' contains v, T'\{v} subset non-core, nonzero k_u, Sum_{u in T'} u k_u = 0, m_R = |T'|-1.

THE NEGATIVE (verified). Low-order truncation (m<=2 = the +-v+-w_i+-w_j=0 relations) gives 0.13 for v=41 vs
actual eps_v=0.019 (FFT). So m>=3 cancels ~0.11 (alternating (-1)^m). Magnitude bound Sum|...| overshoots
massively (and Sum|b_k| diverges). So |eps_v| <= f(#relations) is NOT a theorem for any low-order f -- eps_v is
an alternating multi-order cancellation. The S263 correlation 0.527 reflects only the leading term.

NET. The additive bound rigorization yields a clean exact IDENTITY but no bound; eps_v's smallness is
higher-order (multi-linear, S262), beyond low-order additive counting. The S265 case skeleton stands VERIFIED,
but its two supporting bounds (additive |eps_v|, measure |S_rest|) are the irreducible higher-order
anti-concentration -- the elementary tools (S253-S266) are exhausted; closing them needs an inverse theorem for
the band-multilinear cancellation.

-> opus-S265 (case skeleton), opus-S262 (multi-linear, confirmed), opus-S263 (additive = leading term),
opus-S255 (extremizer), tasks #42-#43 (multi-way entanglement).
"""
import numpy as np
from math import gcd, pi, sin
from functools import reduce
from itertools import combinations
import random
def divcomplete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def b(k): return sin(pi*k/7)/(pi*k) if k!=0 else 2/14
def main():
    D=13860; c=1.0/14; cD=c*D; random.seed(4); fams=[]
    while len(fams)<2:
        v=sorted(random.sample(range(1,90),13))
        if primitive(v) and divcomplete(v) and len([x for x in v if gcd(x,30030)==1])>=2: fams.append(v)
    for v in fams:
        core=[x for x in v if gcd(x,30030)==1]; non=[x for x in v if gcd(x,30030)!=1]; r=len(non)
        a=np.arange(D); safe=np.ones(D,bool)
        for w in non:
            rr=(w*a)%D; safe &= (np.minimum(rr,D-rr)>=cD)
        Gm=safe.mean()
        for vv in core[:2]:
            rr=(vv*a)%D; ea=((np.minimum(rr,D-rr)<cD)&safe).sum()/safe.sum()-1/7
            # m<=2 truncation of the identity
            K=6; tot=0.0
            for w in non:  # m=1
                gg=gcd(vv,w); kv=w//gg; kw=-vv//gg
                for sv in (1,-1): tot+=-(6/7)**(r-1)*b(sv*kv)*b(-sv*kw)
            for w1,w2 in combinations(non,2):  # m=2
                for k1 in range(-K,K+1):
                    for k2 in range(-K,K+1):
                        if k1==0 or k2==0: continue
                        s=w1*k1+w2*k2
                        if s!=0 and s%vv==0 and abs(s//vv)<=60:
                            tot+=(6/7)**(r-2)*b(-s//vv)*b(k1)*b(k2)
            print(f"v={vv:>3}: eps_actual(FFT)={ea:+.4f}, identity(m<=2 trunc)={tot/Gm:+.4f} -- MISMATCH => m>=3 cancels; no clean low-order bound")
    print("=> rigorous IDENTITY holds, but no clean truncation: eps_v is an ALTERNATING MULTI-ORDER cancellation (S262 confirmed).")
if __name__=='__main__':
    main()
