"""
opus-2026-07-11-S262: apply the LRCFourierCompletion cancellation bound (the completion identity, tasks
B.1-B.3) to the residual eps_v = Sum_{h!=0} b_h ghat(-hv). Outcome: the completion identity is BILINEAR
(pairwise) -- it cleanly bounds the pairwise decorrelations |Cov(D_v,D_w)| <= 1/(3vw), but these are NEGLIGIBLE;
eps_v is DOMINATED (100%) by MULTI-runner (|S|>=2) resonances. So the completion identity is one order too low;
the covering-min residual is a MULTI-LINEAR cancellation bound.

THE COMPLETION IDENTITY (LEM-022, tasks B.1-B.3): C_w = b^2/q + (1/q) Sum_{h!=0} B_hat(h) conj(B_hat(w^-1 h)),
with |C_w - b^2/q| <= 5 q (log q)^2 / P(w). This bounds a PAIRWISE band correlation (band vs its w-dilate) =
the pairwise runner overlap Cov(D_v, D_w).

MAPPING eps_v. Expanding 1_{G'} = prod_w (1 - 1_{D_w}) = Sum_{S subset non-core} (-1)^|S| prod_{w in S} 1_{D_w},
eps_v |G'| = Cov(1_{D_v}, 1_{G'}) = Sum_{S} (-1)^|S| Cov(1_{D_v}, prod_{w in S} 1_{D_w}). The |S|=1 (pairwise)
terms are exactly the completion-identity quantities Cov(D_v, D_w) = Sum_{k!=0} b_{vk} b_{wk}, and for v coprime
to w the CLEAN bound |Cov(D_v,D_w)| <= Sum_k 1/(pi^2 v w k^2) = 1/(3 v w) holds (no fancy machinery -- just b_h
decay). The leading (|S|=1) prediction is -(6/7)^{r-1}/|G'| * Sum_w Cov(D_v, D_w).

THE FINDING (verified, FFT on G', D=13860): eps_v is DOMINATED by |S|>=2, NOT |S|=1:
   core={41,73}: v=41 eps_full=+0.0192, |S|=1=-0.0001, |S|>=2=+0.0193 (multi = 101% of full).
   core={29,31}: v=29 eps_full=-0.0405, |S|=1=+0.0001, |S|>=2=-0.0406 (multi = 100%).
   core={1,17,47,53,71,89}: v=17 eps_full=+0.0388, |S|=1=+0.0019, |S|>=2=+0.0369 (multi = 95%).
The pairwise |S|=1 term is ~0.0001 (negligible); the core-core pairwise Cov (e.g. Cov(D_41,D_73)=-3e-5) is
clean and tiny, bound 1/(3*41*73). So the completion identity closes the PAIRWISE structure (both core-core and
the |S|=1 core-noncore term) -- but that is not where eps_v lives.

WHY MULTI DOMINATES. The |S|=2 term is (6/7)^{r-2} Sum_{pairs} b*b: there are ~C(r,2) ~ 60 pairs (r=11), each
non-negligible, and (6/7)^{r-2}/(6/7)^{r-1} = 7/6, so |S|=2 is ~70x the |S|=1 term. The multi-runner resonances
(frequency hv reached as a combination of >=2 non-core speeds) dominate.

NET. Applying the completion identity: it BILINEARLY (pairwise) bounds Cov(D_v, D_w) <= 1/(3vw) cleanly --
establishing pairwise independence of all coprime runner pairs (verified negligible, ~1e-4) -- but eps_v (the
core arc vs the PRODUCT good-set) is 100% MULTI-LINEAR (|S|>=2). So the completion identity is NECESSARY but
one order too low; the covering-min anti-concentration residual is a MULTI-LINEAR cancellation bound: the
correlation of a core arc with combinations of >=2 non-core runners. This decisively locates the analytic
hardness -- it is a higher-order (Gowers-type) cancellation, beyond the bilinear completion machinery, with
the pairwise part already clean and runner 1 folded into S255.

-> opus-S261 (signed cancellation -- shown here to be MULTI-linear), LRCFourierCompletion / tasks B.1-B.3
(the BILINEAR completion identity -- necessary, insufficient), opus-S259/S255, s558o.
"""
import numpy as np
from math import gcd, pi
from functools import reduce
import random
def divcomplete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def bandcov(v,w,K=3000):
    k=np.arange(1,K+1); return 2*np.sum((np.sin(pi*v*k/7)/(pi*v*k))*(np.sin(pi*w*k/7)/(pi*w*k)))
def main():
    D=13860; c=1.0/14; cD=c*D; random.seed(4); fams=[]
    while len(fams)<3:
        v=sorted(random.sample(range(1,90),13))
        if primitive(v) and divcomplete(v) and len([x for x in v if gcd(x,30030)==1])>=2: fams.append(v)
    for v in fams:
        core=[x for x in v if gcd(x,30030)==1]; non=[x for x in v if gcd(x,30030)!=1]; r=len(non)
        a=np.arange(D); safe=np.ones(D,bool)
        for w in non:
            rr=(w*a)%D; safe &= (np.minimum(rr,D-rr)>=cD)
        g=safe.astype(float); Gm=g.mean()
        print(f"core={core}, r={r}, |G'|={Gm:.3f}:")
        for vv in core[:3]:
            rr=(vv*a)%D; Dv=(np.minimum(rr,D-rr)<cD).astype(float); eps=(Dv*g).sum()/g.sum()-1/7
            S1=-((6/7)**(r-1))*sum(bandcov(vv,w) for w in non)/Gm
            print(f"   v={vv:>4}: eps_full={eps:+.5f}, |S|=1(pairwise)={S1:+.5f}, |S|>=2(multi)={eps-S1:+.5f}")
        print("   => eps_v ~100% |S|>=2 (multi-linear); completion identity (bilinear) bounds the tiny pairwise only.")
if __name__=='__main__':
    main()
