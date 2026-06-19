#!/usr/bin/env python3
"""
Focused adversarial hunt: does CONSECUTIVE truly minimize mu_{1/7}(E) for k=8..13,
and is min mu_{1/7} >= thr_k always?  This is the SINGLE surviving gap in the
gp-intersection-uniform angle's k>=8 branch.

We attack with the shape classes most likely to LOWER mu (push points toward a
shifted 1/q-grid so the max circular gap stays small):
  - consecutive 0..k-1
  - perforated-AP (drop one element from 0..k)
  - near-AP with small step d, gcd reduced (scale-invariance => gcd=1 WLOG)
  - sub-torus / common-factor blocks
  - "two-AP" unions
  - random bounded & large spread
mu_{1/7} computed EXACTLY.
"""
from fractions import Fraction as F
from itertools import combinations
import random

src=open('04-computation/lrc14_Bk_verify_gpintersectionun_kps-S5-wf.py').read().split('if __name__')[0]
ns={}; exec(src,ns)
mu_exact=ns['mu_exact']; safe_GP=ns['safe_GP']; meas=ns['meas']; per=None

th17=F(1,7)

# thr_k from min meas(G_P) over |P|=13-k (recompute exactly)
def min_GP(sz, speeds=range(1,14)):
    best=None
    for Pset in combinations(speeds, sz):
        m=meas(safe_GP(list(Pset)))
        if best is None or m<best: best=m
    return best

thr={}
for k in range(8,14):
    thr[k]=1-min_GP(13-k)
print("thr_k:", {k:str(thr[k]) for k in thr})

def gcd_reduce(E):
    from math import gcd
    g=0
    for e in E: g=gcd(g,e)
    if g>1: E=[e//g for e in E]
    return sorted(set(E))

random.seed(20260618)
overall={}
for k in range(8,14):
    cons=list(range(k))
    mu_cons=mu_exact(cons,th17)
    best=mu_cons; bestE=cons; cnt=0
    # structured families
    fams=[]
    fams.append(cons)
    # perforated: drop j from 0..m for m=k..k+3
    for span in range(k, k+4):
        for drop in combinations(range(1,span), span-k):
            E=[0]+[i for i in range(1,span+1) if i not in drop][:k-1]
            E=[0]+sorted(set(range(1,span+1))-set(drop))
            E=[0]+sorted(list(set(range(1,span+1))-set(drop)))
            if len(set(E))>=k:
                fams.append(gcd_reduce(sorted(set(E))[:k]))
    # APs scaled (gcd_reduce -> consecutive, so skip) ; near-AP perturb one point
    for d in [1,2,3]:
        base=[d*i for i in range(k)]
        for idx in range(1,k):
            for delta in [-1,1,-2,2]:
                E=base[:]; E[idx]+=delta
                E=gcd_reduce(sorted(set(E)))
                if len(E)>=k: fams.append(E[:k])
    # two-block (sub-torus): {0..a-1} U {M..M+b-1}
    for a in range(2,k-1):
        b=k-a
        for M in [k, k+1, k+3, 2*k, 5, 7]:
            E=list(range(a))+[M+i for i in range(b)]
            E=gcd_reduce(sorted(set(E)))
            if len(E)>=k: fams.append(E[:k])
    # random bounded + large spread
    for _ in range(6000):
        cap=random.choice([k+1,k+3,k+8,2*k,4*k,8*k])
        if cap < k-1: continue
        rest=random.sample(range(1,cap+1), k-1)
        E=gcd_reduce([0]+sorted(rest))
        if len(E)>=k: fams.append(E)
    for E in fams:
        if len(E)<k: continue
        m=mu_exact(E,th17); cnt+=1
        if m<best:
            best=m; bestE=E
    viol = best < thr[k]
    cons_is_min = (best==mu_cons) or (best>=mu_cons)
    print(f"k={k}: mu_cons={mu_cons}={float(mu_cons):.5f}  min_found={best}={float(best):.5f} at {bestE}  "
          f"[{cnt} shapes]  consec_is_argmin={best>=mu_cons}  thr={thr[k]}={float(thr[k]):.5f}  "
          f"min>=thr={best>=thr[k]}  VIOL={viol}")
    overall[k]=(mu_cons,best,bestE,thr[k])

print("\nSUMMARY")
allok=True
for k in range(8,14):
    mc,bf,be,t=overall[k]
    if bf<mc: print(f"  !! k={k}: found shape BELOW consecutive: {be} mu={bf} < cons {mc}"); allok=False
    if bf<t: print(f"  !! k={k}: found shape BELOW thr: {be} mu={bf} < thr {t}"); allok=False
print("ALL OK (consec minimal & >=thr in this search):", allok)
