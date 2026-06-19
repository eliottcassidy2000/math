#!/usr/bin/env python3
"""k=10..13 non-consecutive hunt, bounded spread, fast, flushing per-k. EXACT."""
from fractions import Fraction as F
from itertools import combinations
import random, sys
src=open('04-computation/lrc14_Bk_verify_gpintersectionun_kps-S5-wf.py').read().split('if __name__')[0]
ns={}; exec(src,ns)
mu_exact=ns['mu_exact']; safe_GP=ns['safe_GP']; meas=ns['meas']
th17=F(1,7)
def minGP(sz):
    b=None
    for Pset in combinations(range(1,14),sz):
        m=meas(safe_GP(list(Pset)))
        if b is None or m<b: b=m
    return b
thr={k:1-minGP(13-k) for k in range(10,14)}
print("thr:", {k:str(thr[k]) for k in thr}, flush=True)
random.seed(2026)
viol=0
for k in range(10,14):
    cons=list(range(k)); mu_cons=mu_exact(cons,th17)
    best=mu_cons; bestE=cons; cnt=0
    shapes=[cons]
    # perforated 0..m
    for span in range(k,k+5):
        for drop in combinations(range(1,span), span-k):
            E=[0]+sorted(set(range(1,span+1))-set(drop))
            if len(E)>=k and E[0]==0: shapes.append(E[:k])
    # near-AP
    for d in [1,2,3]:
        base=[d*i for i in range(k)]
        for idx in range(1,k):
            for delta in [-1,1,-2,2]:
                E=sorted(set([0]+base[1:idx]+[base[idx]+delta]+base[idx+1:]))
                if 0 in E and len(E)>=k: shapes.append(E[:k])
    # two-block sub-torus
    for a in range(2,k-1):
        b=k-a
        for M in [k,k+1,k+3,5,7,11]:
            E=sorted(set(list(range(a))+[M+i for i in range(b)]))
            if 0 in E and len(E)>=k: shapes.append(E[:k])
    # random bounded spread <= 2k (extremal is bounded-spread; large spread RAISES mu)
    for _ in range(600):
        cap=random.choice([k+1,k+2,k+4,k+8, 2*k])
        rest=random.sample(range(1,cap+1), k-1)
        shapes.append([0]+sorted(rest))
    for E in shapes:
        E=sorted(set(E))
        if len(E)<k or E[0]!=0: continue
        m=mu_exact(E[:k],th17); cnt+=1
        if m<best: best=m; bestE=E[:k]
        if m<thr[k]:
            viol+=1; print(f"  VIOLATION k={k} E={E[:k]} mu={m}<thr={thr[k]}",flush=True)
    print(f"k={k}: {cnt} shapes; min mu_1/7={best}={float(best):.5f} at {bestE}; "
          f"consec={mu_cons}; consec_is_argmin={best>=mu_cons}; "
          f"thr={thr[k]}={float(thr[k]):.5f}; min>=thr={best>=thr[k]}", flush=True)
print(f"TOTAL violations: {viol}", flush=True)
print("DONE", flush=True)
