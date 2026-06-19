#!/usr/bin/env python3
"""
Finish claims (3) k=10..13 hunt + (4) R_{1/7}, leanly (bounded spread, fewer randoms).
EXACT.  Prints incrementally and flushes.
"""
from fractions import Fraction as F
from itertools import combinations
import random, sys

src=open('04-computation/lrc14_Bk_verify_gpintersectionun_kps-S5-wf.py').read().split('if __name__')[0]
ns={}; exec(src,ns)
mu_exact=ns['mu_exact']; safe_GP=ns['safe_GP']; meas=ns['meas']
rho_star=ns['rho_star']; good_E=ns['good_E']; intersect_intervals=ns['intersect_intervals']
th17=F(1,7)

def flush(*a):
    print(*a); sys.stdout.flush()

# thr_k
def min_GP(sz):
    best=None
    for Pset in combinations(range(1,14), sz):
        m=meas(safe_GP(list(Pset)))
        if best is None or m<best: best=m
    return best
thr={k:1-min_GP(13-k) for k in range(8,14)}
flush("thr_k:", {k:str(thr[k]) for k in thr})

random.seed(99)
flush("=== k=10..13 NON-consecutive hunt (bounded spread <= 3k, structured + random) ===")
viol=0
for k in range(10,14):
    cons=list(range(k)); mu_cons=mu_exact(cons,th17)
    best=mu_cons; bestE=cons; cnt=0
    shapes=[cons]
    # perforated 0..m drop (m up to k+4)
    for span in range(k,k+5):
        for drop in combinations(range(1,span), span-k):
            E=[0]+sorted(set(range(1,span+1))-set(drop))
            E=sorted(set(E))
            if len(E)>=k: shapes.append(E[:k] if E[0]==0 else [0]+E[:k-1])
    # near-AP perturb
    for d in [1,2,3]:
        base=[d*i for i in range(k)]
        for idx in range(1,k):
            for delta in [-1,1,-2,2]:
                E=sorted(set(base[:idx]+[base[idx]+delta]+base[idx+1:]))
                if 0 in E and len(E)>=k: shapes.append(E[:k])
    # two-block sub-torus
    for a in range(2,k-1):
        b=k-a
        for M in [k,k+1,k+3,5,7,11]:
            E=sorted(set(list(range(a))+[M+i for i in range(b)]))
            if 0 in E and len(E)>=k: shapes.append(E[:k])
    # random bounded spread <= 3k
    for _ in range(1500):
        cap=random.choice([k+1,k+3,2*k,3*k])
        rest=random.sample(range(1,cap+1), k-1)
        shapes.append([0]+sorted(rest))
    for E in shapes:
        E=sorted(set(E))
        if len(E)<k or E[0]!=0: continue
        m=mu_exact(E[:k],th17); cnt+=1
        if m<best: best=m; bestE=E[:k]
        if m<thr[k]:
            viol+=1
            if viol<=5: flush(f"   VIOLATION k={k} E={E[:k]} mu={m} < thr={thr[k]}")
    flush(f"k={k}: {cnt} shapes; min mu_1/7={best}={float(best):.5f} at {bestE}; "
          f"consec={mu_cons}; consec_is_argmin={best>=mu_cons}; thr={thr[k]}={float(thr[k]):.5f}; min>=thr={best>=thr[k]}")
flush(f"TOTAL thr violations: {viol}")

flush("=== CLAIM (4): R_1/7 = rho*/(meas(G_P)*mu) over consecutive E, k=8..13 ===")
Rmin=None; Rarg=None
for k in range(8,14):
    sz=13-k; E=list(range(k)); muc=mu_exact(E,th17)
    GE=good_E(E,th17)
    for Pset in combinations(range(1,14), sz):
        GP=safe_GP(list(Pset)); gp=meas(GP)
        rho=meas(intersect_intervals(GP,GE))
        denom=gp*muc
        if denom>0:
            R=rho/denom
            if Rmin is None or R<Rmin: Rmin=R; Rarg=(k,Pset)
flush(f"min R_1/7={Rmin}={float(Rmin):.6f} at {Rarg}; claimed >= 67053/84241={float(F(67053,84241)):.6f}; "
      f"Rmin>=67053/84241? {Rmin>=F(67053,84241)}; exact match to claim? {Rmin==F(67053,84241)}")
flush("DONE")
