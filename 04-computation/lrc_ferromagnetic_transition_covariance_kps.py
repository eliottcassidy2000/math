"""
lrc_ferromagnetic_transition_covariance_kps.py  (kind-pasteur-2026-06-27-S31aj)

The 1/7 was a near-coincidence (Sigma-k3/S3 crosses ~1/7 at k=8 by the crossing of an
increasing function). The GENUINE structure: the total covariance Sigma-kappa_2 of the
empty-sector indicators changes SIGN -- negative (antiferromagnetic) for k<=5, positive
(ferromagnetic) for k>=6 -- a PHASE TRANSITION at the n=5 compression boundary (mac-mini).

This script:
 (1) the ferromagnetic transition: Sigma-kappa_2(consec_k) for k=3..13, locate the sign change;
 (2) ROBUSTNESS: exhaustive over ALL bounded k=8 clusters (8-subsets of {0..14} with 0=anchor):
     does consec={0..7} maximize Sigma-kappa_2 (covariance)? (the owner's "covariance-max robust?")
 (3) the binding k=8 minimizer {0,1,5,7,8,9,..} check for exact 1/7 (thoroughness).
"""
import sys, itertools
from fractions import Fraction as F
INNER=list(range(1,7))
def sector_of(p): return int((p%1)*7)
def cells_float(E):
    E=sorted(set(E)); b=set([0.0,1.0])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): b.add(mm/(7*e))
    b=sorted(b); out=[]
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1-x0<1e-15: continue
        mid=(x0+x1)/2
        cov=set(int((e*mid)%1*7) for e in E)
        out.append((cov,x1-x0))
    return out
def cov_sum_float(E):
    cl=cells_float(E)
    # p_i = P(sector i empty); P(i,j empty); Sigma Cov = sum_{i<j}[P(ij)-p_i p_j]
    p={i:0.0 for i in INNER}; pij={}
    for cov,w in cl:
        empty=[i for i in INNER if i not in cov]
        for i in empty: p[i]+=w
        for a in range(len(empty)):
            for b2 in range(a+1,len(empty)):
                key=(empty[a],empty[b2]); pij[key]=pij.get(key,0.0)+w
    tot=0.0
    for i,j in itertools.combinations(INNER,2):
        Pij=pij.get((i,j),0.0)
        tot += Pij - p[i]*p[j]
    return tot

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    print("=== (1) FERROMAGNETIC TRANSITION: Sigma-kappa_2(consec_k) ===")
    prev=None
    for k in range(3,14):
        E=tuple(range(k)); s=cov_sum_float(E)
        phase = "ANTIFERRO(<0)" if s<0 else "FERRO(>0)"
        flip = "  <-- SIGN CHANGE" if prev is not None and (prev<0)!=(s<0) else ""
        print(f"  k={k:>2}: Sigma-kappa_2 = {s:+.5f}  {phase}{flip}")
        prev=s
    print("\n=== (2) ROBUSTNESS: exhaustive bounded k=8 -- does consec maximize Sigma-kappa_2? ===")
    consec=tuple(range(8)); cmax=cov_sum_float(consec)
    best=-9; bestE=None; nbeat=0; total=0
    for sub in itertools.combinations(range(1,15),7):   # 7-subsets of {1..14}, plus anchor 0
        E=(0,)+sub; total+=1
        s=cov_sum_float(E)
        if s>best: best=s; bestE=E
        if s>cmax+1e-9: nbeat+=1
    print(f"  consec_8 Sigma-kappa_2 = {cmax:.5f}")
    print(f"  exhaustive over {total} bounded k=8 clusters (span<=14): configs beating consec = {nbeat}")
    print(f"  global argmax = {bestE} val={best:.5f}  (is consec? {bestE==consec})")
    print(f"  --> consec {'MAXIMIZES' if nbeat==0 else 'does NOT maximize'} Sigma-kappa_2 over bounded k=8")
    print("\n=== (3) k=8 minimizer-family exact 1/7 check (exact fractions) ===")
    def cov_exact_ratio_k3(E):
        # exact Sigma-k3/S3 via fractions
        E=sorted(set(E)); b=set([F(0),F(1)])
        for e in E:
            if e==0: continue
            for mm in range(0,7*e+1): b.add(F(mm,7*e))
        b=sorted(b); cl=[]
        for i in range(len(b)-1):
            x0,x1=b[i],b[i+1]
            if x1<=x0: continue
            cov=set(sector_of(e*((x0+x1)/2)) for e in E)
            cl.append((cov,x1-x0))
        def Pe(S):
            S=set(S); t=F(0)
            for cov,w in cl:
                if S.isdisjoint(cov): t+=w
            return t
        Sk3=F(0); S3=F(0)
        for tr in itertools.combinations(INNER,3):
            P3=Pe(tr); S3+=P3
            i,j,k=tr
            k3=P3 - Pe((i,j))*Pe((k,)) - Pe((i,k))*Pe((j,)) - Pe((j,k))*Pe((i,)) + 2*Pe((i,))*Pe((j,))*Pe((k,))
            Sk3+=k3
        return Sk3/S3 if S3 else F(0)
    for E in [(0,1,2,3,4,5,6,7),(0,1,5,7,8,9,11,13),(1,5,7,8,9,11,12,13),(0,2,4,6,8,10,12,14)]:
        r=cov_exact_ratio_k3(E)
        print(f"  E={E}: Sk3/S3={r}={float(r):.7f} (=1/7? {r==F(1,7)})")
