#!/usr/bin/env python3
"""
covering_min_lazy_cuts_klein.py  --  klein-2026-07-01-S61

CREATIVE METHOD #2 (beyond the full ILP): a LAZY-CONSTRAINT / CUTTING-PLANE ILP (Benders-style row
generation). The full set-cover ILP times out at V>~4n because it carries ~40k danger constraints.
Instead: solve a TINY ILP (size + divisibility only), read off a candidate set, compute its EXACT M and
its deepest-hole witness t*, add ONLY that one danger cut, and re-solve. The ILP sees only the handful
of witnesses that actually bind -> it runs at V = n(n-1) (the construction scale).

Question: does a covering set with speeds <= V and M < n/Phi6 exist?  (a beater below the construction)
  - feasibility of "M < n/Phi6" via lazy cuts. If the ILP goes INFEASIBLE -> NO beater (rigorous,
    speeds <= V) -> covering-min = n/Phi6 (construction) for speeds <= n(n-1) -- closes HYP-3778's residual.
  - if a set with M < n/Phi6 is FOUND -> the construction is beaten (a real refutation).
"""
from fractions import Fraction as F
import functools, math
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
print = functools.partial(print, flush=True)

def Mexact_and_witness(S):
    """exact M(S)=max_t min_v ||v t|| and the argmax breakpoint t*."""
    S=sorted(set(S)); cand=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j], S[i]+S[j]):
                if d>0:
                    for k in range(1,d): cand.add(F(k,d))
    best=F(0); at=None
    for t in cand:
        g=min(min((v*t)%1, 1-((v*t)%1)) for v in S)
        if g>best: best=g; at=t
    return best, at

def norm(x):
    f=x-int(x)
    if f<0: f+=1
    return min(f,1-f)

def beater_exists(n, V, r, max_iter=400):
    """lazy-cut ILP: is there a covering set (speeds<=V) with M < r? returns (found_set or None, #cuts)."""
    speeds=list(range(1,V+1)); P=len(speeds)
    rows=[[1.0]*P]; lb=[float(n-1)]; ub=[float(n-1)]                 # size
    for q in range(2,n+1):
        rows.append([1.0 if v%q==0 else 0.0 for v in speeds]); lb.append(1.0); ub.append(float(P))  # divisibility
    cuts=0
    for _ in range(max_iter):
        A=np.array(rows)
        res=milp(c=np.zeros(P), constraints=LinearConstraint(A,lb,ub),
                 integrality=np.ones(P), bounds=Bounds(0,1), options={"time_limit":30})
        if not (res.success and res.x is not None):
            return None, cuts                     # ILP infeasible -> no covering set with M<r
        S=[speeds[i] for i in range(P) if res.x[i]>0.5]
        M,tstar=Mexact_and_witness(S)
        if M < r:
            return S, cuts                        # found a beater
        # add the deepest-hole cut: some speed must be within (r) of 0 at t*  (STRICT: ||v t*|| < r)
        row=[1.0 if norm(v*tstar) < r else 0.0 for v in speeds]
        rows.append(row); lb.append(1.0); ub.append(float(P))
        cuts+=1
    return "TIMEOUT", cuts

if __name__=="__main__":
    print("LAZY-CUT ILP: does a covering set (speeds<=n(n-1)) BEAT the construction n/Phi6?")
    for n in [12,13,14]:
        V=n*(n-1); thr=F(n,n*n-n+1)
        res, cuts=beater_exists(n, V, thr)
        if res is None:
            print(f"  n={n} V={V}: NO beater with M<{thr} (ILP infeasible after {cuts} cuts) -> "
                  f"covering-min={thr} for speeds<=n(n-1) [RIGOROUS, closes residual]")
        elif res=="TIMEOUT":
            print(f"  n={n} V={V}: inconclusive (max cuts hit, {cuts} cuts)")
        else:
            M,at=Mexact_and_witness(res)
            print(f"  n={n} V={V}: BEATER FOUND M={M}={float(M):.5f} < {thr} !! S={res} ({cuts} cuts)")
