#!/usr/bin/env python3
"""
covering_min_lp_lower_bound_klein.py  --  klein-2026-07-01-S61

CREATIVE METHOD #1 (beyond the ILP): an LP-RELAXATION INFEASIBILITY CERTIFICATE for a rigorous
covering-min LOWER bound, at LARGE V where the exact ILP times out (last session's residual).

Set-cover: a covering set with M(S) < r exists iff (size n-1) + (divisibility: mult of every q in
{2..n}) + (danger: every breakpoint t=k/d, d<=2V, has a selected speed with ||v t|| < r) is FEASIBLE.
The LP RELAXATION (x_v in [0,1]) is a relaxation: if the LP is INFEASIBLE, the ILP is infeasible, so
NO covering set (speeds <= V) has M < r -- a RIGOROUS lower bound covering-min >= r. LP (linprog/HiGHS)
scales to V = n(n-1) (the construction scale), unlike the ILP.

Test: r = n/Phi6 (the construction value). If the LP with 'M < n/Phi6' danger constraints is INFEASIBLE
at V = n(n-1), then the construction IS the covering-min (speeds <= n(n-1)) -- closing the (4n, n(n-1)]
residual of HYP-3778. If FEASIBLE, the LP relaxation is too weak (integrality gap) -- report honestly.
"""
from fractions import Fraction as F
import math, functools
import numpy as np
from scipy.optimize import linprog
print = functools.partial(print, flush=True)

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1-f)

def lp_lower_bound(n, V):
    thr = F(n, n*n-n+1)              # construction value n/Phi6
    speeds = list(range(1, V+1)); P = len(speeds)
    univ = sorted({F(k,d) for d in range(2, 2*V+1) for k in range(1,d) if math.gcd(k,d)==1})
    # constraints A_ub x <= b_ub  (we use >= as -A x <= -b) and A_eq x = b_eq
    rows_ge = []; b_ge = []
    # size = n-1 (equality)
    A_eq = [[1.0]*P]; b_eq = [float(n-1)]
    # divisibility: sum_{q|v} x_v >= 1
    for q in range(2, n+1):
        rows_ge.append([1.0 if v % q == 0 else 0.0 for v in speeds]); b_ge.append(1.0)
    # danger: every breakpoint covered by a speed STRICTLY closer than thr  (=> M < thr)
    ncov0 = 0
    for t in univ:
        row = [1.0 if norm(v*t) < thr else 0.0 for v in speeds]
        if sum(row) == 0:
            ncov0 += 1                     # a breakpoint NO speed<=V can cover within thr => M>=thr FORCED
        rows_ge.append(row); b_ge.append(1.0)
    # if any danger row is all-zero, the LP is trivially infeasible (rigorous!)
    A_ub = -np.array(rows_ge); b_ub = -np.array(b_ge)
    res = linprog(c=np.zeros(P), A_ub=A_ub, b_ub=b_ub, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0,1)]*P, method="highs")
    feasible = res.status == 0
    return thr, feasible, len(univ), ncov0

if __name__=="__main__":
    print("LP-RELAXATION LOWER BOUND: is 'M < n/Phi6' infeasible (=> construction is the covering-min)?")
    print(f"{'n':>3} {'V':>5} {'n/Phi6':>10} {'#witnesses':>11} {'all-0 rows':>10} {'LP M<thr':>10} {'verdict':>34}")
    for n, V in [(12,132),(13,156),(14,182)]:
        thr, feas, U, z = lp_lower_bound(n, V)
        verdict = "LP feasible -> too weak (integrality gap)" if feas else "LP INFEASIBLE -> covmin>=n/Phi6 RIGOROUS!"
        print(f"{n:>3} {V:>5} {str(thr):>10} {U:>11} {z:>10} {str(not feas):>10} {verdict:>34}")
    print()
    print("If LP infeasible: NO covering set with speeds<=n(n-1) beats the construction (rigorous, closes")
    print("the (4n,n(n-1)] residual of HYP-3778). If feasible: the LP relaxation has an integrality gap;")
    print("need a tighter relaxation (Lovasz theta / SDP / Lasserre) or the exact ILP.")
