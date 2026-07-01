#!/usr/bin/env python3
"""
covering_min_delsarte_lp_klein.py  --  klein-2026-07-01-S64

Reasoning about the DELSARTE dual, done RIGHT (correcting S63's conflation of the naive Fejer AVERAGE
with the Delsarte DUAL certificate). The Delsarte covering dual is a POINTWISE-valid lower bound, NOT
an average.

STRUCTURED DELSARTE LP (uses the covering structure -- a covering set has one multiple per q in {2..n}):
  find a weight w(t) >= 0, integral 1, minimizing  SUM_q beta_q  where  beta_q = max_{v: q|v} integral w[||vt||<L].
  If a covering set has M < L, its danger arcs cover [0,1), so  1 = int w <= SUM_{v in S} int w[A_v]
  <= SUM_q beta_q.  Hence  SUM_q beta_q^min < 1  =>  covering-min >= L  (RIGOROUS lower bound).
This is an LP over w (the pointwise dual). We scan L and report the certified lower bound, comparing to:
  trivial 1/(2(n-1)),  the LRC floor 1/n,  the covering-min n/Phi6(n).
If it only reaches ~1/(2(n-1)), the LP LEVEL is weak (=> the SDP/Lasserre tightening is genuinely needed).
"""
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog

def dist_grid(v, M):
    # ||v * i/M|| for i=0..M-1
    x = (v*np.arange(M)) % M
    return np.minimum(x, M-x)/M

def delsarte_bound(n, V, M, L):
    """min sum_q beta_q for weight w on M-grid; returns the LP optimum (bound holds iff < 1)."""
    speeds_by_q = {q: [v for v in range(q, V+1, q)] for q in range(2, n+1)}
    qs = list(range(2, n+1)); Q = len(qs)
    # variables: w_0..w_{M-1}, beta_2..beta_n   (total M+Q)
    nvar = M + Q
    c = np.zeros(nvar);
    for j in range(Q): c[M+j] = 1.0     # minimize sum beta_q
    A_ub = []; b_ub = []
    for qi, q in enumerate(qs):
        for v in speeds_by_q[q]:
            ind = (dist_grid(v, M) < L).astype(float)     # danger indicator on grid
            row = np.zeros(nvar); row[:M] = ind; row[M+qi] = -1.0   # sum w*ind - beta_q <= 0
            A_ub.append(row); b_ub.append(0.0)
    A_eq = np.zeros((1, nvar)); A_eq[0,:M] = 1.0; b_eq = [1.0]      # int w = 1
    bounds = [(0,None)]*M + [(0,None)]*Q
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    return res.fun if res.success else None

if __name__=="__main__":
    n=14; V=56; M=900
    triv = 1/(2*(n-1)); floor=1/n; cov=n/(n*n-n+1)
    print(f"STRUCTURED DELSARTE LP (n={n}, V={V}, grid={M}): sum_q beta_q vs 1 (certifies covmin>=L iff <1)")
    print(f"  reference: trivial 1/(2(n-1))={triv:.4f}, LRC floor 1/n={floor:.4f}, covering-min n/Phi6={cov:.4f}")
    print(f"  {'L':>8} {'sum beta_q':>11} {'certifies covmin>=L?':>22}")
    best=0.0
    for L in [0.030,0.0385,0.045,0.050,0.055,0.060,0.065,0.0714,0.0765]:
        s = delsarte_bound(n, V, M, L)
        if s is None: print(f"  {L:>8.4f}  (LP failed)"); continue
        cert = s < 1.0
        if cert and L>best: best=L
        print(f"  {L:>8.4f} {s:>11.4f} {str(cert):>22}")
    print(f"  => Delsarte LP certifies covering-min >= {best:.4f}  (trivial {triv:.4f}, floor {floor:.4f})")
    print("  If ~trivial: the LP LEVEL is weak (grid-approx); the SDP/Lasserre tightening is needed to")
    print("  reach the floor/covering-min -- confirming the certificate must add PSD-moment structure.")
