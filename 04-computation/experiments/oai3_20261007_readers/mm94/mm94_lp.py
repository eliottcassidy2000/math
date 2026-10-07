#!/usr/bin/env python3
"""
mm94 lane: reproduce the optimisation behind Lemma 4.1 of openai/math #107 as a linear program.

Variables P(a,b), 1 <= a <= b <= N (symmetry built in).  Constraints (only those living on the grid):
  P(1,b) = b;  2P(a,b) >= P(a,b+1) + P(a,b-1) (all a, 2 <= b <= N-1, using P(a,c) = P(c,a));
  P(a,3h+a-1) >= 3P(a,h)  (3h+a-1 <= N);   optional extra gadgets (k = 5, 7 alternating blocks).
LP_n := min P(n,n).   Prediction (mm94): LP_n = D_n^min = (3n-1)/2 * prod_{j<n}(1+1/(3j)) once N >= 4n-1.
Also: the largest t for which the rank bound P(a,b) <= (a+b-1)^(1/t) is LP-feasible on the grid.
"""
import numpy as np, scipy.sparse as sp, sys, math
from scipy.optimize import linprog
from fractions import Fraction as Fr

def Dmin(n):
    s = Fr(1)
    for j in range(1, n):
        s *= Fr(3 * j + 1, 3 * j)
    return Fr(3 * n - 1, 2) * s

def build(N, extra_k=()):
    idx = {}
    for a in range(1, N + 1):
        for b in range(a, N + 1):
            idx[(a, b)] = len(idx)
    I = lambda a, b: idx[(min(a, b), max(a, b))]
    rows, cols, vals, rhs = [], [], [], []   # A_ub x <= b_ub
    r = 0
    def add(coefs, bound):
        nonlocal r
        for (a, b), c in coefs:
            rows.append(r); cols.append(I(a, b)); vals.append(c)
        rhs.append(bound); r += 1
    for a in range(1, N + 1):                       # concavity in b (and in a by symmetry)
        for b in range(2, N):
            add([((a, b + 1), 1.0), ((a, b - 1), 1.0), ((a, b), -2.0)], 0.0)
    for a in range(1, N + 1):                       # tripling and optional k-block gadgets
        for k in (3,) + tuple(extra_k):
            for h in range(1, N + 1):
                B = k * h + (k - 1) * (a - 1) // 2
                if B > N: break
                add([((a, h), float(k)), ((a, B), -1.0)], 0.0)
    A_ub = sp.csr_matrix((vals, (rows, cols)), shape=(r, len(idx)))
    b_ub = np.array(rhs)
    eq_rows, eq_cols, eq_vals, eq_rhs = [], [], [], []
    for b in range(1, N + 1):
        eq_rows.append(len(eq_rhs)); eq_cols.append(I(1, b)); eq_vals.append(1.0); eq_rhs.append(float(b))
    A_eq = sp.csr_matrix((eq_vals, (eq_rows, eq_cols)), shape=(len(eq_rhs), len(idx)))
    return idx, I, A_ub, b_ub, A_eq, np.array(eq_rhs)

def lp_min_diag(N, n, extra_k=()):
    idx, I, A_ub, b_ub, A_eq, b_eq = build(N, extra_k)
    c = np.zeros(len(idx)); c[I(n, n)] = 1.0
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    return res.fun if res.status == 0 else None

print("LP_n = min P(n,n) over the paper's hypotheses on an N x N grid, vs closed form D_n^min")
for n in (2, 3, 4, 5, 8, 12, 16, 20, 25):
    N = 4 * n - 1
    v = lp_min_diag(N, n)
    v2 = lp_min_diag(N + 20, n, extra_k=(5, 7))
    d = Dmin(n)
    print(f"  n={n:3d} N={N:4d}: LP = {v:.10f}   with k=5,7 gadgets & N+20: {v2:.10f}   D_n^min = {str(d):>22} = {float(d):.10f}   |diff| = {abs(v-float(d)):.1e}")
# truncation: what if the grid is too small for the paper's chain?
print("Truncation check (grid too small: N < 4n-1 drops the tripling constraint P(n,4n-1) >= 3P(n,n)):")
for n, N in ((5, 12), (5, 16), (5, 19), (8, 24), (8, 31)):
    v = lp_min_diag(N, n)
    print(f"  n={n} N={N}: LP = {v:.6f}  vs D_n^min = {float(Dmin(n)):.6f}")

# feasibility threshold for the rank bound P(a,b) <= (a+b-1)^(1/t)
def feasible(N, t):
    idx, I, A_ub, b_ub, A_eq, b_eq = build(N)
    ub = np.zeros(len(idx))
    for (a, b), k in idx.items():
        ub[k] = (a + b - 1) ** (1.0 / t)
    res = linprog(np.zeros(len(idx)), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=list(zip(np.zeros(len(idx)), ub)), method="highs")
    return res.status == 0
print("Rank-bound feasibility threshold t_N (bisection) vs omega_a/3 at a = (N+1)/4:")
for N in (39, 79, 119):
    lo, hi = 0.70, 1.0
    for _ in range(18):
        mid = (lo + hi) / 2
        if feasible(N, mid): lo = mid
        else: hi = mid
    a = (N + 1) // 4
    om = 3 * math.log(2 * a - 1) / math.log(float(Dmin(a)))
    print(f"  N={N:4d}: t_N = {lo:.6f}  (3 t_N = {3*lo:.6f});  omega_a/3 at a={a}: {om/3:.6f}")
print("At t = 3/4 the grid LP is feasible for every N (witness: P_min).  Check N=119:", feasible(119, 0.75))
