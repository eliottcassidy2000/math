# Independent LP: min P(a,b) at off-diagonal points subject to #107 Lemma 4.1 hypotheses on a finite grid
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction as Fr
def sig(m):
    s = Fr(1)
    for j in range(1, m): s *= Fr(3*j+1, 3*j)
    return s
def Pmin(a, b):
    m, M = min(a, b), max(a, b); return sig(m)*Fr(2*M+m-1, 2)
def lpmin(N, target):
    idx = {}
    for a in range(1, N+1):
        for b in range(a, N+1): idx[(a, b)] = len(idx)
    I = lambda a, b: idx[(min(a, b), max(a, b))]
    rows, cols, vals, rhs = [], [], [], []; r = 0
    def add(co):
        nonlocal r
        for (a, b), c in co: rows.append(r); cols.append(I(a, b)); vals.append(c)
        rhs.append(0.0); r += 1
    for a in range(1, N+1):
        for b in range(2, N):
            add([((a, b+1), 1.0), ((a, b-1), 1.0), ((a, b), -2.0)])
    for a in range(1, N+1):
        for h in range(1, N+1):
            B = 3*h + a - 1
            if B > N: break
            add([((a, h), 3.0), ((a, B), -1.0)])
    A = sp.csr_matrix((vals, (rows, cols)), shape=(r, len(idx)))
    Ae = sp.csr_matrix((np.ones(N), (np.arange(N), [I(1, b) for b in range(1, N+1)])), shape=(N, len(idx)))
    c = np.zeros(len(idx)); c[I(*target)] = 1.0
    res = linprog(c, A_ub=A, b_ub=np.zeros(r), A_eq=Ae, b_eq=np.arange(1, N+1, dtype=float), bounds=(0, None), method="highs")
    return res.fun
for (a, b) in [(2, 5), (3, 7), (3, 10), (4, 9), (5, 6), (6, 15), (2, 30), (7, 8)]:
    N = 3*(b-1) + a - 1 + 2
    N = max(N, 4*a)
    v = lpmin(N, (a, b))
    print(f"(a,b)=({a},{b}) grid N={N}: LP min = {v:.10f}  P_min = {float(Pmin(a,b)):.10f}  diff = {v-float(Pmin(a,b)):.2e}")
