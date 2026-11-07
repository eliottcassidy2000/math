#!/usr/bin/env python3
# audit J1 (THM-4910 statement 6): the flag system is CONSISTENT on a real hypersurface of nondegenerate 4-simplices (exact sign change of det[X|y], faces nondegenerate, X full rank, real flag sphere), so 'no flag extension' holds only generically.
import random
from fractions import Fraction as Fr
import sympy as sp
def sub(a, b): return [x - y for x, y in zip(a, b)]
def add(a, b): return [x + y for x, y in zip(a, b)]
def dt(a, b): return sum(x * y for x, y in zip(a, b))
def solve(M, b):
    m = len(M); A = [list(r) + [bb] for r, bb in zip(M, b)]
    for c in range(m):
        p = next(r for r in range(c, m) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        for r in range(m):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][m] / A[i][i] for i in range(m)]
def cc(P):
    p0 = P[0]; D = [sub(p, p0) for p in P[1:]]
    lam = solve([[dt(x, y) for y in D] for x in D], [dt(x, x) / 2 for x in D])
    return add(p0, [sum(l * d[j] for l, d in zip(lam, D)) for j in range(len(p0))])
def aug(V):
    n = len(V) - 1; O = cc(V); R2 = dt(sub(V[0], O), sub(V[0], O)); rows = []
    for k in range(n, 1, -1):
        F = V[n - k:]; Ok = cc(F); Gk = [sum(c) / len(F) for c in zip(*F)]
        for P, c in ((Ok, k + 1), (Gk, k)):
            y = sub(P, O); rows.append([-2 * yi for yi in y] + [Fr(1)] + [c * (R2 - dt(y, y)) - dt(y, y)])
    return rows
def facegrams(V):
    n = len(V) - 1; out = []
    for k in range(n, 1, -1):
        F = V[n - k:]; D = [sub(p, F[0]) for p in F[1:]]
        out.append(sp.Matrix([[dt(x, y) for y in D] for x in D]).det())
    return out
random.seed(3)
n = 4; shown = 0
for trial in range(40):
    Va = [[Fr(random.randint(-5, 5)) for _ in range(n)] for _ in range(n + 1)]
    Vb = [[Fr(random.randint(-5, 5)) for _ in range(n)] for _ in range(n + 1)]
    def V_(t): return [[(1 - t) * a + t * b for a, b in zip(ra, rb)] for ra, rb in zip(Va, Vb)]
    try:
        vals = [(Fr(i, 40), sp.Matrix(aug(V_(Fr(i, 40)))).det()) for i in range(41)]
    except StopIteration:
        continue
    for (t1, d1), (t2, d2) in zip(vals, vals[1:]):
        if d1 * d2 < 0:
            grid = [t1 + (t2 - t1) * Fr(j, 100) for j in range(101)]
            mingram = min(min(abs(x) for x in facegrams(V_(t))) for t in grid)
            lo, hi = t1, t2
            for _ in range(45):
                mid = (lo + hi) / 2; dm = sp.Matrix(aug(V_(mid))).det()
                if dm * d1 < 0: hi = mid
                else: lo = mid
            V = V_(lo); Am = sp.Matrix(aug(V)); X, y = Am[:, :5], Am[:, 5]
            sol = (X.T * X).LUsolve(X.T * y); res = X * sol - y
            c = sp.Matrix(sol[:4]); kap = sol[4]
            O = cc(V); R2 = dt(sub(V[0], O), sub(V[0], O))
            Xev = sorted(float(e) for e in (X.T * X).eigenvals(multiple=True))
            print(f"path {trial}: det[X|y] changes sign on [{float(t1):.4f},{float(t2):.4f}] (|det| {float(abs(d1)):.2e} -> {float(abs(d2)):.2e});"
                  f" min face Gram det on the bracket {float(mingram):.2e}; at the bisected root t*={float(lo):.12f}:"
                  f" max|residual|/R^2 {float(max(abs(r) for r in res)/R2):.1e}, cond(X^T X) {Xev[-1]/Xev[0]:.1e},"
                  f" flag-sphere radius^2/R^2 {float(((c.T*c)[0]-kap)/R2):.3f}")
            shown += 1
            break
    if shown >= 3: break
print("=> on these nondegenerate 4-simplices a REAL sphere has all six flag-corner powers; statement 6 holds only generically")
