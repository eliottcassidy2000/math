#!/usr/bin/env python3
# audit J1 (THM-4910 statement 6): exact rank test of the flag extension (no sphere with Pow(O_k)=(k+1)(R^2-|OO_k|^2), Pow(G_k)=k(R^2-|OG_k|^2), k=n..2) plus the MINIMUM float residual over random simplices.
import random
from fractions import Fraction as Fr
import numpy as np

def sub(a, b): return [x - y for x, y in zip(a, b)]
def add(a, b): return [x + y for x, y in zip(a, b)]
def dt(a, b): return sum(x * y for x, y in zip(a, b))
def rank(M):
    A = [list(r) for r in M]; r = 0; rows = len(A); cols = len(A[0])
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == rows: break
    return r
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

random.seed(6)
for n in (3, 4, 5, 6):
    incons = 0; tested = 0
    while tested < 10:
        V = [[Fr(random.randint(-6, 6), random.randint(1, 3)) for _ in range(n)] for _ in range(n + 1)]
        try:
            O = cc(V)
        except StopIteration:
            continue
        R2 = dt(sub(V[0], O), sub(V[0], O))
        rows, rhs = [], []
        try:
            for k in range(n, 1, -1):
                F = V[n - k:]; Ok = cc(F); Gk = [sum(c) / len(F) for c in zip(*F)]
                for P, c in ((Ok, k + 1), (Gk, k)):
                    y = sub(P, O)
                    rows.append([-2 * yi for yi in y] + [Fr(1)]); rhs.append(c * (R2 - dt(y, y)) - dt(y, y))
        except StopIteration:
            continue
        tested += 1
        r0 = rank(rows); r1 = rank([r + [b] for r, b in zip(rows, rhs)])
        incons += (r1 > r0)
    print(f"n={n}: {tested} exact rational simplices; flag system (2(n-1) = {2*(n-1)} equations, {n+1} unknowns) inconsistent in {incons}/{tested}")

rng = np.random.default_rng(2026)
def ccf(P):
    p0 = P[0]; D = P[1:] - p0
    return p0 + np.linalg.solve(D @ D.T, 0.5 * (D**2).sum(1)) @ D
for n in (4, 5, 6):
    rels = []
    for _ in range(500):
        V = rng.normal(size=(n + 1, n)); O = ccf(V); R2 = ((V[0] - O)**2).sum()
        pts, vals = [], []
        for k in range(n, 1, -1):
            F = V[n - k:]; Ok = ccf(F); Gk = F.mean(0)
            pts += [Ok, Gk]; vals += [(k + 1) * (R2 - ((Ok - O)**2).sum()), k * (R2 - ((Gk - O)**2).sum())]
        P = np.array(pts) - O; y = np.array(vals) - (P**2).sum(1)
        X = np.hstack([-2 * P, np.ones((len(P), 1))])
        sol, *_ = np.linalg.lstsq(X, y, rcond=None)
        rels.append(np.abs(X @ sol - y).max() / np.abs(y).max())
    rels = np.array(rels)
    print(f"n={n}: 500 random simplices: relative affine-fit residual min {rels.min():.2e}, median {np.median(rels):.2e}, max {rels.max():.2e}")
