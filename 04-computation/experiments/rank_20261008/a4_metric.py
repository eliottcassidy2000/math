#!/usr/bin/env python3
"""Full rank on Z_5 (five multipliers with affinely independent log-vectors, rank 4): in coordinates x (unit at position 0,
non-unit positions = basis vectors) the form Q = 5I - J is the A_4 root-lattice metric (|x|^2 + (sum x)^2 = |y|^2 with y the
S_5-symmetric coordinates).  Check exactly, for every affine coupling j -> a j + b on Z/5:
  margin sign of C_pi under Q (balanced iff tr(C Q^-1) > 2 lambda_max(C Q^-1)): claim < 0 only never, = 0 exactly for the
  involutions a = -1, > 0 otherwise (except the identity, C = 0);
and that two different involutions have no common top eigenvector (so their sum is strictly balanced)."""
from fractions import Fraction as Fr
import itertools
import numpy as np
d = 5
v = [[0, 0, 0, 0]] + [[int(i == j) for i in range(4)] for j in range(4)]
def D(a, b):
    M = [[0] * 4 for _ in range(4)]
    for j in range(d):
        xi = [p - q for p, q in zip(v[(a * j + b) % d], v[j])]
        for x in range(4):
            for y in range(4): M[x][y] += xi[x] * xi[y]
    return M
Q = [[5 * (i == j) - 1 for j in range(4)] for i in range(4)]
Qinv = [[Fr(int(i == j) + 1, 5) for j in range(4)] for i in range(4)]     # (5I - J)^-1 = (I + J)/5
assert all(sum(Q[i][k] * Qinv[k][j] for k in range(4)) == (i == j) for i in range(4) for j in range(4))
def tr_lmax(M):
    W = np.array([[float(sum(Fr(M[x][k]) * Qinv[k][y] for k in range(4))) for y in range(4)] for x in range(4)])
    ev = np.sort(np.linalg.eigvals(W).real)
    return sum(Fr(M[x][k]) * Qinv[k][x] for x in range(4) for k in range(4)), ev[-1]
def det(M):
    M = [list(map(Fr, r)) for r in M]; n = len(M); dv = Fr(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None: return Fr(0)
        if p != c: M[c], M[p] = M[p], M[c]; dv = -dv
        dv *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]; M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return dv
def classify(M):
    t, lm = tr_lmax(M)
    S = [[t / 2 * Q[x][y] - M[x][y] for y in range(4)] for x in range(4)]
    minors = [det([r[:k] for r in S[:k]]) for k in range(1, 5)]
    if all(m > 0 for m in minors): return 'balanced (strict)', float(t / 2), lm
    if det(S) == 0 and all(np.linalg.eigvalsh(np.array([[float(x) for x in r] for r in S])) > -1e-12): return 'EQUALITY (PSD, singular)', float(t / 2), lm
    return 'NOT balanced', float(t / 2), lm
for a in range(1, d):
    for b in range(d):
        if (a, b) == (1, 0): continue
        c, ht, lm = classify(D(a, b))
        print(f"coupling j -> {a} j + {b}: {c}; tr/2 = {ht:.3f}, lambda_max = {lm:.3f}  (both in Q^-1 units, times 1/d omitted)")
# common top eigenvectors of different involutions
tops = {}
for b in range(d):
    M = np.array(D(4, b), dtype=float); W = M @ np.array([[float(x) for x in r] for r in Qinv])
    w, V = np.linalg.eig(W); w = w.real; V = V.real
    tops[b] = V[:, np.isclose(w, w.max())]
for b1, b2 in itertools.combinations(range(d), 2):
    S = np.hstack([tops[b1], tops[b2]])
    print(f"involutions b = {b1}, {b2}: top eigenspaces dims {tops[b1].shape[1]}, {tops[b2].shape[1]}, combined rank {np.linalg.matrix_rank(S, tol=1e-9)} (intersection trivial iff = sum of dims)")
