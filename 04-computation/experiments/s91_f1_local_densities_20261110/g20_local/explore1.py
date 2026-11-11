#!/usr/bin/env python3
"""exploration: bounded functionals (flags) of the minimal representation; orbit counts for short words."""
import itertools
import sys
from fractions import Fraction

import minrep
from lattice_tools import rank_q, solve_q

rep = minrep.build()
R = rep["R"]
A = rep["A"]
L0 = rep["L0"]
FIN = rep["FIN"]


def nullspace_q(mat, ncols):
    """basis of {x : mat x = 0} over Q"""
    Af = [[Fraction(x) for x in row] for row in mat]
    m = len(Af)
    piv = []
    r = 0
    for c in range(ncols):
        pr = None
        for i in range(r, m):
            if Af[i][c]:
                pr = i
                break
        if pr is None:
            continue
        Af[r], Af[pr] = Af[pr], Af[r]
        inv = 1 / Af[r][c]
        Af[r] = [x * inv for x in Af[r]]
        for i in range(m):
            if i != r and Af[i][c]:
                f = Af[i][c]
                Af[i] = [x - f * y for x, y in zip(Af[i], Af[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fc in free:
        x = [Fraction(0)] * ncols
        x[fc] = Fraction(1)
        for i, c in enumerate(piv):
            x[c] = -Af[i][fc]
        basis.append(x)
    return basis


def matmul(X, Y):
    return [[sum(X[i][k] * Y[k][j] for k in range(len(Y))) for j in range(len(Y[0]))] for i in range(len(X))]


def matpow(X, n):
    out = [[int(i == j) for j in range(len(X))] for i in range(len(X))]
    for _ in range(n):
        out = matmul(out, X)
    return out


I = [[int(i == j) for j in range(R)] for i in range(R)]
# functionals lam (columns): l(w c^k).lam = l(w) A_c^k lam bounded for all w  <=>  A_c^N (A_c - I) lam = 0
cons = []
for c in range(4):
    AN = matpow(A[c], R)
    AmI = [[A[c][i][j] - I[i][j] for j in range(R)] for i in range(R)]
    cons += matmul(AN, AmI)
F = nullspace_q(cons, R)
print("dim of functionals bounded along every single-letter power:", len(F))
# largest letter-invariant subspace
while True:
    # lam in span(F) with A_c lam in span(F) for all c
    k = len(F)
    if k == 0:
        break
    # coordinates: lam = sum t_i F_i ; need A_c lam in span F: use a complement test via rank
    # build the projection onto a complement: solve by nullspace of [A_c F | -F] ...
    rows = []
    # constraint: exists s with A_c F t = F s for each c -> unknown (t, s_0..s_3)
    nun = k + 4 * k
    for c in range(4):
        AF = [[sum(A[c][i][j] * F[m][j] for j in range(R)) for m in range(k)] for i in range(R)]  # R x k
        for i in range(R):
            row = [Fraction(0)] * nun
            for m in range(k):
                row[m] = Fraction(AF[i][m])
                row[k + c * k + m] = -F[m][i]
            rows.append(row)
    ns = nullspace_q(rows, nun)
    newF = []
    seen = []
    for v in ns:
        t = v[:k]
        lam = [sum(t[m] * F[m][j] for m in range(k)) for j in range(R)]
        newF.append(lam)
    # reduce to a basis
    basis = []
    for lam in newF:
        if rank_q(basis + [lam]) > len(basis):
            basis.append(lam)
    if len(basis) == k:
        break
    F = basis
print("dim of the largest letter-invariant subspace of such functionals:", len(F))
# values of the flags on words of length <= 7; check finiteness by closure
words = [w for n in range(0, 8) for w in itertools.product(range(4), repeat=n)]
vals = {}
cur = {(): L0}
flagsets = []
allflags = set()
lv_by_len = [{(): L0}]
for n in range(1, 8):
    nxt = {}
    for w, v in lv_by_len[-1].items():
        for c in range(4):
            nxt[w + (c,)] = [sum(v[i] * A[c][i][j] for i in range(R)) for j in range(R)]
    lv_by_len.append(nxt)
for n in range(0, 8):
    for w, v in lv_by_len[n].items():
        fl = tuple(sum(v[j] * lam[j] for j in range(R)) for lam in F)
        allflags.add(fl)
    print(f"  distinct flag vectors over words of length <= {n}: {len(allflags)}")
# distinct left vectors over Z by length
for n in range(0, 8):
    s = set(tuple(v) for v in lv_by_len[n].values())
    print(f"  length {n}: {len(lv_by_len[n])} words, {len(s)} distinct left vectors over Z")
