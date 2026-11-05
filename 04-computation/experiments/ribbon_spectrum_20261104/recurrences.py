#!/usr/bin/env python3
"""recurrences.py -- rigorous linear recurrences for H(a^n) (constant-letter ribbon words):
the Krylov space of the initial state under M_a is finite-dimensional, its minimal polynomial gives a
recurrence satisfied by every coordinate, hence by H.  Exact integer/Fraction arithmetic.  Also
rank(M_g - I) on the reachable space for short blocks g (identity-plus-nilpotent search)."""
import sys, itertools
from fractions import Fraction
from ribbon import transfer_system, letters, count_with_matrices, word_str

def matvec(M, v, d):
    out = [0] * d
    for (j, i), m in M.items():
        if v[i]:
            out[j] += m * v[i]
    return out

def dependency(vectors):
    """return coefficients c (Fractions) with sum c_k vectors[k] = 0 and c[-1] = 1, or None if independent."""
    k = len(vectors); d = len(vectors[0])
    # solve A c = 0 with c_{k-1} = 1: sum_{i<k-1} c_i v_i = -v_{k-1}
    rows = [[Fraction(vectors[i][r]) for i in range(k - 1)] + [Fraction(-vectors[k - 1][r])] for r in range(d)]
    # gaussian elimination
    ncol = k - 1
    piv_rows = []
    r = 0
    pivcols = []
    A = rows
    for c in range(ncol):
        pr = None
        for i in range(r, d):
            if A[i][c] != 0:
                pr = i; break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(d):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        pivcols.append(c); r += 1
    # consistency: rows beyond r must have zero rhs
    for i in range(r, d):
        if A[i][ncol] != 0:
            return None
    if r < ncol:
        return None  # underdetermined: shouldn't happen for a Krylov chain (first dependency)
    c = [Fraction(0)] * ncol
    for i, col in enumerate(pivcols):
        c[col] = A[i][ncol]
    return c + [Fraction(1)]

for w in (2, 3):
    order, states, mats, final = transfer_system(w)
    d = len(order)
    print(f"=== span {w}")
    for a in letters(w):
        if not a: continue
        v = [0] * d; v[0] = 1
        kry = [v]
        while True:
            nv = matvec(mats[a], kry[-1], d)
            kry.append(nv)
            dep = dependency(kry)
            if dep is not None:
                break
        k = len(kry) - 1  # degree of the minimal polynomial of the restriction
        # recurrence: sum_{i=0}^{k} dep[i] M^i v = 0  => H(n+k) = -sum_{i<k} dep[i] H(n+i)
        coeffs = [-dep[i] for i in range(k)]
        assert all(c.denominator == 1 for c in coeffs)
        coeffs = [int(c) for c in coeffs]
        seq = [count_with_matrices(w, [a] * n, (order, states, mats, final)) for n in range(0, 40)]
        ok = all(seq[n + k] == sum(coeffs[i] * seq[n + i] for i in range(k)) for n in range(0, 40 - k))
        # minimal recurrence order for the H sequence itself may be smaller; report the Krylov one
        poly = " + ".join(f"({c})x^{i}" for i, c in enumerate(coeffs))
        print(f"letter {sorted(a)}: Krylov dimension {k}; H(n+{k}) = sum_i c_i H(n+i), c = {coeffs}; verified on 40 terms: {ok}")
        print(f"    first terms: {seq[1:16]}")
print("=== identity-plus-nilpotent search (span 2, blocks of length <= 4 over all letters), rank(M_g - I) on the reachable space")
w = 2
order, states, mats, final = transfer_system(w)
d = len(order)
L = letters(w)
# reachable space basis
def rank_of(vectors):
    rows = [[Fraction(x) for x in v] for v in vectors]
    r = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        pr = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0: pr = i; break
        if pr is None: continue
        rows[r], rows[pr] = rows[pr], rows[r]
        pv = rows[r][c]; rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]; rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        r += 1
    return r
reach = []
v0 = [0] * d; v0[0] = 1
frontier = [v0]; reach = [v0]
while frontier:
    nf = []
    for v in frontier:
        for a in L:
            nv = matvec(mats[a], v, d)
            if rank_of(reach + [nv]) > len(reach):
                reach.append(nv); nf.append(nv)
    frontier = nf
r = len(reach)
print("reachable dimension", r)
best = []
for k in (1, 2, 3, 4):
    for g in itertools.product(L, repeat=k):
        # images of the reachable basis under M_g minus identity
        imgs = []
        for v in reach:
            x = v
            for a in g: x = matvec(mats[a], x, d)
            imgs.append([xi - vi for xi, vi in zip(x, v)])
        # rank of (M_g - I) restricted: rank of the image vectors
        rk = rank_of(imgs)
        best.append((rk, word_str(list(g), w)))
best.sort()
print("smallest ranks of M_g - I on the reachable space:", best[:12])
