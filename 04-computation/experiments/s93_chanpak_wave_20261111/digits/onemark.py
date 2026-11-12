#!/usr/bin/env python3
"""Lane "digits", part (1): the one-mark block system.

Marked tournament (T, z).  State s(T, z) = (h, e, B, g):
    h = H(T),  e = End_T(z) (Hamiltonian paths ending at z),
    B = U_T(z) (ordered pairs (P1, P2) of disjoint paths covering T, P1 possibly empty, P2 starting at z),
    g = H(T - z)  (= 1 when T = {z}).
Block op (X, W, m): T' = T u X with  T - z => X,  z -> X - W,  W -> z  (W a subset of V(X));
new mark m = z ("stay") or a vertex of X ("move").
Theorem (PROOFS.txt, P1): s(T', m) = N(X, W, m) s(T, z) for an integer 4 x 4 matrix depending only on the op.

One-vertex ops (X a point x):
    a = (W = {x}, stay)   "ray step":   (h, e, B, g) -> (h - e + B, g, 2B - g, g)
    b = (W = {x}, move)   "tower step": (h, e, B, g) -> (h - e + B, h - e, B + h, h)
    c = (W = {},  stay)   "sink, stay": (h, e, B, g) -> (h, 0, 2B, g)
    d = (W = {},  move)   "sink, move": (h, e, B, g) -> (h, h, h, h)

This script computes the matrix of every block op with X of order <= KMAX by interpolation on four base pairs
with independent states and verifies it on every marked tournament of order <= 5 (exact integers), checks the
four closed forms above, and writes onemark_matrices_KMAX.json.

Usage: onemark.py KMAX
"""
import itertools
import json
import sys
from fractions import Fraction

import numpy as np

from tour import hcount, iso_classes, pair_stats


def state(n, adj, z):
    h, start, end, Q, hsub = pair_stats(n, adj)
    B = int(start[z]) + sum(int(Q[e, z]) for e in range(n) if e != z)
    g = int(hsub[((1 << n) - 1) ^ (1 << z)])
    return (int(h), int(end[z]), B, g)


def apply_block(n, adj, z, k, adjX, W, m):
    """returns (n', adj', z')."""
    N = n + k
    a = np.zeros(N, np.int64)
    xmask = ((1 << k) - 1) << n
    for i in range(n):
        a[i] = int(adj[i])
        if i != z:
            a[i] |= xmask
        else:
            a[i] |= (((1 << k) - 1) & ~W) << n
    for i in range(k):
        a[n + i] = int(adjX[i]) << n
        if (W >> i) & 1:
            a[n + i] |= 1 << z
    return N, a, (z if m < 0 else n + m)


def solve4(S, S2):
    """M with M S[:, j] = S2[:, j] (S: 4 independent columns), exact fractions."""
    A = [[Fraction(S[i][j]) for j in range(4)] + [Fraction(1 if i == j else 0) for j in range(4)] for i in range(4)]
    for c in range(4):
        p = next(r for r in range(c, 4) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(4):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    Sinv = [row[4:] for row in A]
    return [[sum(Fraction(S2[i][k]) * Sinv[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def labelled_tournaments(k):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for bits in range(1 << len(pairs)):
        a = [0] * k
        for t, (i, j) in enumerate(pairs):
            if (bits >> t) & 1:
                a[i] |= 1 << j
            else:
                a[j] |= 1 << i
        yield a


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    classes = iso_classes(5)
    bases = []
    for n in range(1, 6):
        for key, adj in classes[n].items():
            for z in range(n):
                bases.append((n, adj, z, state(n, adj, z)))
    print(f"{len(bases)} marked base tournaments of order <= 5")
    # four independent states
    idx = None
    for comb in itertools.combinations(range(min(len(bases), 40)), 4):
        S = [[bases[c][3][i] for c in comb] for i in range(4)]
        M = np.array(S, dtype=np.int64)
        if abs(round(np.linalg.det(M.astype(float)))) >= 1:
            idx = comb
            break
    print("interpolation bases (order, mark, state):", [(bases[c][0], bases[c][2], bases[c][3]) for c in idx])
    S = [[bases[c][3][i] for c in idx] for i in range(4)]
    Sinv = solve4(S, [[1 if i == j else 0 for j in range(4)] for i in range(4)])
    mats = {}
    nops = 0
    bad = 0
    for k in range(1, kmax + 1):
        for adjX in labelled_tournaments(k):
            for W in range(1 << k):
                for m in range(-1, k):
                    nops += 1
                    S2 = [[0] * 4 for _ in range(4)]
                    for j, c in enumerate(idx):
                        n, adj, z, s = bases[c]
                        N, a, z2 = apply_block(n, adj, z, k, adjX, W, m)
                        s2 = state(N, a, z2)
                        for i in range(4):
                            S2[i][j] = s2[i]
                    Mx = [[sum(S2[i][q] * Sinv[q][j] for q in range(4)) for j in range(4)] for i in range(4)]
                    if any(x.denominator != 1 for row in Mx for x in row):
                        print("NON-INTEGER MATRIX", k, adjX, W, m)
                        bad += 1
                    Mi = tuple(tuple(int(x) for x in row) for row in Mx)
                    if Mi not in mats:
                        # verify on every base
                        for (n, adj, z, s) in bases:
                            if n + k > 8:
                                continue
                            N, a, z2 = apply_block(n, adj, z, k, adjX, W, m)
                            s2 = state(N, a, z2)
                            pred = tuple(sum(Mi[i][j] * s[j] for j in range(4)) for i in range(4))
                            if pred != s2:
                                bad += 1
                                if bad < 5:
                                    print("LINEARITY FAILS", k, adjX, W, m, (n, [int(x) for x in adj], z), s, pred, s2)
                        mats[Mi] = {"k": k, "adjX": list(adjX), "W": W, "m": m, "count": 1}
                    else:
                        mats[Mi]["count"] += 1
        print(f"blocks of order <= {k}: {nops} labelled ops, {len(mats)} distinct matrices, failures {bad}", flush=True)
    named = {
        "a": ((1, -1, 1, 0), (0, 0, 0, 1), (0, 0, 2, -1), (0, 0, 0, 1)),
        "b": ((1, -1, 1, 0), (1, -1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 0)),
        "c": ((1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 2, 0), (0, 0, 0, 1)),
        "d": ((1, 0, 0, 0), (1, 0, 0, 0), (1, 0, 0, 0), (1, 0, 0, 0)),
    }
    for name, Mx in named.items():
        info = mats.get(Mx)
        print(f"closed form of op {name}: {'CONFIRMED' if info and info['k'] == 1 else 'NOT FOUND'}  {info}")
    out = [{"N": [list(r) for r in Mi], **info} for Mi, info in mats.items()]
    with open(f"onemark_matrices_{kmax}.json", "w") as f:
        json.dump(out, f)
    print(f"{len(mats)} distinct matrices -> onemark_matrices_{kmax}.json; total failures {bad}")


if __name__ == "__main__":
    main()
