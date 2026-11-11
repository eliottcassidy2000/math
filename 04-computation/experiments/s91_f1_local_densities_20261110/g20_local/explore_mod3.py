#!/usr/bin/env python3
"""exploration: the shape of the orbit modulo 3 inside each affine class (3^10 points in an 11-dimensional class)."""
import itertools
import sys

import numpy as np

import minrep
from orbit_modp import rref_mod, in_span, nb_closure_bitmap, R

p = 3
rep = minrep.build()
A = np.array([[[rep["A"][c][i][j] % p for j in range(R)] for i in range(R)] for c in range(4)], dtype=np.int64)
L0 = np.array([x % p for x in rep["L0"]], dtype=np.int64)
FIN = np.array([x % p for x in rep["FIN"]], dtype=np.int64)
codes, size = nb_closure_bitmap(A, L0, p, 40_000_000)
x = codes.copy()
digs = np.zeros((codes.size, R), dtype=np.int64)
for j in range(R):
    digs[:, j] = x % p
    x //= p
base = {0: L0.copy()}
for j in (1, 2, 3):
    base[j] = (L0 @ A[j]) % p
# class by flags: recompute spans from the orbit itself
for j in range(4):
    # points of class j: those reachable as f_j(x) for x in the orbit, plus l(empty) for j = 0
    img = (digs @ A[j]) % p
    pts = np.unique(img, axis=0)
    if j == 0:
        pts = np.unique(np.vstack([pts, L0.reshape(1, -1)]), axis=0)
    d = (pts - base[j]) % p
    rs = np.random.RandomState(7)
    E, piv = rref_mod(d[rs.choice(d.shape[0], 4000, replace=False)], p)
    co = d[:, piv]
    assert np.array_equal((co @ E) % p, d)
    k = len(piv)
    print(f"class {j}: {pts.shape[0]} orbit points, affine span of dimension {k}")
    # coordinates y in F_3^k; stabiliser: v with Y + v = Y
    Y = co
    pw = p ** np.arange(k, dtype=np.int64)
    Yc = set((Y @ pw).tolist())
    stab = []
    for i in range(k):
        e = np.zeros(k, dtype=np.int64)
        e[i] = 1
        if all(int(((y + e) % p) @ pw) in Yc for y in Y[:2000]):
            stab.append(i)
    print(f"   coordinate directions e_i with Y + e_i = Y (tested on 2000 points): {stab}")
    # polynomial relations of degree <= 2 on Y: monomials 1, y_i, y_i y_j (i <= j)
    mons = [np.ones(Y.shape[0], dtype=np.int64)] + [Y[:, i] for i in range(k)]
    names = ["1"] + [f"y{i}" for i in range(k)]
    for i in range(k):
        for jj in range(i, k):
            mons.append((Y[:, i] * Y[:, jj]) % p)
            names.append(f"y{i}*y{jj}")
    Mx = np.stack(mons, axis=1) % p
    # nullspace mod 3 of Mx (columns = monomials): relations c with Mx c = 0
    # compute rref of Mx^T rows? use rank via rref on a subsample of rows then verify on all
    Mu = np.unique(Mx, axis=0)
    Rm, pv = rref_mod(Mu[rs.choice(Mu.shape[0], min(Mu.shape[0], 12000), replace=False)], p)
    rank = len(pv)
    free = [c for c in range(Mx.shape[1]) if c not in pv]
    rels = []
    for fc in free:
        rel = np.zeros(Mx.shape[1], dtype=np.int64)
        rel[fc] = 1
        for i, c in enumerate(pv):
            rel[c] = (-Rm[i, fc]) % p
        rels.append(rel)
    good = [rel for rel in rels if not ((Mx @ rel) % p).any()]
    print(f"   monomials of degree <= 2: {Mx.shape[1]}, rank on a sample: {rank}; candidate relations {len(rels)}, "
          f"valid on all points {len(good)}")
    for rel in good[:6]:
        terms = [f"{int(rel[i])}*{names[i]}" for i in range(len(names)) if rel[i]]
        print("   relation: " + " + ".join(terms) + " = 0")
    # count the F_3 points of the common zero set of the valid relations (brute force over 3^k if small)
    if good and k <= 11:
        allpts = np.array(list(itertools.product(range(p), repeat=k)), dtype=np.int64)
        mons2 = [np.ones(allpts.shape[0], dtype=np.int64)] + [allpts[:, i] for i in range(k)]
        for i in range(k):
            for jj in range(i, k):
                mons2.append((allpts[:, i] * allpts[:, jj]) % p)
        M2 = np.stack(mons2, axis=1) % p
        zero = np.ones(allpts.shape[0], dtype=bool)
        for rel in good:
            zero &= ((M2 @ rel) % p) == 0
        print(f"   points of F_3^{k} satisfying the valid quadratic relations: {int(zero.sum())} (orbit: {Y.shape[0]})")
