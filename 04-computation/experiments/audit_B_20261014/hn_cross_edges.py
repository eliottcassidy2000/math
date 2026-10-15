#!/usr/bin/env python3
"""Audit B: structure of the hinge (cross) edges of W = M3 u psi M3 and of M3 u psi P1, from the independent rebuild
(hn_spindle_free_rebuild_zeta24.py).  A vertex lies in M3 (or P1) iff its sqrt5-part is zero; psi B = rho B has a
nonzero sqrt5-part except at 0.  For every cross edge we report |z|^2, |w|^2 and the angle between z and w, and
compare with the claim "the new edges join radius-2 points at the chord angle theta_4 = arccos(7/8) = 28.955 deg".
Also: the near-miss non-edges found by the 1e-6 float window (how close to 1 they are, in high precision).

Usage: python3 hn_cross_edges.py DIR
"""
import sys
import math
from collections import Counter

import numpy as np
import mpmath as mp

D_ = sys.argv[1]
mp.mp.dps = 50
Z = [mp.exp(1j * mp.pi * j / 12) for j in range(8)]
S5 = mp.sqrt(5)


def hp(row, D):
    z = sum(int(row[j]) * Z[j] for j in range(8)) + S5 * sum(int(row[8 + j]) * Z[j] for j in range(8))
    return z / D


for name in ('W_M3_u_psiM3', 'M3_u_psiP1'):
    f = np.load(f"{D_}/{name}.npz")
    X, D, I, J, xy = f['X'], int(f['D'][0]), f['I'], f['J'], f['xy']
    inA = np.all(X[:, 8:] == 0, axis=1)
    cross = inA[I] != inA[J]
    origin = np.all(X == 0, axis=1)
    # edges at the origin are counted as internal to both pieces; report them separately
    o = origin[I] | origin[J]
    cr = cross & ~o
    print(f"{name}: |V| = {len(X)}, |E| = {len(I)}, cross edges (excluding the origin) = {cr.sum()}, "
          f"edges at origin = {o.sum()}")
    stats = Counter()
    ex = []
    for i, j in zip(I[cr], J[cr]):
        a, b = (i, j) if inA[i] else (j, i)
        za, zb = complex(*xy[a]), complex(*xy[b])
        ra, rb = abs(za) ** 2, abs(zb) ** 2
        ang = math.degrees(math.atan2((zb / za).imag, (zb / za).real))
        stats[(round(ra, 6), round(rb, 6), round(abs(ang), 3))] += 1
        if len(ex) < 3:
            ex.append((za, zb))
    for k, v in sorted(stats.items()):
        print(f"   |z|^2 = {k[0]}, |psi w|^2 = {k[1]}, angle = {k[2]} deg: {v} edges")
    print(f"   arccos(7/8) = {math.degrees(math.acos(7 / 8)):.3f} deg")
    # near misses
    x, y = xy[:, 0], xy[:, 1]
    Eset = set(zip(I.tolist(), J.tolist()))
    near = []
    n = len(X)
    for s in range(0, n, 128):
        e = min(n, s + 128)
        d2 = (x[s:e, None] - x[None, :]) ** 2 + (y[s:e, None] - y[None, :]) ** 2
        ii, jj = np.nonzero(np.abs(d2 - 1) < 1e-6)
        ii = ii + s
        for a, b in zip(ii.tolist(), jj.tolist()):
            if b > a and (a, b) not in Eset:
                near.append((a, b))
    gaps = []
    for a, b in near:
        d2 = abs(hp(X[a], D) - hp(X[b], D)) ** 2
        gaps.append(abs(d2 - 1))
    if gaps:
        print(f"   near-miss non-edges: {len(near)}; |d^2 - 1| in [{mp.nstr(min(gaps), 5)}, {mp.nstr(max(gaps), 5)}] "
              f"(50-digit evaluation)")
        print(f"   distinct values of |d^2-1|: {sorted(set(mp.nstr(g, 8) for g in gaps))[:6]}")
