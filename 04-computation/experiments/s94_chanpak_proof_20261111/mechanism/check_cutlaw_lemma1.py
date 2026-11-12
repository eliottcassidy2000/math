#!/usr/bin/env python3
"""Direct exhaustive checks (exact) of PROOFS.txt (0.1) and Lemma 1, core order c = 3..6:
 (C1) for every labelled core C and every cut B:  table value H(C + v_B)  ==  cut law
        sum_{y not in B} s(y) + sum_{x in B} e(x) + sum_{x in B, y not in B} Q(x, y);
      and == H(C) + sum_{x in B} (U'(x) - s(x)) - sum_{pairs of B} W_xy.
 (C2) for every labelled core and every pair x != y:  W_xy >= 2 H(C - x - y); number of equalities reported."""
import numpy as np
from numba import njit, prange
from arc_cube import table

@njit(cache=True)
def one(c, m, pi, pj, Hbig, cm_scatter, cutbit):
    adj = np.zeros(c, np.int64); radj = np.zeros(c, np.int64)
    na = c * (c - 1) // 2
    for k in range(na):
        if (m >> k) & 1:
            a = pi[k]; b = pj[k]
        else:
            a = pj[k]; b = pi[k]
        adj[a] |= 1 << b; radj[b] |= 1 << a
    full = (1 << c) - 1
    f = np.zeros((1 << c, c), np.int64); g = np.zeros((1 << c, c), np.int64)
    for v in range(c):
        f[1 << v, v] = 1; g[1 << v, v] = 1
    for S in range(1, 1 << c):
        for last in range(c):
            cf = f[S, last]
            if cf != 0:
                nxt = adj[last] & ~S
                for w in range(c):
                    if (nxt >> w) & 1: f[S | (1 << w), w] += cf
            cg = g[S, last]
            if cg != 0:
                nxt = radj[last] & ~S
                for w in range(c):
                    if (nxt >> w) & 1: g[S | (1 << w), w] += cg
    Q = np.zeros((c, c), np.int64)
    for S in range(1, full):
        R = full ^ S
        for x in range(c):
            fx = f[S, x]
            if fx == 0: continue
            for y in range(c):
                if (R >> y) & 1: Q[x, y] += fx * g[R, y]
    h = 0
    for x in range(c): h += f[full, x]
    bad1 = 0; bad2 = 0; bad3 = 0; eq = 0
    for B in range(1 << c):
        t = 0
        for x in range(c):
            if (B >> x) & 1:
                t += f[full, x]
                for y in range(c):
                    if not (B >> y) & 1: t += Q[x, y]
            else:
                t += g[full, x]
        # second form
        t2 = h
        for x in range(c):
            if (B >> x) & 1:
                up = f[full, x]
                for y in range(c):
                    if y != x: up += Q[x, y]
                t2 += up - g[full, x]
                for y in range(x + 1, c):
                    if (B >> y) & 1: t2 -= Q[x, y] + Q[y, x]
        mm = cm_scatter
        for x in range(c):
            if (B >> x) & 1: mm |= cutbit[x]
        if Hbig[mm] != t: bad1 += 1
        if t2 != t: bad2 += 1
    for x in range(c):
        for y in range(x + 1, c):
            rest = full ^ (1 << x) ^ (1 << y)
            hr = 1
            if rest != 0:
                hr = 0
                for z in range(c): hr += f[rest, z]
            w = Q[x, y] + Q[y, x]
            if w < 2 * hr: bad3 += 1
            if w == 2 * hr: eq += 1
    return bad1, bad2, bad3, eq

@njit(parallel=True, cache=True)
def run(c, pi, pj, Hbig, scatter, cutbit):
    M = 1 << (c * (c - 1) // 2)
    out = np.zeros((M, 4), np.int64)
    for m in prange(M):
        b1, b2, b3, eq = one(c, m, pi, pj, Hbig, scatter[m], cutbit)
        out[m, 0] = b1; out[m, 1] = b2; out[m, 2] = b3; out[m, 3] = eq
    return out

for c in (3, 4, 5, 6):
    n = c + 1
    pairs_big = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pib = np.array([p[0] for p in pairs_big], np.int64); pjb = np.array([p[1] for p in pairs_big], np.int64)
    Hbig, _ = table(n, pib, pjb)
    pairs = [(i, j) for i in range(c) for j in range(i + 1, c)]
    pi = np.array([p[0] for p in pairs], np.int64); pj = np.array([p[1] for p in pairs], np.int64)
    core_idx = [k for k, p in enumerate(pairs_big) if p[1] != c]
    M = 1 << len(pairs)
    scatter = np.zeros(M, np.int64)
    for t, k in enumerate(core_idx):
        scatter |= ((np.arange(M) >> t) & 1) << k
    cutbit = np.array([1 << pairs_big.index((x, c)) for x in range(c)], np.int64)
    out = run(c, pi, pj, Hbig, scatter, cutbit)
    print(f"core order {c}: {M} labelled cores x {1 << c} cuts: cut-law mismatches with the table: {int(out[:,0].sum())}; "
          f"mismatches between the two forms: {int(out[:,1].sum())}; pairs with W_xy < 2H(C-x-y): {int(out[:,2].sum())}; "
          f"pairs with equality: {int(out[:,3].sum())} of {M * c * (c - 1) // 2}")
