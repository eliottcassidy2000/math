#!/usr/bin/env python3
"""Direction (i): balanced arcs.  For every labelled tournament T of order n (exhaustive, exact) and every arc x -> y,
mu(x -> y) = number of Hamiltonian paths through the arc.  The arc is BALANCED if 2 mu = H - 1 (then T with this arc
reversed, together with {x, y}, is a seed: S93 seeds/PROOFS.txt Theorem 1) and ANTI-BALANCED if 2 mu = H + 1.
Reported per order: number of tournaments with a balanced / anti-balanced arc, the H-values that admit one, and the
H-values of strong tournaments that admit none.   usage: python balanced_arcs.py n"""
import sys, time
import numpy as np
from numba import njit, prange

@njit(cache=True)
def stats(c, m, pi, pj):
    adj = np.zeros(c, np.int64); radj = np.zeros(c, np.int64)
    na = c * (c - 1) // 2
    for k in range(na):
        if (m >> k) & 1:
            a = pi[k]; b = pj[k]
        else:
            a = pj[k]; b = pi[k]
        adj[a] |= 1 << b
        radj[b] |= 1 << a
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
    h = 0
    for x in range(c): h += f[full, x]
    nb = 0; nab = 0
    for x in range(c):
        for y in range(c):
            if not (adj[x] >> y) & 1: continue
            mu = 0
            for S in range(1, full):
                if (S >> x) & 1 and not (S >> y) & 1:
                    fx = f[S, x]
                    if fx: mu += fx * g[full ^ S, y]
            if 2 * mu == h - 1: nb += 1
            if 2 * mu == h + 1: nab += 1
    # strong?
    sc = np.zeros(c, np.int64)
    for v in range(c):
        t = adj[v]; cnt = 0
        while t:
            cnt += 1; t &= t - 1
        sc[v] = cnt
    sc.sort(); strong = 1; acc = 0
    for k in range(1, c):
        acc += sc[k - 1]
        if acc == k * (k - 1) // 2:
            strong = 0; break
    return h, nb, nab, strong

@njit(parallel=True, cache=True)
def scan(c, pi, pj):
    M = 1 << (c * (c - 1) // 2)
    Hh = np.zeros(M, np.int64); NB = np.zeros(M, np.int32); NAB = np.zeros(M, np.int32); ST = np.zeros(M, np.int8)
    for m in prange(M):
        h, nb, nab, st = stats(c, m, pi, pj)
        Hh[m] = h; NB[m] = nb; NAB[m] = nab; ST[m] = st
    return Hh, NB, NAB, ST

def main():
    n = int(sys.argv[1])
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pi = np.array([p[0] for p in pairs], np.int64); pj = np.array([p[1] for p in pairs], np.int64)
    t0 = time.time()
    Hh, NB, NAB, ST = scan(n, pi, pj)
    M = len(Hh)
    print(f"n = {n}: {M} labelled tournaments ({time.time()-t0:.1f}s); with a balanced arc: {int((NB > 0).sum())} "
          f"({(NB > 0).mean():.4f}); with an anti-balanced arc: {int((NAB > 0).sum())} ({(NAB > 0).mean():.4f}); "
          f"strong: {int(ST.sum())}, strong with a balanced arc: {int(((NB > 0) & (ST == 1)).sum())} "
          f"({((NB > 0) & (ST == 1)).sum() / max(ST.sum(), 1):.4f}), strong with an anti-balanced arc: "
          f"{int(((NAB > 0) & (ST == 1)).sum())}")
    vals = sorted(set(Hh.tolist()))
    hb = sorted(set(Hh[NB > 0].tolist())); hab = sorted(set(Hh[NAB > 0].tolist()))
    print(f"   H-values at this order: {len(vals)}; values admitting a balanced arc: {len(hb)}; an anti-balanced arc: {len(hab)}")
    print(f"   values WITHOUT a balanced arc at this order: {[v for v in vals if v not in set(hb)]}")
    print(f"   values WITHOUT an anti-balanced arc at this order: {[v for v in vals if v not in set(hab)]}")
    print(f"   values with neither: {[v for v in vals if v not in set(hb) and v not in set(hab)]}")
    print(f"   max number of balanced arcs in one tournament: {int(NB.max())}; of anti-balanced: {int(NAB.max())}")

if __name__ == '__main__':
    main()
