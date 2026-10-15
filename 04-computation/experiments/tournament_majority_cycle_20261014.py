#!/usr/bin/env python3
"""Conjecture M: the majority tournament of the uniform Hamiltonian-path measure is transitive only for transitive T.

For a tournament T, M(T) orients each pair x, y from x to y when p(x before y) > 1/2.  H(T) is odd (Redei), so p is
never 1/2 and M(T) is a tournament.  If M(T) has a directed 3-cycle x -> y -> z -> x, then
p(x,y) + p(y,z) + p(z,x) <= 2 (every linear order satisfies at most two of the three cyclic precedences), so one of
them lies in (1/2, 2/3]: a balanced pair.  Hence Conjecture M implies the tournament 1/3-2/3 conjecture (HYP-9277
part 1).  By the component reduction it suffices to treat strong tournaments.
This script: exhaustive over labelled strong tournaments with n <= NMAX, reports
  - tournaments with M(T) transitive (counterexamples to Conjecture M),
  - the number of inverted arcs (arcs x -> y of T with p(x,y) < 1/2) -- maximum and distribution,
  - whether some directed 3-cycle of T is entirely majority-forward (a stronger certificate).
Session opus-2026-10-14-S36.
"""
import sys, time
import numpy as np
import numba as nb

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
NMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 3


@nb.njit(cache=True)
def fb_tables(n, adj):
    N = 1 << n
    F = np.zeros((N, n), dtype=np.int64)
    B = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        F[1 << v, v] = 1
        B[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            if not (S >> v) & 1:
                continue
            c = F[S, v]
            if c:
                for w in range(n):
                    if not (S >> w) & 1 and (adj[v] >> w) & 1:
                        F[S | (1 << w), w] += c
            c = B[S, v]
            if c:
                for w in range(n):
                    if not (S >> w) & 1 and (adj[w] >> v) & 1:
                        B[S | (1 << w), w] += c
    return F, B


@nb.njit(cache=True)
def is_strong(n, adj):
    full = (1 << n) - 1
    seen = 1; fr = 1
    while fr:
        nx = 0
        for v in range(n):
            if (fr >> v) & 1:
                nx |= adj[v]
        nx &= ~seen; seen |= nx; fr = nx
    if seen != full:
        return False
    seen = 1; fr = 1
    while fr:
        nx = 0
        for v in range(n):
            if (fr >> v) & 1:
                for u in range(n):
                    if (adj[u] >> v) & 1:
                        nx |= 1 << u
        nx &= ~seen; seen |= nx; fr = nx
    return seen == full


@nb.njit(cache=True)
def pair_counts(n, adj):
    F, B = fb_tables(n, adj)
    full = (1 << n) - 1
    H = 0
    for v in range(n):
        H += F[full, v]
    C = np.zeros((n, n), dtype=np.int64)
    for x in range(n):
        for y in range(n):
            if x == y:
                continue
            cnt = 0
            for S in range(1, full + 1):
                if not (S >> x) & 1 or (S >> y) & 1:
                    continue
                b = B[full ^ S, y]
                if b == 0:
                    continue
                for u in range(n):
                    if (S >> u) & 1 and (adj[u] >> y) & 1:
                        cnt += F[S, u] * b
            C[x, y] = cnt
    return H, C


@nb.njit(cache=True)
def majority_info(n, adj):
    H, C = pair_counts(n, adj)
    madj = np.zeros(n, dtype=np.int64)
    ninv = 0
    for x in range(n):
        for y in range(n):
            if x != y and 2 * C[x, y] > H:
                madj[x] |= 1 << y
                if not (adj[x] >> y) & 1:
                    ninv += 1      # T has y -> x but x usually first: the arc y -> x is inverted
    # M transitive iff its score sequence is 0..n-1
    seen = 0
    trans = True
    for v in range(n):
        s = 0
        m = madj[v]
        while m:
            m &= m - 1
            s += 1
        if (seen >> s) & 1:
            trans = False
        seen |= 1 << s
    # forward 3-cycle of T (all arcs majority-forward)?
    fwd3 = False
    for x in range(n):
        for y in range(n):
            if y == x or not ((adj[x] >> y) & 1 and (madj[x] >> y) & 1):
                continue
            for z in range(n):
                if z == x or z == y:
                    continue
                if (adj[y] >> z) & 1 and (madj[y] >> z) & 1 and (adj[z] >> x) & 1 and (madj[z] >> x) & 1:
                    fwd3 = True
                    break
            if fwd3:
                break
        if fwd3:
            break
    return trans, ninv, fwd3


@nb.njit(parallel=True, cache=True)
def sweep(n, pu, pv, nchunk):
    E = pu.shape[0]
    total = 1 << E
    nstrong = np.zeros(nchunk, dtype=np.int64)
    ntrans = np.zeros(nchunk, dtype=np.int64)
    nofwd3 = np.zeros(nchunk, dtype=np.int64)
    hist = np.zeros((nchunk, 64), dtype=np.int64)
    witness = np.full(nchunk, -1, dtype=np.int64)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        adj = np.zeros(n, dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for i in range(n):
                adj[i] = 0
            for e in range(E):
                if (code >> e) & 1:
                    adj[pu[e]] |= 1 << pv[e]
                else:
                    adj[pv[e]] |= 1 << pu[e]
            if not is_strong(n, adj):
                continue
            nstrong[c] += 1
            trans, ninv, fwd3 = majority_info(n, adj)
            if trans:
                ntrans[c] += 1
                witness[c] = code
            if not fwd3:
                nofwd3[c] += 1
            hist[c, min(ninv, 63)] += 1
    return nstrong, ntrans, nofwd3, hist, witness


if __name__ == '__main__':
    for n in range(max(3, NMIN), NMAX + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in pairs], dtype=np.int64); pv = np.array([p[1] for p in pairs], dtype=np.int64)
        t0 = time.time()
        ns, nt, nf, hist, wit = sweep(n, pu, pv, 1024)
        h = hist.sum(axis=0)
        dist = {k: int(h[k]) for k in range(64) if h[k]}
        print(f'n = {n}: strong labelled {int(ns.sum())}; M(T) transitive: {int(nt.sum())} '
              f'(witness codes {[int(w) for w in wit if w >= 0][:3]}); no majority-forward 3-cycle of T: {int(nf.sum())}; '
              f'#inverted arcs distribution {dist}  [{time.time() - t0:.1f} s]', flush=True)
