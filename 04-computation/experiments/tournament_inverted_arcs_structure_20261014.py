#!/usr/bin/env python3
"""Inverted arcs and local statistics of uniform random Hamiltonian paths, all labelled strong tournaments n <= NMAX (HYP-9280).

p(x, y) = Pr[x before y] in a uniform random Hamiltonian path; H(T) is odd (Redei), so p != 1/2.  An arc x -> y of T
is forward if p(x, y) > 1/2 and inverted otherwise; F(T) is the digraph of forward arcs.  Exhaustive over labelled
strong tournaments, this script reports:
  (a) the number of T whose forward digraph F(T) is NOT strongly connected (HYP-9280 predicts 0);
  (b) the distribution of the number of inverted arcs; the histogram of (outdeg(head), indeg(tail)) over inverted
      arcs x -> y; the 'squeezed triangle' check for type (1, 1): the unique out-neighbour z of y is the unique
      in-neighbour of x (y -> z -> x); the least p over inverted arcs;
  (c) the false strengthening "max over in-neighbours u of p(u before w) + Pr[w first] >= 1": the number of
      (T, w) pairs violating it and the number of T with at least one violating vertex.
Session opus-2026-10-14-S36.  Shipped after audit item E9 (these checks previously lived in scratch scripts);
the audit's independent code (s36_audit_20261014/) reproduces (a)-(c) for n <= 7.
"""
import sys, os, time
import importlib.util
from fractions import Fraction
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('mc', os.path.join(HERE, 'tournament_majority_cycle_20261014.py'))
mc = importlib.util.module_from_spec(spec)
_argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(mc); sys.argv = _argv
fb_tables, is_strong, pair_counts = mc.fb_tables, mc.is_strong, mc.pair_counts

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7


@nb.njit(cache=True)
def popc(m):
    s = 0
    while m:
        m &= m - 1
        s += 1
    return s


@nb.njit(parallel=True, cache=True)
def sweep(n, pu, pv, nchunk):
    E = pu.shape[0]
    total = 1 << E
    nstrong = np.zeros(nchunk, dtype=np.int64)
    f_not_strong = np.zeros(nchunk, dtype=np.int64)
    ninv_hist = np.zeros((nchunk, 64), dtype=np.int64)
    type_hist = np.zeros((nchunk, n, n), dtype=np.int64)
    sq_ok = np.zeros(nchunk, dtype=np.int64)
    sq_bad = np.zeros(nchunk, dtype=np.int64)
    minnum = np.ones(nchunk, dtype=np.int64)
    minden = np.ones(nchunk, dtype=np.int64)
    viol_pairs = np.zeros(nchunk, dtype=np.int64)
    viol_tourn = np.zeros(nchunk, dtype=np.int64)
    chunk = (total + nchunk - 1) // nchunk
    full = (1 << n) - 1
    for c in nb.prange(nchunk):
        adj = np.zeros(n, dtype=np.int64)
        fadj = np.zeros(n, dtype=np.int64)
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
            H, C = pair_counts(n, adj)
            # (a) forward digraph
            for x in range(n):
                fadj[x] = 0
                for y in range(n):
                    if x != y and (adj[x] >> y) & 1 and 2 * C[x, y] > H:
                        fadj[x] |= 1 << y
            if not is_strong(n, fadj):
                f_not_strong[c] += 1
            # (b) inverted arcs
            ninv = 0
            for x in range(n):
                for y in range(n):
                    if x != y and (adj[x] >> y) & 1 and 2 * C[x, y] < H:
                        ninv += 1
                        od = popc(adj[y])
                        idg = 0
                        for u in range(n):
                            if (adj[u] >> x) & 1:
                                idg += 1
                        type_hist[c, od, idg] += 1
                        if od == 1 and idg == 1:
                            z = 0
                            while not (adj[y] >> z) & 1:
                                z += 1
                            if (adj[z] >> x) & 1:
                                sq_ok[c] += 1
                            else:
                                sq_bad[c] += 1
                        if C[x, y] * minden[c] < minnum[c] * H:
                            minnum[c] = C[x, y]
                            minden[c] = H
            ninv_hist[c, ninv] += 1
            # (c) best in-neighbour precedence + Pr[w first] >= 1 ?
            F, B = fb_tables(n, adj)
            any_bad = False
            for w in range(n):
                best = 0
                for u in range(n):
                    if u != w and (adj[u] >> w) & 1 and C[u, w] > best:
                        best = C[u, w]
                if best + B[full, w] < H:
                    viol_pairs[c] += 1
                    any_bad = True
            if any_bad:
                viol_tourn[c] += 1
    return (nstrong, f_not_strong, ninv_hist, type_hist, sq_ok, sq_bad, minnum, minden, viol_pairs, viol_tourn)


if __name__ == '__main__':
    for n in range(3, NMAX + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in pairs], dtype=np.int64)
        pv = np.array([p[1] for p in pairs], dtype=np.int64)
        t0 = time.time()
        ns, fns, nh, th, sok, sbad, mn, md, vp, vt = sweep(n, pu, pv, 512)
        nh = nh.sum(axis=0)
        th = th.sum(axis=0)
        dist = {k: int(nh[k]) for k in range(64) if nh[k]}
        types = {(a, b): int(th[a, b]) for a in range(n) for b in range(n) if th[a, b]}
        if th.sum():
            k = min(range(len(mn)), key=lambda i: Fraction(int(mn[i]), int(md[i])))
            minp = str(Fraction(int(mn[k]), int(md[k])))
        else:
            minp = '-'
        print(f'n = {n}: strong labelled {int(ns.sum())}; F(T) not strongly connected: {int(fns.sum())}; '
              f'#inverted arcs distribution {dist}; inverted-arc types (outdeg head, indeg tail) {types}; '
              f'squeezed-triangle check for type (1,1): {int(sok.sum())} hold, {int(sbad.sum())} fail; '
              f'least p over inverted arcs {minp}; '
              f'"best in-neighbour p + Pr[w first] >= 1" violated at {int(vp.sum())} (T, w) pairs '
              f'in {int(vt.sum())} tournaments  [{time.time() - t0:.1f} s]', flush=True)
