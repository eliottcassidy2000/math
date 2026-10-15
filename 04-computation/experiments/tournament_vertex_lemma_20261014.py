#!/usr/bin/env python3
"""Vertex-level forms of the majority-cycle conjecture (HYP-9280), exhaustive over labelled strong tournaments n <= NMAX.

For each strong tournament T and vertex w:
  in_max(w)  = max over in-neighbours u of p(u before w)        (vertex lemma: > 1/2?)
  D_in(w)    = sum over in-neighbours u of (2 p(u before w) - 1)  (averaged in-arc forwardness: > 0?)
  D_out(w)   = sum over out-neighbours v of (2 p(w before v) - 1)
  f(w)       = Pr[w first]                                         (PROVED < 1/2 for non-sources)
Reports the minima of in_max and D_in, D_out, and the maximum of f over all strong T and vertices.
Session opus-2026-10-14-S36.
"""
import sys, time, os
import importlib.util
import numpy as np
import numba as nb
from fractions import Fraction

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('mc', os.path.join(HERE, 'tournament_majority_cycle_20261014.py'))
mc = importlib.util.module_from_spec(spec)
_argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(mc); sys.argv = _argv
pair_counts, is_strong, fb_tables = mc.pair_counts, mc.is_strong, mc.fb_tables


@nb.njit(parallel=True, cache=True)
def sweep(n, pu, pv, nchunk):
    E = pu.shape[0]; total = 1 << E
    im_n = np.full(nchunk, 2, dtype=np.int64); im_d = np.ones(nchunk, dtype=np.int64)
    di_n = np.full(nchunk, 10 ** 12, dtype=np.int64); di_d = np.ones(nchunk, dtype=np.int64)
    do_n = np.full(nchunk, 10 ** 12, dtype=np.int64); do_d = np.ones(nchunk, dtype=np.int64)
    f_n = np.zeros(nchunk, dtype=np.int64); f_d = np.ones(nchunk, dtype=np.int64)
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
            F, B = fb_tables(n, adj)
            full = (1 << n) - 1
            H, C = pair_counts(n, adj)
            for w in range(n):
                best = 0; din = 0; dout = 0
                for u in range(n):
                    if u == w:
                        continue
                    if (adj[u] >> w) & 1:
                        if C[u, w] > best:
                            best = C[u, w]
                        din += 2 * C[u, w] - H
                    else:
                        dout += 2 * C[w, u] - H
                if best * im_d[c] < im_n[c] * H:
                    im_n[c] = best; im_d[c] = H
                if din * di_d[c] < di_n[c] * H:
                    di_n[c] = din; di_d[c] = H
                if dout * do_d[c] < do_n[c] * H:
                    do_n[c] = dout; do_d[c] = H
                fw = B[full, w]   # paths starting at w
                if fw * f_d[c] > f_n[c] * H:
                    f_n[c] = fw; f_d[c] = H
    return im_n, im_d, di_n, di_d, do_n, do_d, f_n, f_d


if __name__ == '__main__':
    for n in range(3, NMAX + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in pairs], dtype=np.int64); pv = np.array([p[1] for p in pairs], dtype=np.int64)
        t0 = time.time()
        a, b, c_, d, e_, f_, g, h = sweep(n, pu, pv, 512)
        mn = lambda xs, ys: min(Fraction(int(x), int(y)) for x, y in zip(xs, ys))
        mx = lambda xs, ys: max(Fraction(int(x), int(y)) for x, y in zip(xs, ys))
        print(f'n = {n}: min in_max = {mn(a, b)}; min D_in = {mn(c_, d)}; min D_out = {mn(e_, f_)}; max Pr[w first] = {mx(g, h)}  '
              f'(prediction (2^(n-3)+1)/(2^(n-2)+1) = {Fraction(2 ** (n - 3) + 1, 2 ** (n - 2) + 1) if n >= 3 else None}, 1/H(T_n) = {Fraction(1, 2 ** (n - 2) + 1)})  [{time.time() - t0:.1f} s]', flush=True)
