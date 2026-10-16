#!/usr/bin/env python3
"""The per-in-neighbour refinement of (SS) at k = 1 is false (HYP-9280; session opus-2026-10-14-S36, post-audit).

For a vertex w of in-degree d >= 3 and an in-neighbour u let
    a_u = #{Hamiltonian paths in which u is the ONLY in-neighbour of w before w},
    b_u = #{Hamiltonian paths in which u is the ONLY in-neighbour of w after w}.
(SS) at k = 1 reads  sum_u b_u + #(N = d) >= sum_u a_u + #(N = 0).  The refinement a_u <= b_u for every u would give it
from the proved k = 0 case.  This script counts the triples (T, w, u) with a_u > b_u over all labelled strong
tournaments with n <= NMAX: the refinement fails massively, so (SS) at k = 1 (no counterexample known) needs a global
argument.
"""
import sys, os, time
import importlib.util
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('nts', os.path.join(HERE, 'tournament_near_transitive_sweep_20261014.py'))
nts = importlib.util.module_from_spec(spec)
_argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(nts); sys.argv = _argv
is_strong_bits, popc, tables = nts.is_strong_bits, nts.popc, nts.tables


@nb.njit(cache=True)
def eval_T(n, adj, F, B, A_, B_, res):
    full = (1 << n) - 1
    tables(n, adj, F, B)
    for w in range(n):
        inm = 0
        for u in range(n):
            if u != w and (adj[u] >> w) & 1:
                inm |= 1 << u
        d = popc(inm)
        if d < 3:
            continue
        for u in range(n):
            A_[u] = 0; B_[u] = 0
        for S in range(1, full + 1):
            if (S >> w) & 1:
                continue
            b = B[full ^ S, w]
            if b == 0:
                continue
            si = S & inm
            c = popc(si)
            if c != 1 and c != d - 1:
                continue
            a = 0
            for x in range(n):
                if (S >> x) & 1 and (inm >> x) & 1:
                    a += F[S, x]
            if a == 0:
                continue
            if c == 1:
                u = 0
                while not (si >> u) & 1:
                    u += 1
                A_[u] += a * b
            if c == d - 1:
                rest = inm & ~si
                u = 0
                while not (rest >> u) & 1:
                    u += 1
                B_[u] += a * b
        for u in range(n):
            if (inm >> u) & 1:
                res[0] += 1
                if A_[u] > B_[u]:
                    res[1] += 1


@nb.njit(parallel=True, cache=True)
def sweep_lab(n, pu, pv, nchunk):
    E = pu.shape[0]; total = 1 << E; N = 1 << n
    res = np.zeros((nchunk, 2), dtype=np.int64)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        F = np.zeros((N, n), dtype=np.int64); B = np.zeros((N, n), dtype=np.int64)
        adj = np.zeros(n, dtype=np.int64); A_ = np.zeros(n, dtype=np.int64); B_ = np.zeros(n, dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for i in range(n):
                adj[i] = 0
            for e in range(E):
                if (code >> e) & 1:
                    adj[pu[e]] |= 1 << pv[e]
                else:
                    adj[pv[e]] |= 1 << pu[e]
            if not is_strong_bits(n, adj):
                continue
            eval_T(n, adj, F, B, A_, B_, res[c])
    return res.sum(axis=0)


if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for n in range(4, NMAX + 1):
        prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in prs], dtype=np.int64); pv = np.array([p[1] for p in prs], dtype=np.int64)
        t0 = time.time()
        r = sweep_lab(n, pu, pv, 256)
        frac = r[1] / r[0] if r[0] else 0.0
        print(f'n = {n} labelled strong: (T, w, u) with in-degree d >= 3: {int(r[0])}; refinement a_u <= b_u fails at {int(r[1])} '
              f'({100 * frac:.1f}%)  [{time.time() - t0:.1f} s]', flush=True)
