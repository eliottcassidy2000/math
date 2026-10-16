#!/usr/bin/env python3
"""The vertex lemma at vertices of out-degree 1 (HYP-9280; session opus-2026-10-14-S36, post-audit).

Let w have the unique out-neighbour o, U = V - {w, o} (all in-neighbours of w), U+ = {u in U : o -> u},
U- = {u in U : u -> o}.  This is where the averaged form D_in > 0 first fails (n = 11) and where (SS) fails (n = 9).
For every such (T, w) the script checks:
  (1) PROVED case: if U- is empty, then C(u, w) - C(w, u) = H(T[U]) for EVERY u in U (exact identity), so every
      in-neighbour is a witness;  counted: cases with U- empty, identity failures.
  (2) a witness (p(u before w) > 1/2) exists in U+;
  (3) if U- is nonempty: the average of (2p - 1) over U- is positive;  and a witness exists in U-;
  (4) the average of (2p - 1) over U+ (reported: number of failures, i.e. this average is NOT the mechanism).
Modes:  'lab NMAX'   all labelled strong tournaments n <= NMAX;
        'near n k'   all strong 'transitive + k reversed arcs' tournaments (see tournament_near_transitive_sweep).
"""
import sys, os, time, itertools
import importlib.util
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('nts', os.path.join(HERE, 'tournament_near_transitive_sweep_20261014.py'))
nts = importlib.util.module_from_spec(spec)
_argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(nts); sys.argv = _argv
is_strong_bits, popc, tables = nts.is_strong_bits, nts.popc, nts.tables


@nb.njit(cache=True)
def eval_T(n, adj, F, B, Cw, res):
    full = (1 << n) - 1
    have = False
    for w in range(n):
        if popc(adj[w]) == 1:
            have = True
    if not have:
        return
    tables(n, adj, F, B)
    H = 0
    for v in range(n):
        H += F[full, v]
    for w in range(n):
        if popc(adj[w]) != 1:
            continue
        o = 0
        while not (adj[w] >> o) & 1:
            o += 1
        U = full & ~(1 << w) & ~(1 << o)
        um = 0
        for u in range(n):
            if (U >> u) & 1 and (adj[u] >> o) & 1:
                um |= 1 << u
        for u in range(n):
            Cw[u] = 0
        for S in range(0, full + 1):
            if (S >> w) & 1:
                continue
            b = B[full ^ S, w]
            if b == 0:
                continue
            if S == 0:
                a = 1
            else:
                a = 0
                for x in range(n):
                    if (S >> x) & 1 and (U >> x) & 1:
                        a += F[S, x]
            if a == 0:
                continue
            for u in range(n):
                if (S >> u) & 1:
                    Cw[u] += a * b
        res[0] += 1
        if um == 0:
            res[1] += 1
            HU = 0
            for v in range(n):
                HU += F[U, v]
            for u in range(n):
                if (U >> u) & 1 and (2 * Cw[u] - H) != HU:
                    res[2] += 1
                    break
        wp = False; wm = False; sp = 0; sm = 0
        for u in range(n):
            if (U >> u) & 1:
                mg = 2 * Cw[u] - H
                if (um >> u) & 1:
                    sm += mg
                    if mg > 0:
                        wm = True
                else:
                    sp += mg
                    if mg > 0:
                        wp = True
        if not wp:
            res[3] += 1
        if um != 0:
            res[4] += 1
            if sm <= 0:
                res[5] += 1
            if not wm:
                res[6] += 1
        if sp <= 0:
            res[7] += 1


@nb.njit(parallel=True, cache=True)
def sweep_lab(n, pu, pv, nchunk):
    E = pu.shape[0]; total = 1 << E; N = 1 << n
    res = np.zeros((nchunk, 8), dtype=np.int64)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        F = np.zeros((N, n), dtype=np.int64); B = np.zeros((N, n), dtype=np.int64)
        adj = np.zeros(n, dtype=np.int64); Cw = np.zeros(n, dtype=np.int64)
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
            eval_T(n, adj, F, B, Cw, res[c])
    return res.sum(axis=0)


@nb.njit(parallel=True, cache=True)
def sweep_near(n, arcsets, nthreads):
    m = arcsets.shape[0]; k = arcsets.shape[1]; N = 1 << n
    res = np.zeros((nthreads, 8), dtype=np.int64)
    chunk = (m + nthreads - 1) // nthreads
    for t in nb.prange(nthreads):
        F = np.zeros((N, n), dtype=np.int64); B = np.zeros((N, n), dtype=np.int64)
        adj = np.zeros(n, dtype=np.int64); Cw = np.zeros(n, dtype=np.int64)
        for idx in range(t * chunk, min(m, (t + 1) * chunk)):
            for a in range(n):
                adj[a] = 0
                for b in range(a):
                    adj[a] |= 1 << b
            for e in range(k):
                lo = arcsets[idx, e] // n; hi = arcsets[idx, e] % n
                adj[hi] &= ~(1 << lo); adj[lo] |= 1 << hi
            if not is_strong_bits(n, adj):
                continue
            eval_T(n, adj, F, B, Cw, res[t])
    return res.sum(axis=0)


def report(tag, r, t0):
    print(f'{tag}: out-degree-1 (T, w) pairs {int(r[0])}; U- empty {int(r[1])} (identity C(u,w) - C(w,u) = H(T[U]) fails: {int(r[2])}); '
          f'no witness in U+: {int(r[3])}; U- nonempty {int(r[4])}: average over U- <= 1/2 in {int(r[5])}, no witness in U- in {int(r[6])}; '
          f'average over U+ <= 1/2 in {int(r[7])}  [{time.time() - t0:.1f} s]', flush=True)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'near':
        jobs = [('near', int(sys.argv[2]), int(sys.argv[3]))]
    elif len(sys.argv) > 1 and sys.argv[1] == 'lab':
        jobs = [('lab', n, 0) for n in range(4, int(sys.argv[2]) + 1)]
    else:
        jobs = ([('lab', n, 0) for n in range(4, 8)] + [('near', n, 2) for n in range(10, 17, 2)] +
                [('near', n, 3) for n in (10, 11, 12)] + [('near', 10, 4)])
    for mode, n, k in jobs:
        t0 = time.time()
        if mode == 'lab':
            prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
            pu = np.array([p[0] for p in prs], dtype=np.int64); pv = np.array([p[1] for p in prs], dtype=np.int64)
            report(f'n = {n} labelled', sweep_lab(n, pu, pv, 256), t0)
        else:
            cand = [lo * n + hi for lo in range(n) for hi in range(lo + 2, n)]
            combos = np.array(list(itertools.combinations(cand, k)), dtype=np.int64)
            report(f'n = {n}, transitive + {k} reversed arcs', sweep_near(n, combos, 12), t0)
