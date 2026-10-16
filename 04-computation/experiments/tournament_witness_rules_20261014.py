#!/usr/bin/env python3
"""Which in-neighbour witnesses the vertex lemma?  Natural candidate rules all fail (HYP-9280; session opus-2026-10-14-S36, post-audit).

For every vertex w of every tournament in the sweep, with p(u, w) = Pr[u before w] in a uniform random Hamiltonian
path and L(u) = Pr[u immediately precedes w]:
  W  the predecessor-weighted average  sum_u L(u) (2p(u, w) - 1) > 0 ?
  A  in-neighbour of minimum out-degree;      D  minimum score inside T[U] (U = in-neighbourhood of w);
  E  maximum score inside T[U];               G  beats the fewest out-neighbours of w;
  P  the most frequent immediate predecessor of w.
For the rules, 'some-fail' counts (T, w) where NO tied choice is a witness (p > 1/2) and 'all-fail' where some tied
choice is not.  Modes: 'lab NMAX' (all labelled strong tournaments) and 'near n k' (transitive + k reversed arcs).
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

NR = 5


@nb.njit(cache=True)
def rule_key(r, u, adj, inm, outm, La):
    if r == 0:
        return popc(adj[u])
    if r == 1:
        return popc(adj[u] & inm)
    if r == 2:
        return -popc(adj[u] & inm)
    if r == 3:
        return popc(adj[u] & outm)
    return -La[u]


@nb.njit(cache=True)
def eval_T(n, adj, F, B, Cw, La, res):
    full = (1 << n) - 1
    tables(n, adj, F, B)
    H = 0
    for v in range(n):
        H += F[full, v]
    for w in range(n):
        inm = 0
        for u in range(n):
            if u != w and (adj[u] >> w) & 1:
                inm |= 1 << u
        outm = adj[w]
        for u in range(n):
            Cw[u] = 0; La[u] = 0
        for S in range(1, full + 1):
            if (S >> w) & 1:
                continue
            b = B[full ^ S, w]
            if b == 0:
                continue
            a = 0
            for x in range(n):
                if (S >> x) & 1 and (inm >> x) & 1:
                    a += F[S, x]
                    La[x] += F[S, x] * b
            if a == 0:
                continue
            for u in range(n):
                if (S >> u) & 1:
                    Cw[u] += a * b
        res[NR + 1] += 1
        sw = 0.0
        for u in range(n):
            if (inm >> u) & 1:
                sw += float(La[u]) * float(2 * Cw[u] - H)
        if sw <= 0:
            res[NR] += 1
        for r in range(NR):
            bestkey = 1 << 62
            for u in range(n):
                if (inm >> u) & 1:
                    kk = rule_key(r, u, adj, inm, outm, La)
                    if kk < bestkey:
                        bestkey = kk
            some = False; allw = True
            for u in range(n):
                if (inm >> u) & 1 and rule_key(r, u, adj, inm, outm, La) == bestkey:
                    if 2 * Cw[u] > H:
                        some = True
                    else:
                        allw = False
            if not some:
                res[r] += 1
            if not allw:
                res[NR + 2 + r] += 1


@nb.njit(parallel=True, cache=True)
def sweep_lab(n, pu, pv, nchunk):
    E = pu.shape[0]; total = 1 << E; N = 1 << n
    res = np.zeros((nchunk, 2 * NR + 2), dtype=np.int64)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        F = np.zeros((N, n), dtype=np.int64); B = np.zeros((N, n), dtype=np.int64)
        adj = np.zeros(n, dtype=np.int64); Cw = np.zeros(n, dtype=np.int64); La = np.zeros(n, dtype=np.int64)
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
            eval_T(n, adj, F, B, Cw, La, res[c])
    return res.sum(axis=0)


@nb.njit(parallel=True, cache=True)
def sweep_near(n, arcsets, nthreads):
    m = arcsets.shape[0]; k = arcsets.shape[1]; N = 1 << n
    res = np.zeros((nthreads, 2 * NR + 2), dtype=np.int64)
    chunk = (m + nthreads - 1) // nthreads
    for t in nb.prange(nthreads):
        F = np.zeros((N, n), dtype=np.int64); B = np.zeros((N, n), dtype=np.int64)
        adj = np.zeros(n, dtype=np.int64); Cw = np.zeros(n, dtype=np.int64); La = np.zeros(n, dtype=np.int64)
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
            eval_T(n, adj, F, B, Cw, La, res[t])
    return res.sum(axis=0)


NAMES = ['A min out-degree', 'D min score in T[U]', 'E max score in T[U]', 'G beats fewest out-nbrs of w', 'P most frequent predecessor']


def report(tag, r, t0):
    print(f'{tag}: (T, w) pairs {int(r[NR + 1])}; predecessor-weighted average <= 1/2: {int(r[NR])}; ' +
          '; '.join(f'{NAMES[i]}: some-fail {int(r[i])}, all-fail {int(r[NR + 2 + i])}' for i in range(NR)) +
          f'  [{time.time() - t0:.1f} s]', flush=True)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'near':
        jobs = [('near', int(sys.argv[2]), int(sys.argv[3]))]
    elif len(sys.argv) > 1 and sys.argv[1] == 'lab':
        jobs = [('lab', n, 0) for n in range(4, int(sys.argv[2]) + 1)]
    else:
        jobs = [('lab', n, 0) for n in range(4, 8)] + [('near', n, 2) for n in range(9, 15)] + [('near', n, 3) for n in range(9, 13)]
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
