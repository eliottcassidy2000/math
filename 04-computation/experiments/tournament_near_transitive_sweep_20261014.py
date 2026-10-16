#!/usr/bin/env python3
"""HYP-9280 on all strong tournaments with small feedback arc number, every vertex (session opus-2026-10-14-S36, post-audit).

A tournament with feedback arc number k is a transitive order on positions 0..n-1 (higher beats lower) with k arcs
reversed; for an optimal order no reversed arc joins positions at distance 1.  For given (n, k) this script runs over
all k-sets of reversed arcs {lo -> hi : hi - lo >= 2} (with the k' < k sweeps this covers every tournament of feedback
arc number <= k), keeps the strong ones, and computes exactly (subset DP, integer counts):
  - M(T) transitive?  F(T) strongly connected?                                       (HYP-9280)
  - vertex lemma:      min over w of max over in-neighbours u of (2 p(u before w) - 1), also split by
                       min(indeg w, outdeg w) = 1, 2, >= 3
  - majority source:   min over w of max over all x of (2 p(x before w) - 1)
  - averaged form:     min over w of D_in(w) = sum over in-neighbours of (2p - 1); number of T with some D_in <= 0;
                       min of E[N] - d/2 by in-degree d <= 6
  - (SS) at each k <= 4: number of (T, w) with Pr[N >= d - k] < Pr[N <= k]   (only k with d - k > k)
Usage: python tournament_near_transitive_sweep_20261014.py n k [n k ...]   (default: the session's list)
"""
import sys, time, itertools
import numpy as np
import numba as nb


@nb.njit(cache=True)
def is_strong_bits(n, adj):
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
def popc(m):
    s = 0
    while m:
        m &= m - 1
        s += 1
    return s


@nb.njit(cache=True)
def tables(n, adj, F, B):
    """F[S, v] = # Hamiltonian paths of T[S] ending at v; B[S, v] = # starting at v."""
    N = 1 << n
    for S in range(N):
        for v in range(n):
            F[S, v] = 0; B[S, v] = 0
    for v in range(n):
        F[1 << v, v] = 1; B[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            if not (S >> v) & 1:
                continue
            c = F[S, v]
            if c:
                for x in range(n):
                    if not (S >> x) & 1 and (adj[v] >> x) & 1:
                        F[S | (1 << x), x] += c
            c = B[S, v]
            if c:
                for x in range(n):
                    if not (S >> x) & 1 and (adj[x] >> v) & 1:
                        B[S | (1 << x), x] += c


@nb.njit(cache=True)
def pair_counts(n, adj, F, B, C):
    """C[x, y] = # Hamiltonian paths with x before y.  Returns H."""
    full = (1 << n) - 1
    for x in range(n):
        for y in range(n):
            C[x, y] = 0
    H = 0
    for v in range(n):
        H += F[full, v]
    for S in range(1, full):
        R = full ^ S
        for y in range(n):
            if not (R >> y) & 1:
                continue
            b = B[R, y]
            if b == 0:
                continue
            g = 0
            for v in range(n):
                if (S >> v) & 1 and (adj[v] >> y) & 1:
                    g += F[S, v]
            if g == 0:
                continue
            val = g * b
            for x in range(n):
                if (S >> x) & 1:
                    C[x, y] += val
    return H


@nb.njit(cache=True)
def n_distribution(n, adj, F, B, w, Nd):
    """Nd[k] = # Hamiltonian paths with exactly k in-neighbours of w before w.  Returns in-degree d."""
    full = (1 << n) - 1
    inm = 0
    for u in range(n):
        if u != w and (adj[u] >> w) & 1:
            inm |= 1 << u
    for k in range(n + 1):
        Nd[k] = 0
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
                if (S >> x) & 1 and (inm >> x) & 1:
                    a += F[S, x]
        if a == 0:
            continue
        Nd[popc(S & inm)] += a * b
    return popc(inm)


@nb.njit(parallel=True, cache=True)
def sweep(n, arcsets, nthreads):
    m = arcsets.shape[0]; k = arcsets.shape[1]
    N = 1 << n
    # integer outputs per arc set: strong, H, M_trans, F_strong, vl, ms, dmin, vl_cls1, vl_cls2, vl_cls3
    out = np.zeros((m, 10), dtype=np.int64)
    ssf = np.zeros((nthreads, 5), dtype=np.int64)
    dbyd = np.ones((nthreads, 7)) * 9.0
    chunk = (m + nthreads - 1) // nthreads
    for t in nb.prange(nthreads):
        F = np.zeros((N, n), dtype=np.int64); B = np.zeros((N, n), dtype=np.int64)
        C = np.zeros((n, n), dtype=np.int64)
        adj = np.zeros(n, dtype=np.int64); fadj = np.zeros(n, dtype=np.int64)
        Nd = np.zeros(n + 1, dtype=np.int64)
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
            out[idx, 0] = 1
            tables(n, adj, F, B)
            H = pair_counts(n, adj, F, B, C)
            out[idx, 1] = H
            seen = 0; trans = 1
            for v in range(n):
                s = 0
                for x in range(n):
                    if x != v and 2 * C[v, x] > H:
                        s += 1
                if (seen >> s) & 1:
                    trans = 0
                seen |= 1 << s
            out[idx, 2] = trans
            for v in range(n):
                fadj[v] = 0
                for x in range(n):
                    if x != v and (adj[v] >> x) & 1 and 2 * C[v, x] > H:
                        fadj[v] |= 1 << x
            out[idx, 3] = 1 if is_strong_bits(n, fadj) else 0
            BIG = 1 << 62
            vl = BIG; ms = BIG; dm = BIG; c1 = BIG; c2 = BIG; c3 = BIG
            for w in range(n):
                best = -BIG; bestany = -BIG; dsum = 0
                for u in range(n):
                    if u == w:
                        continue
                    mg = 2 * C[u, w] - H
                    if mg > bestany:
                        bestany = mg
                    if (adj[u] >> w) & 1:
                        dsum += mg
                        if mg > best:
                            best = mg
                od = popc(adj[w]); idg = n - 1 - od
                cls = min(od, idg)
                if best < vl:
                    vl = best
                if cls == 1 and best < c1:
                    c1 = best
                if cls == 2 and best < c2:
                    c2 = best
                if cls >= 3 and best < c3:
                    c3 = best
                if bestany < ms:
                    ms = bestany
                if dsum < dm:
                    dm = dsum
                d = n_distribution(n, adj, F, B, w, Nd)
                for kk in range(5):
                    if 2 * kk + 1 > d:
                        break
                    hi_ = 0; lo_ = 0
                    for j in range(d - kk, d + 1):
                        hi_ += Nd[j]
                    for j in range(0, kk + 1):
                        lo_ += Nd[j]
                    if hi_ < lo_:
                        ssf[t, kk] += 1
                if d <= 6:
                    s1 = 0
                    for j in range(d + 1):
                        s1 += j * Nd[j]
                    val = (2.0 * s1 - d * H) / (2.0 * H)
                    if val < dbyd[t, d]:
                        dbyd[t, d] = val
            out[idx, 4] = vl; out[idx, 5] = ms; out[idx, 6] = dm
            out[idx, 7] = c1; out[idx, 8] = c2; out[idx, 9] = c3
    return out, ssf, dbyd


def run(n, k):
    cand = [lo * n + hi for lo in range(n) for hi in range(lo + 2, n)]
    combos = np.array(list(itertools.combinations(cand, k)), dtype=np.int64)
    t0 = time.time()
    out, ssf, dbyd = sweep(n, combos, 12)
    s = out[:, 0] == 1
    o = out[s]; cs = combos[s]
    arcs = lambda row: [(int(c) // n, int(c) % n) for c in row]
    print(f'n = {n}, k = {k} reversed arcs: {len(combos)} arc sets, {int(s.sum())} strong  [{time.time() - t0:.1f} s]')
    if not s.any():
        return
    H = o[:, 1].astype(float)
    print(f'  M(T) transitive: {int(o[:, 2].sum())};  F(T) not strongly connected: {int((o[:, 3] == 0).sum())}')
    for col, name in ((4, 'vertex-lemma margin'), (5, 'majority-source margin'), (6, 'D_in')):
        rel = o[:, col] / H
        i = int(np.argmin(rel))
        print(f'  min {name} = {int(o[i, col])}/{int(o[i, 1])} = {rel[i]:.6f} at {arcs(cs[i])}')
    print(f'  tournaments with some vertex-lemma margin <= 0: {int((o[:, 4] <= 0).sum())}; with some majority-source margin <= 0: '
          f'{int((o[:, 5] <= 0).sum())}; with some D_in <= 0: {int((o[:, 6] <= 0).sum())}')
    cls = []
    for col, name in ((7, '1'), (8, '2'), (9, '>= 3')):
        mask = o[:, col] < (1 << 61)
        cls.append(f'{name}: {(o[mask, col] / H[mask]).min():.6f}' if mask.any() else f'{name}: -')
    print('  min vertex-lemma margin by min(indeg, outdeg): ' + '; '.join(cls))
    f = ssf.sum(axis=0); dm = dbyd.min(axis=0)
    print('  (SS) failures (T, w) by k: ' + ', '.join(f'k={kk}: {int(f[kk])}' for kk in range(5)) +
          ';  min E[N] - d/2 by d: ' + ', '.join(f'd={d}: {dm[d]:.4f}' for d in range(1, 7) if dm[d] < 8), flush=True)


if __name__ == '__main__':
    args = [int(a) for a in sys.argv[1:]]
    if args:
        jobs = list(zip(args[0::2], args[1::2]))
    else:
        jobs = [(n, 2) for n in range(8, 17)] + [(n, 3) for n in range(8, 15)] + [(n, 4) for n in range(8, 13)] + [(10, 5)]
    for n, k in jobs:
        run(n, k)
