#!/usr/bin/env python3
"""Labelled sweep of the tournament balance constant delta(T) for n <= 8 (numba, parallel).

delta(T) = max over pairs of min(p, 1 - p), p = Pr[x before y] in a uniform random Hamiltonian path of T.
Outputs per n:
  - the number of labelled tournaments with delta = 1/3, against the prediction for "transitive compositions of
    singletons and cyclic triangles": sum over compositions of n into parts 1 and 3 with b >= 1 parts equal to 3 of
    n!/3^b (each 3-block carries 2 labelled cyclic triangles, n!/6^b ordered set partitions);
  - the least delta over strong tournaments and the least delta > 1/3 overall (exact fractions);
  - whether any non-transitive tournament has delta < 1/3.
The prediction is exactly the claim that strong tournaments on >= 4 vertices have delta > 1/3 (delta of a tournament
is the maximum over its strong components; pairs in different components are ordered in every Hamiltonian path).
Session opus-2026-10-14-S35.
"""
import sys, time, math
import numpy as np
import numba as nb
from fractions import Fraction

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
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
    seen = 1; frontier = 1
    while frontier:
        nxt = 0
        for v in range(n):
            if (frontier >> v) & 1:
                nxt |= adj[v]
        nxt &= ~seen
        seen |= nxt; frontier = nxt
    if seen != full:
        return False
    # reverse reachability
    radj = np.zeros(n, dtype=np.int64)
    for u in range(n):
        for v in range(n):
            if (adj[u] >> v) & 1:
                radj[v] |= 1 << u
    seen = 1; frontier = 1
    while frontier:
        nxt = 0
        for v in range(n):
            if (frontier >> v) & 1:
                nxt |= radj[v]
        nxt &= ~seen
        seen |= nxt; frontier = nxt
    return seen == full


@nb.njit(cache=True)
def balance(n, adj):
    F, B = fb_tables(n, adj)
    full = (1 << n) - 1
    H = 0
    for v in range(n):
        H += F[full, v]
    best = 0
    for x in range(n):
        for y in range(x + 1, n):
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
            m = cnt if cnt < H - cnt else H - cnt
            if m > best:
                best = m
    return H, best


@nb.njit(parallel=True, cache=True)
def sweep(n, pu, pv, nchunk):
    E = pu.shape[0]
    total = 1 << E
    at13 = np.zeros(nchunk, dtype=np.int64)
    below = np.zeros(nchunk, dtype=np.int64)
    s_num = np.ones(nchunk, dtype=np.int64); s_den = np.ones(nchunk, dtype=np.int64)   # min over strong
    a_num = np.ones(nchunk, dtype=np.int64); a_den = np.ones(nchunk, dtype=np.int64)   # least > 1/3
    s_code = np.full(nchunk, -1, dtype=np.int64); a_code = np.full(nchunk, -1, dtype=np.int64)
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
            H, m = balance(n, adj)
            if H == 1:
                continue
            if 3 * m == H:
                at13[c] += 1
            elif 3 * m < H:
                below[c] += 1
            else:
                if m * a_den[c] < a_num[c] * H:
                    a_num[c] = m; a_den[c] = H; a_code[c] = code
            if is_strong(n, adj):
                if m * s_den[c] < s_num[c] * H:
                    s_num[c] = m; s_den[c] = H; s_code[c] = code
    return at13, below, s_num, s_den, s_code, a_num, a_den, a_code


def predicted_labelled(n):
    tot = 0
    def rec(rem, b):
        nonlocal tot
        if rem == 0:
            if b >= 1:
                tot += math.factorial(n) // 3 ** b
            return
        rec(rem - 1, b)
        if rem >= 3:
            rec(rem - 3, b + 1)
    rec(n, 0)
    return tot


def main():
    for n in range(NMIN, NMAX + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in pairs], dtype=np.int64)
        pv = np.array([p[1] for p in pairs], dtype=np.int64)
        t0 = time.time()
        at13, below, sn, sd, sc, an, ad, ac = sweep(n, pu, pv, 1024)
        fs = [Fraction(int(sn[k]), int(sd[k])) for k in range(len(sn)) if sc[k] >= 0]
        fa = [Fraction(int(an[k]), int(ad[k])) for k in range(len(an)) if ac[k] >= 0]
        ks = min(range(len(sn)), key=lambda k: Fraction(int(sn[k]), int(sd[k])) if sc[k] >= 0 else Fraction(2))
        ka = min(range(len(an)), key=lambda k: Fraction(int(an[k]), int(ad[k])) if ac[k] >= 0 else Fraction(2))
        print(f'n = {n}: labelled with delta = 1/3: {int(at13.sum())} (prediction {predicted_labelled(n)}); '
              f'delta < 1/3: {int(below.sum())}; least delta over strong: {min(fs) if fs else None} '
              f'(code {int(sc[ks])}); least delta > 1/3: {min(fa) if fa else None} (code {int(ac[ka])})  '
              f'[{time.time() - t0:.1f} s]', flush=True)


if __name__ == '__main__':
    main()
