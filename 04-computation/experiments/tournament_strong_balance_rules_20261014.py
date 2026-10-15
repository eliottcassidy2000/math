#!/usr/bin/env python3
"""Candidate certificates for the tournament 1/3-2/3 conjecture (HYP-9277) on STRONG tournaments, exhaustive n <= NMAX.

p(x, y) = Pr[x before y] in a uniform random Hamiltonian path.  A pair is balanced if 1/3 <= p <= 2/3.
Rules tested on every labelled strong tournament:
  R1  some ARC is balanced;
  R2  some directed 3-cycle has all three arcs with p(tail before head) >= 1/3 (then one of them is <= 2/3, since
      p(x,y) + p(y,z) + p(z,x) <= 2 on a cyclic triple: hence balanced);
  R3  every vertex lies in a balanced pair;
  R4  some vertex x has f(x) + l(x) >= 2/3 where f, l = Pr[x first], Pr[x last]  (diagnostic only);
  R5  the most balanced pair is an arc.
Also records the minimum of delta over strong tournaments, and the minimum over strong tournaments of the best ARC
balance and of the best 3-cycle certificate value (max over 3-cycles of min over its arcs of p(tail first)).
Session opus-2026-10-14-S36.
"""
import sys, time
import numpy as np
import numba as nb
from fractions import Fraction

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
def analyse(n, adj):
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
    # best pair and best arc (as numerator over H of min(c, H - c))
    best_pair = 0; best_arc = 0; best_pair_is_arc = False
    for x in range(n):
        for y in range(x + 1, n):
            m = C[x, y] if C[x, y] < H - C[x, y] else H - C[x, y]
            if m > best_pair:
                best_pair = m
    for x in range(n):
        for y in range(x + 1, n):
            m = C[x, y] if C[x, y] < H - C[x, y] else H - C[x, y]
            if m == best_pair:
                best_pair_is_arc = True  # every pair is an arc in a tournament
    # arcs: every pair is an arc; R1 = R(best pair) trivially in tournaments!  (kept for clarity)
    # R2: 3-cycle certificate value = max over directed 3-cycles of min over its arcs of C[tail, head]
    cyc_best = -1
    for x in range(n):
        for y in range(n):
            if y == x or not (adj[x] >> y) & 1:
                continue
            for z in range(n):
                if z == x or z == y:
                    continue
                if (adj[y] >> z) & 1 and (adj[z] >> x) & 1:
                    m = C[x, y]
                    if C[y, z] < m:
                        m = C[y, z]
                    if C[z, x] < m:
                        m = C[z, x]
                    if m > cyc_best:
                        cyc_best = m
    # R3: every vertex in a balanced pair
    r3 = True
    for x in range(n):
        ok = False
        for y in range(n):
            if y != x:
                m = C[x, y] if C[x, y] < H - C[x, y] else H - C[x, y]
                if 3 * m >= H:
                    ok = True
                    break
        if not ok:
            r3 = False
            break
    # f, l
    fl_best = 0
    for v in range(n):
        f = 0
        for S in range(1):
            pass
        fv = B[full, v]   # paths starting at v
        lv = F[full, v]   # paths ending at v
        if fv + lv > fl_best:
            fl_best = fv + lv
    return H, best_pair, cyc_best, r3, fl_best


@nb.njit(parallel=True, cache=True)
def sweep(n, pu, pv, nchunk):
    E = pu.shape[0]
    total = 1 << E
    nstrong = np.zeros(nchunk, dtype=np.int64)
    fail_r2 = np.zeros(nchunk, dtype=np.int64)
    fail_r3 = np.zeros(nchunk, dtype=np.int64)
    pass_r4 = np.zeros(nchunk, dtype=np.int64)
    # min over strong of delta (best_pair/H) and of the 3-cycle certificate (cyc_best/H)
    dn = np.ones(nchunk, dtype=np.int64); dd = np.ones(nchunk, dtype=np.int64)
    cn = np.ones(nchunk, dtype=np.int64); cd = np.ones(nchunk, dtype=np.int64); cc = np.full(nchunk, -1, dtype=np.int64)
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
            H, bp, cb, r3, fl = analyse(n, adj)
            if 3 * cb < H:
                fail_r2[c] += 1
            if not r3:
                fail_r3[c] += 1
            if 3 * fl >= 2 * H:
                pass_r4[c] += 1
            if bp * dd[c] < dn[c] * H:
                dn[c] = bp; dd[c] = H
            if cb * cd[c] < cn[c] * H:
                cn[c] = cb; cd[c] = H; cc[c] = code
    return nstrong, fail_r2, fail_r3, pass_r4, dn, dd, cn, cd, cc


if __name__ == '__main__':
    for n in range(max(3, NMIN), NMAX + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in pairs], dtype=np.int64); pv = np.array([p[1] for p in pairs], dtype=np.int64)
        t0 = time.time()
        ns, f2, f3, p4, dn, dd, cn, cd, cc = sweep(n, pu, pv, 1024)
        dmin = min(Fraction(int(a), int(b)) for a, b in zip(dn, dd))
        k = min(range(len(cn)), key=lambda i: Fraction(int(cn[i]), int(cd[i])))
        print(f'n = {n}: strong labelled {int(ns.sum())}; min delta {dmin}; R2 failures (no 3-cycle with all arcs '
              f'tail-first >= 1/3): {int(f2.sum())}; min 3-cycle certificate {Fraction(int(cn[k]), int(cd[k]))} (code {int(cc[k])}); '
              f'R3 failures (a vertex in no balanced pair): {int(f3.sum())}; R4 passes: {int(p4.sum())}  [{time.time() - t0:.1f} s]', flush=True)
