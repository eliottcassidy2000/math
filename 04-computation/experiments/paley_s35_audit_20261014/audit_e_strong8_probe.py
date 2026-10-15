#!/usr/bin/env python3
"""Audit (e3): do strong 8-tournaments reach delta = 4/11?  (S35 reports that annealing over strong tournaments found
nothing below 64/161 at n = 8.)  Collect labelled strong 7-tournaments with delta = 4/11, extend each by a vertex in
all 2^7 ways, and record the least delta among strong extensions.  Also checks whether the Claim-A counterexamples
(n = 8, 9, 10) have a balanced arc.
"""
import time
from fractions import Fraction
import numpy as np
import numba as nb
from audit_e_tournaments import stats_dfs, strong, adj_from_code   # this audit's own module


@nb.njit(cache=True)
def collect(n, pu, pv, want_num, want_den, limit):
    E = pu.shape[0]
    out = np.zeros(limit, dtype=np.int64)
    k = 0
    adj = np.zeros(n, dtype=np.int64)
    C = np.zeros((n, n), dtype=np.int64)
    for code in range(1 << E):
        for i in range(n):
            adj[i] = 0
        for e in range(E):
            if (code >> e) & 1:
                adj[pu[e]] |= 1 << pv[e]
            else:
                adj[pv[e]] |= 1 << pu[e]
        if not strong(n, adj):
            continue
        H = stats_dfs(n, adj, C)
        m = 0
        for x in range(n):
            for y in range(x + 1, n):
                v = C[x, y]
                if H - v < v:
                    v = H - v
                if v > m:
                    m = v
        if m * want_den == want_num * H:
            out[k] = code
            k += 1
            if k == limit:
                break
    return out[:k]


def delta_of(n, adj):
    C = np.zeros((n, n), dtype=np.int64)
    H = int(stats_dfs(n, np.array(adj, dtype=np.int64), C))
    d = max(Fraction(min(int(C[x, y]), H - int(C[x, y])), H) for x in range(n) for y in range(x + 1, n))
    bal_arc = any(3 * int(C[x, y]) >= H and 3 * int(C[x, y]) <= 2 * H for x in range(n) for y in range(n) if x != y and adj[x] >> y & 1)
    return H, d, bal_arc


if __name__ == '__main__':
    t0 = time.time()
    n = 7
    pr = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pu = np.array([p[0] for p in pr], dtype=np.int64)
    pv = np.array([p[1] for p in pr], dtype=np.int64)
    codes = collect(n, pu, pv, 4, 11, 400)
    print(f'collected {len(codes)} labelled strong 7-tournaments with delta = 4/11  [{time.time() - t0:.0f} s]')
    best = Fraction(1)
    ex = None
    cnt = 0
    for code in codes:
        a7 = adj_from_code(7, int(code))
        for pat in range(1 << 7):
            adj = [a | ((1 << 7) if (pat >> i) & 1 == 0 else 0) for i, a in enumerate(a7)]   # i -> 7 when bit i is 0
            a8 = sum(1 << i for i in range(7) if (pat >> i) & 1)                           # 7 -> i when bit i is 1
            adj = adj + [a8]
            if not strong(8, np.array(adj, dtype=np.int64)):
                continue
            cnt += 1
            H, d, _ = delta_of(8, adj)
            if d < best:
                best, ex = d, (adj, H)
    print(f'strong one-vertex extensions examined: {cnt}; least delta = {best} = {float(best):.4f} (64/161 = {64 / 161:.4f}); example adj = {ex}')
    for nn, adj in ((8, [190, 176, 250, 242, 64, 144, 163, 16]), (9, [460, 493, 328, 256, 367, 461, 264, 92, 128]),
                    (10, [734, 220, 208, 116, 128, 983, 144, 8, 735, 222])):
        H, d, bal = delta_of(nn, adj)
        print(f'Claim-A counterexample n = {nn}: H = {H}, delta = {d}, has a balanced arc: {bal}')
    print(f'[{time.time() - t0:.0f} s]')
