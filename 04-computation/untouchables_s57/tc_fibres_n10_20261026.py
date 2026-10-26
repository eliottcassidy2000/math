#!/usr/bin/env python3
"""HYP-9352(2) at n = 10: no tournament on 10 vertices has tc = H/|Aut| in {7, 21} (opus-2026-10-26-S57).

Reduction (proved in the note):
  * tc is multiplicative over strong components (H and Aut both factor along the condensation order); components of a
    non-strong 10-tournament have <= 9 vertices, where THM-4741's census found no tc in {7, 21}; a product 21 = 3*7 would
    need a component with tc = 7.  So only STRONG 10-tournaments S matter.
  * rigid S: tc = H, and H is never 7 or 21 (THM-1370).
  * Busch's minimum f(10) = 125 (THM-4744): tc in {7,21} forces H = 7|Aut| or 21|Aut| >= 125, so |Aut| = 3 or 5
    (H in {21, 63, 35, 105}) is impossible.  Every remaining odd |Aut| >= 7 has an element of order 7, of order 5, or of
    order 3 with at least two 3-cycles (a group of order 9 on 10 points cannot act with a single 3-cycle per element).
Hence it suffices to enumerate all tournaments invariant under one permutation of each cycle type:
  (3)(3) + 4 fixed, (3)(3)(3) + 1 fixed, (5) + 5 fixed, (5)(5), (7) + 3 fixed, (9) + 1 fixed,
compute H by the Held-Karp DP, keep the strong ones with H in {7a, 21a} for an odd a >= 7, and compute |Aut| exactly.
Usage: python tc_fibres_n10_20261026.py
"""
import itertools, json, os, sys, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
N = 10


def perm_from_cycles(cycles, n=N):
    p = list(range(n))
    for c in cycles:
        for i in range(len(c)):
            p[c[i]] = c[(i + 1) % len(c)]
    return p


def pair_orbits(p, n=N):
    seen = set()
    orbs = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in seen:
                continue
            orb = []
            a, b = i, j
            while True:
                key = (min(a, b), max(a, b))
                if key in seen:
                    break
                seen.add(key)
                orb.append((a, b))      # orientation representative: a -> b
                a, b = p[a], p[b]
            orbs.append(orb)
    return orbs


@nb.njit(cache=False)
def ham_paths(adj, n):
    # dp[mask][v] = number of Hamiltonian paths of the subtournament on mask ending at v
    full = 1 << n
    dp = np.zeros((full, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, full):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            for w in range(n):
                if (mask >> w) & 1:
                    continue
                if adj[v, w]:
                    dp[mask | (1 << w), w] += c
    s = 0
    for v in range(n):
        s += dp[full - 1, v]
    return s


@nb.njit(cache=False)
def strong(adj, n):
    # strongly connected iff every vertex reaches every vertex: reach from 0 and to 0
    for direction in range(2):
        seen = np.zeros(n, np.uint8)
        stack = np.zeros(n, np.int64)
        top = 0
        stack[top] = 0
        top += 1
        seen[0] = 1
        while top > 0:
            top -= 1
            v = stack[top]
            for w in range(n):
                e = adj[v, w] if direction == 0 else adj[w, v]
                if e and not seen[w]:
                    seen[w] = 1
                    stack[top] = w
                    top += 1
        for w in range(n):
            if not seen[w]:
                return False
    return True


def aut_order(adj, n=N):
    """|Aut| by backtracking with score-sequence pruning."""
    score = adj.sum(axis=1)
    order = sorted(range(n), key=lambda v: score[v])
    cnt = 0

    def bt(k, img, used):
        nonlocal cnt
        if k == n:
            cnt += 1
            return
        v = order[k]
        for w in range(n):
            if used[w] or score[w] != score[v]:
                continue
            ok = True
            for t in range(k):
                u = order[t]
                if adj[u, v] != adj[img[u], w]:
                    ok = False
                    break
            if ok:
                img[v] = w
                used[w] = True
                bt(k + 1, img, used)
                used[w] = False
    bt(0, [0] * n, [False] * n)
    return cnt


if __name__ == '__main__':
    types = {
        '(3)(3)+4': [(0, 1, 2), (3, 4, 5)],
        '(3)(3)(3)+1': [(0, 1, 2), (3, 4, 5), (6, 7, 8)],
        '(5)+5': [(0, 1, 2, 3, 4)],
        '(5)(5)': [(0, 1, 2, 3, 4), (5, 6, 7, 8, 9)],
        '(7)+3': [(0, 1, 2, 3, 4, 5, 6)],
        '(9)+1': [(0, 1, 2, 3, 4, 5, 6, 7, 8)],
    }
    out = {}
    t0 = time.time()
    total_found = []
    for name, cyc in types.items():
        p = perm_from_cycles(cyc)
        orbs = pair_orbits(p)
        m = len(orbs)
        nstrong = 0
        cand = 0
        hits = []
        hvals = {}
        for bits in range(1 << m):
            adj = np.zeros((N, N), np.uint8)
            for t, orb in enumerate(orbs):
                o = (bits >> t) & 1
                for (a, b) in orb:
                    if o:
                        adj[a, b] = 1
                    else:
                        adj[b, a] = 1
            if not strong(adj, N):
                continue
            nstrong += 1
            H = int(ham_paths(adj, N))
            for tc in (7, 21):
                if H % tc == 0:
                    a = H // tc
                    if a >= 7 and a % 2 == 1:
                        cand += 1
                        A = aut_order(adj)
                        if A == a:
                            hits.append((tc, H, A, bits))
        out[name] = {'pair_orbits': m, 'tournaments': 1 << m, 'strong': nstrong, 'candidates_checked': cand,
                     'tc_7_or_21_found': hits}
        total_found += hits
        print(f'{name}: 2^{m} invariant tournaments, {nstrong} strong, {cand} H-candidates needing |Aut|, '
              f'hits {hits} ({time.time() - t0:.0f}s)', flush=True)
    print('TOTAL tc in {7,21} at n = 10:', total_found)
    json.dump(out, open(os.path.join(HERE, 'tc_fibres_n10.json'), 'w'), indent=1)
