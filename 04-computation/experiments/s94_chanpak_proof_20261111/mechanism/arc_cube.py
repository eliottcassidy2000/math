#!/usr/bin/env python3
"""S94 lane "mechanism": additive structure of H on the arc cube of order n (exhaustive, exact).

For n <= 7 the table H[m] of all 2^C(n,2) labelled tournaments is computed by the subset DP (exact integers).
Arc k = pair (i, j), i < j, in lexicographic order; bit k of m set  <=>  i -> j.

(A) adjacent arc pairs (sharing a vertex): mixed second difference
        D = H[m] - H[m ^ e] - H[m ^ f] + H[m ^ e ^ f]
    Theorem A (PROOFS.txt) says |D| = W_xy(T - v) >= 2.  Checked for every m.
(B) disjoint arc pairs: how often is D = 0 (additive squares), with both increments nonzero, all four strong.
(C) additive 3-cubes on matchings of three arcs (n >= 6): all three pairs additive on both layers.
    Reports the multiset of increment triples, and whether some cube has 8 distinct values / 8 consecutive odd values.

usage: python arc_cube.py n
"""
import sys, time
import numpy as np
from numba import njit, prange

@njit(cache=True)
def ham_count(n, m, pi, pj):
    adj = np.zeros(n, np.int64)
    na = n * (n - 1) // 2
    for k in range(na):
        if (m >> k) & 1:
            adj[pi[k]] |= 1 << pj[k]
        else:
            adj[pj[k]] |= 1 << pi[k]
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            c = dp[mask, last]
            if c == 0:
                continue
            nxt = adj[last] & ~mask
            for w in range(n):
                if (nxt >> w) & 1:
                    dp[mask | (1 << w), w] += c
    s = 0
    for v in range(n):
        s += dp[full, v]
    # strong test by scores
    sc = np.zeros(n, np.int64)
    for v in range(n):
        x = adj[v]
        cnt = 0
        while x:
            cnt += 1
            x &= x - 1
        sc[v] = cnt
    sc.sort()
    strong = 1
    acc = 0
    for k in range(1, n):
        acc += sc[k - 1]
        if acc == k * (k - 1) // 2:
            strong = 0
            break
    return s, strong

@njit(parallel=True, cache=True)
def table(n, pi, pj):
    na = n * (n - 1) // 2
    N = 1 << na
    H = np.zeros(N, np.int32)
    S = np.zeros(N, np.int8)
    for m in prange(N):
        h, st = ham_count(n, m, pi, pj)
        H[m] = h
        S[m] = st
    return H, S

@njit(cache=True)
def pair_stats(H, S, e, f):
    """over all m with bits e, f cleared: second difference stats."""
    N = H.shape[0]
    be = 1 << e
    bf = 1 << f
    mn = 1 << 30
    zeros = 0
    zeros_nontriv = 0
    zeros_nontriv_strong = 0
    tot = 0
    for m in range(N):
        if m & be or m & bf:
            continue
        tot += 1
        d = H[m] - H[m | be] - H[m | bf] + H[m | be | bf]
        a = d if d >= 0 else -d
        if a < mn:
            mn = a
        if d == 0:
            zeros += 1
            de = H[m | be] - H[m]
            df = H[m | bf] - H[m]
            if de != 0 and df != 0:
                zeros_nontriv += 1
                if S[m] and S[m | be] and S[m | bf] and S[m | be | bf]:
                    zeros_nontriv_strong += 1
    return tot, mn, zeros, zeros_nontriv, zeros_nontriv_strong

@njit(cache=True)
def cube3(H, S, e, f, g, out_inc, cap):
    """additive 3-cubes on arcs e, f, g: returns count, count nontrivial (all increments nonzero), count nontrivial
    with 8 strong, count with 8 distinct values, count with 8 consecutive odd values; stores increments."""
    N = H.shape[0]
    be = 1 << e; bf = 1 << f; bg = 1 << g
    cnt = 0; nontriv = 0; nts = 0; dist = 0; consec = 0
    k = 0
    vals = np.zeros(8, np.int64)
    for m in range(N):
        if m & (be | bf | bg):
            continue
        h0 = H[m]
        de = H[m | be] - h0
        df = H[m | bf] - h0
        dg = H[m | bg] - h0
        if H[m | be | bf] != h0 + de + df: continue
        if H[m | be | bg] != h0 + de + dg: continue
        if H[m | bf | bg] != h0 + df + dg: continue
        if H[m | be | bf | bg] != h0 + de + df + dg: continue
        cnt += 1
        if de == 0 or df == 0 or dg == 0:
            continue
        nontriv += 1
        allstrong = 1
        idx = 0
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    mm = m | (be * a) | (bf * b) | (bg * c)
                    if S[mm] == 0:
                        allstrong = 0
                    vals[idx] = H[mm]
                    idx += 1
        if allstrong:
            nts += 1
        vs = np.sort(vals)
        d8 = 1
        for i in range(7):
            if vs[i] == vs[i + 1]:
                d8 = 0
        if d8:
            dist += 1
            if vs[7] - vs[0] == 14:
                consec += 1
        if k < cap:
            out_inc[k, 0] = h0; out_inc[k, 1] = de; out_inc[k, 2] = df; out_inc[k, 3] = dg
            out_inc[k, 4] = allstrong; out_inc[k, 5] = m
            k += 1
    return cnt, nontriv, nts, dist, consec, k

def main():
    n = int(sys.argv[1])
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pi = np.array([p[0] for p in pairs], np.int64)
    pj = np.array([p[1] for p in pairs], np.int64)
    t0 = time.time()
    H, S = table(n, pi, pj)
    print(f"n = {n}: table of {H.shape[0]} labelled tournaments in {time.time()-t0:.1f}s; max H = {H.max()}, "
          f"strong: {int(S.sum())}; all H odd: {bool(np.all(H % 2 == 1))}")
    vals = np.unique(H)
    print(f"  distinct values: {len(vals)}; odd numbers <= max missing: "
          f"{[v for v in range(1, int(H.max()) + 1, 2) if v not in set(vals.tolist())][:40]}")
    na = len(pairs)
    # (A), (B)
    adj_min = 1 << 30; adj_zero = 0; adj_tot = 0
    dis_tot = 0; dis_zero = 0; dis_nt = 0; dis_nts = 0; dis_min = 1 << 30
    for e in range(na):
        for f in range(e + 1, na):
            shared = len(set(pairs[e]) & set(pairs[f]))
            tot, mn, z, znt, znts = pair_stats(H, S, e, f)
            if shared:
                adj_tot += tot; adj_zero += z; adj_min = min(adj_min, mn)
            else:
                dis_tot += tot; dis_zero += z; dis_nt += znt; dis_nts += znts; dis_min = min(dis_min, mn)
    print(f"(A) adjacent arc pairs: {adj_tot} squares, zero second differences: {adj_zero}, min |D| = {adj_min}")
    print(f"(B) disjoint arc pairs: {dis_tot} squares, D = 0: {dis_zero} ({dis_zero/max(dis_tot,1):.4f}), "
          f"with both increments nonzero: {dis_nt}, of these with 4 strong tournaments: {dis_nts}")
    # (C)
    if n >= 6:
        from itertools import combinations
        from collections import Counter
        tot = Counter()
        incs = Counter()
        cap = 200000
        buf = np.zeros((cap, 6), np.int64)
        examples = []
        for e, f, g in combinations(range(na), 3):
            vs = set(pairs[e]) | set(pairs[f]) | set(pairs[g])
            if len(vs) != 6:
                continue
            c, nt, nts, dist, consec, k = cube3(H, S, e, f, g, buf, cap)
            tot['cubes'] += c; tot['nontrivial'] += nt; tot['nontrivial_strong'] += nts
            tot['distinct8'] += dist; tot['consecutive8'] += consec
            for r in range(k):
                key = tuple(sorted(abs(int(x)) for x in buf[r, 1:4]))
                incs[(key, int(buf[r, 4]))] += 1
                if len(examples) < 12 and buf[r, 4] == 1:
                    examples.append((int(buf[r, 5]), (e, f, g), [int(x) for x in buf[r, :4]]))
        print(f"(C) additive 3-cubes on 3-matchings: {dict(tot)}")
        top = sorted(incs.items(), key=lambda kv: -kv[1])[:40]
        print("    |increment| triples (sorted), all-strong flag : count")
        for (key, st), c in top:
            print(f"      {key} strong={st} : {c}")
        print("    examples (mask, arcs, [h0, de, df, dg]):")
        for ex in examples:
            print("     ", ex, [pairs[a] for a in ex[1]])

if __name__ == '__main__':
    main()
