#!/usr/bin/env python3
"""S93, evaluation of the pre-registered W1 and W2 (PREREGISTRATION.md, pushed in commit 3a675e5b7a before this run).

Two-apex family G(2,0): spine s_1 -> ... -> s_n transitive, apexes u -> v, letter (a_i, b_i): a_i = 1 iff u -> s_i,
b_i = 1 iff v -> s_i.  H = number of Hamiltonian paths, by the 36-state transfer of seed/two_apex_probe.py (checked
there against brute force for n <= 5; re-checked here for n <= 4).

W1: the words of length at most NMAX = 16 attain every odd N in [2^21, 2^23).
W2: over the odd N of [2^21, 2^22) the minimum multiplicity (words of length at most 16) is at least 100.
Usage: python3 eval_w1_w2_two_apex.py [NMAX]
"""
import itertools
import sys

import numpy as np
from numba import njit, prange

BOUND = 1 << 23


@njit(cache=True)
def step(vec, a, b, out):
    for k in range(36):
        out[k] = 0
    for arr in range(2):
        bx = a if arr == 0 else b
        by = b if arr == 0 else a
        for r1 in range(3):
            for r2 in range(3):
                for r3 in range(2):
                    c = vec[arr * 18 + r1 * 6 + r2 * 2 + r3]
                    if c == 0:
                        continue
                    n1 = 1 if bx == 0 else 2
                    out[arr * 18 + n1 * 6 + r2 * 2 + r3] += c
                    if r2 != 0 or bx == 1:
                        n2 = 1 if by == 0 else 2
                        out[arr * 18 + r1 * 6 + n2 * 2 + r3] += c
                    if r3 != 0 or by == 1:
                        out[arr * 18 + r1 * 6 + r2 * 2 + 1] += c


@njit(cache=True)
def final(vec):
    tot = 0
    for arr in range(2):
        for r1 in range(3):
            if r1 == 2:
                continue
            for r2 in range(3):
                if r2 == 2:
                    continue
                if r2 == 0 and arr == 1:
                    continue
                for r3 in range(2):
                    tot += vec[arr * 18 + r1 * 6 + r2 * 2 + r3]
    return tot


@njit(cache=True)
def count_word(n, a, b):
    vec = np.zeros(36, np.int64)
    out = np.zeros(36, np.int64)
    vec[0] = 1
    vec[18] = 1
    for i in range(n):
        step(vec, a[i], b[i], out)
        for k in range(36):
            vec[k] = out[k]
    return final(vec)


@njit(cache=True)
def dfs(start, n, d0, mult, bound):
    """all words of length n extending the state `start` reached after d0 letters."""
    stack = np.zeros((n + 1, 36), np.int64)
    for k in range(36):
        stack[d0, k] = start[k]
    letters = np.zeros(n + 1, np.int64)
    out = np.zeros(36, np.int64)
    big = 0
    if d0 == n:
        h = final(stack[n])
        if h < bound:
            mult[h >> 1] += 1
        else:
            big += 1
        return big
    depth = d0
    letters[depth] = -1
    while depth >= d0:
        letters[depth] += 1
        if letters[depth] == 4:
            depth -= 1
            continue
        l = letters[depth]
        step(stack[depth], l & 1, (l >> 1) & 1, out)
        for k in range(36):
            stack[depth + 1, k] = out[k]
        if depth + 1 == n:
            h = final(stack[n])
            if h < bound:
                mult[h >> 1] += 1
            else:
                big += 1
        else:
            depth += 1
            letters[depth] = -1
    return big


@njit(parallel=True, cache=True)
def census(n, bound, mult_all):
    """words of length n >= 2, split over the 16 two-letter prefixes."""
    bigs = np.zeros(16, np.int64)
    for pre in prange(16):
        vec = np.zeros(36, np.int64)
        out = np.zeros(36, np.int64)
        vec[0] = 1
        vec[18] = 1
        l0, l1 = pre & 3, pre >> 2
        step(vec, l0 & 1, (l0 >> 1) & 1, out)
        for k in range(36):
            vec[k] = out[k]
        step(vec, l1 & 1, (l1 >> 1) & 1, out)
        for k in range(36):
            vec[k] = out[k]
        bigs[pre] = dfs(vec, n, 2, mult_all[pre], bound)
    return bigs.sum()


def brute(n, a, b):
    V = n + 2
    u, v = n, n + 1
    adj = [[0] * V for _ in range(V)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = 1
        if a[i]:
            adj[u][i] = 1
        else:
            adj[i][u] = 1
        if b[i]:
            adj[v][i] = 1
        else:
            adj[i][v] = 1
    adj[u][v] = 1
    dp = [[0] * V for _ in range(1 << V)]
    for x in range(V):
        dp[1 << x][x] = 1
    for S in range(1, 1 << V):
        for x in range(V):
            c = dp[S][x]
            if c:
                for y in range(V):
                    if adj[x][y] and not (S >> y) & 1:
                        dp[S | (1 << y)][y] += c
    return sum(dp[(1 << V) - 1])


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    ok = all(brute(n, a, b) == count_word(n, np.array(a, np.int64), np.array(b, np.int64))
             for n in (1, 2, 3, 4) for a in itertools.product((0, 1), repeat=n) for b in itertools.product((0, 1), repeat=n))
    print("transfer count = brute force for n <= 4:", ok)
    total = np.zeros(BOUND >> 1, np.int64)
    for n in (1,):
        for l in range(4):
            h = count_word(1, np.array([l & 1], np.int64), np.array([(l >> 1) & 1], np.int64))
            total[h >> 1] += 1
    for n in range(2, NMAX + 1):
        mult_all = np.zeros((16, BOUND >> 1), np.int32)
        big = census(n, BOUND, mult_all)
        m = mult_all.sum(axis=0, dtype=np.int64)
        total += m
        print(f"n = {n:2d}: {4 ** n} words, {int(big)} with H >= 2^23, distinct odd values below 2^23: {int((m > 0).sum())}", flush=True)
    miss = np.flatnonzero(total == 0)
    miss_all = [int(2 * i + 1) for i in miss]
    print(f"odd N < 2^23 not attained by words of length <= {NMAX}: {miss_all[:20]}{' ...' if len(miss_all) > 20 else ''} ({len(miss_all)} in all)")
    w1_range = [N for N in miss_all if (1 << 21) <= N < (1 << 23)]
    print(f"W1: odd N in [2^21, 2^23) not attained: {len(w1_range)} -> {'CONFIRMED' if not w1_range else 'REFUTED'}")
    for k in range(3, 23):
        lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
        seg = total[lo:hi]
        print(f"  odd N in [2^{k}, 2^{k+1}): min multiplicity {int(seg.min())}, median {int(np.median(seg))}, mean {seg.mean():.1f}, zeros {int((seg == 0).sum())}")
    seg = total[(1 << 21) >> 1:(1 << 22) >> 1]
    print(f"W2: minimum multiplicity over the odd N of [2^21, 2^22): {int(seg.min())} (registered: at least 100) -> "
          f"{'CONFIRMED' if seg.min() >= 100 else 'REFUTED'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
