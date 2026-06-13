#!/usr/bin/env python3
"""
H=21 forbidden? Test directly by enumerating which Omega structures arise
in tournaments and which independence vectors can equal (4,3,0,...).

If we can show no tournament achieves alpha=(4,3,0,...) or alpha=(6,2,0,...) etc.,
we get partial progress on the H=21 obstruction.

opus-2026-05-28-S5
"""
import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import random
from itertools import combinations, permutations
from collections import Counter

random.seed(0xc0ffee)


def is_strongly_connected(T):
    n = len(T)
    if n <= 1: return True
    def reach(s, m):
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for v in range(n):
                if m[u][v] and v not in seen:
                    seen.add(v); st.append(v)
        return seen
    if len(reach(0, T)) != n: return False
    Trev = [[T[j][i] for j in range(n)] for i in range(n)]
    return len(reach(0, Trev)) == n


def enumerate_odd_cycles(T):
    n = len(T)
    seen = set()
    for L in range(3, n+1, 2):
        for subset in combinations(range(n), L):
            v0 = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(L):
                    a,b = cycle[i], cycle[(i+1)%L]
                    if not T[a][b]:
                        ok = False; break
                if ok:
                    seen.add(cycle)
    return list(seen)


def independence_poly(T):
    cycles = enumerate_odd_cycles(T)
    nc = len(cycles)
    Vsets = [frozenset(c) for c in cycles]
    n = len(T)
    counts = [0] * (n+1)
    counts[0] = 1
    def recurse(start, blocked, k):
        for i in range(start, nc):
            if not (Vsets[i] & blocked):
                counts[k+1] += 1
                recurse(i+1, blocked | Vsets[i], k+1)
    recurse(0, frozenset(), 0)
    # truncate trailing zeros
    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()
    return tuple(counts), cycles


def count_HP(T, n):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for S in range(1<<n):
        for v in range(n):
            if dp[S][v] == 0: continue
            if not ((S >> v) & 1): continue
            for u in range(n):
                if (S >> u) & 1: continue
                if T[v][u]:
                    dp[S | (1<<u)][u] += dp[S][v]
    full = (1<<n) - 1
    return sum(dp[full][v] for v in range(n))


def gen_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(1 << m):
        T = [[False]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = True
            else:
                T[j][i] = True
        yield T


def search_h21_n7(samples):
    """Sample n=7 tournaments and look for H=21."""
    n = 7
    found = []
    h21_count = 0
    for trial in range(samples):
        T = [[False]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                if random.random() < 0.5:
                    T[i][j] = True
                else:
                    T[j][i] = True
        H = count_HP(T, n)
        if H == 21:
            h21_count += 1
            ips, _ = independence_poly(T)
            found.append((ips, T))
    return h21_count, found


def main():
    # Enumerate all n=6 tournaments and look at the alpha vectors
    n = 6
    print(f"Enumerating n={n} tournaments — distinct alpha vectors:")
    av = Counter()
    av_to_H = {}
    for T in gen_tournaments(n):
        ips, _ = independence_poly(T)
        av[ips] += 1
        H = sum(c * (2**i) for i, c in enumerate(ips))
        av_to_H[ips] = H
    # sort by H
    by_H = sorted(av.items(), key=lambda x: av_to_H[x[0]])
    print(f"\n  H | alpha vector | count")
    for ips, cnt in by_H:
        H = av_to_H[ips]
        print(f"  H={H:3d}: {ips} -> {cnt}")
    # check what alpha vectors give H=21
    print(f"\nH=21 decompositions among n=6: any present?")
    h21_present = [ips for ips in av if av_to_H[ips] == 21]
    print(f"  H=21 alpha vectors at n=6: {h21_present}")

    # Now sample n=7 and look for H=21
    print(f"\nSampling n=7 for H=21:")
    samples = 1000000
    found_count, found = search_h21_n7(samples)
    print(f"  H=21 hits in {samples} samples at n=7: {found_count}")
    if found:
        from collections import Counter as C
        ips_distrib = C(f[0] for f in found)
        for ips, cnt in ips_distrib.most_common():
            print(f"    alpha = {ips}: {cnt}")


if __name__ == "__main__":
    main()
