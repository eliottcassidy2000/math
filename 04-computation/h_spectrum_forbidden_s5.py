#!/usr/bin/env python3
"""
H-spectrum analysis: which H values are achievable at n=3..7,
and which (alpha_1, alpha_2, ...) decompositions are forbidden?

The goal is to identify forbidden H values beyond H=7 (THM-343).

opus-2026-05-28-S5
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import permutations, combinations
from collections import Counter


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


def count_HP(T):
    n = len(T)
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


def is_strongly_connected(T):
    n = len(T)
    if n <= 1: return True
    def reach(start, mat):
        seen = {start}; stack = [start]
        while stack:
            u = stack.pop()
            for v in range(n):
                if mat[u][v] and v not in seen:
                    seen.add(v); stack.append(v)
        return seen
    if len(reach(0, T)) != n: return False
    Trev = [[T[j][i] for j in range(n)] for i in range(n)]
    return len(reach(0, Trev)) == n


def enumerate_odd_cycles(T):
    """Return list of distinct directed odd cycles, each as a frozen tuple
    canonicalized by starting at smallest vertex."""
    n = len(T)
    seen = set()
    for L in range(3, n+1, 2):
        for subset in combinations(range(n), L):
            v0 = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(L):
                    a, b = cycle[i], cycle[(i+1)%L]
                    if not T[a][b]:
                        ok = False; break
                if ok:
                    seen.add(cycle)
    return list(seen)


def independence_poly_omega(T):
    """Returns the independence polynomial coefficients alpha_k of Omega(T).
    Omega has vertices = odd cycles, edges = pairs sharing >= 1 vertex.
    alpha_k = # of size-k independent sets = # of k-tuples of pairwise vertex-disjoint odd cycles.
    """
    cycles = enumerate_odd_cycles(T)
    nc = len(cycles)
    Vsets = [frozenset(c) for c in cycles]
    n = len(T)
    coeffs = []
    # recursive enumeration
    counts = [0] * (n+1)  # alpha_k
    counts[0] = 1
    def recurse(start, blocked, k):
        for i in range(start, nc):
            if not (Vsets[i] & blocked):
                counts[k+1] += 1
                recurse(i+1, blocked | Vsets[i], k+1)
    recurse(0, frozenset(), 0)
    return counts


def main():
    for n in range(3, 7):
        print(f"\n=== n = {n} ===")
        spectrum = Counter()
        alpha_vecs = Counter()
        sc_spectrum = Counter()
        for T in gen_tournaments(n):
            H = count_HP(T)
            spectrum[H] += 1
            ips = independence_poly_omega(T)
            # truncate trailing zeros
            while len(ips) > 1 and ips[-1] == 0:
                ips.pop()
            alpha_vecs[tuple(ips)] += 1
            if is_strongly_connected(T):
                sc_spectrum[H] += 1
        print(f"  H spectrum: {sorted(spectrum.items())}")
        H_values = sorted(spectrum.keys())
        H_min, H_max = H_values[0], H_values[-1]
        missing = [h for h in range(H_min, H_max+1, 2) if h not in spectrum]
        print(f"  Missing odd H in [{H_min},{H_max}]: {missing}")
        print(f"  SC-only H spectrum: {sorted(sc_spectrum.items())}")
        print(f"  Distinct independence polynomials: {len(alpha_vecs)}")
        # show top 5 alpha vectors
        top = alpha_vecs.most_common(8)
        for vec, cnt in top:
            H = sum(c * (2**i) for i, c in enumerate(vec))
            print(f"    {vec} (H={H}): {cnt} tournaments")


if __name__ == "__main__":
    main()
