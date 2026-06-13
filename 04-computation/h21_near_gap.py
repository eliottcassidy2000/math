#!/usr/bin/env python3
"""
Study the gap near H=21 at n=7 (exhaustive) and n=8 (sampling).

What H values are near 21? What score sequences produce them?
Why does the gap persist?

kind-pasteur-2026-03-07-S31
"""

from collections import defaultdict
import random
import time


def held_karp(n, adj):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            c = dp[S][v]
            if c == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v] & (1 << u):
                    dp[S | (1 << u)][u] += c
    return sum(dp[(1 << n) - 1])


def score_sequence(n, adj):
    scores = []
    for v in range(n):
        scores.append(bin(adj[v]).count('1'))
    return tuple(sorted(scores))


def n7_exhaustive():
    """Exhaustive H analysis at n=7, focusing on H in [15..27]."""
    n = 7
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    H_count = defaultdict(int)
    H_scores = defaultdict(set)

    t0 = time.time()
    for bits in range(1 << m):
        adj = [0]*n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[j] |= (1 << i)
            else:
                adj[i] |= (1 << j)

        H = held_karp(n, adj)
        H_count[H] += 1

        if 13 <= H <= 29:
            ss = score_sequence(n, adj)
            H_scores[H].add(ss)

    elapsed = time.time() - t0
    print(f"n=7 exhaustive ({elapsed:.0f}s):")
    print(f"\nH values near 21:")
    for h in range(13, 30, 2):
        cnt = H_count.get(h, 0)
        scores = H_scores.get(h, set())
        print(f"  H={h}: {cnt:>7} tournaments, {len(scores)} score sequences")
        if scores:
            for ss in sorted(scores)[:3]:
                print(f"         score={ss}")

    print(f"\nGaps below 50:")
    for h in range(1, 50, 2):
        if h not in H_count:
            print(f"  H={h}: GAP")


def n8_sampling():
    """Sample n=8 tournaments, focus on H near 21."""
    n = 8
    random.seed(42)
    num_samples = 500000

    H_count = defaultdict(int)
    H_scores = defaultdict(set)

    t0 = time.time()
    for _ in range(num_samples):
        adj = [0]*n
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i] |= (1 << j)
                else:
                    adj[j] |= (1 << i)

        H = held_karp(n, adj)
        H_count[H] += 1

        if 13 <= H <= 29:
            ss = score_sequence(n, adj)
            H_scores[H].add(ss)

    elapsed = time.time() - t0
    print(f"\nn=8 sampling ({num_samples} samples, {elapsed:.0f}s):")
    print(f"\nH values near 21:")
    for h in range(13, 30, 2):
        cnt = H_count.get(h, 0)
        scores = H_scores.get(h, set())
        print(f"  H={h}: {cnt:>7} tournaments, {len(scores)} score sequences")
        if scores:
            for ss in sorted(scores)[:3]:
                print(f"         score={ss}")

    print(f"\nGaps below 50:")
    for h in range(1, 50, 2):
        if h not in H_count:
            print(f"  H={h}: GAP")


if __name__ == '__main__':
    n7_exhaustive()
    n8_sampling()
