#!/usr/bin/env python3
"""paley11_subtournament_uniformity.py -- Does Paley T_11's advantage come from
sub-tournament uniformity, like BIBD at n=7?

Session: kind-pasteur-2026-03-20-S1

At n=7, the BIBD (Paley T_7) has ALL 21 five-vertex sub-tournaments identical.
Does the same hold for Paley T_11 on its sub-tournaments?

We check: for each size k=5,7,9, do the C(11,k) sub-tournaments of Paley T_11
all have the same H? And how does this compare with Interval T_11?
"""

from itertools import combinations
from collections import Counter
from math import comb

def paley_adj(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in qr:
            j = (i + s) % p
            adj[i][j] = 1
    return adj

def interval_adj(p):
    S = set(range(1, (p+1)//2))
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S:
            j = (i + s) % p
            adj[i][j] = 1
    return adj

def held_karp_H(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def sub_H(adj_global, verts):
    """H of sub-tournament on given vertices."""
    k = len(verts)
    sub = [[adj_global[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
    return held_karp_H(sub, k)

def sub_score(adj_global, verts):
    """Score sequence of sub-tournament."""
    k = len(verts)
    sub = [[adj_global[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
    return score_seq(sub, k)


def analyze(adj, p, name):
    print(f"\n{'='*60}")
    print(f"  {name} (p={p})")
    print(f"{'='*60}")

    for k in [5, 7]:
        n_subsets = comb(p, k)
        H_dist = Counter()
        score_dist = Counter()

        for verts in combinations(range(p), k):
            H = sub_H(adj, verts)
            H_dist[H] += 1
            sc = sub_score(adj, verts)
            score_dist[sc] += 1

        avg_H = sum(h*c for h,c in H_dist.items()) / n_subsets
        num_distinct_H = len(H_dist)
        num_distinct_scores = len(score_dist)

        print(f"\n  {k}-vertex sub-tournaments ({n_subsets} total):")
        print(f"    Distinct H values: {num_distinct_H}")
        print(f"    Distinct score sequences: {num_distinct_scores}")
        print(f"    Average H: {avg_H:.2f}")

        if num_distinct_H <= 10:
            print(f"    H distribution: {dict(sorted(H_dist.items()))}")
        else:
            print(f"    H range: {min(H_dist.keys())} to {max(H_dist.keys())}")
            top5 = H_dist.most_common(5)
            print(f"    Top 5 H values: {top5}")

        # Variance as uniformity measure
        var_H = sum((h - avg_H)**2 * c for h,c in H_dist.items()) / n_subsets
        print(f"    Variance(H): {var_H:.2f}")
        print(f"    StdDev(H): {var_H**0.5:.2f}")
        print(f"    CV(H) = StdDev/Mean: {(var_H**0.5/avg_H)*100:.1f}%")

    # For 9-vertex sub-tournaments, only do a sample (C(11,9)=55)
    k = 9
    n_subsets = comb(p, k)
    if n_subsets <= 100:
        H_dist = Counter()
        for verts in combinations(range(p), k):
            H = sub_H(adj, verts)
            H_dist[H] += 1

        avg_H = sum(h*c for h,c in H_dist.items()) / n_subsets
        var_H = sum((h - avg_H)**2 * c for h,c in H_dist.items()) / n_subsets

        print(f"\n  {k}-vertex sub-tournaments ({n_subsets} total):")
        print(f"    Distinct H values: {len(H_dist)}")
        print(f"    Average H: {avg_H:.2f}")
        print(f"    Variance(H): {var_H:.2f}, CV: {(var_H**0.5/avg_H)*100:.1f}%")
        if len(H_dist) <= 10:
            print(f"    H distribution: {dict(sorted(H_dist.items()))}")


def main():
    p = 11
    adj_p = paley_adj(p)
    adj_i = interval_adj(p)

    analyze(adj_p, p, "Paley T_11")
    analyze(adj_i, p, "Interval T_11")

    # Summary
    print(f"\n{'='*60}")
    print(f"  UNIFORMITY COMPARISON SUMMARY")
    print(f"{'='*60}")
    print(f"  If Paley's sub-tournaments are MORE uniform (lower CV)")
    print(f"  than Interval's, this confirms the BIBD uniformity mechanism.")


if __name__ == "__main__":
    main()
