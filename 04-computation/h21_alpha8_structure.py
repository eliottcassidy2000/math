#!/usr/bin/env python3
"""
Deep structural analysis of alpha_1=8 tournaments.

Key question: WHY does alpha_1=8 force i_2 in {0, 7} but never 1?

At n=7: alpha_1=8 always gives i_2=0 (all 8 cycles pairwise conflict).
At n=8: alpha_1=8 gives i_2=0 or 7 (never 1,2,3,4,5,6).

Strategy:
1. Enumerate all alpha_1=8 tournaments at n=7 exhaustively
2. Study their score sequences, cycle compositions, Omega graph structure
3. At n=8: check whether ALL alpha_1=8 have source/sink (reducing to n=7)
4. For those without source/sink at n=8: analyze the jump i_2=0 -> i_2=7

Instance: kind-pasteur-2026-03-07-S33
"""

import os, sys
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from collections import Counter, defaultdict

def find_directed_cycles_dp(adj, n, k):
    result = []
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        full = (1 << k) - 1
        for S in range(1, full + 1):
            for i in range(k):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(k):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, k):
            if (full, j) in dp and adj[v[j]][v[0]]:
                count += dp[(full, j)]
        if count > 0:
            result.append((frozenset(verts), count))
    return result


def get_all_cycles(adj, n):
    all_cycles = []
    max_k = n if n % 2 == 1 else n - 1
    for k in range(3, max_k + 1, 2):
        for vs, d in find_directed_cycles_dp(adj, n, k):
            for _ in range(d):
                all_cycles.append(vs)
    return all_cycles


def compute_i2(all_cycles):
    i2 = 0
    for a in range(len(all_cycles)):
        for b in range(a+1, len(all_cycles)):
            if not (all_cycles[a] & all_cycles[b]):
                i2 += 1
    return i2


def omega_graph(all_cycles):
    """Return adjacency list of Omega(T) (conflict = edge)."""
    n = len(all_cycles)
    adj = [[] for _ in range(n)]
    for a in range(n):
        for b in range(a+1, n):
            if all_cycles[a] & all_cycles[b]:
                adj[a].append(b)
                adj[b].append(a)
    return adj


def main():
    # n=7 exhaustive: enumerate ALL tournaments with alpha_1=8
    n = 7
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    num_edges = len(edges)  # 21

    print(f"=== n={n} EXHAUSTIVE: alpha_1=8 analysis ===")
    print(f"Total tournaments: {2**num_edges}")

    results = []

    for bits in range(2**num_edges):
        adj = [[0]*n for _ in range(n)]
        for k_idx, (i, j) in enumerate(edges):
            if (bits >> k_idx) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1

        all_cyc = get_all_cycles(adj, n)
        alpha1 = len(all_cyc)

        if alpha1 != 8:
            continue

        i2 = compute_i2(all_cyc)
        scores = tuple(sorted([sum(adj[i]) for i in range(n)]))

        # Cycle composition
        comp = {}
        for c in all_cyc:
            k = len(c)
            comp[k] = comp.get(k, 0) + 1
        comp_key = tuple(sorted(comp.items()))

        # Omega structure: degree sequence
        omega_adj = omega_graph(all_cyc)
        omega_degs = tuple(sorted([len(omega_adj[i]) for i in range(8)]))

        # Check for source/sink
        has_source_sink = 0 in [sum(adj[i]) for i in range(n)] or (n-1) in [sum(adj[i]) for i in range(n)]

        results.append({
            'i2': i2, 'scores': scores, 'comp': comp_key,
            'omega_degs': omega_degs, 'source_sink': has_source_sink,
            'bits': bits
        })

    print(f"Found {len(results)} tournaments with alpha_1=8")

    # i_2 distribution
    i2_dist = Counter(r['i2'] for r in results)
    print(f"i_2 distribution: {dict(sorted(i2_dist.items()))}")

    # Score distribution
    score_dist = Counter(r['scores'] for r in results)
    print(f"Score sequences ({len(score_dist)} distinct):")
    for s, c in sorted(score_dist.items(), key=lambda x: -x[1])[:10]:
        print(f"  {s}: {c}")

    # Composition distribution
    comp_dist = Counter(r['comp'] for r in results)
    print(f"Cycle compositions ({len(comp_dist)} distinct):")
    for comp, c in sorted(comp_dist.items(), key=lambda x: -x[1]):
        print(f"  {dict(comp)}: {c}")

    # Omega degree sequence distribution
    omega_dist = Counter(r['omega_degs'] for r in results)
    print(f"Omega degree sequences ({len(omega_dist)} distinct):")
    for od, c in sorted(omega_dist.items(), key=lambda x: -x[1])[:10]:
        print(f"  {od}: {c}")

    # Source/sink
    ss_dist = Counter(r['source_sink'] for r in results)
    print(f"Has source/sink: {dict(ss_dist)}")

    # Cross-tabulate: for each composition, what i_2 values occur?
    print(f"\n=== Cross-tabulation: composition -> i_2 ===")
    by_comp = defaultdict(list)
    for r in results:
        by_comp[r['comp']].append(r['i2'])
    for comp, i2s in sorted(by_comp.items()):
        i2_c = Counter(i2s)
        print(f"  {dict(comp)}: i_2 = {dict(sorted(i2_c.items()))}")

    # For i_2=0 cases: what's the Omega structure?
    print(f"\n=== i_2=0 case: Omega is complete graph K_8 ===")
    i2_0 = [r for r in results if r['i2'] == 0]
    if i2_0:
        omega_dist_0 = Counter(r['omega_degs'] for r in i2_0)
        print(f"  Omega degs: {dict(omega_dist_0)}")
        score_dist_0 = Counter(r['scores'] for r in i2_0)
        print(f"  Scores: {dict(sorted(score_dist_0.items(), key=lambda x: -x[1]))}")

    # For non-zero i_2 cases
    print(f"\n=== i_2 > 0 cases ===")
    i2_pos = [r for r in results if r['i2'] > 0]
    if i2_pos:
        omega_dist_p = Counter(r['omega_degs'] for r in i2_pos)
        print(f"  Omega degs: {dict(omega_dist_p)}")
        score_dist_p = Counter(r['scores'] for r in i2_pos)
        print(f"  Scores: {dict(sorted(score_dist_p.items(), key=lambda x: -x[1]))}")
        # Example
        ex = i2_pos[0]
        print(f"  Example: bits={ex['bits']}, i_2={ex['i2']}, scores={ex['scores']}")

    # KEY CHECK: vertex coverage of the 8 cycles
    print(f"\n=== Vertex coverage analysis ===")
    coverage_dist = Counter()
    for r in results[:200]:  # sample for speed
        bits = r['bits']
        adj = [[0]*n for _ in range(n)]
        for k_idx, (i, j) in enumerate(edges):
            if (bits >> k_idx) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1
        all_cyc = get_all_cycles(adj, n)
        verts = set()
        for c in all_cyc:
            verts |= c
        coverage_dist[len(verts)] += 1
    print(f"  Coverage (first 200): {dict(sorted(coverage_dist.items()))}")


if __name__ == "__main__":
    main()
