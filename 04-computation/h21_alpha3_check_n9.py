#!/usr/bin/env python3
"""
At n=9, cycle-rich, t3 <= 10: does alpha_3 >= 1 always hold?

If YES: Part C immediately blocks H=21 at n=9.
If NO: need alternative argument.

Also check: what is the minimum of alpha_1 + 2*alpha_2 + 4*alpha_3?
For H=21: need this to equal 10.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from collections import Counter
from math import comb
import random

def count_3cycles(adj, n):
    scores = [sum(adj[i]) for i in range(n)]
    return comb(n, 3) - sum(comb(s, 2) for s in scores)

def vertex_in_3cycle(adj, n, v):
    out_v = [j for j in range(n) if adj[v][j]]
    in_v = [j for j in range(n) if adj[j][v]]
    for u in out_v:
        for w in in_v:
            if u != w and adj[u][w]:
                return True
    return False

def find_3cycle_sets(adj, n):
    """Find all 3-cycle vertex sets (as frozensets)."""
    cycle_sets = []
    for vs in combinations(range(n), 3):
        a, b, c = vs
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            cycle_sets.append(frozenset(vs))
    return cycle_sets

def find_all_cycles(adj, n):
    """Find all directed odd cycles as frozensets (with multiplicity)."""
    all_cyc = []
    max_k = n if n % 2 == 1 else n - 1
    for k in range(3, max_k + 1, 2):
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
            for _ in range(count):
                all_cyc.append(frozenset(verts))
    return all_cyc

def max_independent_set_3cycles(cycle_sets):
    """Find max pairwise-disjoint set among 3-cycle vertex sets."""
    n = len(cycle_sets)
    if n == 0:
        return 0
    best = 1
    # Try size 3 first
    for a in range(n):
        for b in range(a+1, n):
            if cycle_sets[a] & cycle_sets[b]:
                continue
            for c in range(b+1, n):
                if cycle_sets[c] & cycle_sets[a] or cycle_sets[c] & cycle_sets[b]:
                    continue
                return 3  # Found 3 pairwise-disjoint
            best = max(best, 2)
    return best

def compute_alpha_k(all_cyc, k):
    """Count independent sets of size k in Omega (disjoint cycle tuples)."""
    n = len(all_cyc)
    if k == 1:
        return n
    if k == 2:
        count = 0
        for a in range(n):
            for b in range(a+1, n):
                if not (all_cyc[a] & all_cyc[b]):
                    count += 1
        return count
    if k == 3:
        count = 0
        for a in range(n):
            for b in range(a+1, n):
                if all_cyc[a] & all_cyc[b]:
                    continue
                for c in range(b+1, n):
                    if (not (all_cyc[c] & all_cyc[a])) and (not (all_cyc[c] & all_cyc[b])):
                        count += 1
        return count
    return 0  # Skip k >= 4

def main():
    n = 9
    random.seed(42)

    print(f"=== n={n}: alpha_3 check for cycle-rich tournaments with t3 <= 10 ===")

    alpha3_dist = Counter()
    found = 0
    alpha3_zero_examples = []
    decomp_vals = []

    for trial in range(2000000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        scores = [sum(adj[i]) for i in range(n)]
        if 0 in scores or (n-1) in scores:
            continue

        t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)
        if t3 > 10 or t3 < 3:
            continue

        all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
        if not all_in:
            continue

        found += 1

        # Check alpha_3 for 3-cycle vertex sets
        c3_sets = find_3cycle_sets(adj, n)
        max_disj = max_independent_set_3cycles(c3_sets)

        alpha3_dist[max_disj >= 3] += 1

        if max_disj < 3 and len(alpha3_zero_examples) < 5:
            alpha3_zero_examples.append({
                'trial': trial,
                't3': t3,
                'scores': tuple(sorted(scores)),
                'c3_sets': c3_sets,
                'max_disj': max_disj
            })
            print(f"  alpha_3<1 example! trial={trial}, t3={t3}, scores={tuple(sorted(scores))}, max_disj={max_disj}")
            print(f"    3-cycle sets: {c3_sets}")

        # For found examples, compute full decomposition
        if found <= 50 or (max_disj < 3):
            all_cyc = find_all_cycles(adj, n)
            alpha1 = len(all_cyc)
            alpha2 = compute_alpha_k(all_cyc, 2)
            alpha3 = compute_alpha_k(all_cyc, 3)
            decomp = alpha1 + 2*alpha2 + 4*alpha3
            decomp_vals.append(decomp)

            if found <= 5:
                print(f"  #{found}: t3={t3}, alpha_1={alpha1}, alpha_2={alpha2}, alpha_3={alpha3}, sum={decomp}")

        if found % 500 == 0:
            print(f"  ... {found} cycle-rich found, alpha3_dist: {dict(alpha3_dist)}")

    print(f"\nTotal trials: {trial+1}, cycle-rich found: {found}")
    print(f"alpha_3 >= 1 (via 3-cycle sets): {dict(alpha3_dist)}")
    if decomp_vals:
        print(f"Decomposition alpha_1+2*alpha_2+4*alpha_3 range: [{min(decomp_vals)}, {max(decomp_vals)}]")
        if min(decomp_vals) > 10:
            print(f"*** Minimum decomposition = {min(decomp_vals)} > 10 => H=21 impossible! ***")

    if alpha3_zero_examples:
        print(f"\nExamples with alpha_3=0 (via 3-cycle sets only):")
        for ex in alpha3_zero_examples:
            print(f"  t3={ex['t3']}, scores={ex['scores']}, max_disj={ex['max_disj']}")
    else:
        print(f"\n*** ALL cycle-rich n=9 tournaments have 3 disjoint 3-cycles! ***")
        print(f"*** Part C blocks H=21 at n=9. ***")


if __name__ == "__main__":
    main()
