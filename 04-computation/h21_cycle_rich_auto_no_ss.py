#!/usr/bin/env python3
"""
CRITICAL CHECK: Does "every vertex in 3-cycle" automatically imply
"no source/sink"?

If YES: the source/sink avoidance in the dichotomy proof is UNNECESSARY.
The poisoning graph DAG argument alone completes the proof!

Proof sketch:
- Source (score n-1): beats everyone. Can't be in 3-cycle (no in-neighbor).
- Sink (score 0): beats nobody. Can't be in 3-cycle (no out-neighbor).

This script verifies exhaustively at small n.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from math import comb
import random

def vertex_in_3cycle(adj, n, v):
    out_v = [j for j in range(n) if adj[v][j]]
    in_v = [j for j in range(n) if adj[j][v]]
    for u in out_v:
        for w in in_v:
            if u != w and adj[u][w]:
                return True
    return False

def main():
    print("=== Verifying: 'every vertex in 3-cycle' implies 'no source/sink' ===")
    print()

    # Exhaustive check at n=4,5,6,7
    for n in range(4, 8):
        total = 0
        all_in_3cycle = 0
        all_in_and_has_ss = 0

        num_edges = n * (n - 1) // 2
        for bits in range(1 << num_edges):
            adj = [[0]*n for _ in range(n)]
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if bits & (1 << idx):
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1
                    idx += 1

            total += 1
            if all(vertex_in_3cycle(adj, n, v) for v in range(n)):
                all_in_3cycle += 1
                scores = [sum(adj[i]) for i in range(n)]
                if 0 in scores or (n-1) in scores:
                    all_in_and_has_ss += 1
                    print(f"  COUNTEREXAMPLE at n={n}: scores={scores}")

        print(f"n={n}: {total} tournaments, {all_in_3cycle} have all-in-3-cycle, "
              f"{all_in_and_has_ss} of those have source/sink")

    # Random check at n=8,9
    for n in [8, 9]:
        random.seed(42)
        checked = 0
        all_in_3cycle = 0
        counterexamples = 0

        for trial in range(500000):
            adj = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1

            checked += 1
            if all(vertex_in_3cycle(adj, n, v) for v in range(n)):
                all_in_3cycle += 1
                scores = [sum(adj[i]) for i in range(n)]
                if 0 in scores or (n-1) in scores:
                    counterexamples += 1

        print(f"n={n}: {checked} checked, {all_in_3cycle} all-in-3-cycle, "
              f"{counterexamples} with source/sink")

    print()
    print("=== PROOF ===")
    print("Source (score n-1): v -> u for ALL u. For v in 3-cycle v->a->b->v:")
    print("  Need b->v. But v beats b (score n-1). Contradiction.")
    print("Sink (score 0): u -> v for ALL u. For v in 3-cycle v->a->b->v:")
    print("  Need v->a. But v beats nobody (score 0). Contradiction.")
    print("Therefore: 'every vertex in 3-cycle' => 'no source/sink'. QED")


if __name__ == "__main__":
    main()
