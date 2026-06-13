#!/usr/bin/env python3
"""Analyze the 70 oriented graphs on n=5 with β_2 > 0.
What structural property do they share that tournaments lack?"""
import numpy as np
from itertools import combinations
import sys
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import path_betti_numbers

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
n_pairs = len(pairs)

print("=" * 70)
print("ANALYSIS OF ORIENTED GRAPHS WITH β_2 > 0 (n=5)")
print("=" * 70)

counterexamples = []
for mask in range(3**n_pairs):
    A = [[0]*n for _ in range(n)]
    m = mask
    for i, (a, b) in enumerate(pairs):
        choice = m % 3
        m //= 3
        if choice == 1:
            A[a][b] = 1
        elif choice == 2:
            A[b][a] = 1

    betti = path_betti_numbers(A, n, max_dim=3)
    if betti[2] > 0:
        counterexamples.append((A, betti))

print(f"Total counterexamples: {len(counterexamples)}")

# Analyze structure
edge_counts = Counter()
missing_edge_counts = Counter()
out_deg_seqs = Counter()
in_deg_seqs = Counter()

for A, betti in counterexamples:
    total_edges = sum(A[i][j] for i in range(n) for j in range(n))
    missing = n*(n-1)//2 - total_edges  # pairs with NO edge
    edge_counts[total_edges] += 1
    missing_edge_counts[missing] += 1

    out_degs = sorted([sum(A[i]) for i in range(n)])
    in_degs = sorted([sum(A[j][i] for j in range(n)) for i in range(n)])
    out_deg_seqs[tuple(out_degs)] += 1

print(f"\nEdge count distribution:")
for k in sorted(edge_counts): print(f"  {k} edges: {edge_counts[k]}")

print(f"\nMissing-pair count (pairs with no edge in either direction):")
for k in sorted(missing_edge_counts): print(f"  {k} missing: {missing_edge_counts[k]}")

print(f"\nOut-degree sequences:")
for k in sorted(out_deg_seqs): print(f"  {list(k)}: {out_deg_seqs[k]}")

# Show a few examples in detail
print(f"\n--- First 5 counterexamples ---")
for idx, (A, betti) in enumerate(counterexamples[:5]):
    print(f"\nExample {idx+1}: β = {betti}")
    for i in range(n):
        out = [j for j in range(n) if A[i][j]]
        print(f"  {i} → {out}")
    # Which pairs have no edge?
    no_edge = []
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 0 and A[j][i] == 0:
                no_edge.append((i,j))
    print(f"  Missing pairs: {no_edge}")

    # Check: does adding a directed edge for each missing pair always kill β_2?
    # (i.e., completing to a tournament)
    killed = 0
    survived = 0
    for bits in range(2**len(no_edge)):
        B = [row[:] for row in A]
        for k, (a, b) in enumerate(no_edge):
            if (bits >> k) & 1:
                B[a][b] = 1
            else:
                B[b][a] = 1
        betti_B = path_betti_numbers(B, n, max_dim=3)
        if betti_B[2] == 0:
            killed += 1
        else:
            survived += 1
    print(f"  Completing to tournament: {killed} kill β_2, {survived} keep β_2>0")

# Key question: what is the minimum number of missing edges needed for β_2 > 0?
print(f"\n\n{'='*70}")
print("MINIMUM MISSING EDGES FOR β_2 > 0")
print("="*70)
min_missing = min(missing_edge_counts.keys())
print(f"Minimum missing edges: {min_missing}")
examples = [(A, b) for A, b in counterexamples
            if sum(1 for i in range(n) for j in range(i+1,n)
                   if A[i][j]==0 and A[j][i]==0) == min_missing]
print(f"Number with {min_missing} missing: {len(examples)}")

# Check n=4 too — are there oriented graphs close to tournaments with β_2>0?
print(f"\n\n{'='*70}")
print("n=4: ORIENTED GRAPHS BY MISSING EDGES")
print("="*70)
n = 4
pairs4 = [(i,j) for i in range(n) for j in range(i+1, n)]
n_p4 = len(pairs4)

by_missing = Counter()
for mask in range(3**n_p4):
    A = [[0]*n for _ in range(n)]
    m = mask
    for i, (a, b) in enumerate(pairs4):
        choice = m % 3
        m //= 3
        if choice == 1:
            A[a][b] = 1
        elif choice == 2:
            A[b][a] = 1

    missing = sum(1 for i in range(n) for j in range(i+1,n)
                  if A[i][j]==0 and A[j][i]==0)
    betti = path_betti_numbers(A, n, max_dim=3)
    by_missing[(missing, betti[2])] += 1

print(f"  (missing_pairs, β_2): count")
for k in sorted(by_missing): print(f"    {k}: {by_missing[k]}")

print("\nDone.")
