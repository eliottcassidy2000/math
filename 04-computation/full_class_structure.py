#!/usr/bin/env python3
"""
Analysis of the full tiling class: why |class(T_full)| = 1 + 2^{n-2} for n >= 4.
Instance: opus-2026-03-06-S11

The "full" tiling (all bits = 1) corresponds to the tournament where:
  - Path edges: i -> i+1 (forward)
  - All non-adjacent edges: j -> i for j - i >= 2 (backward)

This means vertex i beats only i+1 (forward) and all j with j <= i-2 (backward).
The score of vertex i is: 1_{i<n-1} + max(0, i-1).

H(T_full) = number of Hamiltonian paths of this tournament = |class(T_full)|.

Conjecture: H(T_full) = 1 + 2^{n-2} for n >= 4.
Known: Tribonacci for T_full was disproved (see tiling-class-structure.md).
Need to distinguish between T_full H-value and class size.
"""
from itertools import permutations
from collections import defaultdict
import sys

def full_tournament(n):
    """Build the full tiling tournament."""
    # Path edges: i -> i+1
    # Non-adjacent: j -> i for j - i >= 2
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1  # path edge
    for i in range(n):
        for j in range(i+2, n):
            A[j][i] = 1  # backward non-adjacent
    return A

def count_hamiltonian_paths(A, n):
    """Count directed Hamiltonian paths."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def enumerate_ham_paths(A, n):
    """Enumerate all Hamiltonian paths (as vertex orderings)."""
    paths = []
    def backtrack(path, visited):
        if len(path) == n:
            paths.append(tuple(path))
            return
        last = path[-1] if path else -1
        for v in range(n):
            if v in visited:
                continue
            if last == -1 or A[last][v]:
                backtrack(path + [v], visited | {v})
    backtrack([], set())
    return paths

def analyze_full_class(n):
    """Analyze the structure of T_full and its Hamiltonian paths."""
    A = full_tournament(n)

    print(f"\n{'='*60}")
    print(f"n = {n}: Full tiling tournament")
    print(f"{'='*60}")

    scores = sorted([sum(A[i]) for i in range(n)], reverse=True)
    print(f"Scores: {scores}")

    H = count_hamiltonian_paths(A, n)
    formula = 1 + 2**(n-2) if n >= 4 else None
    print(f"H(T_full) = {H}")
    if formula:
        print(f"1 + 2^(n-2) = {formula}, match: {H == formula}")

    # Enumerate paths for small n
    if n <= 8:
        paths = enumerate_ham_paths(A, n)
        print(f"Enumerated {len(paths)} Hamiltonian paths:")
        for p in paths:
            print(f"  {p}")

        # Look for patterns in the path structure
        # The canonical path is 0,1,2,...,n-1
        # What are the other paths?
        if len(paths) > 1:
            print(f"\nPath structure analysis:")
            for p in paths:
                # Where does each path deviate from 0,1,...,n-1?
                deviations = []
                for i in range(n):
                    if p[i] != i:
                        deviations.append(i)
                if deviations:
                    print(f"  {p}: deviates at positions {deviations}")
                else:
                    print(f"  {p}: canonical path")

    # Adjacency matrix display
    print(f"\nAdjacency matrix (A[i][j]=1 means i->j):")
    for i in range(n):
        row = ''.join(str(A[i][j]) for j in range(n))
        print(f"  {i}: {row}")

    return H

def main():
    print("Full tiling tournament: H(T_full) analysis")
    print("="*60)

    H_values = []
    for n in range(3, 10):
        H = analyze_full_class(n)
        H_values.append((n, H))

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'n':>3} {'H(T_full)':>10} {'1+2^(n-2)':>10} {'match':>6}")
    for n, H in H_values:
        f = 1 + 2**(n-2) if n >= 2 else '?'
        m = '✓' if H == f else '✗'
        print(f"{n:>3} {H:>10} {f:>10} {m:>6}")

if __name__ == '__main__':
    main()
