#!/usr/bin/env python3
"""beta2_algebraic_proof.py — Analyze when cancel pairs are in Omega3

Author: opus-2026-03-08-S43
"""
import sys, numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

print("="*70)
print("CANCEL PAIR OMEGA3 MEMBERSHIP + DETAILED STRUCTURE")
print("="*70)

n = 5

# For each tournament, examine non-DT Omega3 elements
print(f"\nDetailed Omega3 non-DT structure at n={n}:")
examined = 0

for A in all_tournaments(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a2_set = set(a2)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d3 = om3.shape[1] if om3.ndim == 2 else 0

    has_nondt = False
    for j in range(d3):
        col = om3[:, j]
        support = [(a3[i], col[i]) for i in range(len(col)) if abs(col[i]) > 1e-10]
        if len(support) > 1:
            has_nondt = True
            break

    if not has_nondt: continue
    examined += 1
    if examined > 5: break

    t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]))

    print(f"\nTournament t3={t3}:")
    for j in range(d3):
        col = om3[:, j]
        support = [(a3[i], col[i]) for i in range(len(col)) if abs(col[i]) > 1e-10]
        if len(support) <= 1:
            p = support[0][0]
            print(f"  Omega3[{j}]: DT {p}")
            continue

        print(f"  Omega3[{j}]: {len(support)} terms:")
        for p, c in support:
            bad = []
            if not A[p[0]][p[2]]: bad.append(f"({p[0]},{p[2]},{p[3]})")
            if not A[p[1]][p[3]]: bad.append(f"({p[0]},{p[1]},{p[3]})")
            n_bad = len(bad)
            print(f"    {c:+.4f} * {p}  [{n_bad} bad: {bad}]")

# Count: how many Omega3 non-DT elements have 2 terms vs 3+ terms?
print(f"\n{'='*70}")
print("OMEGA3 NON-DT ELEMENT SIZES")
print("="*70)

size_counts = defaultdict(int)
bad_pattern_counts = defaultdict(int)

for A in all_tournaments(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d3 = om3.shape[1] if om3.ndim == 2 else 0

    for j in range(d3):
        col = om3[:, j]
        support = [(a3[i], col[i]) for i in range(len(col)) if abs(col[i]) > 1e-10]
        if len(support) > 1:
            size_counts[len(support)] += 1

            # Bad face pattern
            bads_per_path = []
            for p, c in support:
                nb = (0 if A[p[0]][p[2]] else 1) + (0 if A[p[1]][p[3]] else 1)
                bads_per_path.append(nb)
            pattern = tuple(sorted(bads_per_path))
            bad_pattern_counts[pattern] += 1

print(f"\nSupport sizes:")
for sz, cnt in sorted(size_counts.items()):
    print(f"  {sz} terms: {cnt}")

print(f"\nBad face patterns (sorted tuple of #bad per path):")
for pattern, cnt in sorted(bad_pattern_counts.items()):
    print(f"  {pattern}: {cnt}")

# KEY QUESTION: do ALL 2-term non-DT elements share exactly one bad face?
print(f"\n{'='*70}")
print("2-TERM ELEMENTS: DO THEY SHARE ONE BAD FACE?")
print("="*70)

shares_one = 0
shares_both = 0
shares_none = 0

for A in all_tournaments(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a2_set = set(a2)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d3 = om3.shape[1] if om3.ndim == 2 else 0

    for j in range(d3):
        col = om3[:, j]
        support = [(a3[i], col[i]) for i in range(len(col)) if abs(col[i]) > 1e-10]
        if len(support) != 2: continue

        p1, c1 = support[0]
        p2, c2 = support[1]

        # Get bad faces for each path
        bads1 = set()
        if not A[p1[0]][p1[2]]: bads1.add((p1[0], p1[2], p1[3]))
        if not A[p1[1]][p1[3]]: bads1.add((p1[0], p1[1], p1[3]))

        bads2 = set()
        if not A[p2[0]][p2[2]]: bads2.add((p2[0], p2[2], p2[3]))
        if not A[p2[1]][p2[3]]: bads2.add((p2[0], p2[1], p2[3]))

        shared = bads1 & bads2
        if len(shared) == 1:
            shares_one += 1
        elif len(shared) == 2:
            shares_both += 1
        else:
            shares_none += 1
            print(f"  NO SHARED BAD FACE: {p1} bads={bads1}, {p2} bads={bads2}")

print(f"\n2-term elements: shares_one={shares_one}, shares_both={shares_both}, shares_none={shares_none}")

print("\nDONE")
