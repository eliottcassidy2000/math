#!/usr/bin/env python3
"""Quick test: is β_2 = 0 for ALL oriented graphs (no mutual edges)?
If so, the result is much stronger than just tournaments."""
import numpy as np
from itertools import combinations
import sys, random
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import path_betti_numbers

# All oriented graphs on 4 vertices
print("=" * 70)
print("β_2 FOR ALL ORIENTED GRAPHS ON n=4")
print("=" * 70)

n = 4
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
n_pairs = len(pairs)

beta2_nonzero = 0
total = 0
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

    total += 1
    betti = path_betti_numbers(A, n, max_dim=3)
    if betti[2] > 0:
        beta2_nonzero += 1
        if beta2_nonzero <= 5:
            print(f"  β_2 = {betti[2]}: β = {betti}")
            for i in range(n):
                nbrs = [j for j in range(n) if A[i][j]]
                in_nbrs = [j for j in range(n) if A[j][i]]
                print(f"    {i} → {nbrs}, ← {in_nbrs}")

print(f"\nn=4: {beta2_nonzero}/{total} oriented graphs have β_2 > 0")

# All oriented graphs on 5 vertices (3^10 = 59049)
print(f"\n{'='*70}")
print("β_2 FOR ALL ORIENTED GRAPHS ON n=5")
print("="*70)

n = 5
pairs5 = [(i,j) for i in range(n) for j in range(i+1, n)]
n_pairs5 = len(pairs5)
total5 = 3**n_pairs5

beta2_nonzero5 = 0
for mask in range(total5):
    A = [[0]*n for _ in range(n)]
    m = mask
    for i, (a, b) in enumerate(pairs5):
        choice = m % 3
        m //= 3
        if choice == 1:
            A[a][b] = 1
        elif choice == 2:
            A[b][a] = 1

    betti = path_betti_numbers(A, n, max_dim=3)
    if betti[2] > 0:
        beta2_nonzero5 += 1
        if beta2_nonzero5 <= 5:
            print(f"  β_2 = {betti[2]}: β = {betti}")
            for i in range(n):
                nbrs = [j for j in range(n) if A[i][j]]
                print(f"    {i} → {nbrs}")

    if (mask+1) % 10000 == 0:
        print(f"  ... {mask+1}/{total5}", flush=True)

print(f"\nn=5: {beta2_nonzero5}/{total5} oriented graphs have β_2 > 0")

# If β_2 > 0 found, also check: do they all have mutual edges?
if beta2_nonzero > 0 or beta2_nonzero5 > 0:
    print(f"\n{'='*70}")
    print("CHECKING: Do β_2>0 graphs require mutual edges?")
    print("="*70)

    for n in [4, 5]:
        pairs_n = [(i,j) for i in range(n) for j in range(i+1, n)]
        n_p = len(pairs_n)
        mutual_required = 0
        no_mutual_b2 = 0

        for mask in range(3**n_p):
            A = [[0]*n for _ in range(n)]
            m = mask
            for i, (a, b) in enumerate(pairs_n):
                choice = m % 3
                m //= 3
                if choice == 1:
                    A[a][b] = 1
                elif choice == 2:
                    A[b][a] = 1

            betti = path_betti_numbers(A, n, max_dim=3)
            if betti[2] > 0:
                has_mutual = any(A[i][j] and A[j][i]
                                for i in range(n) for j in range(i+1, n))
                if has_mutual:
                    mutual_required += 1
                else:
                    no_mutual_b2 += 1

        print(f"  n={n}: β_2>0 with mutual edges: {mutual_required}")
        print(f"  n={n}: β_2>0 WITHOUT mutual edges: {no_mutual_b2}")

print("\nDone.")
