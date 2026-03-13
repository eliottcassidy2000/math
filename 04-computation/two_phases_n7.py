#!/usr/bin/env python3
"""
two_phases_n7.py — opus-2026-03-13-S67k
Check if the two H-maximization phases persist at n=7.

At n=7, the regular score (3,3,3,3,3,3,3) has 3 iso classes:
  H=171 (Type C), H=175 (Type B), H=189 (Paley P_7)

Question: Do these use different (α₁,α₂,α₃) channel balances?

At n=7, we can have:
  - c3 (3-cycles): up to C(7,3) = 35 vertex sets
  - c5 (5-cycles): C(7,5) = 21 vertex sets
  - c7 (7-cycles): 1 vertex set (all 7 vertices)
  - Disjoint pairs: (3,3) needs 6≤7 ✓, (3,5) needs 8>7 ✗
  - So α₂ = number of disjoint (3,3) cycle pairs (using 6 of 7 vertices)
  - α₃ requires 9 vertices → α₃=0

Strategy: construct the three regular n=7 tournaments explicitly.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def find_all_directed_odd_cycles(A, n):
    cycles = []
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = True
                for k in range(length):
                    if not A[perm[k]][perm[(k+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = perm.index(min(perm))
                    normalized = tuple(perm[min_idx:] + perm[:min_idx])
                    cycles.append((normalized, frozenset(combo)))
    seen = set()
    unique = []
    for c, vs in cycles:
        if c not in seen:
            seen.add(c)
            unique.append((c, vs))
    return unique

def compute_alphas(cycles, max_k=3):
    nc = len(cycles)
    vertex_sets = [vs for c, vs in cycles]
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if vertex_sets[i] & vertex_sets[j]:
                adj[i][j] = True
                adj[j][i] = True

    alpha = [0] * (max_k + 1)
    alpha[0] = 1

    # For large nc, enumerate only up to pairs and triples
    for i in range(nc):
        alpha[1] += 1  # Each cycle is an independent set of size 1

    for i in range(nc):
        for j in range(i+1, nc):
            if not adj[i][j]:
                alpha[2] += 1

    if max_k >= 3:
        for i in range(nc):
            for j in range(i+1, nc):
                if adj[i][j]:
                    continue
                for k in range(j+1, nc):
                    if not adj[i][k] and not adj[j][k]:
                        alpha[3] += 1

    return alpha

def eigenvalue_stats(A, n):
    M = np.array(A, dtype=complex)
    eigs = np.linalg.eigvals(M)
    mags = np.abs(eigs)
    trivial = (n-1)/2
    nt_mags = [abs(e) for e in eigs if abs(abs(e) - trivial) > 0.1]
    if len(nt_mags) > 0:
        uniformity = np.std(nt_mags) / np.mean(nt_mags)
    else:
        uniformity = 0
    return {'uniformity': uniformity, 'nt_mags': nt_mags}

# Construct the Paley tournament P_7 (QR construction)
# In Z_7, the quadratic residues are {1, 2, 4} (1²=1, 2²=4, 3²=2)
# i→j if j-i is a QR mod 7
n = 7
QR = {1, 2, 4}
paley = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j - i) % n in QR:
            paley[i][j] = 1

print("=" * 70)
print("TWO PHASES AT n=7?")
print("=" * 70)

print(f"\nPaley P_7 adjacency:")
for row in paley:
    print("  ", row)

H_paley = count_ham_paths(paley, n)
print(f"\nH(Paley P_7) = {H_paley}")

# Find all 3 regular tournaments at n=7
# We'll enumerate all tournaments with score (3,3,3,3,3,3,3) by random sampling
# since full enumeration is too slow

import random

def random_regular_tournament(n, max_tries=100000):
    """Try to generate a random regular tournament (all scores = (n-1)/2)."""
    target = (n-1) // 2
    for _ in range(max_tries):
        A = [[0]*n for _ in range(n)]
        # Random tournament
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        scores = [sum(A[i]) for i in range(n)]
        if all(s == target for s in scores):
            return A
    return None

# Instead, let's construct the known regular tournaments on 7 vertices.
# There are exactly 3 iso classes.
# Paley P_7 is one. The other two can be found by systematic construction.

# Method: enumerate regular tournaments using a structured search
print("\nSearching for all 3 regular tournament classes at n=7...")
regular_classes = {}

# The Paley is class 1
def canonical_form_n7(A):
    """Canonical form using only a subset of permutations for speed."""
    n = 7
    best = None
    # Use all 5040 permutations of 7 elements
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or form < best:
            best = form
    return best

cf_paley = canonical_form_n7(paley)
regular_classes[cf_paley] = paley
print(f"  Found Paley class (H={H_paley})")

# Search by random sampling
random.seed(42)
attempts = 0
while len(regular_classes) < 3 and attempts < 500000:
    T = random_regular_tournament(7, max_tries=1)
    if T is not None:
        cf = canonical_form_n7(T)
        if cf not in regular_classes:
            H = count_ham_paths(T, 7)
            regular_classes[cf] = T
            print(f"  Found new class! H={H} (total classes: {len(regular_classes)})")
    attempts += 1
    if attempts % 100000 == 0:
        print(f"  ... {attempts} attempts, {len(regular_classes)} classes found")

print(f"\nFound {len(regular_classes)} regular tournament classes at n=7")

# Analyze each class
print(f"\n{'='*70}")
print(f"ANALYSIS OF REGULAR n=7 TOURNAMENTS")
print(f"{'='*70}")

for cf, A in sorted(regular_classes.items(), key=lambda x: count_ham_paths(x[1], 7)):
    H = count_ham_paths(A, 7)
    cycles = find_all_directed_odd_cycles(A, 7)
    alphas = compute_alphas(cycles, max_k=3)
    stats = eigenvalue_stats(A, 7)

    # Count by length
    by_length = defaultdict(int)
    for c, vs in cycles:
        by_length[len(c)] += 1

    print(f"\n  H = {H}")
    print(f"  Total directed cycles: {len(cycles)}")
    print(f"  By length: {dict(sorted(by_length.items()))}")
    print(f"  α₀={alphas[0]}, α₁={alphas[1]}, α₂={alphas[2]}, α₃={alphas[3]}")
    print(f"  H check: 1 + 2*{alphas[1]} + 4*{alphas[2]} + 8*{alphas[3]} = {1 + 2*alphas[1] + 4*alphas[2] + 8*alphas[3]}")
    print(f"  Match: {H == 1 + 2*alphas[1] + 4*alphas[2] + 8*alphas[3]}")
    print(f"  2α₁/H = {2*alphas[1]/H:.3f}, 4α₂/H = {4*alphas[2]/H:.3f}, 8α₃/H = {8*alphas[3]/H:.3f}")
    print(f"  Spectral uniformity: {stats['uniformity']:.4f}")
    print(f"  Nontrivial |λ|: {[f'{x:.4f}' for x in sorted(stats['nt_mags'])]}")

print(f"\n{'='*70}")
print("TWO-PHASE ANALYSIS")
print(f"{'='*70}")
print("""
If the two-phase phenomenon persists at n=7:
- The Paley class (H=189) should have the most uniform eigenvalues
  but possibly NOT the highest α₂.
- Another class might achieve high H through different channel balance.

At n=7, the possible disjoint pairs are (3,3)-pairs using 6 of 7 vertices.
For each of the C(7,6)=7 choices of 6 vertices, we can have at most
C(6,3)/2 = 10 ways to partition into two 3-sets.
But only some of these will both form directed 3-cycles.

The α₂ maximum at n=7 tells us how many such disjoint 3-cycle pairs exist.
""")
