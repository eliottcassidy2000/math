#!/usr/bin/env python3
"""
PHASE 2: Deep dive into path homology findings

Key discoveries from Phase 1:
1. Tournaments n<=5: only β=(1,0,...) or β=(1,1,0,...)
2. n=7 sampling: β=(1,0,0,1,0,0) appeared! First higher Betti number for tournaments!
3. Circulant digraphs show very rich topology

This script investigates:
A. The β_3=1 tournaments at n=7 — what structure causes them?
B. When do tournaments first get β_2>0?
C. Systematic circulant patterns (Euler characteristic, Poincaré duality)
D. Connection set symmetry and topology
"""
import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
import random
import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, circulant_digraph, count_3cycles, ham_path_count
)

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# ===== SECTION A: Hunt for β_3=1 tournaments at n=7 =====
print("=" * 70)
print("SECTION A: HUNTING β_3>0 TOURNAMENTS AT n=7")
print("=" * 70)

found_b3 = []
found_b2 = []
n = 7
total = 0
for trial in range(500):
    A = random_tournament(n)
    total += 1
    betti = path_betti_numbers(A, n, max_dim=6)

    if betti[3] > 0:
        t3 = count_3cycles(A, n)
        H = ham_path_count(A, n)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        found_b3.append((A, betti, t3, H, scores))
        if len(found_b3) <= 5:
            print(f"\n  β_3>0 tournament #{len(found_b3)}:")
            print(f"    β = {betti}, t3={t3}, H={H}")
            print(f"    scores = {scores}")
            # Show adjacency
            for i in range(n):
                print(f"    {i}: → {[j for j in range(n) if A[i][j]]}")

    if betti[2] > 0:
        t3 = count_3cycles(A, n)
        H = ham_path_count(A, n)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        found_b2.append((betti, t3, H, scores))

print(f"\n  Summary after {total} trials at n={n}:")
print(f"    β_3>0: {len(found_b3)} ({100*len(found_b3)/total:.1f}%)")
print(f"    β_2>0: {len(found_b2)} ({100*len(found_b2)/total:.1f}%)")

if found_b3:
    print(f"\n  β_3 values: {Counter([b[1][3] for b in found_b3])}")
    print(f"  t3 range: {min(b[2] for b in found_b3)} to {max(b[2] for b in found_b3)}")
    print(f"  H range: {min(b[3] for b in found_b3)} to {max(b[3] for b in found_b3)}")
    print(f"  Score sequences: {Counter([b[4] for b in found_b3])}")

# ===== SECTION B: Euler characteristic =====
print("\n\n" + "=" * 70)
print("SECTION B: EULER CHARACTERISTIC χ = Σ(-1)^p β_p")
print("=" * 70)

# For tournaments
for n in [4, 5]:
    print(f"\nn={n}:")
    chi_dist = Counter()
    for mask in range(1 << (n*(n-1)//2)):
        edges = [(i,j) for i in range(n) for j in range(i+1,n)]
        m = len(edges)
        if mask >= (1 << m):
            break
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        betti = path_betti_numbers(A, n, max_dim=n-1)
        chi = sum((-1)**p * betti[p] for p in range(len(betti)))
        chi_dist[chi] += 1

    print(f"  Euler char distribution: {dict(sorted(chi_dist.items()))}")

# For circulant digraphs
print(f"\nCirculant Euler characteristics:")
for n in [5, 6, 7]:
    print(f"\n  n={n}:")
    for size in range(1, n):
        for S in combinations(range(1, n), size):
            A = circulant_digraph(n, list(S))
            betti = path_betti_numbers(A, n, max_dim=min(n-1, 6))
            chi = sum((-1)**p * betti[p] for p in range(len(betti)))
            if chi != 1:  # Only print non-trivial
                print(f"    C_{n}^{set(S)}: β={betti}, χ={chi}")

# ===== SECTION C: Analyze β_3=1 tournament structure =====
print("\n\n" + "=" * 70)
print("SECTION C: STRUCTURE OF β_3>0 TOURNAMENTS")
print("=" * 70)

if found_b3:
    # Take first β_3>0 tournament and analyze its cycle structure
    A, betti, t3, H, scores = found_b3[0]
    n = 7

    # Count k-cycles for k=3,4,5
    def count_kcycles(A, n, k):
        """Count directed k-cycles."""
        count = 0
        for combo in combinations(range(n), k):
            from itertools import permutations
            for perm in permutations(combo):
                is_cycle = True
                for i in range(k):
                    if A[perm[i]][perm[(i+1)%k]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    count += 1
        return count // k  # each k-cycle counted k times (rotations)

    t4 = count_kcycles(A, n, 4)
    t5 = count_kcycles(A, n, 5)
    print(f"  β_3>0 tournament: t3={t3}, t4={t4}, t5={t5}, H={H}")
    print(f"  scores={scores}")

    # Compare with a random β_1>0, β_3=0 tournament
    for trial in range(200):
        A2 = random_tournament(n)
        b2 = path_betti_numbers(A2, n, max_dim=6)
        if b2[1] > 0 and b2[3] == 0:
            t3_2 = count_3cycles(A2, n)
            t4_2 = count_kcycles(A2, n, 4)
            t5_2 = count_kcycles(A2, n, 5)
            H_2 = ham_path_count(A2, n)
            print(f"  β_1>0 comparison: t3={t3_2}, t4={t4_2}, t5={t5_2}, H={H_2}")
            break

# ===== SECTION D: Circulant symmetry patterns =====
print("\n\n" + "=" * 70)
print("SECTION D: CIRCULANT SYMMETRY AND COMPLEMENT PATTERNS")
print("=" * 70)

print("\nKey patterns in circulant path homology:")
print("\n  1. Complement pairs C_n^S vs C_n^{[n-1]\\S}:")

for n in [5, 6, 7]:
    full = set(range(1, n))
    seen = set()
    for size in range(1, n):
        for S in combinations(range(1, n), size):
            S_set = frozenset(S)
            comp = frozenset(full - S_set)
            if comp in seen or S_set in seen:
                continue
            seen.add(S_set)

            A1 = circulant_digraph(n, list(S))
            A2 = circulant_digraph(n, list(comp)) if comp else None
            b1 = path_betti_numbers(A1, n, max_dim=min(n-1, 6))
            b2 = path_betti_numbers(A2, n, max_dim=min(n-1, 6)) if A2 else [1]

            if any(b > 0 for b in b1[1:]) or any(b > 0 for b in b2[1:]):
                print(f"    C_{n}^{set(S)}: β={b1}")
                if comp:
                    print(f"    C_{n}^{set(comp)}: β={b2}")
                    # Check Poincaré-like duality
                    if len(b1) == len(b2):
                        b2_rev = list(reversed(b2))
                        if b1 == b2_rev:
                            print(f"      → POINCARÉ DUAL!")
                print()

# ===== SECTION E: n=6 exhaustive for tournaments =====
print("\n\n" + "=" * 70)
print("SECTION E: TOURNAMENTS n=6 — LARGER SAMPLE")
print("=" * 70)

n = 6
betti_dist = Counter()
b1_by_t3 = defaultdict(list)
total_6 = 0
for trial in range(1000):
    A = random_tournament(n)
    total_6 += 1
    betti = path_betti_numbers(A, n, max_dim=5)
    bt = tuple(betti)
    betti_dist[bt] += 1
    t3 = count_3cycles(A, n)
    b1_by_t3[t3].append(betti[1])

print(f"\nBetti distribution (n={n}, {total_6} samples):")
for bt in sorted(betti_dist.keys()):
    print(f"  β={list(bt)}: {betti_dist[bt]} ({100*betti_dist[bt]/total_6:.1f}%)")

print(f"\nβ_1 vs t3:")
for t3 in sorted(b1_by_t3.keys()):
    vals = Counter(b1_by_t3[t3])
    print(f"  t3={t3}: β_1 = {dict(vals)}")

# ===== SECTION F: Regular tournaments and β =====
print("\n\n" + "=" * 70)
print("SECTION F: REGULAR TOURNAMENTS (all scores = (n-1)/2)")
print("=" * 70)

for n in [5, 7]:
    print(f"\nn={n}: Sampling regular tournaments")
    found_regular = []
    for trial in range(10000):
        A = random_tournament(n)
        scores = sorted([sum(A[i]) for i in range(n)])
        if all(s == (n-1)//2 for s in scores):
            betti = path_betti_numbers(A, n, max_dim=min(n-1, 6))
            t3 = count_3cycles(A, n)
            found_regular.append((tuple(betti), t3))

    if found_regular:
        print(f"  Found {len(found_regular)} regular tournaments")
        bt_dist = Counter([r[0] for r in found_regular])
        for bt, cnt in bt_dist.most_common():
            print(f"    β={list(bt)}: {cnt}")
    else:
        print(f"  No regular tournaments found in sample")

print("\n\nDone.")
