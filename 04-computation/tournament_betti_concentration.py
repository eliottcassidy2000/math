#!/usr/bin/env python3
"""
TOURNAMENT BETTI CONCENTRATION CONJECTURE
==========================================
Conjecture: For any tournament T, at most one β_p (p ≥ 1) is nonzero.
Combined with β₂ = 0 (THM-108), this means:
  β(T) = (1, 0, 0, ...) [point-like]
  or β(T) = (1, 1, 0, ...) [circle-like]
  or β(T) = (1, 0, 0, 1, 0, ...) [S³-like]
  or β(T) = (1, 0, 0, 0, 1, 0, ...) [S⁴-like]
  etc.

This would mean tournaments are "homologically spherical" — each looks
like a single sphere S^k or a point.

opus-2026-03-15-S72d
"""
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
import random
import sys
import time

sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')
from path_homology_v2 import path_betti_numbers

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                wins = [A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]]
                if max(wins) < 2:
                    c3 += 1
    return c3

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

print("=" * 70)
print("TOURNAMENT BETTI CONCENTRATION")
print("=" * 70)

# =================================================================
# n = 5: EXHAUSTIVE — check full Betti numbers
# =================================================================
print("\n--- n=5 EXHAUSTIVE ---")
n = 5
E = n*(n-1)//2
betti_dist = Counter()
nonzero_count = Counter()  # how many β_p > 0 for p >= 1
multi_nonzero = []

for bits in range(1 << E):
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_dist[bt] += 1

    nz = sum(1 for p in range(1, len(betti)) if betti[p] > 0)
    nonzero_count[nz] += 1
    if nz > 1:
        multi_nonzero.append((bits, betti))

print(f"  Total: {1 << E}")
print(f"  Betti distribution:")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")
print(f"  # nonzero β_p (p≥1): {dict(nonzero_count)}")
if multi_nonzero:
    print(f"  *** {len(multi_nonzero)} VIOLATIONS with multiple nonzero β_p ***")
    for bits, betti in multi_nonzero[:5]:
        print(f"    bits={bits}: β={betti}")
else:
    print(f"  ✓ AT MOST 1 nonzero β_p for all {1<<E} tournaments")

# =================================================================
# n = 6: EXHAUSTIVE
# =================================================================
print("\n--- n=6 EXHAUSTIVE ---")
n = 6
E = n*(n-1)//2
betti_dist = Counter()
nonzero_count = Counter()
multi_nonzero = []
t0 = time.time()

for bits in range(1 << E):
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_dist[bt] += 1

    nz = sum(1 for p in range(1, len(betti)) if betti[p] > 0)
    nonzero_count[nz] += 1
    if nz > 1:
        multi_nonzero.append((bits, betti))

print(f"  Total: {1 << E}, time: {time.time()-t0:.1f}s")
print(f"  Betti distribution:")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")
print(f"  # nonzero β_p (p≥1): {dict(nonzero_count)}")
if multi_nonzero:
    print(f"  *** {len(multi_nonzero)} VIOLATIONS ***")
    for bits, betti in multi_nonzero[:5]:
        print(f"    bits={bits}: β={betti}")
else:
    print(f"  ✓ AT MOST 1 nonzero β_p for all {1<<E} tournaments")

# =================================================================
# n = 7: SAMPLED
# =================================================================
print("\n--- n=7 SAMPLED ---")
n = 7
E = n*(n-1)//2
random.seed(42)
N_SAMPLES = 3000
betti_dist = Counter()
nonzero_count = Counter()
multi_nonzero = []
by_nonzero_dim = defaultdict(list)
t0 = time.time()

for trial in range(N_SAMPLES):
    bits = random.randint(0, (1 << E) - 1)
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_dist[bt] += 1

    nz_dims = [p for p in range(1, len(betti)) if betti[p] > 0]
    nonzero_count[len(nz_dims)] += 1
    if len(nz_dims) > 1:
        multi_nonzero.append((bits, betti))
    if nz_dims:
        by_nonzero_dim[tuple(nz_dims)].append(bt)

print(f"  Sampled: {N_SAMPLES}, time: {time.time()-t0:.1f}s")
print(f"  Betti distribution:")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")
print(f"  # nonzero β_p (p≥1): {dict(nonzero_count)}")
print(f"  Nonzero dimensions breakdown:")
for dims, bts in sorted(by_nonzero_dim.items()):
    print(f"    p∈{dims}: {len(bts)} tournaments, Betti types: {Counter(bts).most_common(5)}")
if multi_nonzero:
    print(f"  *** {len(multi_nonzero)} VIOLATIONS ***")
    for bits, betti in multi_nonzero[:10]:
        c3 = count_3cycles(tournament_from_bits(n, bits), n)
        sc = score_seq(tournament_from_bits(n, bits), n)
        print(f"    c₃={c3}, score={sc}: β={betti}")
else:
    print(f"  ✓ AT MOST 1 nonzero β_p in {N_SAMPLES} samples")

# =================================================================
# n = 8: SAMPLED (smaller sample, n=8 is expensive)
# =================================================================
print("\n--- n=8 SAMPLED ---")
n = 8
E = n*(n-1)//2
random.seed(123)
N_SAMPLES = 500
betti_dist = Counter()
nonzero_count = Counter()
multi_nonzero = []
by_nonzero_dim = defaultdict(int)
t0 = time.time()

for trial in range(N_SAMPLES):
    bits = random.randint(0, (1 << E) - 1)
    A = tournament_from_bits(n, bits)
    try:
        betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
        bt = tuple(betti)
        betti_dist[bt] += 1

        nz_dims = [p for p in range(1, len(betti)) if betti[p] > 0]
        nonzero_count[len(nz_dims)] += 1
        if len(nz_dims) > 1:
            multi_nonzero.append((bits, betti))
        if nz_dims:
            by_nonzero_dim[tuple(nz_dims)] += 1
    except Exception as e:
        print(f"  Error at trial {trial}: {e}")
        continue

print(f"  Sampled: {N_SAMPLES}, time: {time.time()-t0:.1f}s")
print(f"  Betti distribution:")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")
print(f"  # nonzero β_p (p≥1): {dict(nonzero_count)}")
print(f"  Nonzero dimensions: {dict(by_nonzero_dim)}")
if multi_nonzero:
    print(f"  *** {len(multi_nonzero)} VIOLATIONS ***")
    for bits, betti in multi_nonzero[:10]:
        print(f"    β={betti}")
else:
    print(f"  ✓ AT MOST 1 nonzero β_p in {N_SAMPLES} samples")

# =================================================================
# SUMMARY
# =================================================================
print("\n" + "=" * 70)
print("SUMMARY: BETTI CONCENTRATION CONJECTURE")
print("=" * 70)
print("""
If the conjecture holds, every tournament T is "homologically spherical":
  β(T) = (1, 0, ..., 0, β_k, 0, ..., 0)  for some k ≥ 0, with β_k ∈ {0, 1}.

This means the path homology "topological type" of a tournament is
determined by a single integer k ∈ {0, 1, 3, 4, ...}:
  k = 0: point (contractible)
  k = 1: circle S¹
  k = 3: sphere S³
  k = 4: sphere S⁴  (if β₂=0 by THM-108)
  etc.

The tournament has "one hole or none" — it's as simple as a sphere.
""")
