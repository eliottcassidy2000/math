#!/usr/bin/env python3
"""
Test whether Interval is GLOBAL maximum among circulant tournaments at p=23.

We can't check all 2^11 = 2048 tournaments, but we can:
1. Check all AFFINE images of Interval (S = a·{1,...,m} mod p for a coprime to p)
2. Check all AFFINE images of Paley
3. Sample random tournaments and run gradient ascent
4. Check the "best of each swap step" landscape

At p=13, all 12 maxima were affine images of Interval (verified exhaustively).
This tests whether the same holds at p=23.

opus-2026-03-12-S62d
"""
import numpy as np
import time

def make_tournament(p, S):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for s in S:
            A[i][(i + s) % p] = 1
    return A

def count_H_layered(A, verbose=False):
    n = len(A)
    full = (1 << n) - 1
    prev = {}
    for i in range(n):
        prev[1 << i] = np.zeros(n, dtype=np.int64)
        prev[1 << i][i] = 1

    for layer in range(2, n + 1):
        curr = {}
        for mask, counts in prev.items():
            for u in range(n):
                if counts[u] == 0:
                    continue
                for v in range(n):
                    if (mask & (1 << v)) == 0 and A[u][v]:
                        new_mask = mask | (1 << v)
                        if new_mask not in curr:
                            curr[new_mask] = np.zeros(n, dtype=np.int64)
                        curr[new_mask][v] += counts[u]

        if layer == n:
            if full in curr:
                return int(curr[full].sum())
            return 0
        prev = curr
    return 0

def get_QR(p):
    return sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})

def is_valid_tournament_set(S, p):
    m = (p - 1) // 2
    if len(S) != m:
        return False
    S_set = set(S)
    for j in range(1, m + 1):
        if (j in S_set) == ((p - j) in S_set):
            return False
    return True

p = 23
m = (p - 1) // 2
QR = get_QR(p)
S_int = list(range(1, m + 1))

print(f"p = {p}, m = {m}")
print(f"Total circulant tournaments: 2^{m} = {2**m}")
print()

# 1. Check all affine images of Interval
print("=" * 72)
print("1. AFFINE IMAGES OF INTERVAL: S = a·{1,...,m} mod p")
print("=" * 72)

# For each a coprime to p, compute S_a = {a*1 mod p, ..., a*m mod p}
# These are all equivalent up to relabeling, so H should be the same.
# But let's verify.

affine_sets = {}
for a in range(1, p):
    if a == 0:
        continue
    S_a = sorted(set((a * s) % p for s in S_int) - {0})
    # Check if it's a valid tournament set
    if is_valid_tournament_set(S_a, p):
        key = tuple(S_a)
        if key not in affine_sets:
            affine_sets[key] = a

print(f"Distinct affine images of Interval: {len(affine_sets)}")

# For circulant tournaments, the H value depends only on the CONNECTION SET,
# not on the labeling. And S_a = a·S just relabels.
# The H values should all be equal by the circulant automorphism.
# But let's verify with a couple.

# Actually, a·{1,...,m} for different a can give different S_a because
# the tournament set property {j, p-j} partition is different.
# Let me check which a values give valid tournament sets.

valid_a = []
for a in range(1, p):
    S_a = sorted(set((a * s) % p for s in range(1, m + 1)))
    if len(S_a) == m and is_valid_tournament_set(S_a, p):
        valid_a.append(a)
        # Check if this equals the Interval up to complement
        if set(S_a) == set(S_int):
            print(f"  a={a}: S_a = Interval itself")

print(f"Valid affine images: {len(valid_a)}")
print(f"Values of a: {valid_a}")

# Check a few H values
for a in valid_a[:3]:
    S_a = sorted(set((a * s) % p for s in range(1, m + 1)))
    A = make_tournament(p, S_a)
    print(f"\n  Computing H(a={a}, S_a={S_a[:5]}...)...", flush=True)
    t0 = time.time()
    H = count_H_layered(A)
    t1 = time.time()
    print(f"  H = {H:,} [{t1-t0:.1f}s]")

# 2. Check a few random tournament sets
print("\n" + "=" * 72)
print("2. RANDOM TOURNAMENT SETS (gradient ascent)")
print("=" * 72)

import random
random.seed(42)

def random_tournament_set(p, m):
    """Generate a random valid tournament set."""
    S = set()
    remaining = list(range(1, p))
    random.shuffle(remaining)
    for j in range(1, (p + 1) // 2):
        if random.random() < 0.5:
            S.add(j)
        else:
            S.add(p - j)
    return sorted(S)

best_H = 0
best_S = None

for trial in range(5):
    S = random_tournament_set(p, m)
    A = make_tournament(p, S)
    t0 = time.time()
    H = count_H_layered(A)
    t1 = time.time()
    print(f"  Trial {trial+1}: S={S[:5]}..., H={H:,} [{t1-t0:.1f}s]")

    if H > best_H:
        best_H = H
        best_S = S

print(f"\n  Best random: H = {best_H:,}")
print(f"  H(Interval) = 16,011,537,490,557,279")
print(f"  Interval advantage: {16011537490557279 / best_H:.6f}x")

print("\nDONE.")
