#!/usr/bin/env python3
"""
p17_max_orbit_analysis.py — Analyze the tied H-maximizers at p=17

At p=17 (≡1 mod 8), multiple orientations achieve the same maximum H.
What algebraic structure explains these ties?

For p≡1 mod 4, the QR group C_m acts on the orientation cube.
Since -I ∈ C_m, the orbits come in complementary pairs: σ and -σ.
The tie structure should correspond to C_m orbits.

Author: opus-2026-03-12-S63
"""

import numpy as np
from itertools import product
import math

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

def chord_type(a, p):
    m = (p-1)//2
    a = a % p
    if a == 0: return 0
    return a if a <= m else p - a

def make_signed_permutation(a, p):
    m = (p-1)//2
    P = np.zeros((m, m))
    for k in range(1, m+1):
        ak = (a * k) % p
        if ak <= m:
            P[ak-1, k-1] = 1
        else:
            P[(p-ak)-1, k-1] = -1
    return P

def held_karp(A):
    n = len(A)
    dp = {}
    for start in range(n):
        dp[(1 << start, start)] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            if (mask, last) not in dp:
                continue
            count = dp[(mask, last)]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if A[last][nxt]:
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def tournament_from_sigma(sigma, p):
    m = (p-1)//2
    n = p
    A = np.zeros((n, n), dtype=int)
    for k in range(1, m+1):
        for i in range(n):
            j = (i + k) % n
            if sigma[k-1] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


p = 17
m = (p-1)//2  # = 8

print("=" * 70)
print(f"ORBIT ANALYSIS OF H-MAXIMIZERS AT p={p}")
print("=" * 70)

# QR elements
qr = sorted(a for a in range(1, p) if legendre(a, p) == 1)
print(f"QR = {qr}")
print(f"|QR| = m = {m}")

# Build QR signed permutations
P_matrices = {}
for a in qr:
    P_matrices[a] = make_signed_permutation(a, p)

# Compute all H values (from previous run, but let's get the full distribution)
print(f"\nComputing all {2**m} orientations...")
all_sigmas = list(product([1, -1], repeat=m))
H_dict = {}
for sigma in all_sigmas:
    A = tournament_from_sigma(np.array(sigma), p)
    H_dict[sigma] = held_karp(A)

# Group by H value
H_groups = {}
for sigma, h in H_dict.items():
    if h not in H_groups:
        H_groups[h] = []
    H_groups[h].append(sigma)

print(f"\nH value distribution ({len(H_groups)} distinct values):")
for h in sorted(H_groups.keys(), reverse=True):
    members = H_groups[h]
    print(f"  H = {h}: {len(members)} orientations")

# Analyze the maximum group
H_max = max(H_groups.keys())
max_group = H_groups[H_max]
print(f"\n--- MAXIMUM H = {H_max}: {len(max_group)} orientations ---")
for sigma in max_group:
    A_val = sum(legendre(k, p) * sigma[k-1] for k in range(1, m+1))
    print(f"  σ = {sigma}, A(σ) = {A_val}")

# Check: are these all in the same QR orbit?
print(f"\n--- QR ORBIT STRUCTURE ---")
# Find QR orbits
visited = set()
orbits = []
for sigma in all_sigmas:
    if sigma in visited:
        continue
    orbit = set()
    # Generate orbit under all P_a
    queue = [sigma]
    while queue:
        s = queue.pop()
        if s in orbit:
            continue
        orbit.add(s)
        for a in qr:
            Pa = P_matrices[a]
            s_arr = np.array(s, dtype=float)
            new_s = Pa @ s_arr
            new_sigma = tuple(int(x) for x in new_s)
            if new_sigma not in orbit:
                queue.append(new_sigma)
    visited.update(orbit)
    orbits.append(frozenset(orbit))

print(f"Total QR orbits: {len(orbits)}")
# Check that H is constant on each orbit
for idx, orbit in enumerate(orbits):
    h_values = set(H_dict[s] for s in orbit)
    if len(h_values) > 1:
        print(f"  WARNING: orbit {idx} has multiple H values: {h_values}")

# Find which orbit(s) contain the maximum
max_orbits = []
for idx, orbit in enumerate(orbits):
    if any(H_dict[s] == H_max for s in orbit):
        max_orbits.append((idx, orbit))

print(f"\nOrbits at H_max:")
for idx, orbit in max_orbits:
    print(f"  Orbit {idx}: {len(orbit)} members")
    for s in sorted(orbit):
        A_val = sum(legendre(k, p) * s[k-1] for k in range(1, m+1))
        print(f"    {s}  A={A_val}")

# How many orbits per H value?
print(f"\n--- ORBITS PER H VALUE ---")
for h in sorted(H_groups.keys(), reverse=True):
    members = set(H_groups[h])
    orbit_ids = set()
    for idx, orbit in enumerate(orbits):
        if orbit & members:
            orbit_ids.add(idx)
    print(f"  H = {h}: {len(members)} orientations in {len(orbit_ids)} QR orbit(s)")

# The key question: do the max orientations form a SINGLE orbit?
print(f"\n--- KEY QUESTION ---")
if len(max_orbits) == 1:
    print(f"All {len(max_group)} H-maximizers are in a SINGLE QR orbit.")
    print(f"This means they are ALL algebraically equivalent!")
    print(f"The orbit structure REPLACES the eigenvector as the H-maximization mechanism.")
else:
    print(f"H-maximizers span {len(max_orbits)} QR orbits.")
    print(f"Multiple algebraically distinct tournaments achieve the maximum.")

# Connection set analysis
print(f"\n--- CONNECTION SET ANALYSIS ---")
for sigma in max_group[:5]:
    S = [k for k in range(1, m+1) if sigma[k-1] == 1]
    print(f"  σ = {sigma}")
    print(f"  S = {S} (connection set within {{1,...,m}})")
    # Full connection set including negatives
    S_full = []
    for k in S:
        S_full.append(k)
        S_full.append(p - k)
    S_full.sort()
    print(f"  S_full = {S_full} (all elements)")
    # Is it an interval (mod p)?
    # Check if S_full is a contiguous set mod p
    is_interval = all((S_full[i+1] - S_full[i]) == 1 for i in range(len(S_full)-1))
    # Or check arithmetic structure
    diffs = [S_full[i+1] - S_full[i] for i in range(len(S_full)-1)]
    print(f"  Consecutive differences: {diffs}")
    # Additive energy
    sums = {}
    for a in S_full:
        for b in S_full:
            s = (a + b) % p
            sums[s] = sums.get(s, 0) + 1
    E = sum(v*v for v in sums.values())
    print(f"  Additive energy E(S) = {E}")

# Now analyze p=13 similarly for comparison
print(f"\n\n{'='*70}")
print(f"COMPARISON: p=13 ORBIT STRUCTURE")
print(f"{'='*70}")

p13 = 13
m13 = 6
qr13 = sorted(a for a in range(1, p13) if legendre(a, p13) == 1)

P13_matrices = {}
for a in qr13:
    P13_matrices[a] = make_signed_permutation(a, p13)

all_sigmas13 = list(product([1, -1], repeat=m13))
H13_dict = {}
for sigma in all_sigmas13:
    A = tournament_from_sigma(np.array(sigma), p13)
    H13_dict[sigma] = held_karp(A)

# Find orbits
visited13 = set()
orbits13 = []
for sigma in all_sigmas13:
    if sigma in visited13:
        continue
    orbit = set()
    queue = [sigma]
    while queue:
        s = queue.pop()
        if s in orbit:
            continue
        orbit.add(s)
        for a in qr13:
            Pa = P13_matrices[a]
            s_arr = np.array(s, dtype=float)
            new_s = Pa @ s_arr
            new_sigma = tuple(int(x) for x in new_s)
            if new_sigma not in orbit:
                queue.append(new_sigma)
    visited13.update(orbit)
    orbits13.append(frozenset(orbit))

H13_max = max(H13_dict.values())
print(f"QR = {qr13}")
print(f"Total QR orbits: {len(orbits13)}")
print(f"Maximum H = {H13_max}")

max_group13 = [s for s in all_sigmas13 if H13_dict[s] == H13_max]
print(f"Number of H-maximizers: {len(max_group13)}")

max_orbits13 = []
for idx, orbit in enumerate(orbits13):
    if any(H13_dict[s] == H13_max for s in orbit):
        max_orbits13.append((idx, orbit))

print(f"Number of orbits at max: {len(max_orbits13)}")

for idx, orbit in max_orbits13:
    print(f"  Orbit: {len(orbit)} members")
    for s in sorted(orbit):
        A_val = sum(legendre(k, p13) * s[k-1] for k in range(1, m13+1))
        S = [k for k in range(1, m13+1) if s[k-1] == 1]
        print(f"    {s}  A={A_val}  S={S}")

# Orbits per H value
print(f"\nOrbits per H value:")
for h in sorted(set(H13_dict.values()), reverse=True):
    members = set(s for s in all_sigmas13 if H13_dict[s] == h)
    orbit_ids = set()
    for idx, orbit in enumerate(orbits13):
        if orbit & members:
            orbit_ids.add(idx)
    print(f"  H = {h}: {len(members)} orientations in {len(orbit_ids)} orbit(s)")


print("\n\nDONE.")
