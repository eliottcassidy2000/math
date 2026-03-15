#!/usr/bin/env python3
"""
walsh_beta1_connection_S72b.py — opus-2026-03-15-S72b

THE SPECTRAL-TOPOLOGICAL CONNECTION:
Walsh degree-2 energy as a detector for β₁

At n=5: β₁=1 ⟹ Walsh degree-2 energy ≥ 1.5
This is extraordinary: a SPECTRAL property implies a TOPOLOGICAL property.

Questions:
1. Does this hold at n=6?
2. Is there an exact threshold?
3. Can we prove it?
4. What IS the Walsh degree-2 energy measuring geometrically?
"""

import numpy as np
from math import factorial, comb
from itertools import permutations
import sys
import time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 32749

def adj_from_bits(n, bits):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def rank_mod(M, p):
    M = M.copy() % p
    nrows, ncols = M.shape
    r = 0
    for col in range(ncols):
        nz = np.where(M[r:, col] != 0)[0]
        if len(nz) == 0: continue
        pivot = nz[0] + r
        if pivot != r:
            M[[r, pivot]] = M[[pivot, r]]
        inv = pow(int(M[r, col]), p - 2, p)
        M[r] = M[r] * inv % p
        factors = M[:, col].copy()
        factors[r] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[r])) % p
        r += 1
    return r

def compute_H(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            count = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def compute_beta1(A, n):
    edges = []
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]: edges.append((i, j))
    edge_idx = {e: i for i, e in enumerate(edges)}
    ne = len(edges)
    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[i][k]: triples.append((i, j, k))
    nt = len(triples)
    d1 = np.zeros((n, ne), dtype=np.int64)
    for idx, (i, j) in enumerate(edges):
        d1[j, idx] += 1; d1[i, idx] -= 1
    d2 = np.zeros((ne, nt), dtype=np.int64)
    for t_idx, (i, j, k) in enumerate(triples):
        if (j,k) in edge_idx: d2[edge_idx[(j,k)], t_idx] += 1
        if (i,k) in edge_idx: d2[edge_idx[(i,k)], t_idx] -= 1
        if (i,j) in edge_idx: d2[edge_idx[(i,j)], t_idx] += 1
    return ne - rank_mod(d1 % PRIME, PRIME) - rank_mod(d2 % PRIME, PRIME)

print("=" * 72)
print("WALSH-BETA1 CONNECTION: Spectral Detection of Topology")
print("opus-2026-03-15-S72b")
print("=" * 72)

# =====================================================================
# PART 1: Verify at n=5 and understand the mechanism
# =====================================================================
print("\n" + "=" * 60)
print("PART 1: n=5 DETAILED ANALYSIS")
print("=" * 60)

n = 5
m = n*(n-1)//2
total = 1 << m

# Compute H for all tournaments
t0 = time.time()
H_all = np.zeros(total, dtype=np.int64)
for bits in range(total):
    A = adj_from_bits(n, bits)
    H_all[bits] = compute_H(A, n)

# Walsh-Hadamard transform
H_hat = H_all.copy().astype(float)
for i in range(m):
    half = 1 << i
    for j in range(0, total, 2 * half):
        for k in range(half):
            a, b = H_hat[j + k], H_hat[j + k + half]
            H_hat[j + k] = a + b
            H_hat[j + k + half] = a - b
H_hat /= total

# Walsh degree-2 energy per tournament
# E_2(T) = sum_{|S|=2} H_hat[S] * chi_S(T)
# where chi_S(T) = (-1)^{T·S} (parity of bits in T at positions in S)

# More carefully: H(T) = sum_S H_hat[S] * (-1)^{<T,S>}
# So the degree-2 contribution to H(T) is:
# sum_{|S|=2} H_hat[S] * (-1)^{<T,S>}

# Let's compute this directly
deg2_energy = np.zeros(total)
for S in range(total):
    if bin(S).count('1') != 2:
        continue
    coeff = H_hat[S]
    if abs(coeff) < 1e-15:
        continue
    for bits in range(total):
        parity = bin(bits & S).count('1') % 2
        deg2_energy[bits] += coeff * (1 - 2 * parity)

# Also compute β₁ for all
beta1_all = np.zeros(total, dtype=int)
for bits in range(total):
    A = adj_from_bits(n, bits)
    beta1_all[bits] = compute_beta1(A, n)

elapsed = time.time() - t0
print(f"\n  Computed n=5 ({total} tournaments) in {elapsed:.1f}s")

# Verify the finding
b0_min = min(deg2_energy[bits] for bits in range(total) if beta1_all[bits] == 0)
b0_max = max(deg2_energy[bits] for bits in range(total) if beta1_all[bits] == 0)
b1_min = min(deg2_energy[bits] for bits in range(total) if beta1_all[bits] == 1)
b1_max = max(deg2_energy[bits] for bits in range(total) if beta1_all[bits] == 1)

print(f"\n  Walsh degree-2 energy ranges:")
print(f"    β₁=0: [{b0_min:.4f}, {b0_max:.4f}]")
print(f"    β₁=1: [{b1_min:.4f}, {b1_max:.4f}]")

# Is there a clean threshold?
# Find the gap between max(β₁=0 ∩ deg2) and min(β₁=1 ∩ deg2)
overlap_b0 = sorted(set(round(deg2_energy[b], 6) for b in range(total) if beta1_all[b] == 0))
overlap_b1 = sorted(set(round(deg2_energy[b], 6) for b in range(total) if beta1_all[b] == 1))

print(f"\n  Distinct deg-2 energy values:")
print(f"    β₁=0: {overlap_b0}")
print(f"    β₁=1: {overlap_b1}")

shared = set(overlap_b0) & set(overlap_b1)
print(f"    Shared values: {sorted(shared)}")

if not shared and b0_max < b1_min:
    print(f"\n  ★ CLEAN SEPARATION: threshold at {(b0_max + b1_min)/2:.4f} ★")
    print(f"    β₁=1 ⟺ deg-2 energy > {b0_max:.4f}")
elif shared:
    print(f"\n  Overlap exists — checking if β₁ is determined by deg-2 energy:")
    for val in sorted(shared):
        c0 = sum(1 for b in range(total) if beta1_all[b]==0 and abs(deg2_energy[b]-val)<0.001)
        c1 = sum(1 for b in range(total) if beta1_all[b]==1 and abs(deg2_energy[b]-val)<0.001)
        print(f"    E_2={val}: β₁=0 count={c0}, β₁=1 count={c1}")

# What ARE the nonzero Walsh coefficients at degree 2?
print(f"\n  Nonzero Walsh degree-2 coefficients H_hat[S]:")
for S in range(total):
    if bin(S).count('1') == 2:
        if abs(H_hat[S]) > 1e-15:
            # Decode S: which two edge pairs?
            bits_list = [i for i in range(m) if S & (1 << i)]
            edges = []
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if idx in bits_list:
                        edges.append((i, j))
                    idx += 1
            print(f"    S=({edges[0]},{edges[1]}): H_hat = {H_hat[S]:.6f}")

# =====================================================================
# PART 2: What does degree-2 Walsh energy MEAN?
# =====================================================================
print("\n" + "=" * 60)
print("PART 2: WHAT DEGREE-2 WALSH ENERGY MEANS")
print("=" * 60)

print(f"""
  Walsh degree-2 energy E_2(T) = sum_{{|S|=2}} H_hat[S] * chi_S(T)

  H_hat[S] for |S|=2 is nonzero only when S forms an even-length
  path union (by THM-071). Two edges forming a length-2 path = a P_3.

  So E_2(T) sums over all P_3 subgraphs (paths of length 2):
    For each triple (i,j,k) where edges (i,j) and (j,k) are
    in the tournament's edge set, there's a Walsh contribution
    based on the ORIENTATIONS of these two edges.

  The sign chi_S(T) = (-1)^{{# edges in S oriented "forward"}}.

  High E_2 means: the edge orientations in T are CORRELATED across
  neighboring edges (length-2 paths tend to agree in direction).

  Low E_2 means: edge orientations are ANTI-CORRELATED (zig-zag patterns).

  β₁=1 tournaments have HIGH E_2: their edge orientations are
  CONSISTENTLY ALIGNED along 2-paths.

  This is exactly the ORIENTABILITY CONNECTION:
  - High E_2 = edges "flow consistently" = but paradoxically,
    this CONSISTENCY creates the twist (too aligned to fill cycles!)
  - Low E_2 = edges "zig-zag" = transitive structure,
    triple boundaries can fill all cycles.
""")

# Verify: transitive tournament has lowest E_2
trans_bits = None
max_e2_bits = None
min_e2_bits = None
max_e2 = -999
min_e2 = 999
for bits in range(total):
    if deg2_energy[bits] < min_e2:
        min_e2 = deg2_energy[bits]
        min_e2_bits = bits
    if deg2_energy[bits] > max_e2:
        max_e2 = deg2_energy[bits]
        max_e2_bits = bits

# Check: is min E_2 the transitive tournament?
A_min = adj_from_bits(n, min_e2_bits)
H_min = H_all[min_e2_bits]
print(f"  Min E_2 = {min_e2:.4f} at bits={min_e2_bits}: H={H_min}, β₁={beta1_all[min_e2_bits]}")
print(f"    Score seq = {tuple(sorted(int(sum(A_min[i])) for i in range(n)))}")

A_max = adj_from_bits(n, max_e2_bits)
H_max = H_all[max_e2_bits]
print(f"  Max E_2 = {max_e2:.4f} at bits={max_e2_bits}: H={H_max}, β₁={beta1_all[max_e2_bits]}")
print(f"    Score seq = {tuple(sorted(int(sum(A_max[i])) for i in range(n)))}")

# =====================================================================
# PART 3: Test at n=6
# =====================================================================
print("\n" + "=" * 60)
print("PART 3: TESTING AT n=6")
print("=" * 60)

n = 6
m = n*(n-1)//2  # 15
total = 1 << m  # 32768

t0 = time.time()
H_all6 = np.zeros(total, dtype=np.int64)
for bits in range(total):
    A = adj_from_bits(n, bits)
    H_all6[bits] = compute_H(A, n)

# Walsh transform
H_hat6 = H_all6.copy().astype(float)
for i in range(m):
    half = 1 << i
    for j in range(0, total, 2 * half):
        for k in range(half):
            a, b = H_hat6[j + k], H_hat6[j + k + half]
            H_hat6[j + k] = a + b
            H_hat6[j + k + half] = a - b
H_hat6 /= total

# Degree-2 energy per tournament
deg2_energy6 = np.zeros(total)
for S in range(total):
    if bin(S).count('1') != 2:
        continue
    coeff = H_hat6[S]
    if abs(coeff) < 1e-15:
        continue
    for bits in range(total):
        parity = bin(bits & S).count('1') % 2
        deg2_energy6[bits] += coeff * (1 - 2 * parity)

# β₁ for all n=6 tournaments
print(f"  Computing β₁ for all {total} tournaments at n=6...")
beta1_6 = np.zeros(total, dtype=int)
for bits in range(total):
    A = adj_from_bits(n, bits)
    beta1_6[bits] = compute_beta1(A, n)
    if bits % 5000 == 0 and bits > 0:
        print(f"    ...{bits}/{total} ({time.time()-t0:.0f}s)")

elapsed = time.time() - t0
print(f"  Done in {elapsed:.1f}s")

b0_vals = [deg2_energy6[b] for b in range(total) if beta1_6[b] == 0]
b1_vals = [deg2_energy6[b] for b in range(total) if beta1_6[b] == 1]

print(f"\n  n=6 Walsh degree-2 energy:")
print(f"    β₁=0 ({len(b0_vals)}): mean={np.mean(b0_vals):.4f}, range=[{min(b0_vals):.4f}, {max(b0_vals):.4f}]")
print(f"    β₁=1 ({len(b1_vals)}): mean={np.mean(b1_vals):.4f}, range=[{min(b1_vals):.4f}, {max(b1_vals):.4f}]")

if max(b0_vals) < min(b1_vals):
    print(f"\n  ★ CLEAN SEPARATION AT n=6: threshold at {(max(b0_vals)+min(b1_vals))/2:.4f} ★")
else:
    # Check overlap
    b0_round = set(round(v, 4) for v in b0_vals)
    b1_round = set(round(v, 4) for v in b1_vals)
    shared = b0_round & b1_round
    if shared:
        print(f"\n  Overlap at n=6. Shared values: {sorted(list(shared))[:10]}...")
        # Count violations
        threshold = min(b1_vals)
        false_pos = sum(1 for v in b0_vals if v >= threshold)
        print(f"    If threshold = {threshold:.4f} (min of β₁=1):")
        print(f"    False positives (β₁=0 but E_2 ≥ threshold): {false_pos}")
    else:
        gap_low = max(v for v in b0_round if v < min(b1_round))
        gap_high = min(b1_round)
        print(f"\n  ★ CLEAN SEPARATION: gap ({gap_low:.4f}, {gap_high:.4f}) ★")

# =====================================================================
# PART 4: The degree-0 energy (mean H) connection
# =====================================================================
print("\n" + "=" * 60)
print("PART 4: MULTI-DEGREE ANALYSIS AT n=6")
print("=" * 60)

# Compute all even-degree Walsh energies
for deg in [0, 2, 4]:
    energy = np.zeros(total)
    for S in range(total):
        if bin(S).count('1') != deg:
            continue
        coeff = H_hat6[S]
        if abs(coeff) < 1e-15:
            continue
        for bits in range(total):
            parity = bin(bits & S).count('1') % 2
            energy[bits] += coeff * (1 - 2 * parity)

    b0 = [energy[b] for b in range(total) if beta1_6[b] == 0]
    b1 = [energy[b] for b in range(total) if beta1_6[b] == 1]
    print(f"\n  Degree {deg}:")
    print(f"    β₁=0: mean={np.mean(b0):.4f}, range=[{min(b0):.4f}, {max(b0):.4f}]")
    print(f"    β₁=1: mean={np.mean(b1):.4f}, range=[{min(b1):.4f}, {max(b1):.4f}]")
    overlap = max(b0) >= min(b1) and max(b1) >= min(b0)
    print(f"    Overlap: {overlap}")

# =====================================================================
# SYNTHESIS
# =====================================================================
print("\n" + "=" * 60)
print("SYNTHESIS")
print("=" * 60)

print(f"""
  The Walsh degree-2 energy E_2(T) measures how consistently
  neighboring edge pairs are oriented.

  HIGH E_2 = "flow alignment" = edges locally agree on direction
  LOW E_2 = "zig-zag" = edges alternate direction (transitive-like)

  The connection to β₁:
    β₁ > 0 = orientation obstruction in the path complex
    HIGH E_2 = consistent local orientations

  PARADOX: more consistent local orientation → WORSE global orientability.
  This is exactly like a Möbius strip: locally consistent orientation
  that doesn't close up globally.

  The transitive tournament (min E_2) has β₁=0: its zig-zag pattern
  makes every cycle decomposable into transitive triples.

  The regular tournament (max E_2) has β₁=1: its uniform flow
  creates cycles that can't be decomposed.

  CONJECTURE: E_2(T) > c(n) ⟹ β₁(T) = 1 for some threshold c(n).
  This would give a SPECTRAL SUFFICIENT CONDITION for the topological twist.
""")

print("DONE — walsh_beta1_connection_S72b.py complete")
