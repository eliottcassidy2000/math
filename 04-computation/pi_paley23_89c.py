#!/usr/bin/env python3
"""
pi_paley23_89c.py — Compute H(P_23) using optimized bitmask DP
opus-2026-03-14-S89c

P_23 has 23 vertices, so 2^23 = 8,388,608 states.
Each state needs at most 23 endpoint counts.
Using arrays instead of dicts for speed.

Memory estimate: 8M × 23 × 8 bytes = ~1.5 GB for int64 arrays.
That's borderline. Use a rolling approach: process masks by popcount.
"""

import time
import numpy as np

def paley_tournament(p):
    """Build adjacency for Paley tournament on Z/pZ."""
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    adj = [[] for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                adj[i].append(j)
    return adj, qr

def count_ham_paths_optimized(adj, n):
    """
    Count Hamiltonian paths using bitmask DP with dict-of-dicts.
    Only keep current popcount layer + next layer.
    """
    # dp[mask] = {endpoint: count}
    # Process by popcount to save memory

    # For n=23, even with rolling, we need masks with popcount ~12
    # which is C(23,12) ≈ 1.3M masks, each with up to 12 endpoints.
    # At popcount 11-12 is the peak: C(23,11) = 1352078, C(23,12) = 1352078

    # Use numpy arrays indexed by mask for each popcount layer
    # Actually let's just use plain Python dicts for flexibility

    print(f"  Starting DP for n={n}...")
    t0 = time.time()

    # Layer 0: single vertices
    current = {}
    for v in range(n):
        mask = 1 << v
        current[mask] = {v: 1}

    for step in range(1, n):
        t1 = time.time()
        next_layer = {}

        for mask, endpoints in current.items():
            for v, count in endpoints.items():
                for u in adj[v]:
                    if mask & (1 << u):
                        continue
                    new_mask = mask | (1 << u)
                    if new_mask not in next_layer:
                        next_layer[new_mask] = {}
                    if u in next_layer[new_mask]:
                        next_layer[new_mask][u] += count
                    else:
                        next_layer[new_mask][u] = count

        t2 = time.time()
        num_masks = len(next_layer)
        total_entries = sum(len(v) for v in next_layer.values())
        print(f"    Step {step}/{n-1}: {num_masks} masks, {total_entries} entries, {t2-t1:.1f}s")

        current = next_layer

    # Sum over full mask
    full = (1 << n) - 1
    if full in current:
        H = sum(current[full].values())
    else:
        H = 0

    total_time = time.time() - t0
    print(f"  Total time: {total_time:.1f}s")
    return H

# First verify on known cases
print("=" * 70)
print("Verification on known cases")
print("=" * 70)

for p in [7, 11]:
    adj, qr = paley_tournament(p)
    H = count_ham_paths_optimized(adj, p)
    print(f"  H(P_{p}) = {H}")

print()
print("=" * 70)
print("Computing H(P_19)")
print("=" * 70)

adj19, qr19 = paley_tournament(19)
H19 = count_ham_paths_optimized(adj19, 19)
print(f"\n  H(P_19) = {H19}")
print(f"  Expected: 1172695746915")
print(f"  Match: {'✓' if H19 == 1172695746915 else '✗'}")

print()
print("=" * 70)
print("Computing H(P_23) — this is the big one")
print("=" * 70)

adj23, qr23 = paley_tournament(23)
print(f"  QR mod 23: {sorted(qr23)}")
print(f"  Number of QR: {len(qr23)} (should be 11)")

H23 = count_ham_paths_optimized(adj23, 23)
print(f"\n  H(P_23) = {H23}")
print(f"  H mod 23 = {H23 % 23}")
print(f"  H/23 = {H23 // 23}")
print(f"  H/23 mod 23 = {(H23 // 23) % 23}")

import sympy
print(f"  Factorization: {sympy.factorint(H23)}")

# Also compute Hamiltonian cycles for P_23
# Fix vertex 0, count paths from 0 back to 0
print(f"\n  Computing directed Hamiltonian cycles for P_23...")

# Modified DP: start only at vertex 0
current_cyc = {}
mask0 = 1
current_cyc[mask0] = {0: 1}

t0 = time.time()
for step in range(1, 23):
    next_layer = {}
    for mask, endpoints in current_cyc.items():
        for v, count in endpoints.items():
            for u in adj23[v]:
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                if new_mask not in next_layer:
                    next_layer[new_mask] = {}
                if u in next_layer[new_mask]:
                    next_layer[new_mask][u] += count
                else:
                    next_layer[new_mask][u] = count
    current_cyc = next_layer
    if step % 5 == 0:
        print(f"    Cycle step {step}/22: {len(current_cyc)} masks")

full23 = (1 << 23) - 1
hc23 = 0
if full23 in current_cyc:
    for v, count in current_cyc[full23].items():
        if 0 in adj23[v]:  # edge back to 0
            hc23 += count

# Multiply by number of starting vertices (23) then divide by cycle length (23)
# Actually: fixing start at 0 counts each cycle exactly once (the cycle visits 0 once)
# No wait: a directed cycle on p vertices starting at 0 has one unique representation.
# So hc23 = number of directed Ham cycles.
hc23_total = hc23  # This counts cycles starting at 0. Multiply by p and divide by p (cancel) = hc23.
# Actually this is wrong. Let me think again.
# A directed Ham cycle visits all p vertices. If we fix the starting vertex as 0,
# we get each cycle exactly once (since a directed cycle has p rotations, and exactly
# one rotation starts at 0). Wait no — a directed Hamiltonian cycle has p distinct
# rotations, each starting at a different vertex. So fixing start at 0 gives
# total_directed_cycles / 1 (each cycle counted once starting at 0).
# Hmm, no. For the cycle 0→1→2→...→22→0, fixing start at 0 gives us exactly the
# path 0→1→2→...→22 that ends at 22 with edge to 0.
# For the cycle 1→2→...→22→0→1, fixing start at 0 gives us 0→1→... wait, starting
# at 0 we'd have 0→{something that goes to 1 eventually}...
# Actually: the way I computed it, I start at vertex 0 and find all Hamiltonian paths
# from 0 to any vertex v where v→0 exists. Each such path corresponds to a unique
# directed Hamiltonian cycle. So hc23 IS the total number of directed Ham cycles.

print(f"\n  Directed Ham cycles in P_23: {hc23}")
print(f"  hc mod 23 = {hc23 % 23}")
print(f"  Expected: (p-1)/2 = 11")
print(f"  Match: {'✓' if hc23 % 23 == 11 else '✗'}")

t1 = time.time()
print(f"  Cycle computation time: {t1-t0:.1f}s")

print("\nDone!")
