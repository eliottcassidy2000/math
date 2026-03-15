#!/usr/bin/env python3
"""
wn_extend_89c.py — Extend W(n) computation to n=22,23,24
opus-2026-03-15-S89c

W(n) = Σ_{σ∈NUD(n)} 2^{adj1(σ)} where NUD = no unit descent.
adj1(σ) = #{i : σ(i+1) = σ(i)+1} (number of adjacencies/ascents by 1).

Uses bitmask DP on the Hamiltonian path counting interpretation:
dp[mask][v] = number of paths (Hamiltonian) through vertices in mask, ending at v,
weighted by 2^(adj1).

This is O(n²·2^n) time, O(n·2^n) space.
For n=22: ~2 billion ops. In pure Python this is too slow.
Use numpy arrays and bitwise tricks.

Actually, let's use a different formulation. We want:
W(n) = Σ_σ (1 if NUD(σ)) · 2^{adj1(σ)}

A permutation σ has a unit descent at position i if σ(i+1) = σ(i)-1.
adj1 counts positions where σ(i+1) = σ(i)+1 (unit ascent).

Equivalently, we look at consecutive pairs (σ(i), σ(i+1)):
- If σ(i+1) = σ(i)+1: unit ascent, contributes factor 2
- If σ(i+1) = σ(i)-1: unit descent, σ is excluded (NUD condition)
- Otherwise: neutral, contributes factor 1

So W(n) = Σ_σ ∏_i w(σ(i), σ(i+1))
where w(a,b) = 2 if b=a+1, 0 if b=a-1, 1 otherwise.

This is exactly a weighted Hamiltonian path count on K_n with edge weights w(a,b).
Bitmask DP: dp[mask][v] = sum of weighted paths through mask ending at v.

For n up to ~22, we can use arrays. Let me try to optimize with arrays.
"""

import sys
import time
from fractions import Fraction
from math import factorial

def compute_W(n, verbose=True):
    """Compute W(n) using bitmask DP. Returns integer."""
    if verbose:
        print(f"Computing W({n})... ", end="", flush=True)
    t0 = time.time()

    full = (1 << n) - 1
    # dp[v] for current mask (we process masks in order of popcount)
    # Use dict of {mask: array[n]} to save memory? No, need all masks.
    # For n=22: 2^22 = 4M masks × 22 vertices = 88M entries. As int64 array: 704 MB.
    # Too much for numpy. Use a dict approach: only keep masks of current popcount.

    # Actually for weighted paths, we process by popcount.
    # At popcount p, we have C(n,p) masks, each with p possible endpoints.
    # Total entries at popcount p: C(n,p) × n (but only p are nonzero per mask).

    # For n=22: max popcount ~11 has C(22,11) ≈ 700K masks. Storage: 700K × 22 ≈ 15M. OK.

    # Use dict: mask -> list of (v, weight) pairs, or mask -> {v: weight}
    # Better: use arrays indexed by mask, but only for one popcount at a time.

    # Actually, let me use a flat dict: (mask, v) -> weight. Prune zeros.
    # At popcount 1: n entries.
    # At popcount p: each entry branches to at most n-p new vertices.

    # Let's use dict of mask -> array.
    # Better yet: use two dicts (current and next popcount).

    # Initialize: paths of length 1 (single vertex)
    current = {}
    for v in range(n):
        mask = 1 << v
        if mask not in current:
            current[mask] = [0] * n
        current[mask][v] = 1

    for step in range(1, n):
        if verbose and step % 5 == 0:
            elapsed = time.time() - t0
            print(f"step {step}/{n-1} ({elapsed:.1f}s, {len(current)} masks)... ", end="", flush=True)

        next_level = {}
        for mask, endpoints in current.items():
            for v in range(n):
                wt = endpoints[v]
                if wt == 0:
                    continue
                # Extend path from v to u (not in mask)
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    # Edge weight
                    if u == v - 1:
                        # unit descent: weight 0 (exclude)
                        continue
                    elif u == v + 1:
                        # unit ascent: weight 2
                        edge_wt = 2
                    else:
                        edge_wt = 1

                    new_mask = mask | (1 << u)
                    if new_mask not in next_level:
                        next_level[new_mask] = [0] * n
                    next_level[new_mask][u] += wt * edge_wt

        current = next_level

    # Sum over all endpoints for the full mask
    W = sum(current.get(full, [0]))
    elapsed = time.time() - t0
    if verbose:
        print(f"done in {elapsed:.1f}s")
        print(f"W({n}) = {W}")
    return W

# Known values for verification
known_W = {
    3: 8, 4: 36, 5: 176, 6: 1080, 7: 8268, 8: 76104, 9: 822880,
    10: 10262160, 11: 145655148, 12: 2326258872, 13: 41383027792,
    14: 814116023280, 15: 17574474321708, 16: 413321769498648,
    17: 10545197020978752, 18: 290565166853498160,
}

# First verify a small case
W8 = compute_W(8)
print(f"W(8) check: computed={W8}, known={known_W[8]}, match={W8 == known_W[8]}")

# Now try n=14 as a timing benchmark
W14 = compute_W(14)
print(f"W(14) check: computed={W14}, known={known_W[14]}, match={W14 == known_W[14]}")

# Estimate time for n=22
# The time scales as n² × 2^n roughly.
# n=14: t14 seconds. n=22: t14 × (22/14)² × 2^8 ≈ t14 × 630. Probably too slow in Python.
print(f"\nEstimated time for n=22: too slow in pure Python.")
print("Need C implementation or numpy optimization.")
print("For now, try n=16,17 to see timing trend.")
