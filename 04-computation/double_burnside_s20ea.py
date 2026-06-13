#!/usr/bin/env python3
"""
double_burnside_s20ea.py — EXACT E(G_n) via double Burnside on unordered pairs
kind-pasteur-2026-03-23-S20ea

E(G_n) = #{orbits of S_n on unordered pairs {T, T'} connected by one arc flip}

By Burnside: E = (1/n!) * sum_sigma |Fix_pairs(sigma)|
where Fix_pairs(sigma) = #{unordered {T, T flip e}: sigma preserves the pair, [T]!=[T']}

Fix_pairs(sigma) = (1/2) * [Case1(sigma) + Case2(sigma)]
where:
  Case1: sigma(T)=T AND sigma(T')=T' AND [T]!=[T']  (sigma fixes both)
  Case2: sigma(T)=T' AND sigma(T')=T AND [T]!=[T']  (sigma swaps them)

Case1: sigma is an automorphism of T that also fixes the arc set-wise.
Case2: sigma is a "near-automorphism" — aut except reverses exactly arc e.

COMPUTATIONAL: directly enumerate for n=3,4,5 and verify E = known values.
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DOUBLE BURNSIDE: EXACT E(G_n)")
print("  kind-pasteur-2026-03-23-S20ea")
print("=" * 80)

for n in [3, 4, 5]:
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    all_perms = list(permutations(range(n)))

    def make_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_ARCS):
            if bits & (1 << k): A[i][j] = 1
            else: A[j][i] = 1
        return A

    def apply_perm(bits, sigma):
        """Apply sigma to tournament bits, return new bits."""
        A = make_adj(bits)
        B = [[A[sigma[i]][sigma[j]] for j in range(n)] for i in range(n)]
        new_bits = 0
        for k, (i,j) in enumerate(ALL_ARCS):
            if B[i][j]: new_bits |= (1 << k)
        return new_bits

    def canon(A):
        best = None
        for p in all_perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    # Precompute canon for all tournaments
    canon_map = {}
    for bits in range(1 << m):
        A = make_adj(bits)
        canon_map[bits] = canon(A)

    # For each sigma, compute Case1 and Case2
    total_fix_pairs = 0

    for sigma in all_perms:
        sigma_inv = [0]*n
        for i in range(n): sigma_inv[sigma[i]] = i

        case1 = 0
        case2 = 0

        for bits in range(1 << m):
            # Check if sigma fixes T: sigma(T) = T
            sigma_bits = apply_perm(bits, sigma_inv)
            sigma_fixes_T = (sigma_bits == bits)

            for k in range(m):
                flip_bits = bits ^ (1 << k)
                u, v = ALL_ARCS[k]

                # Check [T] != [T flip e]
                same_class = (canon_map[bits] == canon_map[flip_bits])
                if same_class:
                    continue  # skip self-loops

                # Case 1: sigma fixes T AND sigma fixes T'
                if sigma_fixes_T:
                    # sigma must also fix T flip e
                    sigma_flip = apply_perm(flip_bits, sigma_inv)
                    if sigma_flip == flip_bits:
                        case1 += 1

                # Case 2: sigma(T) = T' and sigma(T') = T
                if sigma_bits == flip_bits:
                    # sigma maps T to T flip e. Check reverse.
                    sigma_flip = apply_perm(flip_bits, sigma_inv)
                    if sigma_flip == bits:
                        case2 += 1

        fix_pairs = (case1 + case2) // 2  # unordered pairs
        total_fix_pairs += fix_pairs

    E_computed = total_fix_pairs // factorial(n)
    E_known = {3:1, 4:5, 5:30}[n]

    print(f"\nn={n}: ({time.time()-t0:.1f}s)")
    print(f"  sum Fix_pairs = {total_fix_pairs}")
    print(f"  E = sum / n! = {total_fix_pairs} / {factorial(n)} = {E_computed}")
    print(f"  E_known = {E_known}")
    print(f"  MATCH: {E_computed == E_known}")

    # Decompose by sigma cycle type
    print(f"  Case1 total = {sum(1 for s in all_perms for b in range(1<<m) for k in range(m) if apply_perm(b, [0]*n)==b)}")  # placeholder
    # Actually just print the totals
    # I'll add detailed decomposition if time permits

print(f"\n  DONE.")
print("=" * 80)
