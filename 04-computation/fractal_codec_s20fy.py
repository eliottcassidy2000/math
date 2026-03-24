#!/usr/bin/env python3
"""
fractal_codec_s20fy.py — The fractal tournament codec
kind-pasteur-2026-03-24-S20fy

THE BIG IDEA: A tournament on n vertices = a tournament on n-1 vertices
+ n-1 new bits (how the new vertex connects to the existing ones).

This gives a RECURSIVE encoding:
  T(n) = T(n-1) + residual(n-1 bits)
  T(n-1) = T(n-2) + residual(n-2 bits)
  ...
  T(3) = T(2) + residual(1 bit)
  T(2) = 1 bit

Total bits: 1 + 2 + 3 + ... + (n-1) = C(n,2). No savings yet!

BUT: the residual is NOT random — it's PREDICTABLE from the base.
If we know the iso class of T(n-1), the residual has only k possible
values (where k << 2^{n-1}). This is the CONDITIONAL entropy:

  H(residual | class(T(n-1))) < n-1 bits

The savings come from the PREDICTABILITY of the residual given the base class.

MEASURE THIS: for each iso class C at n-1, how many distinct extensions
to n exist? The "extension entropy" H(ext | C) measures the compression.

Also: the TILING MODEL already does this! The base path gives
T(n-1) from T(n) by deleting vertex n, and the n-1 "new" tiles
(those involving vertex n) are the residual. But the tiling model
stores ALL tiles, not just the residual. The fractal codec stores
ONLY the residual at each level.
"""

import sys
from math import factorial, comb, log2
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  FRACTAL TOURNAMENT CODEC")
print("  kind-pasteur-2026-03-24-S20fy")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6, 7]:
    t0 = time.time()
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms_n = list(permutations(range(N)))
    all_perms_nm1 = list(permutations(range(N-1)))

    # Build all tournaments on n vertices (as adjacency matrices)
    # and compute their class + the class of the sub-tournament on {0,...,n-2}
    def canonicalize_n(A):
        best = None
        for p in all_perms_n:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def canonicalize_nm1(A):
        nn = N - 1
        best = None
        for p in all_perms_nm1:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(nn) for j in range(nn))
            if best is None or s < best: best = s
        return best

    # Use tiling model for enumeration
    TILES = get_tiles(n)
    m = len(TILES)
    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]

    def b2a(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tv[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}")
    print(f"{'#'*60}")

    # For each tournament: compute its n-class and the (n-1)-class of the
    # sub-tournament obtained by deleting vertex 0 (the "source" in base path)
    print(f"  Building classes...", end=" ", flush=True)

    # Class of full tournament
    canon_n = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_n[mask] = canonicalize_n(A)

    # Sub-tournament: delete vertex 0 (index 0 in adjacency)
    # This gives an (n-1)-tournament on vertices {1,...,n-1}
    def sub_tournament(A):
        nn = N - 1
        return [[A[i+1][j+1] for j in range(nn)] for i in range(nn)]

    # The "residual": how vertex 0 connects to {1,...,n-1}
    # = the first row of A (excluding diagonal): A[0][1], A[0][2], ..., A[0][n-1]
    def residual(A):
        return tuple(A[0][j] for j in range(1, N))

    canon_sub = {}
    resid_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        sub = sub_tournament(A)
        canon_sub[mask] = canonicalize_nm1(sub)
        resid_map[mask] = residual(A)

    classes_n = sorted(set(canon_n.values()))
    classes_nm1 = sorted(set(canon_sub.values()))
    V_n = len(classes_n)
    V_nm1 = len(classes_nm1)

    print(f"V(n)={V_n}, V(n-1)={V_nm1}")

    # ================================================================
    # CONDITIONAL ENTROPY: H(residual | class(n-1))
    # ================================================================
    # For each (n-1)-class: what residuals appear?
    class_to_residuals = defaultdict(Counter)
    for mask in range(1 << m):
        cn_sub = canon_sub[mask]
        r = resid_map[mask]
        class_to_residuals[cn_sub][r] += 1

    total_tilings = 2 ** m
    cond_entropy = 0
    for cn_sub, resid_counter in class_to_residuals.items():
        total_in_class = sum(resid_counter.values())
        p_class = total_in_class / total_tilings
        # Entropy of residual given this class
        class_entropy = 0
        for r, count in resid_counter.items():
            p = count / total_in_class
            if p > 0:
                class_entropy -= p * log2(p)
        cond_entropy += p_class * class_entropy

    naive_bits = N - 1  # n-1 bits for the residual
    print(f"\n  CONDITIONAL ENTROPY:")
    print(f"    Naive residual: {naive_bits} bits")
    print(f"    H(residual | class(n-1)): {cond_entropy:.4f} bits")
    print(f"    Savings: {naive_bits - cond_entropy:.4f} bits = {(1 - cond_entropy/naive_bits)*100:.1f}%")

    # How many distinct residuals per (n-1)-class?
    resid_counts = [len(resids) for resids in class_to_residuals.values()]
    print(f"    Distinct residuals per class: min={min(resid_counts)}, max={max(resid_counts)}, avg={sum(resid_counts)/len(resid_counts):.1f}")
    print(f"    Max possible: 2^{naive_bits} = {2**naive_bits}")
    print(f"    Residual compression: 2^{naive_bits}/{max(resid_counts)} = {2**naive_bits/max(resid_counts):.1f}x for best class")

    # ================================================================
    # FULL FRACTAL CODEC: recursive conditional entropy
    # ================================================================
    # Total bits = sum over levels k=2..n of H(residual_k | class(k-1))
    # Compare to naive = sum k-1 for k=2..n = C(n,2)

    print(f"\n  FRACTAL COMPRESSION (recursive):")
    print(f"    Naive total: C({n},2) = {comb(n,2)} bits")
    print(f"    Level {n}: H(resid|class({n-1})) = {cond_entropy:.4f} bits")
    # Estimate lower levels from the data we have
    # At each level k: residual is k-1 bits naive, conditional entropy is lower
    # For simplicity: estimate total as sum of conditional entropies

    # Actually compute for all delete-vertex-0 chains
    # But we only have the n -> n-1 step. Extrapolate.
    estimated_total = cond_entropy
    for k in range(n-1, 2, -1):
        # Rough estimate: conditional entropy at level k ~ (k-1) * (1 - savings_fraction)
        # Use the same savings fraction as level n
        frac = cond_entropy / naive_bits
        estimated_total += (k-1) * frac
    estimated_total += 1  # Level 2: 1 bit (no savings possible)

    print(f"    Estimated total (fractal): {estimated_total:.2f} bits")
    print(f"    Fractal compression ratio: {comb(n,2)/estimated_total:.2f}x")
    print(f"    vs Level 3 (class): {comb(n,2)/log2(max(V_n,2)):.2f}x")

    # ================================================================
    # PER-CLASS DETAIL
    # ================================================================
    print(f"\n  TOP (n-1)-CLASSES by extension count:")
    for cn_sub in sorted(class_to_residuals.keys(),
                         key=lambda c: -len(class_to_residuals[c]))[:5]:
        resids = class_to_residuals[cn_sub]
        total = sum(resids.values())
        n_distinct = len(resids)
        ent = sum(-count/total * log2(count/total) for count in resids.values() if count > 0)
        print(f"    class: {n_distinct} extensions ({total} tilings), entropy={ent:.2f} bits (naive={naive_bits})")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
