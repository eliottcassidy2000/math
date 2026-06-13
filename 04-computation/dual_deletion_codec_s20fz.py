#!/usr/bin/env python3
"""
dual_deletion_codec_s20fz.py — Mode B dual deletion for fractal compression
kind-pasteur-2026-03-24-S20fz

Mode A: delete ONE vertex (n-1 residual bits). Score-based saves 35-47%.
Mode B: delete source AND sink (2n-5 residual bits). How much savings?

The dual deletion removes both endpoints of the base path simultaneously.
The residual encodes: bottom wiring (n-2 bits) + top wiring (n-2 bits) + apex (1 bit).
The overlap tournament on {1,...,n-2} is the base.

Question: is the dual residual MORE compressible than two single residuals?
If source and sink wirings are CORRELATED given the overlap class,
then dual deletion captures this correlation.
"""

import sys
from math import factorial, comb, log2
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DUAL DELETION CODEC (Mode B)")
print("  kind-pasteur-2026-03-24-S20fz")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6]:
    t0 = time.time()
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms_n = list(permutations(range(N)))
    all_perms_inner = list(permutations(range(1, N-1)))  # inner vertices

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

    def canonicalize(A, perms):
        nn = len(A)
        best = None
        for p in perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(nn) for j in range(nn))
            if best is None or s < best: best = s
        return best

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}")
    print(f"  Overlap: vertices {{1,...,{n-2}}}, {comb(n-2,2)} arcs")
    print(f"  Dual residual: bottom({n-2}) + top({n-2}) + apex(1) = {2*(n-2)+1} bits")
    print(f"{'#'*60}")

    # Build all tournaments
    canon_n = {}
    adj_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        adj_map[mask] = A
        canon_n[mask] = canonicalize(A, all_perms_n)

    # Overlap tournament: sub-tournament on vertices {1,...,n-2}
    # In adjacency matrix: rows/cols 1..n-2
    def overlap_tournament(A):
        inner = list(range(1, N-1))
        nn = len(inner)
        sub = [[A[inner[i]][inner[j]] for j in range(nn)] for i in range(nn)]
        return sub

    def overlap_canon(A):
        sub = overlap_tournament(A)
        # Canonicalize with perms on inner vertices
        nn = len(sub)
        perms = list(permutations(range(nn)))
        return canonicalize(sub, perms)

    # Dual residual: (bottom_wiring, top_wiring, apex)
    # Bottom: A[0][1], A[0][2], ..., A[0][n-2] (arcs from vertex 0 to inner)
    # Top: A[n-1][1], A[n-1][2], ..., A[n-1][n-2] (arcs from vertex n-1 to inner)
    # Apex: A[0][n-1] (arc from vertex 0 to vertex n-1)
    def dual_residual(A):
        bottom = tuple(A[0][j] for j in range(1, N-1))
        top = tuple(A[N-1][j] for j in range(1, N-1))
        apex = A[0][N-1]
        return (bottom, top, apex)

    # Single residual (Mode A, delete vertex 0)
    def single_residual_v0(A):
        return tuple(A[0][j] for j in range(1, N))

    # Score of vertex 0 and vertex n-1
    def vertex_scores(A):
        s0 = sum(A[0][j] for j in range(N) if j != 0)
        sn = sum(A[N-1][j] for j in range(N) if j != N-1)
        return s0, sn

    print(f"  Computing...", end=" ", flush=True)

    # For each overlap class: what dual residuals appear?
    overlap_to_dual = defaultdict(Counter)
    overlap_to_single = defaultdict(Counter)
    overlap_to_scores = defaultdict(Counter)

    for mask in range(1 << m):
        A = adj_map[mask]
        oc = overlap_canon(A)
        dr = dual_residual(A)
        sr = single_residual_v0(A)
        s0, sn = vertex_scores(A)

        overlap_to_dual[oc][dr] += 1
        overlap_to_single[oc][sr] += 1
        overlap_to_scores[oc][(s0, sn)] += 1

    total = 2 ** m

    # Dual residual conditional entropy
    dual_naive = 2*(n-2) + 1  # total dual residual bits
    single_naive = n - 1

    dual_cond_ent = 0
    single_cond_ent = 0
    score_cond_ent = 0

    for oc in overlap_to_dual:
        p_class = sum(overlap_to_dual[oc].values()) / total

        # Dual entropy
        dual_total = sum(overlap_to_dual[oc].values())
        ent_d = sum(-c/dual_total * log2(c/dual_total) for c in overlap_to_dual[oc].values() if c > 0)
        dual_cond_ent += p_class * ent_d

        # Single entropy
        single_total = sum(overlap_to_single[oc].values())
        ent_s = sum(-c/single_total * log2(c/single_total) for c in overlap_to_single[oc].values() if c > 0)
        single_cond_ent += p_class * ent_s

        # Score-pair entropy
        score_total = sum(overlap_to_scores[oc].values())
        ent_sc = sum(-c/score_total * log2(c/score_total) for c in overlap_to_scores[oc].values() if c > 0)
        score_cond_ent += p_class * ent_sc

    print(f"done ({time.time()-t0:.1f}s)")

    print(f"\n  MODE A (single vertex deletion):")
    print(f"    Naive: {single_naive} bits")
    print(f"    H(single_resid | overlap_class): {single_cond_ent:.4f} bits")
    print(f"    Savings: {single_naive - single_cond_ent:.4f} bits ({(1 - single_cond_ent/single_naive)*100:.1f}%)")

    print(f"\n  MODE B (dual source-sink deletion):")
    print(f"    Naive: {dual_naive} bits")
    print(f"    H(dual_resid | overlap_class): {dual_cond_ent:.4f} bits")
    print(f"    Savings: {dual_naive - dual_cond_ent:.4f} bits ({(1 - dual_cond_ent/dual_naive)*100:.1f}%)")

    print(f"\n  SCORE-PAIR entropy:")
    print(f"    H(score_pair | overlap_class): {score_cond_ent:.4f} bits")

    # Bits per arc comparison
    mode_a_per_arc = single_cond_ent / single_naive
    mode_b_per_arc = dual_cond_ent / dual_naive
    print(f"\n  EFFICIENCY (conditional bits / naive bits):")
    print(f"    Mode A: {mode_a_per_arc:.4f} bits/bit")
    print(f"    Mode B: {mode_b_per_arc:.4f} bits/bit")
    print(f"    Mode B {'BETTER' if mode_b_per_arc < mode_a_per_arc else 'WORSE'} than Mode A")

    # Distinct residuals per overlap class
    dual_counts = [len(resids) for resids in overlap_to_dual.values()]
    single_counts = [len(resids) for resids in overlap_to_single.values()]
    print(f"\n  DISTINCT RESIDUALS per overlap class:")
    print(f"    Dual:   min={min(dual_counts)}, max={max(dual_counts)}, avg={sum(dual_counts)/len(dual_counts):.1f} / {2**dual_naive}")
    print(f"    Single: min={min(single_counts)}, max={max(single_counts)}, avg={sum(single_counts)/len(single_counts):.1f} / {2**single_naive}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
