#!/usr/bin/env python3
"""
beta2_delta_ratio_n6.py - Verify delta_|A_3| = 2*delta_|A_2| at n=6

Key identity discovered at n=5:
  For any tournament T and any arc flip (u->v) -> (v->u):
  |A_3(T')| - |A_3(T)| = 2 * (|A_2(T')| - |A_2(T)|)

where A_p = allowed p-paths.

This is a COUNTING IDENTITY about allowed paths.

Also check if similar ratios hold for higher path dimensions.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import enumerate_allowed_paths
sys.stdout = _saved

def build_adj(n, bits):
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

# n=6 exhaustive
for n in [6]:
    print(f"{'='*70}")
    print(f"DELTA RATIO ANALYSIS AT n={n}")
    print(f"{'='*70}")

    n_arcs = n*(n-1)//2
    total = 1 << n_arcs
    t0 = time.time()

    # For efficiency, precompute |A_p| for all tournaments
    print(f"  Precomputing allowed path counts for {total} tournaments...")
    path_counts = {}  # bits -> {p: count}
    for bits in range(total):
        if bits % 10000 == 0 and bits > 0:
            dt = time.time() - t0
            print(f"    ... {bits}/{total} ({dt:.0f}s)")
        A = build_adj(n, bits)
        counts = {}
        for p in range(5):  # 0,1,2,3,4
            paths = enumerate_allowed_paths(A, n, p)
            counts[p] = len(paths)
            if not paths:
                break
        path_counts[bits] = counts

    print(f"  Precompute done. ({time.time()-t0:.0f}s)")

    # Now check ratios for all arc flips
    # For efficiency, only sample if too many
    violations_32 = 0
    violations_42 = 0
    ratio_32 = Counter()
    ratio_42 = Counter()
    total_flips = 0
    zero_d2 = 0

    for bits in range(total):
        if bits % 10000 == 0 and bits > 0 and n >= 6:
            print(f"    ... checking flips for {bits}/{total}")
        A = build_adj(n, bits)
        counts_before = path_counts[bits]

        for u in range(n):
            for v in range(n):
                if u == v or A[u][v] == 0:
                    continue

                # Flip
                B = [row[:] for row in A]
                B[u][v] = 0
                B[v][u] = 1
                bits_flip = 0
                idx = 0
                for i in range(n):
                    for j in range(i+1, n):
                        if B[i][j] == 1:
                            bits_flip |= (1 << idx)
                        idx += 1

                counts_after = path_counts[bits_flip]

                d2 = counts_after.get(2, 0) - counts_before.get(2, 0)
                d3 = counts_after.get(3, 0) - counts_before.get(3, 0)
                d4 = counts_after.get(4, 0) - counts_before.get(4, 0)

                total_flips += 1

                if d2 == 0:
                    zero_d2 += 1
                    if d3 != 0:
                        violations_32 += 1
                    if d4 != 0:
                        violations_42 += 1
                else:
                    r32 = d3 / d2
                    ratio_32[round(r32, 6)] += 1
                    if abs(r32 - 2.0) > 1e-6:
                        violations_32 += 1

                    r42 = d4 / d2
                    ratio_42[round(r42, 6)] += 1

    dt = time.time() - t0
    print(f"\n  Total flips: {total_flips}")
    print(f"  Time: {dt:.0f}s")
    print(f"  d2=0 cases: {zero_d2}")
    print(f"\n  d3/d2 ratio distribution:")
    for r in sorted(ratio_32.keys()):
        print(f"    {r}: {ratio_32[r]}")
    print(f"  d3=2*d2 violations: {violations_32}")

    print(f"\n  d4/d2 ratio distribution:")
    for r in sorted(ratio_42.keys()):
        print(f"    {r}: {ratio_42[r]}")
    if violations_42 > 0:
        print(f"  d4 != 0 when d2=0: {violations_42}")

print("\nDone.")
