#!/usr/bin/env python3
"""
beta2_delta_ratio_general.py - Verify delta_|A_3| = (n-3)*delta_|A_2| for general n

Discovered pattern:
  n=5: delta_|A_3| = 2*delta_|A_2|  (2 = n-3)
  n=6: delta_|A_3| = 3*delta_|A_2|  (3 = n-3)

Conjecture: delta_|A_3| = (n-3)*delta_|A_2| for all tournaments, all arcs.

Also check: what is the general formula delta_|A_p| = f(n,p) * delta_|A_2|?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import enumerate_allowed_paths
sys.stdout = _saved

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_allowed(A, n, max_p=5):
    counts = {}
    for p in range(max_p + 1):
        paths = enumerate_allowed_paths(A, n, p)
        counts[p] = len(paths)
        if not paths:
            break
    return counts

N_SAMPLES = 200

for n in [4, 5, 6, 7, 8]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    max_p = min(n, 6)
    ratios_p = {p: Counter() for p in range(3, max_p + 1)}
    zero_d2 = 0
    total_flips = 0
    violations = {p: 0 for p in range(3, max_p + 1)}

    t0 = time.time()

    for trial in range(N_SAMPLES):
        A = random_tournament(n)
        counts_before = count_allowed(A, n, max_p)

        # Try a few random arc flips
        for _ in range(5):
            u = random.randint(0, n-1)
            v = random.randint(0, n-1)
            while u == v or A[u][v] == 0:
                u = random.randint(0, n-1)
                v = random.randint(0, n-1)

            B = [row[:] for row in A]
            B[u][v] = 0
            B[v][u] = 1

            counts_after = count_allowed(B, n, max_p)
            total_flips += 1

            d2 = counts_after.get(2, 0) - counts_before.get(2, 0)

            if d2 == 0:
                zero_d2 += 1
                for p in range(3, max_p + 1):
                    dp = counts_after.get(p, 0) - counts_before.get(p, 0)
                    if dp != 0:
                        violations[p] += 1
            else:
                for p in range(3, max_p + 1):
                    dp = counts_after.get(p, 0) - counts_before.get(p, 0)
                    r = dp / d2
                    ratios_p[p][round(r, 6)] += 1

    dt = time.time() - t0
    print(f"  Total flips: {total_flips}, d2=0: {zero_d2}, Time: {dt:.1f}s")

    for p in range(3, max_p + 1):
        print(f"\n  d{p}/d2 ratio distribution:")
        for r in sorted(ratios_p[p].keys()):
            if ratios_p[p][r] >= 3:  # Only show non-rare values
                print(f"    {r}: {ratios_p[p][r]}")
        if violations[p] > 0:
            print(f"    violations (dp!=0 when d2=0): {violations[p]}")
        # What's the dominant ratio?
        if ratios_p[p]:
            dominant = max(ratios_p[p], key=ratios_p[p].get)
            print(f"    DOMINANT ratio = {dominant} (expected n-3 = {n-3})")

print("\nDone.")
