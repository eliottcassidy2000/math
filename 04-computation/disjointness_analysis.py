#!/usr/bin/env python3
"""
disjointness_analysis.py -- Why does Interval have more disjoint cycle pairs?

At p=7: Paley has 80 cycles, 7 disjoint pairs (0.22% of all pairs)
        Interval has 59 cycles, 14 disjoint pairs (0.82% of all pairs)

At p=11: Paley has 21169 cycles, 10879 disjoint pairs (0.00487% of all pairs)
         Interval has 18397 cycles, 11110 disjoint pairs (0.00657% of all pairs)

The DISJOINTNESS RATE (alpha_2 / C(alpha_1, 2)) is always higher for Interval.
Why? This script investigates the structural reason.

Hypothesis: Interval cycles cluster in "local neighborhoods" on Z_p,
leaving space for non-overlapping cycles elsewhere. Paley cycles are
more "spread out" (quasi-random property distributes them evenly),
creating more vertex conflicts.

Author: kind-pasteur-2026-03-12-S58
"""

import time
from itertools import combinations
import sys
sys.path.insert(0, '.')
from alpha_decomp_p19 import build_adj, count_ham_cycles_k3, count_ham_cycles_k5, count_ham_cycles_on_subset


def enumerate_all_cycles(A, p, max_k=None):
    """Enumerate all directed odd cycles, returning (frozenset, length) pairs."""
    if max_k is None:
        max_k = p
    cycles = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(p), k):
            if k == 3:
                n_cyc = count_ham_cycles_k3(A, p, subset)
            elif k == 5:
                n_cyc = count_ham_cycles_k5(A, subset)
            else:
                n_cyc = count_ham_cycles_on_subset(A, subset)
            for _ in range(n_cyc):
                cycles.append((frozenset(subset), k))
    return cycles


def vertex_usage_profile(cycles, p):
    """How many cycles use each vertex?"""
    count = [0] * p
    for fs, k in cycles:
        for v in fs:
            count[v] += 1
    return count


def pair_overlap_distribution(cycles):
    """Distribution of |C_i intersect C_j| over all pairs."""
    n = len(cycles)
    overlap_dist = {}
    for i in range(n):
        for j in range(i+1, n):
            ov = len(cycles[i][0] & cycles[j][0])
            overlap_dist[ov] = overlap_dist.get(ov, 0) + 1
    return overlap_dist


def length_pair_disjointness(cycles):
    """Disjointness rate by (length_i, length_j) pair type."""
    n = len(cycles)
    pair_count = {}
    disjoint_count = {}

    for i in range(n):
        for j in range(i+1, n):
            k1 = min(cycles[i][1], cycles[j][1])
            k2 = max(cycles[i][1], cycles[j][1])
            key = (k1, k2)
            pair_count[key] = pair_count.get(key, 0) + 1
            if not (cycles[i][0] & cycles[j][0]):
                disjoint_count[key] = disjoint_count.get(key, 0) + 1

    return pair_count, disjoint_count


def cycle_spread(cycles, p):
    """Measure how "spread out" cycles are on Z_p.

    For each cycle C, compute its "diameter" = max(C) - min(C) mod p.
    Also compute the "center of mass" on the circle.
    """
    import cmath
    diameters = []
    for fs, k in cycles:
        verts = sorted(fs)
        # Diameter on circular Z_p
        gaps = [(verts[(i+1) % k] - verts[i]) % p for i in range(k)]
        # Min gap tells us how "spread" the cycle is
        max_gap = max(gaps)
        diameter = p - max_gap  # arc length of the smallest containing arc
        diameters.append(diameter)
    return diameters


def main():
    print("=" * 70)
    print("DISJOINTNESS ANALYSIS: WHY INTERVAL HAS MORE DISJOINT PAIRS")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            A = build_adj(p, S)
            cycles = enumerate_all_cycles(A, p)
            n = len(cycles)

            print(f"\n{'='*60}")
            print(f"p={p}, {name}: {n} cycles")
            print(f"{'='*60}")

            # 1. Vertex usage profile
            usage = vertex_usage_profile(cycles, p)
            avg_usage = sum(usage) / p
            var_usage = sum((u - avg_usage)**2 for u in usage) / p
            print(f"\n  Vertex usage: avg={avg_usage:.1f}, var={var_usage:.1f}, std={var_usage**0.5:.1f}")
            print(f"  Per vertex: {usage}")

            # For circulant tournaments, usage should be uniform
            # (each vertex appears in the same number of cycles)
            if len(set(usage)) == 1:
                print(f"  UNIFORM: all vertices in {usage[0]} cycles")

            # 2. Overlap distribution
            overlap_dist = pair_overlap_distribution(cycles)
            total_pairs = n * (n - 1) // 2
            disjoint = overlap_dist.get(0, 0)
            rate = disjoint / total_pairs if total_pairs > 0 else 0
            print(f"\n  Overlap distribution (|C_i cap C_j|):")
            for ov in sorted(overlap_dist.keys()):
                pct = 100 * overlap_dist[ov] / total_pairs
                print(f"    overlap={ov}: {overlap_dist[ov]:>8d} pairs ({pct:.2f}%)")
            print(f"  Disjointness rate: {disjoint}/{total_pairs} = {100*rate:.4f}%")

            # 3. Disjointness by cycle length pair
            pair_ct, disj_ct = length_pair_disjointness(cycles)
            print(f"\n  Disjointness by (len_i, len_j):")
            for key in sorted(pair_ct.keys()):
                total = pair_ct[key]
                disj = disj_ct.get(key, 0)
                rate_k = disj / total if total > 0 else 0
                print(f"    ({key[0]},{key[1]}): {disj}/{total} = {100*rate_k:.2f}%")

            # 4. Cycle spread
            diameters = cycle_spread(cycles, p)
            avg_diam = sum(diameters) / len(diameters) if diameters else 0
            print(f"\n  Cycle diameter (arc length): avg={avg_diam:.2f}")
            diam_dist = {}
            for d in diameters:
                diam_dist[d] = diam_dist.get(d, 0) + 1
            for d in sorted(diam_dist.keys()):
                print(f"    diameter={d}: {diam_dist[d]} cycles")


if __name__ == '__main__':
    main()
