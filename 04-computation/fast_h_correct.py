"""
fast_h_correct.py — opus-2026-04-04-S16

CORRECTED fast H computation via OCF truncation.

For n ≤ 8: H = 1 + 2*alpha_1 + 4*alpha_2 EXACTLY.
alpha_1 = number of directed odd cycles (3-cycles + 5-cycles + 7-cycles)
alpha_2 = number of vertex-disjoint pairs

KEY FIX: At n=7, 7-cycles exist. But max independent set = 2 (since
3 disjoint 3-cycles need 9 > 7 vertices). So alpha_3 = 0 always for n ≤ 8.

For the CYCLE COUNT: we need ALL directed odd cycles.
- 3-cycles: 2 per vertex triple that forms a cycle
- 5-cycles: many per vertex quintuple
- 7-cycles: many per vertex septuple (only at n=7)

The 5-cycle count needs all cyclic orderings checked.
"""
import numpy as np
from itertools import combinations, permutations
from math import factorial
import time
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import make_tiles, tiling_to_tournament, count_hamiltonian_paths

def count_directed_cycles(T, n, k):
    """Count all directed k-cycles in tournament T on n vertices.
    A directed k-cycle: v_0→v_1→...→v_{k-1}→v_0.
    Two cycles are the same iff one is a rotation of the other.
    Returns the count (each cycle counted once, regardless of starting vertex)."""
    count = 0
    for verts in combinations(range(n), k):
        # Try all cyclic orderings of these k vertices
        # Fix the smallest vertex as start, try all (k-1)! orderings of the rest
        start = verts[0]
        rest = list(verts[1:])
        for perm in permutations(rest):
            cycle = (start,) + perm
            # Check if this is a valid directed cycle
            valid = True
            for i in range(k):
                if T[cycle[i]][cycle[(i+1)%k]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
                break  # Found at least one valid cycle for this vertex set
        # Wait: there can be MULTIPLE directed cycles on the same vertex set!
        # Actually for a tournament on k vertices, there are (k-1)!/2 directed
        # Hamiltonian cycles if the tournament is "strongly connected enough".
        # We need to count ALL of them.
        # Let me rewrite to count all distinct directed cycles.

    # REWRITE: count ALL directed k-cycles
    count = 0
    for verts in combinations(range(n), k):
        start = verts[0]
        rest = list(verts[1:])
        for perm in permutations(rest):
            cycle = (start,) + perm
            valid = True
            for i in range(k):
                if T[cycle[i]][cycle[(i+1)%k]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
    # Each cycle counted once (fixed start = smallest vertex, tried all orders of rest)
    return count

def fast_H_v2(T, n):
    """
    Compute H(T) exactly for n ≤ 8 using OCF truncation.
    H = 1 + 2*(c3 + c5 + c7) + 4*alpha_2
    where c_k = number of directed k-cycles, alpha_2 = number of disjoint pairs.
    """
    # Count all odd directed cycles
    cycles_by_vset = []  # list of (vertex_set, count_of_cycles_on_that_set)

    # 3-cycles
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                # Check both orientations
                has = 0
                if T[a][b] and T[b][c] and T[c][a]: has += 1
                if T[a][c] and T[c][b] and T[b][a]: has += 1
                if has > 0:
                    cycles_by_vset.append(frozenset({a,b,c}))
                    # Note: a tournament on 3 vertices has exactly 1 directed 3-cycle
                    # (if it has any). So has is always 0 or 1.

    # 5-cycles
    if n >= 5:
        for verts in combinations(range(n), 5):
            # Count directed 5-cycles on these 5 vertices
            start = verts[0]
            rest = list(verts[1:])
            cycle_count = 0
            for perm in permutations(rest):
                cycle = (start,) + perm
                valid = True
                for i in range(5):
                    if T[cycle[i]][cycle[(i+1)%5]] != 1:
                        valid = False
                        break
                if valid:
                    cycle_count += 1
            # Each directed 5-cycle on these 5 vertices is counted once
            # (with start fixed at smallest vertex)
            for _ in range(cycle_count):
                cycles_by_vset.append(frozenset(verts))

    # 7-cycles (only at n=7)
    if n == 7:
        verts = tuple(range(7))
        start = 0
        rest = list(range(1, 7))
        for perm in permutations(rest):
            cycle = (0,) + perm
            valid = True
            for i in range(7):
                if T[cycle[i]][cycle[(i+1)%7]] != 1:
                    valid = False
                    break
            if valid:
                cycles_by_vset.append(frozenset(verts))

    alpha_1 = len(cycles_by_vset)

    # Count vertex-disjoint pairs
    alpha_2 = 0
    for i in range(len(cycles_by_vset)):
        for j in range(i+1, len(cycles_by_vset)):
            if not (cycles_by_vset[i] & cycles_by_vset[j]):
                alpha_2 += 1

    return 1 + 2*alpha_1 + 4*alpha_2

# ==============================================================
# BENCHMARK: Correct version
# ==============================================================
def benchmark(n, sample_size=None):
    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m
    if sample_size is None:
        sample_size = min(N, 100)

    np.random.seed(42)
    sample = np.random.choice(N, sample_size, replace=False) if N > sample_size else range(N)

    print(f"\n  n={n}, sample={sample_size}:")

    # DP
    t0 = time.time()
    h_dp = []
    for t in sample:
        T = tiling_to_tournament(n, tiles, t)
        h_dp.append(int(count_hamiltonian_paths(T, n)))
    t_dp = time.time() - t0

    # OCF v2
    t0 = time.time()
    h_ocf = []
    for t in sample:
        T = tiling_to_tournament(n, tiles, t)
        h_ocf.append(fast_H_v2(T, n))
    t_ocf = time.time() - t0

    matches = sum(1 for a, b in zip(h_dp, h_ocf) if a == b)
    print(f"    DP: {t_dp:.3f}s ({t_dp/sample_size*1000:.2f}ms each)")
    print(f"    OCF: {t_ocf:.3f}s ({t_ocf/sample_size*1000:.2f}ms each)")
    print(f"    Matches: {matches}/{sample_size}")

    if matches < sample_size:
        # Show mismatches
        for i, (a, b) in enumerate(zip(h_dp, h_ocf)):
            if a != b and i < 5:
                t = list(sample)[i]
                print(f"      MISMATCH: tiling {t}, DP={a}, OCF={b}, diff={b-a}")

if __name__ == "__main__":
    print("CORRECTED FAST H BENCHMARK")
    print("=" * 60)
    for n in [5, 6, 7]:
        benchmark(n)
