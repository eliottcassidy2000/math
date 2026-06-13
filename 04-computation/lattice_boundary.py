#!/usr/bin/env python3
"""
lattice_boundary.py -- The boundary of the achievable (a1, a2) lattice

Key question: What are the constraints on (a1, a2) beyond a1 >= 0, a2 >= 0?

From exhaustive n=5: ONLY a2=0 appears. Every tournament has all cycles
pairwise conflicting.

From exhaustive n=6: a2 in {0,1,2,3,4}. But not all (a1,a2) pairs are
achievable. For instance:
  - a1=3 is IMPOSSIBLE (confirms H=7 gap)
  - a1=10, a2=0 is IMPOSSIBLE (part of H=21 gap)

This script maps the boundary precisely at n=5,6 (exhaustive) and
characterizes the constraints that create gaps.

Also explores: the relationship between a1 (total cycles) and a2
(disjoint pairs) as constrained by tournament geometry.

Author: kind-pasteur-2026-03-14-S63
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict


def bits_to_adj(bits, n):
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


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1 << n):
        bits_set = bin(mask).count('1')
        if bits_set < 2:
            continue
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            prev = mask ^ (1 << v)
            s = 0
            for u in range(n):
                if (prev & (1 << u)) and A[u][v]:
                    s += dp.get((prev, u), 0)
            if s:
                dp[(mask, v)] = s
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp or dp[(mask, v)] == 0:
                continue
            c = dp[(mask, v)]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + c
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and dp[(full, v)] > 0:
            if A[verts[v]][verts[0]]:
                total += dp[(full, v)]
    return total


def get_cycle_data(A, n):
    """Get directed odd cycles with vertex sets."""
    cycles = []
    for k in range(3, n+1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            nc = count_directed_ham_cycles(A, verts)
            for _ in range(nc):
                cycles.append(frozenset(subset))
    a1 = len(cycles)
    a2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                a2 += 1
    return a1, a2, cycles


def score_sequence(A, n):
    """Get score sequence (sorted out-degrees)."""
    return tuple(sorted(sum(A[i]) for i in range(n)))


def main():
    # ============================================================
    print("=" * 70)
    print("PART 1: Exhaustive (a1, a2) lattice at n=5")
    print("=" * 70)

    n = 5
    total_bits = n*(n-1)//2
    lattice_5 = defaultdict(int)
    score_to_a1 = defaultdict(set)

    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        a1, a2, _ = get_cycle_data(A, n)
        lattice_5[(a1, a2)] += 1
        sc = score_sequence(A, n)
        score_to_a1[sc].add(a1)

    print(f"\n  Lattice points: {len(lattice_5)}")
    print(f"\n  a1 | a2 | count | H")
    for (a1, a2) in sorted(lattice_5):
        H = 1 + 2*a1 + 4*a2
        print(f"  {a1:2d} | {a2:2d} | {lattice_5[(a1,a2)]:5d} | {H:3d}")

    # Why is a1=3 impossible?
    # a1 = number of directed odd cycles = c3 + c5 (with multiplicity)
    # At n=5: c5 = 0 or 2 (undirected 5-cycle has 2 directed cycles)
    # c3 counts directed 3-cycles
    print(f"\n  c3 and c5 at n=5:")
    c3c5_dist = defaultdict(int)
    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        # c3
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
                c3 += 1  # count vertex sets
        # directed c3
        dc3 = 0
        for a, b, c in combinations(range(n), 3):
            dc3 += A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]
        # c5: directed Hamiltonian cycles on all 5 vertices
        dc5 = count_directed_ham_cycles(A, list(range(n)))
        c3c5_dist[(dc3, dc5)] += 1

    print(f"  (dc3, dc5) -> count:")
    for key in sorted(c3c5_dist):
        a1 = key[0] + key[1]
        print(f"    dc3={key[0]}, dc5={key[1]}: a1={a1}, count={c3c5_dist[key]}")

    # The key insight: at n=5, c3 and c5 are NOT independent
    # c3=1 forces c5=0 or c5=2?
    print(f"\n  c5 constrained by c3:")
    c3_to_c5 = defaultdict(set)
    for (dc3, dc5), cnt in c3c5_dist.items():
        c3_to_c5[dc3].add(dc5)
    for dc3 in sorted(c3_to_c5):
        print(f"    dc3={dc3}: possible dc5 = {sorted(c3_to_c5[dc3])}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 2: Exhaustive (a1, a2) lattice at n=6")
    print("=" * 70)

    n = 6
    total_bits = n*(n-1)//2
    lattice_6 = defaultdict(int)

    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        a1, a2, _ = get_cycle_data(A, n)
        lattice_6[(a1, a2)] += 1

    print(f"\n  Lattice points: {len(lattice_6)}")
    print(f"\n  The (a1, a2) grid:")

    a2_max = max(a2 for a1, a2 in lattice_6)
    a1_max = max(a1 for a1, a2 in lattice_6)

    # Print as grid
    print(f"\n  a2\\a1", end="")
    for a1 in range(a1_max + 1):
        print(f" {a1:3d}", end="")
    print()

    for a2 in range(a2_max + 1):
        print(f"  {a2:3d} ", end="")
        for a1 in range(a1_max + 1):
            cnt = lattice_6.get((a1, a2), 0)
            if cnt > 0:
                print(f" {cnt:3d}", end="")
            else:
                print("   .", end="")
        print()

    # a1 achievability
    a1_vals_6 = sorted(set(a1 for a1, a2 in lattice_6))
    a1_gaps_6 = [a for a in range(max(a1_vals_6)+1) if a not in set(a1_vals_6)]
    print(f"\n  Missing a1 values: {a1_gaps_6}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 3: The boundary curve a2_max(a1)")
    print("=" * 70)

    # For each a1, what is the max achievable a2?
    print(f"\n  n=6: a2_max(a1)")
    for a1 in a1_vals_6:
        a2s = [a2 for (a1_, a2) in lattice_6 if a1_ == a1]
        a2_mx = max(a2s)
        print(f"    a1={a1:2d}: a2_max={a2_mx}")

    # Also compute a1/a2 ratio at the boundary
    print(f"\n  The boundary is determined by tournament geometry:")
    print(f"  a2 <= floor(a1 * something)")
    print(f"  But the exact bound is NOT linear.")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 4: Score sequence determines cycle counts?")
    print("=" * 70)

    # Does the score sequence determine (a1, a2)?
    n = 5
    total_bits = n*(n-1)//2
    score_to_a1a2 = defaultdict(set)

    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        a1, a2, _ = get_cycle_data(A, n)
        sc = score_sequence(A, n)
        score_to_a1a2[sc].add((a1, a2))

    print(f"\n  n=5: score -> (a1, a2) map:")
    for sc in sorted(score_to_a1a2):
        pairs = sorted(score_to_a1a2[sc])
        if len(pairs) > 1:
            print(f"    score={sc}: {pairs}  [MULTIPLE]")
        else:
            print(f"    score={sc}: {pairs[0]}")

    n = 6
    total_bits = n*(n-1)//2
    score_to_a1a2 = defaultdict(set)

    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        a1, a2, _ = get_cycle_data(A, n)
        sc = score_sequence(A, n)
        score_to_a1a2[sc].add((a1, a2))

    print(f"\n  n=6: score -> (a1, a2) map:")
    for sc in sorted(score_to_a1a2):
        pairs = sorted(score_to_a1a2[sc])
        if len(pairs) > 1:
            print(f"    score={sc}: {pairs}  [MULTIPLE]")
        else:
            print(f"    score={sc}: {pairs[0]}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 5: Why is a1=3 impossible? (n=5 proof)")
    print("=" * 70)

    print("""
  At n=5: a1 = dc3 + dc5 (total directed odd cycles with multiplicity)
  From Part 1 data:
    dc3=0: dc5 in {0, 2}. a1 in {0, 2}.
    dc3=1: dc5=0. a1=1.
    dc3=2: dc5=0. a1=2.
    dc3=4: dc5=0. a1=4.
    dc3=5: dc5=0. a1=5.
    dc3=6: dc5=0. a1=6.
    dc3=5: dc5=2. a1=7.

  So a1=3 requires dc3=3 (since dc5 is 0 or 2, and 3-2=1 needs dc3=1 but
  then dc5=2 is impossible when dc3=1).

  Wait, let me re-examine: dc3=1, dc5=2 gives a1=3. Is this possible?
  dc3=1 means exactly one 3-vertex set forming a directed 3-cycle.
  dc5=2 means there IS a Hamiltonian 5-cycle.

  But at n=5: if there is a Hamiltonian cycle, then at least 2 of the
  C(5,2)=10 vertex triples form 3-cycles. So dc3=1 with dc5=2 is impossible!

  Also: dc3=3 means exactly 3 directed 3-cycles. This means 1 or 2
  undirected 3-cycle vertex sets. Each undirected set gives exactly 2
  directed cycles, so dc3 is always EVEN for 3-cycles.

  WAIT: Each 3-vertex tournament is either transitive (0 directed 3-cycles)
  or cyclic (2 directed 3-cycles: fwd and bwd).
  So dc3 is ALWAYS even!

  dc3 in {0, 2, 4, 6, 8, 10, ...} only.
  At n=5: C(5,3)=10 vertex triples, each contributes 0 or 2 directed cycles.
  So dc3 is even, max = 20 (impossible).
  Observed: dc3 in {0, 2, 4, 6, 8, 10}... let me check.
""")

    # Verify dc3 is always even
    n = 5
    total_bits = n*(n-1)//2
    dc3_vals = set()
    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        dc3 = sum(
            A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]
            for a, b, c in combinations(range(n), 3)
        )
        dc3_vals.add(dc3)
    print(f"  n=5: dc3 values = {sorted(dc3_vals)}")
    print(f"  All even: {all(v % 2 == 0 for v in dc3_vals)}")

    # Same for n=6
    n = 6
    total_bits = n*(n-1)//2
    dc3_vals = set()
    dc5_vals = set()
    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        dc3 = sum(
            A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]
            for a, b, c in combinations(range(n), 3)
        )
        dc3_vals.add(dc3)

        # dc5
        dc5 = 0
        for subset in combinations(range(n), 5):
            dc5 += count_directed_ham_cycles(A, list(subset))
        dc5_vals.add(dc5)
    print(f"\n  n=6: dc3 values = {sorted(dc3_vals)}")
    print(f"  n=6: dc5 values = {sorted(dc5_vals)}")

    # So dc3 is always EVEN, and dc5 is also always EVEN.
    # Therefore a1 = dc3 + dc5 is always EVEN!
    # But wait, a1=1 and a1=5 appear at n=5...

    # Let me recount. At n=5, the "directed 3-cycles" count both fwd and bwd?
    # A 3-cycle on {a,b,c} that is cyclic gives fwd AND bwd, so 2 directed.
    # But in our cycle enumeration for the independence polynomial,
    # each DIRECTED cycle is a separate vertex in Omega.
    # So alpha_1 counts DIRECTED cycles.

    # But we showed dc3 is always even. And dc5 is even (2 per undirected).
    # So alpha_1 should always be even!
    # But our exhaustive data shows a1=1 at n=5 (120 tournaments). Contradiction?

    # Let me recheck a1=1 at n=5
    print(f"\n  Rechecking a1=1 at n=5:")
    n = 5
    total_bits = n*(n-1)//2
    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        a1, a2, cycles = get_cycle_data(A, n)
        if a1 == 1:
            print(f"    bits={bits}: a1=1, cycles={[list(c) for c in cycles]}")
            # Print adjacency
            for i in range(n):
                row = [A[i][j] for j in range(n)]
                print(f"      {row}")
            break

    # If we have a1=1, that means a single directed cycle.
    # But each undirected cycle gives 2 directed...
    # Unless our counting is wrong!
    # Let me check: for a 3-cycle {0,1,2} with 0->1->2->0,
    # count_directed_ham_cycles returns 1 (only fwd, not bwd).
    # NO: we count BOTH directions. So it should return 2.

    print(f"\n  Test: 3-cycle 0->1->2->0:")
    A_test = [[0]*3 for _ in range(3)]
    A_test[0][1] = 1
    A_test[1][2] = 1
    A_test[2][0] = 1
    result = count_directed_ham_cycles(A_test, [0,1,2])
    print(f"    count_directed_ham_cycles = {result}")
    # The backward cycle is 0->2->1->0 which requires A[0][2]=1, but A[2][0]=1, A[0][2]=0.
    # So for a tournament, each 3-vertex cyclic set has EXACTLY 1 directed 3-cycle!
    # NOT 2! Because the backward direction is the REVERSAL tournament.
    print(f"    In a tournament, each cyclic triple has exactly 1 directed 3-cycle")
    print(f"    (the other direction uses the opposite arcs)")

    # So dc3 = number of cyclic triples = c3 (vertex sets)
    # This is NOT always even!
    # c3 at n=5: check
    c3_vals_5 = set()
    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
                c3 += 1
        c3_vals_5.add(c3)
    print(f"\n  n=5: c3 (vertex sets) = {sorted(c3_vals_5)}")

    # And c5 at n=5?
    c5_dir_vals = set()
    for bits in range(2**total_bits):
        A = bits_to_adj(bits, n)
        dc5 = count_directed_ham_cycles(A, list(range(n)))
        c5_dir_vals.add(dc5)
    print(f"  n=5: dc5 (directed 5-cycles) = {sorted(c5_dir_vals)}")

    # Each Hamiltonian cycle on 5 vertices: the tournament determines
    # exactly 1 direction (the one compatible with all arcs).
    # Wait no: a Hamiltonian cycle 0-1-2-3-4-0 requires specific arc directions.
    # The reverse direction 0-4-3-2-1-0 is a DIFFERENT Hamiltonian cycle
    # and may or may not exist.

    # Actually for 5-vertex tournaments: an undirected Hamiltonian cycle
    # has exactly 2 orientations, but typically only 0 or 1 match the tournament.
    # Unless the tournament is a directed 5-cycle in BOTH directions...
    # which is impossible since a tournament has exactly one arc direction per pair.

    # So dc5 counts the actual number of directed Hamiltonian 5-cycles.
    # This CAN be odd (e.g., 1 or 3).

    print(f"\n  CORRECTED understanding:")
    print(f"  Each cyclic triple -> exactly 1 directed 3-cycle")
    print(f"  Each directed 5-cycle is counted individually")
    print(f"  a1 = c3 + dc5 can be any non-negative integer")
    print(f"  The constraint a1 != 3 is genuinely non-trivial!")

    print("\nDone.")


if __name__ == '__main__':
    main()
