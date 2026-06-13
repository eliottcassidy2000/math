#!/usr/bin/env python3
"""
Test whether Omega(T) can have a C5 (5-hole) at n=7 or n=8.

Vertex counting for C5 in the 3-cycle conflict graph:
  5 three-cycles C1,...,C5 forming a 5-hole.
  Consecutive pairs share >= 1 vertex, non-consecutive share 0.

  If each consecutive pair shares exactly 1 vertex: need 5 shared + 5 private = 10 vertices.
  If ADJACENT consecutive pairs share 2 vertices (say C1,C2 share {a,b}):
    We can potentially reduce to 9 or even 8 vertices.

This script constructs explicit tournaments with C5 in Omega at small n.

Instance: opus-2026-03-05-S7
"""

import random
from itertools import combinations


def make_tournament_from_arcs(n, forced_arcs, seed=42):
    """Build a tournament with specified forced arcs."""
    T = [[0]*n for _ in range(n)]

    for a, b in forced_arcs:
        T[a][b] = 1
        T[b][a] = 0

    # Fill remaining arcs randomly
    random.seed(seed)
    for i in range(n):
        for j in range(i+1, n):
            if T[i][j] == 0 and T[j][i] == 0:
                if random.random() < 0.5:
                    T[i][j] = 1
                else:
                    T[j][i] = 1

    # Verify
    for i in range(n):
        for j in range(n):
            if i != j:
                assert T[i][j] + T[j][i] == 1

    return T


def verify_directed_3cycle(T, a, b, c):
    """Check a->b->c->a."""
    return T[a][b] and T[b][c] and T[c][a]


def count_ham_paths(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def find_all_3cycles(T, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if T[a][b] and T[b][c] and T[c][a]:
            cycles.append((a, b, c))
        elif T[a][c] and T[c][b] and T[b][a]:
            cycles.append((a, c, b))
    return cycles


def find_odd_hole_in_3cycles(cycles_3):
    """Find a C5 (5-hole) in the 3-cycle conflict graph."""
    m = len(cycles_3)
    cycle_sets = [set(c) for c in cycles_3]

    # Build adjacency
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = True

    # Find C5
    for combo in combinations(range(m), 5):
        # Check if these 5 form a C5
        # Need: each has exactly 2 neighbors in combo, and it's connected
        v_list = list(combo)
        deg_ok = True
        nbrs = {}
        for v in v_list:
            n_v = [u for u in v_list if u != v and adj[v][u]]
            if len(n_v) != 2:
                deg_ok = False
                break
            nbrs[v] = n_v
        if not deg_ok:
            continue
        # Check connected
        vis = {v_list[0]}
        stk = [v_list[0]]
        while stk:
            u = stk.pop()
            for w in nbrs.get(u, []):
                if w not in vis:
                    vis.add(w)
                    stk.append(w)
        if len(vis) == 5:
            return combo
    return None


def test_c5_at_n8():
    """Construct n=8 tournament with C5 in Omega(T)."""
    print("=" * 60)
    print("TEST: Explicit C5 construction at n=8")
    print("=" * 60)

    # 5 directed 3-cycles forming a C5:
    # C1 = (0,1,5): 0->1->5->0
    # C2 = (1,5,2): 1->5->2->1  (wait, this requires 5->2 AND 1->5 AND 2->1)
    #   C1 requires 1->5 ✓. C2 requires 1->5 ✓ (consistent).
    # But C2 = {1,5,2}: to be a 3-cycle we need a cyclic order.
    # Say 2->1->5->2. Then: 2->1, 1->5, 5->2.
    # C1 has 5->0. C2 has 5->2. OK (different targets).
    # C1 has 0->1. C2 has 2->1. OK.
    # C3 = (2,3,6): 2->3->6->2
    # C4 = (3,4,7): 3->4->7->3
    # C5 = (0,4,7): 0->4->7->0

    # Check vertex sets for C5 conditions:
    # C1={0,1,5}, C2={1,2,5}, C3={2,3,6}, C4={3,4,7}, C5={0,4,7}
    # Adjacent (share vertex):
    #   C1∩C2 = {1,5} ✓
    #   C2∩C3 = {2} ✓
    #   C3∩C4 = {3} ✓
    #   C4∩C5 = {4,7} ✓
    #   C5∩C1 = {0} ✓
    # Non-adjacent (disjoint):
    #   C1∩C3 = ∅ ✓
    #   C1∩C4 = ∅ ✓
    #   C2∩C4 = ∅ ✓
    #   C2∩C5 = ∅ ✓
    #   C3∩C5 = ∅ ✓

    print("Target C5:")
    print("  C1 = (0,1,5): 0->1->5->0")
    print("  C2 = (2,1,5): 2->1->5->2")
    print("  C3 = (2,3,6): 2->3->6->2")
    print("  C4 = (3,4,7): 3->4->7->3")
    print("  C5 = (0,4,7): 0->4->7->0")
    print()

    # Forced arcs
    forced = [
        (0, 1), (1, 5), (5, 0),  # C1
        (2, 1), (5, 2),          # C2 (1->5 already from C1)
        (2, 3), (3, 6), (6, 2),  # C3
        (3, 4), (4, 7), (7, 3),  # C4
        (0, 4), (7, 0),          # C5 (4->7 already from C4)
    ]

    # Check for contradictions
    arc_set = {}
    for a, b in forced:
        if (a, b) in arc_set:
            continue
        if (b, a) in arc_set:
            print(f"  CONTRADICTION: both {a}->{b} and {b}->{a} forced!")
            return False
        arc_set[(a, b)] = True

    print(f"  Forced arcs: {len(arc_set)}, no contradictions")

    for seed in range(100):
        T = make_tournament_from_arcs(8, forced, seed=seed)

        # Verify all 5 cycles exist
        ok = True
        ok &= verify_directed_3cycle(T, 0, 1, 5)
        ok &= verify_directed_3cycle(T, 2, 1, 5)
        ok &= verify_directed_3cycle(T, 2, 3, 6)
        ok &= verify_directed_3cycle(T, 3, 4, 7)
        ok &= verify_directed_3cycle(T, 0, 4, 7)

        if not ok:
            continue

        print(f"\n  Seed {seed}: all 5 cycles verified!")

        # Now check the 3-cycle conflict graph for the C5
        cycles_3 = find_all_3cycles(T, 8)
        print(f"  Total 3-cycles: {len(cycles_3)}")

        hole = find_odd_hole_in_3cycles(cycles_3)
        if hole:
            print(f"  C5 FOUND! Cycles at indices {hole}:")
            for idx in hole:
                c = cycles_3[idx]
                print(f"    {c} -> vertices {set(c)}")

            H = count_ham_paths(T, 8)
            print(f"  H(T) = {H}")

            # Verify OCF still holds
            from itertools import combinations as comb
            # Simple I(Omega,2) via inclusion-exclusion on 3-cycles
            # Actually need ALL odd cycles, not just 3. Use tournament_lib.
            import sys, os
            sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            '..', '03-artifacts', 'code'))
            from tournament_lib import find_odd_cycles, independence_poly_at_fast
            all_cycles = find_odd_cycles(T)
            I2 = independence_poly_at_fast(all_cycles, 2)
            print(f"  I(Omega,2) = {I2}")
            print(f"  OCF holds: {H == I2}")

            return True

    print("  No seed produced all 5 valid cycles?!")
    return False


def test_c5_random_n7_n8():
    """Search for C5 in random n=7, n=8 tournaments."""
    for n in [7, 8]:
        print(f"\n{'='*60}")
        print(f"TEST: Random n={n} tournaments — C5 search")
        print(f"{'='*60}")

        c5_count = 0
        for trial in range(1000):
            T = [[0]*n for _ in range(n)]
            random.seed(trial)
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        T[i][j] = 1
                    else:
                        T[j][i] = 1

            cycles_3 = find_all_3cycles(T, n)
            if len(cycles_3) < 5:
                continue

            hole = find_odd_hole_in_3cycles(cycles_3)
            if hole:
                c5_count += 1
                if c5_count <= 3:
                    print(f"  Trial {trial}: C5 found! ({len(cycles_3)} 3-cycles)")
                    for idx in hole:
                        print(f"    {cycles_3[idx]}")

        print(f"\nC5 found: {c5_count}/1000 ({c5_count/10:.1f}%)")


if __name__ == "__main__":
    test_c5_at_n8()
    test_c5_random_n7_n8()
