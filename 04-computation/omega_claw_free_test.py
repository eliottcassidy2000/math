#!/usr/bin/env python3
"""
Test claw-freeness of Omega(T) at n=9 and beyond.

KEY OBSERVATION: A claw (K_{1,3}) in Omega(T) requires 4 odd cycles C0,C1,C2,C3
where C1,C2,C3 are pairwise vertex-disjoint and each shares a vertex with C0.
Since each Ci has >= 3 vertices and C1,C2,C3 are pairwise v.d.,
|V(C1 union C2 union C3)| >= 9. The total distinct vertices >= 9.

Therefore Omega(T) is TRIVIALLY claw-free for n <= 8.

At n=9, we can construct an explicit claw with all 3-cycles:
  C0 = (0,1,2), C1 = (0,3,4), C2 = (1,5,6), C3 = (2,7,8)

This script verifies:
  1. The construction works at n=9 (Omega(T) has a claw)
  2. Omega(T) may still be PERFECT despite having claws (SPGT doesn't need claw-free)

Instance: opus-2026-03-05-S7
"""

import sys
import os
import random
from itertools import combinations

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '03-artifacts', 'code'))
from tournament_lib import (
    hamiltonian_path_count, find_odd_cycles, conflict_graph,
    independence_poly_at_fast
)


def make_tournament_with_claw_9():
    """Construct a 9-vertex tournament with an explicit claw in Omega(T).

    Cycles:
      C0 = 0->1->2->0
      C1 = 0->3->4->0
      C2 = 1->5->6->1
      C3 = 2->7->8->2
    """
    n = 9
    T = [[0]*n for _ in range(n)]

    def set_arc(a, b):
        T[a][b] = 1
        T[b][a] = 0

    # C0: 0->1->2->0
    set_arc(0, 1)
    set_arc(1, 2)
    set_arc(2, 0)

    # C1: 0->3->4->0
    set_arc(0, 3)
    set_arc(3, 4)
    set_arc(4, 0)

    # C2: 1->5->6->1
    set_arc(1, 5)
    set_arc(5, 6)
    set_arc(6, 1)

    # C3: 2->7->8->2
    set_arc(2, 7)
    set_arc(7, 8)
    set_arc(8, 2)

    # Fill remaining arcs arbitrarily (just need a valid tournament)
    remaining_pairs = []
    for i in range(n):
        for j in range(i+1, n):
            if T[i][j] == 0 and T[j][i] == 0:
                remaining_pairs.append((i, j))

    random.seed(42)
    for i, j in remaining_pairs:
        if random.random() < 0.5:
            set_arc(i, j)
        else:
            set_arc(j, i)

    # Verify it's a valid tournament
    for i in range(n):
        for j in range(n):
            if i != j:
                assert T[i][j] + T[j][i] == 1, f"Not a tournament at ({i},{j})"

    return T


def has_claw(cycles, cg_adj, m):
    """Check if conflict graph has an induced K_{1,3} (claw)."""
    adj_sets = [set(j for j in range(m) if cg_adj[i][j]) for i in range(m)]

    for center in range(m):
        # Find non-neighbors of center
        non_nbrs = [v for v in range(m) if v != center and v not in adj_sets[center]]
        # Neighbors of center
        nbrs = list(adj_sets[center])

        # Need 3 neighbors of center that are pairwise non-adjacent
        if len(nbrs) < 3:
            continue

        for triple in combinations(nbrs, 3):
            a, b, c = triple
            if b not in adj_sets[a] and c not in adj_sets[a] and c not in adj_sets[b]:
                return True, center, triple

    return False, None, None


def cycle_vertices(cycle_tuple):
    """Extract vertex set from a cycle tuple."""
    return set(cycle_tuple)


def test_claw_n9():
    """Test claw existence at n=9."""
    print("=" * 60)
    print("TEST 1: Explicit construction at n=9")
    print("=" * 60)

    T = make_tournament_with_claw_9()
    n = len(T)
    H = hamiltonian_path_count(T)
    cycles = find_odd_cycles(T)
    m = len(cycles)

    print(f"n={n}, H={H}, #odd_cycles={m}")

    if m == 0:
        print("No odd cycles — no claw possible")
        return

    cg = conflict_graph(cycles)
    m = len(cg)

    # Check claw
    found, center, triple = has_claw(cycles, cg, m)
    if found:
        print(f"\nCLAW FOUND!")
        print(f"  Center cycle: {cycles[center]}")
        print(f"  Leaf cycles: {[cycles[t] for t in triple]}")
        # Verify vertex-disjointness of leaves
        leaf_verts = [cycle_vertices(cycles[t]) for t in triple]
        for i in range(3):
            for j in range(i+1, 3):
                inter = leaf_verts[i] & leaf_verts[j]
                assert len(inter) == 0, f"Leaves {i},{j} share vertices {inter}!"
        # Verify each leaf touches center
        center_verts = cycle_vertices(cycles[center])
        for i in range(3):
            inter = leaf_verts[i] & center_verts
            assert len(inter) > 0, f"Leaf {i} doesn't touch center!"
            print(f"    Leaf {i} {cycles[triple[i]]} shares {inter} with center")
        print("  All conditions verified!")
    else:
        print("No claw found in this particular tournament.")
        print("Trying more random completions...")
        # Try many random fill-ins
        found_any = False
        for seed in range(200):
            random.seed(seed)
            T2 = make_tournament_with_claw_9_seed(seed)
            cycles2 = find_odd_cycles(T2)
            if len(cycles2) == 0:
                continue
            cg2 = conflict_graph(cycles2)
            m2 = len(cg2)
            f2, c2, t2 = has_claw(cycles2, cg2, m2)
            if f2:
                print(f"\n  CLAW FOUND with seed={seed}!")
                print(f"  Center: {cycles2[c2]}, Leaves: {[cycles2[t] for t in t2]}")
                found_any = True
                break
        if not found_any:
            print("  No claw found in 200 attempts?!")


def make_tournament_with_claw_9_seed(seed):
    """Like make_tournament_with_claw_9 but with variable seed."""
    n = 9
    T = [[0]*n for _ in range(n)]

    def set_arc(a, b):
        T[a][b] = 1
        T[b][a] = 0

    set_arc(0, 1); set_arc(1, 2); set_arc(2, 0)
    set_arc(0, 3); set_arc(3, 4); set_arc(4, 0)
    set_arc(1, 5); set_arc(5, 6); set_arc(6, 1)
    set_arc(2, 7); set_arc(7, 8); set_arc(8, 2)

    random.seed(seed)
    for i in range(n):
        for j in range(i+1, n):
            if T[i][j] == 0 and T[j][i] == 0:
                if random.random() < 0.5:
                    set_arc(i, j)
                else:
                    set_arc(j, i)
    return T


def test_claw_random_n9(num_trials=500):
    """Test claw-freeness on random n=9 tournaments."""
    print("\n" + "=" * 60)
    print("TEST 2: Random n=9 tournaments")
    print("=" * 60)

    claw_count = 0
    for trial in range(num_trials):
        random.seed(trial)
        n = 9
        T = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    T[i][j] = 1
                else:
                    T[j][i] = 1

        cycles = find_odd_cycles(T)
        if len(cycles) == 0:
            continue
        cg = conflict_graph(cycles)
        m = len(cg)

        found, center, triple = has_claw(cycles, cg, m)
        if found:
            claw_count += 1
            if claw_count <= 3:
                print(f"  Trial {trial}: CLAW! center={cycles[center]}, "
                      f"leaves={[cycles[t] for t in triple]}")

    print(f"\nClaws found: {claw_count}/{num_trials} ({100*claw_count/num_trials:.1f}%)")


def test_perfectness_with_claws(num_trials=200):
    """Check if Omega(T) remains perfect even when it has claws (n=9)."""
    print("\n" + "=" * 60)
    print("TEST 3: Perfectness despite claws at n=9")
    print("=" * 60)

    perf_fail = 0
    claw_and_perfect = 0
    claw_count = 0
    ocf_fail = 0

    for trial in range(num_trials):
        random.seed(trial + 1000)
        n = 9
        T = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    T[i][j] = 1
                else:
                    T[j][i] = 1

        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        if len(cycles) == 0:
            continue

        I2 = independence_poly_at_fast(cycles, 2)
        if H != I2:
            ocf_fail += 1

        cg = conflict_graph(cycles)
        m = len(cg)

        found_claw, _, _ = has_claw(cycles, cg, m)
        if found_claw:
            claw_count += 1

        # Check perfectness (expensive — only for small m)
        if m <= 15:
            adj_sets = [set(j for j in range(m) if cg[i][j]) for i in range(m)]
            perf = is_perfect_spgt(adj_sets, m)
            if not perf:
                perf_fail += 1
                print(f"  Trial {trial}: NOT PERFECT! m={m}")
            elif found_claw:
                claw_and_perfect += 1

    print(f"Tested: {num_trials}")
    print(f"OCF failures: {ocf_fail}")
    print(f"Claws found: {claw_count}")
    print(f"Non-perfect: {perf_fail}")
    print(f"Has claw AND perfect: {claw_and_perfect}")


def has_odd_hole(adj_sets, m):
    """Check for induced odd cycle of length >= 5."""
    if m < 5:
        return False
    for length in range(5, min(m + 1, 12), 2):
        for verts in combinations(range(m), length):
            # Check each vertex has exactly 2 neighbors in the subset
            ok = True
            for v in verts:
                nbr_count = sum(1 for u in verts if u != v and u in adj_sets[v])
                if nbr_count != 2:
                    ok = False
                    break
            if not ok:
                continue
            # Check connectivity
            vis = {verts[0]}
            stk = [verts[0]]
            while stk:
                u = stk.pop()
                for w in verts:
                    if w not in vis and w in adj_sets[u]:
                        vis.add(w)
                        stk.append(w)
            if len(vis) == length:
                return True
    return False


def is_perfect_spgt(adj_sets, m):
    """Check perfectness via Strong Perfect Graph Theorem."""
    if m <= 4:
        return True
    if has_odd_hole(adj_sets, m):
        return False
    # Check complement for odd holes (= odd antiholes in G)
    comp = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if j not in adj_sets[i]:
                comp[i].add(j)
                comp[j].add(i)
    return not has_odd_hole(comp, m)


def test_vertex_count_argument():
    """Verify the vertex-counting argument for claw impossibility."""
    print("=" * 60)
    print("VERTEX COUNT ARGUMENT")
    print("=" * 60)
    print()
    print("A claw in Omega(T) requires cycles C0, C1, C2, C3 where:")
    print("  - C1, C2, C3 are pairwise vertex-disjoint")
    print("  - Each Ci shares >= 1 vertex with C0")
    print()
    print("Since each Ci has >= 3 vertices (odd cycle) and C1,C2,C3")
    print("are pairwise disjoint: |V(C1) u V(C2) u V(C3)| >= 9.")
    print()
    print("Let ki = |V(Ci) n V(C0)|. Disjointness of C1,C2,C3 within V(C0)")
    print("requires k1+k2+k3 <= |V(C0)|.")
    print()
    print("Total vertices = |V(C0)| + sum(|V(Ci)| - ki)")
    print("              >= |V(C0)| + 9 - |V(C0)| = 9")
    print()
    print("Therefore: Omega(T) is TRIVIALLY claw-free for n <= 8.")
    print("The computational verification at n<=8 confirms a trivial fact!")
    print()
    print("At n=9, claws become possible. Testing...")


if __name__ == "__main__":
    test_vertex_count_argument()
    test_claw_n9()
    test_claw_random_n9(300)
    test_perfectness_with_claws(100)
