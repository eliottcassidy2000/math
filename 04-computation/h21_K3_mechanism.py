#!/usr/bin/env python3
"""
Mechanism: 3 pairwise-intersecting 3-cycles force 5-cycles.

All 2880 such tournaments at n=6 have common vertex (0 without common vertex).
All have H=9.

Structural analysis: with common vertex v, 3 cycles {v,a_i,b_i},
the arcs between the a_i, b_i create forced 5-cycles.

THM-029 proof: v has out-degree d among the 6 cycle vertices.
Each 3-cycle through v: {v, out-neighbor, in-neighbor} or the reverse.
With 3 cycles: at least 3 "crossing" arcs between out and in neighbors.
These crossings force 5-cycles through v.

Let me trace the exact mechanism for one example.

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations, permutations


def analyze_example():
    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]

    # bits=8 from earlier output:
    # cycles [{0,1,4}, {0,2,4}, {0,3,4}], common={0,4}
    # Tournament: 0->[1,2,3,5], 1->[2,3,4,5], 2->[3,4,5], 3->[4,5], 4->[0,5], 5->[]
    bits = 8
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(edges):
        if bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1

    print("Example: bits=8")
    print("Tournament arcs:")
    for i in range(n):
        outs = [j for j in range(n) if j != i and A[i][j]]
        print(f"  {i} -> {outs}")

    # 3-cycles
    cycles3 = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    cycles3.append(frozenset({a, b, c}))
                    # Show direction
                    if A[a][b] and A[b][c] and A[c][a]:
                        print(f"  3-cycle: {a}->{b}->{c}->{a}")
                    else:
                        print(f"  3-cycle: {a}->{c}->{b}->{a}")

    # 5-cycles
    print("\n5-cycles:")
    cycles5 = []
    for verts5 in combinations(range(n), 5):
        v0 = verts5[0]
        rest = verts5[1:]
        for p in permutations(rest):
            seq = (v0,) + p
            ok = all(A[seq[q]][seq[(q+1)%5]] for q in range(5))
            if ok:
                cycle_str = "->".join(str(s) for s in seq) + f"->{seq[0]}"
                print(f"  5-cycle: {cycle_str}")
                cycles5.append(frozenset(verts5))
                break

    print(f"\n  {len(cycles5)} vertex sets with 5-cycles")

    # Full Omega
    all_cycles = cycles3.copy()
    # Count directed 5-cycles properly
    for verts5 in combinations(range(n), 5):
        v0 = verts5[0]
        rest = verts5[1:]
        count = 0
        for p in permutations(rest):
            seq = (v0,) + p
            ok = all(A[seq[q]][seq[(q+1)%5]] for q in range(5))
            if ok:
                count += 1
        for _ in range(count):
            all_cycles.append(frozenset(verts5))

    print(f"\nFull Omega: {len(all_cycles)} vertices")
    print(f"  3-cycle vertices: {len(cycles3)}")
    print(f"  5-cycle vertices: {len(all_cycles) - len(cycles3)}")

    # The common vertices are {0, 4}. Every 5-cycle uses 5 of 6 vertices,
    # so it always includes at least one of {0, 4}.
    # Therefore every 5-cycle shares a vertex with the 3-cycles!
    print("\n  All 5-cycles share vertices with 3-cycles:")
    for c5 in cycles5:
        shared = c5 & (frozenset({0, 4}))
        print(f"    5-cycle {set(c5)}: shares {set(shared)} with common verts")


def prove_general_mechanism():
    """
    THEOREM: If tournament T has exactly 3 pairwise-intersecting 3-cycles,
    all sharing a common vertex v, then T has at least one 5-cycle through v.

    Proof:
    Let the 3 cycles be C1={v,a1,b1}, C2={v,a2,b2}, C3={v,a3,b3}.
    Each Ci: v->ai->bi->v or v->bi->ai->v (two orientations).

    Case 1: All same orientation, e.g., v->ai->bi->v for all i.
    Then {a1,a2,a3} are out-neighbors of v, {b1,b2,b3} are in-neighbors.
    The arcs ai->bj (j != i) create potential 5-cycles.

    Actually, the THM-029 argument is simpler:
    v has 3 out-neighbors among cycle vertices and 3 in-neighbors.
    Among the out-neighbors, their mutual arcs and arcs to in-neighbors
    are constrained. Since there are exactly 3 cycles (not 4+), the
    constraint is tight.

    BUT: the key point is that 5-cycles ALWAYS exist.
    At n=6 with 6 vertices and 3 three-cycles: C(5,4) = 5 possible
    5-vertex subsets (each omitting one of 6 vertices). Each subset
    is a 5-vertex tournament and may have Hamiltonian cycles.

    The near-transitive structure (only 3 cycles) limits the tournament
    to being "almost transitive" but the 3 crossings create 5-cycles.
    """
    print("\n=== General mechanism ===")
    print("Common vertex v. 3 cycles through v using vertices a1..a3, b1..b3.")
    print("Without loss: score sequence is nearly transitive.")
    print("The 3 backward arcs that create the 3-cycles also create 5-cycle paths.")
    print()
    print("KEY STRUCTURAL FACT (proved exhaustively at n=6):")
    print("3 pairwise-intersecting 3-cycles => at least 1 five-cycle.")
    print("The 5-cycle shares vertices with the 3-cycles, so in full Omega,")
    print("the component is strictly larger than K3.")
    print("Therefore I(component, 2) > I(K3, 2) = 7.")
    print()
    print("At n=6: all such tournaments have H = I(Omega, 2) = 9.")
    print("Decomposition: single component with I = 9 (not 7).")


if __name__ == '__main__':
    analyze_example()
    prove_general_mechanism()
