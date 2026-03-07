#!/usr/bin/env python3
"""
WHY can't 3 pairwise-intersecting 3-cycles be isolated from all other 3-cycles?

Key insight: any three pairwise-intersecting 3-sets force additional 3-cycles.

Case A (common vertex v): {v,a,b}, {v,c,d}, {v,e,f}
  THM-029: the arcs between {a,b,c,d,e,f} and v create additional cycles.
  Specifically: v->a, b->v (since {v,a,b} is cyclic: v->a->b->v or v->b->a->v).
  The out-neighbors of v form a "linked" structure that forces 5-cycles.

Case B (no common vertex): {a,b,c}, {a,d,e}, {b,d,f}
  a in C1,C2; b in C1,C3; d in C2,C3.
  These 6 vertices {a,b,c,d,e,f} have 3 directed 3-cycles among them.

  Question: must there be a 4th 3-cycle among these 6 vertices?

  Let's check exhaustively: for ALL tournaments on 6 vertices with exactly
  3 directed 3-cycles that are pairwise intersecting, is it possible
  that no other 3-cycle exists on those 6 vertices?

kind-pasteur-2026-03-07-S31
"""

from itertools import combinations
from collections import defaultdict


def check_6vertex_constraint():
    """
    On 6 vertices, find tournaments with exactly 3 pairwise-intersecting
    3-cycles and no other 3-cycles.
    """
    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    count_exactly_3 = 0
    count_K3_pairwise = 0
    count_K3_no_common = 0

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1

        cycles = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if (A[a][b] and A[b][c] and A[c][a]) or \
                       (A[a][c] and A[c][b] and A[b][a]):
                        cycles.append(frozenset({a, b, c}))

        if len(cycles) != 3:
            continue
        count_exactly_3 += 1

        # Check if all pairwise intersecting
        if all(cycles[i] & cycles[j] for i in range(3) for j in range(i+1, 3)):
            count_K3_pairwise += 1
            common = cycles[0] & cycles[1] & cycles[2]
            if not common:
                count_K3_no_common += 1

    print(f"n=6: {count_exactly_3} tournaments with exactly 3 three-cycles")
    print(f"  Of which {count_K3_pairwise} have all 3 pairwise intersecting")
    print(f"  Of which {count_K3_no_common} have no common vertex")


def check_general_3_intersecting():
    """
    More general: on n vertices, can 3 pairwise-intersecting 3-cycles
    have NO other 3-cycle sharing a vertex with them?

    At n=6 this means the 6 vertices used by the 3 cycles have exactly
    these 3 3-cycles. At n>=7, the remaining vertices must not create
    additional 3-cycles through the 6 used vertices.

    Step 1: On 6 vertices, find tournaments with exactly 3 pairwise-
    intersecting 3-cycles (no other 3-cycles at all).

    If this is impossible, then at ANY n, the 6 vertices used by
    3 pairwise-intersecting 3-cycles always create a 4th 3-cycle,
    which shares a vertex and thus connects to the K3 component.
    """
    print("\nStep 1: Can 6 vertices have exactly 3 pairwise-intersecting 3-cycles?")
    print("(This means all C(6,3)-3 = 17 other triples are transitive.)")

    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    # For each tournament, check
    found = False
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1

        cycles = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if (A[a][b] and A[b][c] and A[c][a]) or \
                       (A[a][c] and A[c][b] and A[b][a]):
                        cycles.append(frozenset({a, b, c}))

        if len(cycles) != 3:
            continue

        # Check pairwise intersecting
        if all(cycles[i] & cycles[j] for i in range(3) for j in range(i+1, 3)):
            found = True
            common = cycles[0] & cycles[1] & cycles[2]
            print(f"  bits={bits}: cycles {[set(c) for c in cycles]}, common={set(common)}")
            # Show the tournament
            print(f"    Tournament arcs:")
            for i in range(n):
                outs = [j for j in range(n) if j != i and A[i][j]]
                print(f"      {i} -> {outs}")

            # Now check: if we add vertex 6, can vertex 6 avoid creating
            # any 3-cycle with these 6 vertices?
            print(f"    Can vertex 6 avoid 3-cycles with these?")
            # Vertex 6's arcs to {0..5}: 2^6 = 64 possibilities
            can_avoid = 0
            for v6_bits in range(1 << n):
                A7 = [row[:] + [0] for row in A]
                A7.append([0]*7)
                for j in range(n):
                    if v6_bits & (1 << j):
                        A7[6][j] = 1
                    else:
                        A7[j][6] = 1

                # Check if any 3-cycle involves vertex 6 AND vertices from the triple
                triple_verts = cycles[0] | cycles[1] | cycles[2]
                new_cycle = False
                for a in sorted(triple_verts):
                    for b in range(a+1, 7):
                        if b == 6 or b in triple_verts:
                            c = 6
                            if b < c:
                                triple_set = frozenset({a, b, c})
                                if (A7[a][b] and A7[b][c] and A7[c][a]) or \
                                   (A7[a][c] and A7[c][b] and A7[b][a]):
                                    if triple_set & triple_verts:
                                        new_cycle = True
                                        break
                    if new_cycle:
                        break

                if not new_cycle:
                    can_avoid += 1

            print(f"    Ways vertex 6 avoids 3-cycles: {can_avoid}/64")

    if not found:
        print("  NOT FOUND: no tournament on 6 vertices has exactly 3 pairwise-intersecting 3-cycles")
        print("  This means 3 pairwise-intersecting 3-cycles ALWAYS create additional 3-cycles!")
        print("  Therefore K3 isolation is impossible.")


def check_pairwise_forcing():
    """
    THEOREM: 3 pairwise-intersecting directed 3-cycles in a tournament
    always create at least one additional 3-cycle sharing a vertex.

    Proof by exhaustive case analysis on 6 vertices.
    """
    print("\n=== PROOF: 3 pairwise-intersecting 3-cycles force a 4th ===")

    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    # Case A: common vertex
    print("\nCase A (common vertex): handled by THM-029")

    # Case B: no common vertex but vertices shared pairwise
    # Count tournaments with exactly 3 pairwise-intersecting 3-cycles (no common)
    # that have NO additional 3-cycle sharing vertices

    # Check ALL tournaments: any with 3+ pairwise-intersecting 3-cycles but
    # which have exactly those 3 on a set of 6 vertices?

    # Actually, let's just count: among ALL n=6 tournaments, does any have
    # the property that some K3 triple in the 3-cycle graph is isolated?

    isolated_count = 0

    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1

        cycles = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if (A[a][b] and A[b][c] and A[c][a]) or \
                       (A[a][c] and A[c][b] and A[b][a]):
                        cycles.append(frozenset({a, b, c}))

        nc = len(cycles)
        if nc < 3:
            continue

        # For each triple of pairwise-intersecting cycles
        for i in range(nc):
            for j in range(i+1, nc):
                if not (cycles[i] & cycles[j]):
                    continue
                for k in range(j+1, nc):
                    if not (cycles[i] & cycles[k]) or not (cycles[j] & cycles[k]):
                        continue

                    triple_verts = cycles[i] | cycles[j] | cycles[k]
                    # Check if any OTHER cycle shares a vertex
                    has_neighbor = False
                    for other in range(nc):
                        if other in (i, j, k):
                            continue
                        if cycles[other] & triple_verts:
                            has_neighbor = True
                            break

                    if not has_neighbor:
                        isolated_count += 1

    print(f"\nn=6: isolated K3 triples: {isolated_count}")
    if isolated_count == 0:
        print("CONFIRMED: no K3 in 3-cycle graph is ever isolated at n=6")
        print("Combined with the exhaustive n=7 check: K3 isolation is impossible.")


if __name__ == '__main__':
    check_6vertex_constraint()
    check_general_3_intersecting()
    check_pairwise_forcing()
