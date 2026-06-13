#!/usr/bin/env python3
"""
Investigate WHY f(C) = 2*mu(C).

f(C)/2 = I(Omega(T) \ N[C], 2)
mu(C) = I(Omega(T-v) |_{avoid C\{v}}, 2)

So the identity says these two subgraphs have the same I(G, 2).

Are they actually the SAME graph (isomorphic)?
Or just the same I(G, 2)?

Graph 1: Omega(T) \ N[C]
  = cycles of T that share NO vertex with C

Graph 2: Omega(T-v) |_{avoid C\{v}}
  = cycles of T-v that share NO vertex with C\{v}
  = cycles of T not through v AND not sharing vertex with C\{v}
  = cycles of T not through v AND not sharing vertex with any vertex of C except v
  = cycles of T sharing no vertex with C (since they already don't use v)

Wait! If a cycle D doesn't use v and doesn't share any vertex with C\{v},
then D shares no vertex with C (since D already doesn't use v).

And conversely, if D shares no vertex with C, then D doesn't use v (since v is in C)
and D shares no vertex with C\{v}.

So Graph 1 and Graph 2 have the SAME vertex set:
  {cycles D in T : D shares no vertex with C}

But wait: Graph 1 uses cycles of T, and Graph 2 uses cycles of T-v.
A cycle D in T that doesn't use v is also a cycle of T-v. Conversely,
every cycle of T-v is a cycle of T (since removing v doesn't create new cycles).

SO THEY ARE LITERALLY THE SAME GRAPH!

The graphs Omega(T)\\\\N[C] and Omega(T-v)|_{avoid C\{v}} are IDENTICAL
because:
1. The vertex sets are the same: cycles not sharing any vertex with C
2. These cycles exist in both T and T-v (since they don't use v)
3. The conflict edges are the same (conflict = sharing a vertex, unchanged)

THIS IS WHY f(C) = 2*mu(C) ALWAYS HOLDS!

Let's verify this graph equality computationally.

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
import random

def tournament_from_seed(n, seed):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_directed_odd_cycles(T, n, vertices=None):
    if vertices is None:
        vertices = list(range(n))
    cycles = set()
    nv = len(vertices)
    for length in range(3, nv+1, 2):
        for subset in combinations(vertices, length):
            first = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (first,) + perm
                if all(T[cycle[i]][cycle[(i+1) % length]] for i in range(length)):
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    cycles.add(rotated)
    return list(cycles)

def main():
    print("=" * 70)
    print("PROVING f(C) = 2*mu(C): GRAPH EQUALITY")
    print("=" * 70)

    print("\nCLAIM: For any cycle C through v in tournament T:")
    print("  Omega(T) \\ N[C] = Omega(T-v)|_{avoid C\\{v}}")
    print("  as subgraphs (same vertices, same edges)")
    print()

    # Verify at n=7
    n = 7
    all_equal = True
    total_checks = 0

    for seed in range(20):
        A = tournament_from_seed(n, seed)
        cycles_T = find_directed_odd_cycles(A, n)
        nc = len(cycles_T)

        for v in range(n):
            remaining = [u for u in range(n) if u != v]
            cycles_Tv = find_directed_odd_cycles(A, n, remaining)

            for ci in range(nc):
                if v not in cycles_T[ci]:
                    continue
                C = frozenset(cycles_T[ci])
                total_checks += 1

                # Graph 1: cycles in T not sharing any vertex with C
                graph1_cycles = [c for c in cycles_T if not (frozenset(c) & C)]
                graph1_set = set(c for c in graph1_cycles)

                # Graph 2: cycles in T-v not sharing any vertex with C\{v}
                C_minus_v = C - {v}
                graph2_cycles = [c for c in cycles_Tv if not (frozenset(c) & C_minus_v)]
                graph2_set = set(c for c in graph2_cycles)

                if graph1_set != graph2_set:
                    all_equal = False
                    only1 = graph1_set - graph2_set
                    only2 = graph2_set - graph1_set
                    print(f"  DIFFER seed={seed}, v={v}, C={cycles_T[ci]}")
                    if only1:
                        print(f"    Only in graph1: {[set(c) for c in only1]}")
                    if only2:
                        print(f"    Only in graph2: {[set(c) for c in only2]}")

    print(f"\nChecked {total_checks} (C, v) pairs at n={n}")
    if all_equal:
        print("ALL GRAPHS IDENTICAL!")
        print()
        print("PROOF:")
        print("  Let C be a directed odd cycle through v in T.")
        print("  Let D be any directed odd cycle in T.")
        print()
        print(r"  D in Omega(T) \ N[C]")
        print("    <=> D shares no vertex with C")
        print("    <=> D doesn't use v, AND D shares no vertex with C\\{v}")
        print("    <=> D is a cycle of T-v not sharing vertex with C\\{v}")
        print("    <=> D in Omega(T-v)|_{avoid C\\{v}}")
        print()
        print("  The key step: D shares no vertex with C => D doesn't use v")
        print("  (because v is in C). And D doesn't use v => D is a cycle of T-v")
        print("  (removing v doesn't affect cycles not using v).")
        print()
        print("  Since the vertex sets are identical and the adjacency relation")
        print("  (sharing a vertex) is unchanged, the subgraphs are identical.")
        print()
        print("  Therefore f(C) = 2*I(graph1, 2) = 2*I(graph2, 2) = 2*mu(C). QED")
        print()
        print("  CONSEQUENCE: Claim A follows from Claim B!")
        print("  Claim B: I(Omega(T),2) - I(Omega(T-v),2) = 2*sum_C mu(C)")
        print("  is equivalent to:")
        print("  sum_{S ni some C ni v} 2^|S| = 2*sum_C mu(C)")
        print("  which is: sum_C f(C) - (higher order IE) = 2*sum_C mu(C)")
        print("  Since f(C) = 2*mu(C), this becomes:")
        print("  2*sum_C mu(C) - (higher order IE) = 2*sum_C mu(C)")
        print("  So the higher-order inclusion-exclusion terms MUST vanish!")

if __name__ == "__main__":
    main()
