#!/usr/bin/env python3
"""
Search for a bijection between Ham paths and 2-colored independent cycle sets.

H(T) = I(Omega(T), 2) means there's a bijection:
  {Ham paths of T} <-> {(U, f) : U ⊂ V(Omega) independent, f: U -> {0,1}}

The RHS is the set of vertex-disjoint odd cycle collections U in T,
each cycle colored 0 or 1. The empty collection gives 1 element (the "base path").

Approach: for each tournament, find all Ham paths and all colored cycle sets,
then look for patterns in the correspondence.

Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import hamiltonian_path_count, find_odd_cycles, all_tournaments
from itertools import permutations, combinations


def get_ham_paths(T):
    """Return all Hamiltonian paths of T as tuples."""
    n = len(T)
    paths = []
    for perm in permutations(range(n)):
        if all(T[perm[k]][perm[k+1]] for k in range(n-1)):
            paths.append(perm)
    return paths


def get_colored_cycle_sets(T):
    """Return all 2-colored independent cycle sets of T.
    Each element is (frozenset of cycles, coloring dict).
    Return as list of (cycle_set, coloring) tuples."""
    cycles = find_odd_cycles(T)
    n_cycles = len(cycles)

    # Find all independent sets in the conflict graph
    # Two cycles conflict iff they share a vertex
    cycle_vertices = [set(c) for c in cycles]

    result = []
    # Empty set
    result.append((frozenset(), {}))

    # Enumerate independent sets by inclusion
    def find_independent_sets(idx, current_set, used_vertices):
        """Find all independent sets starting from index idx."""
        for i in range(idx, n_cycles):
            if cycle_vertices[i].isdisjoint(used_vertices):
                new_set = current_set | {i}
                new_used = used_vertices | cycle_vertices[i]
                # For each independent set, 2^|set| colorings
                yield new_set
                yield from find_independent_sets(i+1, new_set, new_used)

    indep_sets = [frozenset()]  # empty set
    for iset in find_independent_sets(0, frozenset(), set()):
        indep_sets.append(iset)

    # Generate all colored versions
    colored = []
    for iset in indep_sets:
        k = len(iset)
        cycle_list = sorted(iset)
        for mask in range(1 << k):
            coloring = {}
            for idx, ci in enumerate(cycle_list):
                coloring[ci] = (mask >> idx) & 1
            colored.append((iset, coloring))

    return colored, cycles


def analyze_tournament(T):
    """Analyze the bijection structure for a single tournament."""
    n = len(T)
    paths = get_ham_paths(T)
    colored_sets, cycles = get_colored_cycle_sets(T)

    H = len(paths)
    I = len(colored_sets)

    if H != I:
        print(f"  MISMATCH: H={H}, I={I}")
        return

    print(f"  H = I = {H}, #cycles = {len(cycles)}")

    if H <= 1:
        return  # transitive, trivial

    # Analyze: which cycles does each path "interact" with?
    # For each path, compute which odd cycles it "uses" in some sense.

    # One natural encoding: for a path pi, look at which consecutive triples
    # form 3-cycles. More generally, which consecutive (2k+1)-tuples form cycles.

    # Actually, let's look at the TYPE II positions from the OCF recursion.
    # Fix vertex v=0. A path has v at some position p.
    # The path contributes to H(T)-H(T-v) iff... well, the Rédei recursion says
    # H(T) - H(T-v) = 2 * sum_C mu(C) by Claim A.

    # But the bijection should be more direct.
    # Let's look at what happens when we remove cycles from the tournament.

    # For the 3-cycle {a,b,c} (say a->b->c->a), the corresponding
    # 2 colored elements should map to 2 specific Ham paths.

    # Let's categorize paths by their "signature" relative to each cycle.
    for pi in paths:
        # For each cycle, determine the path's relationship to it
        sig = []
        for ci, cycle in enumerate(cycles):
            cv = set(cycle)
            # Find where cycle vertices appear in the path
            positions = {v: k for k, v in enumerate(pi)}
            cycle_positions = sorted(positions[v] for v in cv)
            # Are cycle vertices consecutive in the path?
            consec = all(cycle_positions[k+1] - cycle_positions[k] == 1
                        for k in range(len(cycle_positions)-1))
            sig.append('C' if consec else '.')
        print(f"    path {pi}: cycle_sig={''.join(sig)}")

    if len(cycles) > 0 and H <= 20:
        # Show the colored sets too
        for iset, coloring in colored_sets:
            if len(iset) == 0:
                print(f"    colored: empty set (base)")
            else:
                desc = []
                for ci in sorted(iset):
                    desc.append(f"C{ci}({'R' if coloring[ci] else 'B'})")
                print(f"    colored: {' + '.join(desc)}")


def main():
    print("Bijection analysis for small tournaments")
    print("="*60)

    for n in [3, 4, 5]:
        print(f"\nn = {n}:")
        count = 0
        for T in all_tournaments(n):
            H = hamiltonian_path_count(T)
            if H > 1 and H <= 15:  # interesting cases
                count += 1
                if count > 3:
                    break
                print(f"\n  Tournament (H={H}):")
                # Print adjacency
                for i in range(n):
                    row = ''.join(str(T[i][j]) for j in range(n))
                    print(f"    {row}")
                analyze_tournament(T)


if __name__ == "__main__":
    main()
