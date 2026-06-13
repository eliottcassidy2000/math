#!/usr/bin/env python3
"""
Test f(C) = 2*mu(C) at n=7 using FAST cycle finding.

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
    """Fast cycle finding: fix first vertex, permute rest."""
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

def I_at_2_fast(cycles_list):
    """Compute I(conflict graph, 2) using divide-and-conquer."""
    nc = len(cycles_list)
    if nc == 0:
        return 1
    if nc > 20:
        # Brute force for small
        pass

    # Build adjacency
    cvsets = [frozenset(c) for c in cycles_list]
    adj_dict = {}
    for i in range(nc):
        adj_dict[i] = frozenset(j for j in range(nc) if j != i and cvsets[i] & cvsets[j])

    # Divide and conquer independence polynomial
    memo = {}
    def solve(verts):
        if verts in memo:
            return memo[verts]
        if not verts:
            return 1
        v = max(verts, key=lambda u: len(adj_dict[u] & verts))
        p1 = solve(verts - {v})
        p2 = solve(verts - (adj_dict[v] & verts) - {v})
        result = p1 + 2 * p2
        memo[verts] = result
        return result

    return solve(frozenset(range(nc)))

def main():
    n = 7
    print(f"--- f(C) = 2*mu(C) test at n={n} ---\n")

    match_count = 0
    fail_count = 0
    total_cycles = 0

    for seed in range(20):
        A = tournament_from_seed(n, seed)
        cycles_T = find_directed_odd_cycles(A, n)
        nc = len(cycles_T)
        cvsets_T = [frozenset(c) for c in cycles_T]

        H = I_at_2_fast(cycles_T)

        for v in range(n):
            remaining = [u for u in range(n) if u != v]
            cycles_Tv = find_directed_odd_cycles(A, n, remaining)

            for ci in range(nc):
                if v not in cvsets_T[ci]:
                    continue
                C = cycles_T[ci]
                total_cycles += 1

                # f(C) = 2 * I(Omega(T) - N[C], 2)
                non_conflict = [cycles_T[j] for j in range(nc)
                                if j != ci and not (cvsets_T[ci] & cvsets_T[j])]
                f_val = 2 * I_at_2_fast(non_conflict)

                # mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)
                forbidden = set(C) - {v}
                allowed = [c for c in cycles_Tv if not (frozenset(c) & forbidden)]
                mu_val = I_at_2_fast(allowed)

                if f_val == 2 * mu_val:
                    match_count += 1
                else:
                    fail_count += 1
                    if fail_count <= 10:
                        print(f"  FAIL seed={seed}, v={v}, C={C} (len={len(C)})")
                        print(f"    f(C)={f_val}, 2*mu(C)={2*mu_val}, diff={f_val - 2*mu_val}")
                        # Debug: what cycles are in the two graphs?
                        print(f"    Omega-N[C]: {len(non_conflict)} cycles, I={I_at_2_fast(non_conflict)}")
                        print(f"    Omega(T-v)|avoid: {len(allowed)} cycles, I={mu_val}")

                        # Are they different graphs?
                        nc_set = set(frozenset(c) for c in non_conflict)
                        allow_set = set(frozenset(c) for c in allowed)
                        only_nc = nc_set - allow_set
                        only_allow = allow_set - nc_set
                        if only_nc:
                            print(f"    Only in Omega-N[C]: {[set(c) for c in only_nc]}")
                        if only_allow:
                            print(f"    Only in Omega(T-v)|avoid: {[set(c) for c in only_allow]}")

        if seed < 5 or fail_count > 0:
            print(f"  seed={seed}: H={H}, #cycles={nc}, match={match_count}, fail={fail_count}")

    print(f"\n--- TOTAL ---")
    print(f"  Tested: {total_cycles} (cycle, vertex) pairs")
    print(f"  Match: {match_count}")
    print(f"  Fail:  {fail_count}")

    if fail_count == 0:
        print("\n  f(C) = 2*mu(C) HOLDS AT n=7!")
    else:
        print(f"\n  f(C) != 2*mu(C) in {fail_count}/{total_cycles} cases")

if __name__ == "__main__":
    main()
