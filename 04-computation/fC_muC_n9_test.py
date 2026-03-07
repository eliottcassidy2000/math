#!/usr/bin/env python3
"""
Test f(C) = 2*mu(C) at n=9 for random tournaments.

For a tournament T on n vertices with conflict graph Omega(T):
  - f(C) = 2 * I(Omega(T) - N[C], 2)
  - mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)
  - I(G, 2) = independence polynomial of G evaluated at x=2

Previously verified at n=5 (exhaustive) and n=7 (20 random tournaments).
"""
from itertools import permutations, combinations
import random
import time
import sys

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
    """Find all directed odd cycles. Fix first vertex, permute rest."""
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
    """Compute I(conflict graph, 2) via branch-and-bound on max-degree vertex."""
    nc = len(cycles_list)
    if nc == 0:
        return 1
    cvsets = [frozenset(c) for c in cycles_list]
    adj_dict = {}
    for i in range(nc):
        adj_dict[i] = frozenset(j for j in range(nc) if j != i and cvsets[i] & cvsets[j])
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
    n = 9
    seeds = list(range(10))
    print(f"--- f(C) = 2*mu(C) test at n={n} ---")
    print(f"Seeds: {seeds}\n")

    total_match = 0
    total_fail = 0
    total_pairs = 0

    for seed in seeds:
        t0 = time.time()
        A = tournament_from_seed(n, seed)

        print(f"Seed {seed}: finding cycles in T...", flush=True)
        t1 = time.time()
        cycles_T = find_directed_odd_cycles(A, n)
        nc = len(cycles_T)
        cvsets_T = [frozenset(c) for c in cycles_T]
        t2 = time.time()
        print(f"  Found {nc} odd cycles in T ({t2-t1:.1f}s)", flush=True)

        # Precompute N[C] structure: for each cycle, which other cycles share a vertex
        adj_full = {}
        for i in range(nc):
            adj_full[i] = frozenset(j for j in range(nc) if j != i and cvsets_T[i] & cvsets_T[j])

        # Precompute cycles in T-v for each vertex v
        print(f"  Finding cycles in T-v for each v...", flush=True)
        t3 = time.time()
        cycles_Tv = {}
        for v in range(n):
            remaining = [u for u in range(n) if u != v]
            cycles_Tv[v] = find_directed_odd_cycles(A, n, remaining)
        t4 = time.time()
        print(f"  Done ({t4-t3:.1f}s)", flush=True)

        match_count = 0
        fail_count = 0
        pairs_this = 0

        for v in range(n):
            for ci in range(nc):
                if v not in cvsets_T[ci]:
                    continue
                C = cycles_T[ci]
                pairs_this += 1

                # f(C) = 2 * I(Omega(T) - N[C], 2)
                # N[C] = C itself + neighbors in conflict graph
                non_neighbor_indices = [j for j in range(nc)
                                        if j != ci and j not in adj_full[ci]]
                non_conflict = [cycles_T[j] for j in non_neighbor_indices]
                f_val = 2 * I_at_2_fast(non_conflict)

                # mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)
                forbidden = set(C) - {v}
                allowed = [c for c in cycles_Tv[v] if not (frozenset(c) & forbidden)]
                mu_val = I_at_2_fast(allowed)

                if f_val == 2 * mu_val:
                    match_count += 1
                else:
                    fail_count += 1
                    print(f"  FAIL v={v}, C={C} (len={len(C)})")
                    print(f"    f(C)={f_val}, 2*mu(C)={2*mu_val}, diff={f_val - 2*mu_val}")
                    if fail_count >= 5:
                        print("  Too many failures, stopping this seed.")
                        break
            if fail_count >= 5:
                break

        t5 = time.time()
        total_match += match_count
        total_fail += fail_count
        total_pairs += pairs_this
        status = "PASS" if fail_count == 0 else "FAIL"
        print(f"  Seed {seed}: {status} — {pairs_this} (C,v) pairs, "
              f"{match_count} match, {fail_count} fail, total {t5-t0:.1f}s\n", flush=True)

    print(f"=== SUMMARY ===")
    print(f"  n = {n}")
    print(f"  Seeds tested: {seeds}")
    print(f"  Total (C,v) pairs: {total_pairs}")
    print(f"  Match: {total_match}")
    print(f"  Fail:  {total_fail}")
    if total_fail == 0:
        print(f"\n  f(C) = 2*mu(C) HOLDS at n={n} for all tested tournaments!")
    else:
        print(f"\n  FAILURES FOUND: {total_fail}/{total_pairs}")

if __name__ == "__main__":
    main()
