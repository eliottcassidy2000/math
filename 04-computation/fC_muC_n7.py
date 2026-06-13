#!/usr/bin/env python3
"""
Test f(C) = 2*mu(C) at n=7.

f(C) = sum_{S: C in S, S indep in Omega(T)} 2^|S|
     = 2 * I(Omega(T) - N[C], 2)

mu(C) = I(Omega(T-v) |_{avoid C\{v}}, 2)

If f(C) = 2*mu(C), then Claim A follows trivially from the
partition of independent sets into A (not through v) and B (through v).

At n=5, this holds everywhere. Does it hold at n=7?

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import defaultdict
import random

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n, seed):
    rng = random.Random(seed)
    m = n*(n-1)//2
    bits = [rng.randint(0,1) for _ in range(m)]
    return tournament_from_bits(n, bits)

def find_directed_odd_cycles(A, n, vertices=None):
    """Find all directed odd cycles using the faster method."""
    if vertices is None:
        vertices = list(range(n))
    cycles = set()
    nv = len(vertices)
    for length in range(3, nv+1, 2):
        for subset in combinations(vertices, length):
            # Fix first vertex, permute rest
            first = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (first,) + perm
                if all(A[cycle[i]][cycle[(i+1)%length]] for i in range(length)):
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    cycles.add(rotated)
    return list(cycles)

def I_at_2(cycles, vertices_used=None):
    """Compute I(conflict graph, 2) = sum_{S indep} 2^|S|."""
    nc = len(cycles)
    if nc == 0:
        return 1  # just the empty set

    # Build adjacency
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = True

    total = 0
    for mask in range(2**nc):
        nodes = [i for i in range(nc) if (mask >> i) & 1]
        ok = True
        for a in range(len(nodes)):
            for b in range(a+1, len(nodes)):
                if adj[nodes[a]][nodes[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += 2**len(nodes)
    return total

def main():
    n = 7
    print(f"--- Testing f(C) = 2*mu(C) at n={n} ---\n")

    match_count = 0
    fail_count = 0

    for seed in range(30):
        A = random_tournament(n, seed)
        cycles_T = find_directed_odd_cycles(A, n)
        nc = len(cycles_T)

        # Build conflict graph
        adj = [[False]*nc for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if set(cycles_T[i]) & set(cycles_T[j]):
                    adj[i][j] = adj[j][i] = True

        H = I_at_2(cycles_T)

        for v in range(n):
            remaining = [u for u in range(n) if u != v]
            cycles_Tv = find_directed_odd_cycles(A, n, remaining)

            cycles_v_idx = [i for i in range(nc) if v in cycles_T[i]]

            for ci in cycles_v_idx:
                C = cycles_T[ci]

                # f(C) = 2 * I(Omega(T) - N[C], 2)
                # N[C] = {C} union {cycles sharing vertex with C}
                non_conflict = [j for j in range(nc) if j != ci and not adj[ci][j]]
                f_val = 2 * I_at_2([cycles_T[j] for j in non_conflict])

                # mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)
                forbidden = set(C) - {v}
                allowed_Tv = [c for c in cycles_Tv if not (set(c) & forbidden)]
                mu_val = I_at_2(allowed_Tv)

                if f_val == 2 * mu_val:
                    match_count += 1
                else:
                    fail_count += 1
                    if fail_count <= 5:
                        print(f"  FAIL seed={seed}, v={v}, C={C}")
                        print(f"    f(C)={f_val}, 2*mu(C)={2*mu_val}")
                        print(f"    Omega-N[C] has {len(non_conflict)} cycles")
                        print(f"    Omega(T-v)|_avoid has {len(allowed_Tv)} cycles")

                        # What is the discrepancy?
                        nc_omega = [cycles_T[j] for j in non_conflict]
                        print(f"    I(Omega-N[C], 2) = {I_at_2(nc_omega)}")
                        print(f"    I(Omega(T-v)|_avoid, 2) = {mu_val}")

        if seed < 5 or fail_count > 0:
            print(f"  seed={seed}: H={H}, #cycles={nc}, "
                  f"match={match_count}, fail={fail_count}")

    print(f"\n--- TOTAL ---")
    print(f"  Match: {match_count}")
    print(f"  Fail:  {fail_count}")

    if fail_count == 0:
        print("\n  f(C) = 2*mu(C) HOLDS AT n=7!")
        print("  This means Claim A follows from the simple partition argument!")
    else:
        print(f"\n  f(C) != 2*mu(C) in general at n=7")
        print("  The inclusion-exclusion does NOT trivially collapse.")
        print("  Need to understand the higher-order terms.")

if __name__ == "__main__":
    main()
