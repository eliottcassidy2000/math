#!/usr/bin/env python3
"""
Higher-order IE terms at n=5 (exhaustive).

At n=5, through-v cycles are all 3-cycles and possibly 5-cycles.
Two 3-cycles through v can share at most 1 other vertex.
A 3-cycle and 5-cycle through v always share v, so adjacent in Omega.
Two 3-cycles sharing v are adjacent if they share another vertex.

Key question: are there pairs of non-adjacent through-v cycles at n=5?

opus-2026-03-07-S36
"""
from itertools import permutations, combinations

def tournament_from_bits(n, bits_int):
    m = n*(n-1)//2
    b = [(bits_int >> k) & 1 for k in range(m)]
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if b[idx] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def find_directed_odd_cycles(A, n, vertices=None):
    if vertices is None:
        vertices = list(range(n))
    cycles = set()
    for length in range(3, len(vertices)+1, 2):
        for subset in combinations(vertices, length):
            first = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (first,) + perm
                if all(A[cycle[i]][cycle[(i+1) % length]] for i in range(length)):
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    cycles.add(rotated)
    return list(cycles)

def main():
    n = 5
    m = n*(n-1)//2
    print(f"=== Higher-order IE at n={n} (exhaustive) ===\n")

    nonzero_higher = 0
    zero_higher = 0
    total = 0

    for bits_int in range(2**m):
        A = tournament_from_bits(n, bits_int)
        cycles = find_directed_odd_cycles(A, n)
        nc = len(cycles)
        cvsets = [frozenset(c) for c in cycles]

        adj = [[False]*nc for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if cvsets[i] & cvsets[j]:
                    adj[i][j] = adj[j][i] = True

        for v in range(n):
            cycles_v = [i for i in range(nc) if v in cvsets[i]]
            nv = len(cycles_v)
            if nv <= 1:
                total += 1
                zero_higher += 1
                continue

            # Check: are any PAIRS of through-v cycles non-adjacent?
            has_nonadj = False
            for a in range(len(cycles_v)):
                for b in range(a+1, len(cycles_v)):
                    if not adj[cycles_v[a]][cycles_v[b]]:
                        has_nonadj = True
                        break
                if has_nonadj:
                    break

            if not has_nonadj:
                # All through-v cycles are pairwise adjacent -> form a clique
                # Higher-order IE terms are all zero (can't pick 2+ non-adjacent)
                total += 1
                zero_higher += 1
            else:
                total += 1
                nonzero_higher += 1
                if nonzero_higher <= 5:
                    print(f"  bits={bits_int}, v={v}: non-adjacent through-v pair found!")
                    for ci in cycles_v:
                        print(f"    C: {cycles[ci]} (verts={cvsets[ci]})")
                    for a in range(len(cycles_v)):
                        for b in range(a+1, len(cycles_v)):
                            if not adj[cycles_v[a]][cycles_v[b]]:
                                print(f"    Non-adj: {cycles[cycles_v[a]]} -- {cycles[cycles_v[b]]}")

    print(f"\nTotal (T,v) pairs: {total}")
    print(f"Zero higher-order IE: {zero_higher}")
    print(f"Nonzero higher-order IE: {nonzero_higher}")

    if nonzero_higher == 0:
        print("\n*** ALL through-v cycles form a CLIQUE in Omega(T) at n=5 ***")
        print("*** This means higher-order IE terms are TRIVIALLY zero ***")
        print("*** (No pair of non-adjacent through-v cycles exists) ***")
    else:
        print(f"\n{nonzero_higher} cases with non-adjacent through-v pairs")
        print("Need to check if their IE contributions still cancel")

if __name__ == "__main__":
    main()
