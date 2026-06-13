#!/usr/bin/env python3
"""
Alpha_0 = 1 and Claim A: the empty-set contribution cancels in deletion.

Claim A: H(T) - H(T-v) = 2 * sum_{C ni v} mu(C)

Using OCF: H(T) = I(Omega(T), 2) = sum_{S indep} 2^|S|

The empty set contributes 1 to both H(T) and H(T-v).
So the difference H(T) - H(T-v) comes entirely from CYCLES.

More precisely:
  H(T) - H(T-v) = [sum_{S indep in Omega(T)} 2^|S|] - [sum_{S indep in Omega(T-v)} 2^|S|]

Since alpha_0 = 1 in both:
  = sum_{S non-empty indep in Omega(T)} 2^|S| - sum_{S non-empty indep in Omega(T-v)} 2^|S|

Every non-empty independent set S in Omega(T) either:
  (a) contains no cycle through v -> S is also in Omega(T-v) [if the cycles survive]
  (b) contains at least one cycle through v -> S is NOT in Omega(T-v)

Wait, this is NOT quite right because:
- Cycles in Omega(T) that don't involve v might not exist in Omega(T-v)
  Actually NO: if a cycle in T doesn't involve v, it's still a cycle in T-v.
  So cycles NOT through v are preserved.
- The CONFLICT graph might change: two cycles might become
  non-conflicting after removing v. But conflict means sharing a vertex,
  and if neither cycle involves v, their vertex sets are unchanged.

So cycles not through v and their conflict relationships are UNCHANGED.

This means:
  I(Omega(T), x) - I(Omega(T-v), x) = [contribution of independent sets
  containing at least one cycle through v] - [contribution of independent sets
  of T-v that are NOT independent sets of T]

The second term is 0: every indep set of T-v is also an indep set of T
(fewer cycles means fewer conflicts, so more independent sets if anything...
 wait, T-v has FEWER cycles, so Omega(T-v) is a SUBGRAPH of Omega(T).
 Every indep set of Omega(T-v) corresponds to a set of cycles in T-v,
 which are also cycles in T, and they're still vertex-disjoint, so still
 independent in Omega(T). So indep sets of Omega(T-v) INJECT into Omega(T).)

Actually wait: Omega(T-v) has cycles of T-v, which are cycles in T that
DON'T pass through v. So Omega(T-v) is the INDUCED subgraph of Omega(T)
on cycles not through v. An independent set in this induced subgraph is
also independent in Omega(T).

Therefore:
  I(Omega(T), x) = I(Omega(T-v), x)  [as subgraph contribution]
                 + sum_{S: S ni some cycle through v} x^|S| * [something]

More precisely: partition independent sets of Omega(T) into:
  A = those using NO cycle through v
  B = those using >= 1 cycle through v

Then sum_A 2^|S| = I(Omega(T-v), 2) = H(T-v)  [since A = indep sets of T-v]
And sum_B 2^|S| = H(T) - H(T-v)

Claim A says this equals 2 * sum_{C ni v} mu(C).

VERIFY THIS AT n=5.

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb

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

def find_ALL_directed_odd_cycles(A, n, vertices=None):
    if vertices is None:
        vertices = list(range(n))
    all_cycles = []
    nv = len(vertices)
    for length in range(3, nv+1, 2):
        for subset in combinations(vertices, length):
            for perm in permutations(subset):
                is_cycle = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1)%length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = list(perm).index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    all_cycles.append(canon)
    return list(set(all_cycles))

def main():
    n = 5
    m = n*(n-1)//2
    print(f"--- Verifying alpha_0 cancellation for Claim A at n={n} ---\n")

    total_checks = 0
    all_ok = True

    for bits_int in range(2**m):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)

        # Get cycles of T
        cycles_T = find_ALL_directed_odd_cycles(A, n)
        nc = len(cycles_T)

        # Build conflict graph
        adj = [[False]*nc for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if set(cycles_T[i]) & set(cycles_T[j]):
                    adj[i][j] = adj[j][i] = True

        # Find all independent sets
        indep_T = []
        for mask in range(2**nc):
            nodes = [i for i in range(nc) if (mask >> i) & 1]
            ok = True
            for a in range(len(nodes)):
                for b_idx in range(a+1, len(nodes)):
                    if adj[nodes[a]][nodes[b_idx]]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                indep_T.append(frozenset(nodes))

        H_T = sum(2**len(s) for s in indep_T)

        for v in range(n):
            # Partition independent sets: A (no cycle through v), B (has cycle through v)
            remaining = [u for u in range(n) if u != v]
            cycles_through_v = {i for i in range(nc) if v in cycles_T[i]}
            cycles_not_v = {i for i in range(nc) if v not in cycles_T[i]}

            A_sets = [s for s in indep_T if not (s & cycles_through_v)]
            B_sets = [s for s in indep_T if (s & cycles_through_v)]

            sum_A = sum(2**len(s) for s in A_sets)
            sum_B = sum(2**len(s) for s in B_sets)

            # Verify sum_A = H(T-v)
            cycles_Tv = find_ALL_directed_odd_cycles(A, n, remaining)
            nc_v = len(cycles_Tv)
            if nc_v == 0:
                H_Tv = 1
            else:
                adj_v = [[False]*nc_v for _ in range(nc_v)]
                for i in range(nc_v):
                    for j in range(i+1, nc_v):
                        if set(cycles_Tv[i]) & set(cycles_Tv[j]):
                            adj_v[i][j] = adj_v[j][i] = True
                indep_v = []
                for mask in range(2**nc_v):
                    nodes = [i for i in range(nc_v) if (mask >> i) & 1]
                    ok = True
                    for a in range(len(nodes)):
                        for b_idx in range(a+1, len(nodes)):
                            if adj_v[nodes[a]][nodes[b_idx]]:
                                ok = False
                                break
                        if not ok:
                            break
                    if ok:
                        indep_v.append(frozenset(nodes))
                H_Tv = sum(2**len(s) for s in indep_v)

            if sum_A != H_Tv:
                print(f"FAIL: bits={bits_int}, v={v}: sum_A={sum_A}, H(T-v)={H_Tv}")
                all_ok = False

            if sum_B != H_T - H_Tv:
                print(f"FAIL: bits={bits_int}, v={v}: sum_B={sum_B}, H-H_v={H_T-H_Tv}")
                all_ok = False

            total_checks += 1

    print(f"Total checks: {total_checks}")
    if all_ok:
        print("ALL PASS!")
        print()
        print("CONFIRMED: The empty set cancellation is exact.")
        print("  sum_{S not through v} 2^|S| = H(T-v)")
        print("  sum_{S through v} 2^|S| = H(T) - H(T-v)")
        print()
        print("This means Claim A reduces to:")
        print("  sum_{S with >= 1 cycle through v} 2^|S| = 2 * sum_{C ni v} mu(C)")
        print()
        print("The alpha_0 = 1 contribution CANCELS in the deletion.")
        print("Only cycles through v matter for the difference.")

if __name__ == "__main__":
    main()
