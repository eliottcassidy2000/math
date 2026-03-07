#!/usr/bin/env python3
"""
Investigate why higher-order inclusion-exclusion terms vanish in Claim A.

We know:
  H(T) - H(T-v) = sum_{S indep, S ni some C ni v} 2^|S|  (from OCF)

By inclusion-exclusion over cycles through v:
  LHS = sum_C f(C) - sum_{C1<C2} f(C1,C2) + sum_{C1<C2<C3} f(C1,C2,C3) - ...

where f(C1,...,Cr) = sum_{S indep: C1,...,Cr all in S} 2^|S|
     = 2^r * I(Omega - N[C1] - ... - N[Cr], 2)
     (if C1,...,Cr are mutually non-adjacent; 0 otherwise)

And we know f(C) = 2*mu(C) for all C (THM-069).

So: sum_C f(C) = 2 * sum_C mu(C) = RHS of Claim A.

This means the HIGHER-ORDER IE terms MUST sum to zero:
  sum_{C1<C2} f(C1,C2) - sum_{C1<C2<C3} f(C1,C2,C3) + ... = 0

WHY? Let's investigate at n=7.

opus-2026-03-07-S36
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

def I_at_2(cycles_list):
    """Compute I(conflict graph, 2) via brute force."""
    nc = len(cycles_list)
    if nc == 0:
        return 1
    cvsets = [frozenset(c) for c in cycles_list]
    total = 0
    for mask in range(2**nc):
        selected = [i for i in range(nc) if (mask >> i) & 1]
        ok = True
        for a in range(len(selected)):
            for b in range(a+1, len(selected)):
                if cvsets[selected[a]] & cvsets[selected[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += 2**len(selected)
    return total

def main():
    n = 7
    print(f"=== Higher-order IE terms at n={n} ===\n")

    for seed in range(10):
        A = tournament_from_seed(n, seed)
        cycles = find_directed_odd_cycles(A, n)
        nc = len(cycles)
        cvsets = [frozenset(c) for c in cycles]

        # Build adjacency
        adj = [[False]*nc for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if cvsets[i] & cvsets[j]:
                    adj[i][j] = adj[j][i] = True

        H = I_at_2(cycles)

        for v in range(min(n, 2)):  # just check v=0,1 for speed
            cycles_v = [i for i in range(nc) if v in cvsets[i]]
            nv = len(cycles_v)

            if nv == 0:
                continue

            # Compute f(C) for each cycle through v
            f1_sum = 0
            for ci in cycles_v:
                non_conflict = [cycles[j] for j in range(nc) if j != ci and not adj[ci][j]]
                f1_sum += 2 * I_at_2(non_conflict)

            # Compute higher-order IE terms
            ie_terms = {}  # ie_terms[r] = sum over r-subsets
            for r in range(2, nv+1):
                ie_r = 0
                for subset in combinations(cycles_v, r):
                    # Check mutual non-adjacency
                    mutual_ok = True
                    for a in range(len(subset)):
                        for b in range(a+1, len(subset)):
                            if adj[subset[a]][subset[b]]:
                                mutual_ok = False
                                break
                        if not mutual_ok:
                            break
                    if not mutual_ok:
                        continue

                    # Non-conflicting with all cycles in subset
                    non_conflict = [cycles[j] for j in range(nc)
                                    if j not in subset and
                                    all(not adj[j][si] for si in subset)]
                    ie_r += (2**r) * I_at_2(non_conflict)
                ie_terms[r] = ie_r

            # Alternating sum of higher orders
            higher_sum = sum((-1)**(r+1) * ie_terms.get(r, 0) for r in range(2, nv+1))

            # Compute actual Claim A quantities
            remaining = [u for u in range(n) if u != v]
            cycles_Tv = find_directed_odd_cycles(A, n, remaining)
            H_Tv = I_at_2(cycles_Tv)
            diff = H - H_Tv

            print(f"seed={seed}, v={v}: #cycles_v={nv}, diff={diff}")
            print(f"  f1_sum={f1_sum} (should = diff = {diff})")
            print(f"  higher_order_IE = {higher_sum}")
            for r in range(2, min(nv+1, 6)):
                if ie_terms.get(r, 0) != 0:
                    print(f"    r={r}: {(-1)**(r+1)} * {ie_terms[r]} = {(-1)**(r+1) * ie_terms[r]}")

            # Key check: does f1_sum = diff + higher_sum?
            # IE says: diff = f1_sum - higher_sum (alternating)
            ie_lhs = f1_sum + sum((-1)**(r) * ie_terms.get(r, 0) for r in range(2, nv+1))
            print(f"  IE reconstruction: {ie_lhs} (should = {diff})")

            if f1_sum == diff:
                print(f"  => f1_sum = diff: higher-order terms EXACTLY CANCEL")
            elif ie_lhs == diff:
                print(f"  => IE reconstruction correct, but f1_sum != diff")
                print(f"  => higher_sum = {f1_sum - diff} (nonzero!)")
            else:
                print(f"  => IE ERROR")
            print()

if __name__ == "__main__":
    main()
