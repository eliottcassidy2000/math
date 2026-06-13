#!/usr/bin/env python3
"""
GS bridge at n=7: typed independence polynomial with disjoint cycle pairs.

At n=7, the cycle types in independent sets are:
  () - empty set
  (3,) - single 3-cycle
  (5,) - single 5-cycle
  (7,) - single 7-cycle (Hamiltonian cycle)
  (3,3) - disjoint pair of 3-cycles (using 6 of 7 vertices)

The typed IP has 5 terms. Our OCF has 4 invariants (t3, t5, t7, bc33).
The typed IP is FINER: it distinguishes 3-cycles from 5-cycles from 7-cycles,
while our alpha_1 merges them.

Question: does U_T at n=7 exactly match the typed IP?
We can't do exhaustive n=7 (over 2M tournaments), but we sample.

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import factorial
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
    return tournament_from_bits(n, bits), bits

def find_ALL_directed_odd_cycles(A, n):
    all_cycles = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            for perm in permutations(subset):
                is_cycle = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1)%length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = list(perm).index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    all_cycles.append((canon, length))
    seen = set()
    result = []
    for c, l in all_cycles:
        if c not in seen:
            seen.add(c)
            result.append((c, l))
    return result

def conflict_graph_adj(cycles):
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(cycles[i][0]) & set(cycles[j][0]):
                adj[i][j] = adj[j][i] = True
    return adj

def enumerate_independent_sets(adj, nc):
    result = []
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
            result.append(nodes)
    return result

def typed_indep_poly(cycles):
    nc = len(cycles)
    if nc == 0:
        return {(): 1}
    adj = conflict_graph_adj(cycles)
    indep = enumerate_independent_sets(adj, nc)
    result = defaultdict(int)
    for s in indep:
        lengths = tuple(sorted([cycles[i][1] for i in s], reverse=True))
        result[lengths] += 1
    return dict(result)

def compute_UT(A, n):
    """Compute U_T expansion — SLOW for n=7."""
    T_bar = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                T_bar[i][j] = 1 - A[i][j]

    result = defaultdict(int)
    for perm in permutations(range(n)):
        visited = [False]*n
        cycles_perm = []
        for i in range(n):
            if not visited[i]:
                cycle = []
                j = i
                while not visited[j]:
                    visited[j] = True
                    cycle.append(j)
                    j = perm[j]
                cycles_perm.append(tuple(cycle))

        valid = True
        phi = 0
        for cyc in cycles_perm:
            k = len(cyc)
            if k == 1:
                continue
            is_T = all(A[cyc[i]][cyc[(i+1)%k]] == 1 for i in range(k))
            is_Tbar = all(T_bar[cyc[i]][cyc[(i+1)%k]] == 1 for i in range(k))
            if is_T:
                phi += k - 1
            elif is_Tbar:
                pass
            else:
                valid = False
                break

        if valid:
            cycle_type = tuple(sorted([len(c) for c in cycles_perm], reverse=True))
            result[cycle_type] += (-1)**phi

    return dict(result)

def reconstruct_UT_from_typed(tip, n):
    """Reconstruct U_T from typed independence polynomial."""
    result = defaultdict(int)
    for cycle_lengths, count in tip.items():
        s = len(cycle_lengths)
        remaining = n - sum(cycle_lengths)
        partition = tuple(sorted(list(cycle_lengths) + [1]*remaining, reverse=True))
        coeff = count * (2**s)
        result[partition] += coeff
    return dict(result)

def main():
    print("=" * 70)
    print("GS BRIDGE AT n=7: TYPED INDEPENDENCE POLYNOMIAL")
    print("=" * 70)

    n = 7
    print(f"\nSampling {20} random tournaments at n={n}...")
    print(f"(U_T computation takes ~{factorial(7)} permutations per tournament)")

    all_match = True
    for seed in range(20):
        A, bits = random_tournament(n, seed)

        # Compute typed IP
        cycles = find_ALL_directed_odd_cycles(A, n)
        tip = typed_indep_poly(cycles)
        H = sum(2**len(ls) * c for ls, c in tip.items())

        # Reconstruct expected U_T
        expected_ut = reconstruct_UT_from_typed(tip, n)

        # Compute actual U_T
        actual_ut = compute_UT(A, n)

        # Compare
        match = True
        all_parts = set(list(expected_ut.keys()) + list(actual_ut.keys()))
        for part in all_parts:
            if expected_ut.get(part, 0) != actual_ut.get(part, 0):
                match = False
                break

        if not match:
            all_match = False
            print(f"\n  MISMATCH seed={seed}")
            print(f"    Typed IP: {tip}")
            print(f"    Expected: {dict(sorted(expected_ut.items()))}")
            print(f"    Actual:   {dict(sorted(actual_ut.items()))}")
        else:
            # Show type composition
            t3 = sum(1 for c, l in cycles if l == 3)
            t5 = sum(1 for c, l in cycles if l == 5)
            t7 = sum(1 for c, l in cycles if l == 7)
            bc33 = tip.get((3,3), 0)
            print(f"  seed={seed:2d}: H={H:3d}, t3={t3}, t5={t5}, t7={t7}, bc33={bc33}, "
                  f"tip_types={list(tip.keys())} OK")

    if all_match:
        print(f"\n  ALL {20} MATCH! GS bridge verified at n=7.")
        print(f"\n  THEOREM (confirmed n=3,5,7):")
        print(f"    U_T = sum_{{S indep in Omega}} 2^|S| * prod_{{c in S}} p_{{len(c)}} * p_1^{{n-sum len}}")
        print(f"\n  This means:")
        print(f"    - coeff of p_3*p_1^4 in U_T = 2 * #{'{'}3-cycles{'}'}")
        print(f"    - coeff of p_5*p_1^2 in U_T = 2 * #{'{'}5-cycles{'}'}")
        print(f"    - coeff of p_7 in U_T = 2 * #{'{'}7-cycles (Hamiltonian){'}'}")
        print(f"    - coeff of p_3^2*p_1 in U_T = 4 * #{'{'}disjoint 3,3-pairs{'}'}")
    else:
        print(f"\n  MISMATCHES FOUND")

    # Analysis: what GS adds beyond our OCF
    print(f"\n\n--- WHAT GS ADDS BEYOND OCF ---")
    print(f"Our OCF uses alpha_k = #{'{'}ind. sets of size k in Omega{'}'}")
    print(f"GS uses I_typed which separates by cycle lengths.")
    print(f"At n=7:")
    print(f"  alpha_1 = t3 + t5 + t7  (all single cycles)")
    print(f"  alpha_2 = bc33           (only (3,3) pairs at n=7)")
    print(f"  alpha_0 = 1              (empty set)")
    print(f"")
    print(f"So I(Omega, x) = 1 + (t3+t5+t7)*x + bc33*x^2")
    print(f"But I_typed(y3,y5,y7) = 1 + t3*y3 + t5*y5 + t7*y7 + bc33*y3^2")
    print(f"")
    print(f"I(Omega, x) = I_typed(x, x, x)")
    print(f"")
    print(f"KEY: at n>=9, we get (3,5) pairs, (3,7) pairs, (5,5) pairs, etc.")
    print(f"The typed IP always refines the standard IP.")
    print(f"But our G_T(t,x) uses the STANDARD IP (all cycles weighted equally by x).")
    print(f"Could we define a typed G_T?")
    print(f"  G_T(t, y3, y5, y7, ...) = A_n(t) + sum_I prod y_{'{'}len{'}'} * I(T) * A_f(t) * (t-1)^d")
    print(f"This would be an INFINITE-variable generating function!")

if __name__ == "__main__":
    main()
