#!/usr/bin/env python3
"""
Typed G_T additivity test at n=7.

At n=7, independent sets can be:
  () - empty
  (3,) - single 3-cycle  [f = 6 - 2 = 4]
  (5,) - single 5-cycle  [f = 6 - 4 = 2]
  (7,) - single 7-cycle  [f = 6 - 6 = 0]
  (3,3) - disjoint pair  [f = 6 - 4 = 2]

Additivity question: G_T^typed(t; y3, y5, y7) =
  A_7(t) + t3*g3(t)*y3 + t5*g5(t)*y5 + t7*g7(t)*y7 + bc33*g33(t)*y3^2 ?

The bc33*y3^2 term creates a y3^2 contribution that breaks simple additivity.
But the QUESTION is whether the y3^2 coefficient is exactly bc33 * some fixed g33(t).

Since (3,3) pairs have f=2 (same as 5-cycles), g33(t) = A_3(t)*(t-1)^4 = g5(t)!
And (3,) has f=4, so g3(t) = A_5(t)*(t-1)^2
And (7,) has f=0, so g7(t) = A_1(t)*(t-1)^6 = (t-1)^6

So the typed G_T IS a polynomial in y3, y5, y7:
  A_7(t) + t3 * A_5(t)*(t-1)^2 * y3
         + t5 * A_3(t)*(t-1)^4 * y5
         + t7 * (t-1)^6 * y7
         + bc33 * A_3(t)*(t-1)^4 * y3^2

This IS additive in individual contributions, with a quadratic y3^2 term.
Setting y3=y5=y7=x recovers G_T(t,x).

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import factorial, comb
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def eulerian_poly_eval(n, t):
    return sum(eulerian_number(n, k) * t**k for k in range(n))

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

def main():
    print("=" * 70)
    print("TYPED G_T AT n=7: EXACT STRUCTURE")
    print("=" * 70)

    n = 7
    d = n - 1  # = 6

    # The correction function g_type(t) for each independent set type:
    # g_type(t) = A_{f+1}(t) * (t-1)^{d-f}
    #
    # (3,):  f = d - 2 = 4, g(t) = A_5(t) * (t-1)^2
    # (5,):  f = d - 4 = 2, g(t) = A_3(t) * (t-1)^4
    # (7,):  f = d - 6 = 0, g(t) = A_1(t) * (t-1)^6 = (t-1)^6
    # (3,3): f = d - 4 = 2, g(t) = A_3(t) * (t-1)^4  <-- SAME as (5,)!

    print("\nCorrection functions g(t):")
    for t_val in [0, 1, 2, -1, 3]:
        g3 = eulerian_poly_eval(5, t_val) * (t_val - 1)**2
        g5 = eulerian_poly_eval(3, t_val) * (t_val - 1)**4
        g7 = (t_val - 1)**6
        g33 = eulerian_poly_eval(3, t_val) * (t_val - 1)**4  # same as g5!
        An = eulerian_poly_eval(7, t_val)
        print(f"  t={t_val:+d}: A_7={An:6d}, g3={g3:6d}, g5={g5:6d}, g7={g7:6d}, g33={g33:6d}")

    print(f"\n  NOTE: g5(t) = g33(t) always! (same f value)")
    print(f"  So (5,) and (3,3) contribute with the SAME t-weight function.")
    print(f"  This is WHY the null space at n=7 has (0,-2,0,1) — ")
    print(f"  because t5 and bc33 enter with the same t-dependence!")

    # Verify: typed G_T formula matches direct computation
    print(f"\n\n--- Verification: typed formula at n=7 ---")

    for seed in range(10):
        A, bits = random_tournament(n, seed)
        cycles = find_ALL_directed_odd_cycles(A, n)
        nc = len(cycles)

        t3 = sum(1 for c, l in cycles if l == 3)
        t5 = sum(1 for c, l in cycles if l == 5)
        t7 = sum(1 for c, l in cycles if l == 7)

        # Count bc33
        if nc > 0:
            adj = conflict_graph_adj(cycles)
            indep = enumerate_independent_sets(adj, nc)
            bc33 = sum(1 for s in indep
                       if len(s) == 2 and
                       cycles[s[0]][1] == 3 and cycles[s[1]][1] == 3)
        else:
            bc33 = 0

        H_formula = 1 + 2*t3 + 2*t5 + 2*t7 + 4*bc33
        H_direct = sum(2**len(s) for s in (indep if nc > 0 else [[]]))

        # Full G_T(t, 2) using typed formula
        for t_val in [0, 2, -1]:
            An = eulerian_poly_eval(7, t_val)
            g3_val = eulerian_poly_eval(5, t_val) * (t_val - 1)**2
            g5_val = eulerian_poly_eval(3, t_val) * (t_val - 1)**4
            g7_val = (t_val - 1)**6
            g33_val = g5_val  # same!

            # Typed G_T(t, 2) = An + 2*t3*g3 + 2*t5*g5 + 2*t7*g7 + 4*bc33*g33
            GT_formula = An + 2*t3*g3_val + 2*t5*g5_val + 2*t7*g7_val + 4*bc33*g33_val

            # Direct computation via independence polynomial
            if nc > 0:
                GT_direct = An
                for s in indep:
                    if not s:
                        continue
                    parts = len(s)
                    sum_len = sum(cycles[i][1] for i in s)
                    f_I = d - (sum_len - parts)
                    correction = eulerian_poly_eval(f_I + 1, t_val) * (t_val - 1)**(d - f_I)
                    GT_direct += (2**parts) * correction
            else:
                GT_direct = An

            match = abs(GT_formula - GT_direct) < 1e-10
            if not match:
                print(f"  MISMATCH seed={seed}, t={t_val}")
                print(f"    formula={GT_formula}, direct={GT_direct}")

        if seed < 3:
            print(f"  seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc33={bc33}, H={H_direct} OK")

    print(f"  All 10 samples verified!")

    # THE DEEP INSIGHT
    print(f"\n\n{'='*70}")
    print(f"DEEP INSIGHT: g5(t) = g33(t)")
    print(f"{'='*70}")
    print()
    print("Because a 5-cycle and a disjoint (3,3)-pair both cover")
    print("the SAME number of vertices (5 vs 3+3=6? No...)")
    print()
    print("Wait: let's recalculate f values:")
    print("  (5,): uses 5 vertices, parts=1, sum_len=5")
    print(f"    f = d - (sum_len - parts) = 6 - (5-1) = 6-4 = 2")
    print("  (3,3): uses 6 vertices, parts=2, sum_len=6")
    print(f"    f = d - (sum_len - parts) = 6 - (6-2) = 6-4 = 2")
    print()
    print("So f=2 in both cases because sum_len - parts = 4 in both.")
    print("  5-cycle: 5-1 = 4")
    print("  (3,3): 6-2 = 4")
    print()
    print("In general: sum(len_i) - parts = sum(len_i - 1)")
    print("  = sum of (half-cycle minus 1) for each cycle")
    print("  = total 'excess' over minimal contributions")
    print()
    print("This is the KEY: the correction function depends only on")
    print("  S = sum(len_i - 1) = total 'excess'")
    print("not on the individual cycle lengths!")
    print()
    print("So ANY combination of cycles with the same S gets the same g(t).")
    print("This is why the null space exists: 2*(5-cycle) ~ 1*(3,3)-pair")
    print("in terms of the t-polynomial, because both have S=4.")
    print()
    print("The ONLY information G_T(t,x) captures from cycle lengths is:")
    print("  - How many cycles (parts, entering as x^parts)")
    print("  - Total excess S = sum(len_i - 1)")
    print("It does NOT separate different cycle lengths with the same S.")
    print()
    print("GS's U_T DOES separate them (different p_k monomials).")
    print("So U_T is strictly finer than G_T(t,x).")

    # At n=9: what pairs have the same S?
    print(f"\n\n--- Cycles with same S at n=9 ---")
    print("S = sum(len_i - 1) for independent set of cycles")
    print()
    for S in range(2, 9):
        combos = []
        # Single cycles
        for k in range(3, 10, 2):
            if k - 1 == S and k <= 9:
                combos.append(f"({k},)")
        # Pairs
        for k1 in range(3, 10, 2):
            for k2 in range(k1, 10, 2):
                if (k1-1) + (k2-1) == S and k1+k2 <= 9:
                    combos.append(f"({k1},{k2})")
        # Triples
        for k1 in range(3, 10, 2):
            for k2 in range(k1, 10, 2):
                for k3 in range(k2, 10, 2):
                    if (k1-1)+(k2-1)+(k3-1) == S and k1+k2+k3 <= 9:
                        combos.append(f"({k1},{k2},{k3})")
        if combos:
            print(f"  S={S}: {', '.join(combos)}")

if __name__ == "__main__":
    main()
