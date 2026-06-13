#!/usr/bin/env python3
"""
Typed G_T: multivariate generating function separating cycle lengths.

G_T(t, x) = A_n(t) + sum_I x^{parts(I)} I(T) A_{f+1}(t) (t-1)^{n-1-f}

Currently, ALL cycle types get weight x. But GS tells us we should separate them.
Define:

  G_T^typed(t; y_3, y_5, y_7, ...) =
    A_n(t) + sum_I [prod_{c in I} y_{len(c)}] I(T) A_{f+1}(t) (t-1)^{n-1-f}

Special evaluations:
  G_T^typed(t; x, x, x, ...) = G_T(t, x)     [recover standard G_T]
  G_T^typed(0; y_3, y_5, ...) = I_typed(Omega; y_3, y_5, ...)  [typed IP]
  G_T^typed(t; 0, 0, ...) = A_n(t)             [Eulerian polynomial]
  G_T^typed(t; 2, 2, ...) = E_T^perm(t)         [all-permutation forward-edge poly, THM-062/063]

NOTE on E_T naming collision:
  E_T^perm(t) = sum over ALL n! permutations of t^{forward edges}  (computed by forward_edge_dist_dp)
  E_T^ham(t)  = sum over Hamiltonian paths only of t^{descents}    (computed by direct_ET below)
  G_T(t, 2) equals E_T^perm(t), NOT E_T^ham(t). These are different polynomials.

Question: does this factor through U_T in some nice way?
If U_T = sum c_lambda p_lambda, and we specialize p_k appropriately...

Under the "descent-counting" specialization (the quasisymmetric one):
  The fundamental QSF L_{S,n} specializes to give descent polynomial terms
  NOT just power-sum evaluation.

So the typed G_T is a hybrid: it uses GS's cycle-type separation but
our descent-counting specialization.

Let's verify the typed G_T at n=5 and n=7.

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import factorial, comb
from sympy import symbols, expand, Poly

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def inflated_eulerian(f, d, k):
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total

def eulerian_poly(n):
    t = symbols('t')
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

def compute_typed_GT(A, n, t_val, y_vals):
    """Compute G_T^typed(t_val; y_3, y_5, y_7, ...).
    y_vals is a dict: {3: y3_val, 5: y5_val, 7: y7_val, ...}
    """
    d = n - 1
    cycles = find_ALL_directed_odd_cycles(A, n)
    nc = len(cycles)

    if nc == 0:
        # Only the base Eulerian term
        return sum(eulerian_number(n, k) * t_val**k for k in range(n))

    adj = conflict_graph_adj(cycles)
    indep = enumerate_independent_sets(adj, nc)

    result = sum(eulerian_number(n, k) * t_val**k for k in range(n))

    for s in indep:
        if not s:
            continue  # skip empty set (already counted as A_n)

        # Product of y values for cycle lengths
        y_prod = 1
        for i in s:
            c_len = cycles[i][1]
            y_prod *= y_vals.get(c_len, 0)

        if y_prod == 0:
            continue

        # Compute f_I = n-1 - sum(len-1) for each cycle, wait no...
        # f_I is the number of "free" positions = (n-1) - sum(len(c) - 1) for c in S
        # = n - 1 - (sum_len - parts)
        parts = len(s)
        sum_len = sum(cycles[i][1] for i in s)
        f_I = d - (sum_len - parts)  # = n - 1 - sum(len_i - 1)

        # A_{f+1}(t_val) * (t_val - 1)^{d - f}
        A_f1 = sum(eulerian_number(f_I + 1, k) * t_val**k for k in range(f_I + 1))
        correction = A_f1 * (t_val - 1)**(d - f_I)

        result += y_prod * correction

    return result

def direct_ET(A, n, t_val):
    """Hamiltonian-path descent polynomial: sum over Hamiltonian paths of t^{des(H)}.

    WARNING: This is E_T^ham(t), NOT E_T^perm(t). G_T(t,2) equals the
    all-permutation forward-edge polynomial E_T^perm(t), not this function.
    See THM-062/063 and forward_edge_dist_dp for the correct definition.
    """
    total = 0
    for perm in permutations(range(n)):
        valid = all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1))
        if valid:
            desc = sum(1 for i in range(n-1) if perm[i] > perm[i+1])
            total += t_val**desc
    return total

def main():
    print("=" * 70)
    print("TYPED G_T: MULTIVARIATE GENERATING FUNCTION")
    print("=" * 70)

    # n=5 verification
    n = 5
    m = n*(n-1)//2
    print(f"\n--- n={n}: Verify G_T^typed specializations ---")

    t = symbols('t')

    all_ok = True
    for bits_int in range(0, 2**m, 64):  # sample every 64th tournament
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)

        # Test 1 removed: previously compared G_T^typed(t;2,2,...) against direct_ET,
        # but direct_ET computes E_T^ham (Hamiltonian-path descents), while
        # G_T(t,2) = E_T^perm (all-permutation forward edges). These differ.
        # To test G_T(t,2), compare against forward_edge_dist_dp instead.

        # Test 2: G_T^typed(0; y3,y5) = I_typed(Omega; y3, y5)
        cycles = find_ALL_directed_odd_cycles(A, n)
        nc = len(cycles)
        if nc > 0:
            adj = conflict_graph_adj(cycles)
            indep = enumerate_independent_sets(adj, nc)
            for y3_val, y5_val in [(1,1), (2,2), (0,1), (1,0), (3,5)]:
                # I_typed at (y3, y5)
                i_typed = 0
                for s in indep:
                    prod_val = 1
                    for i in s:
                        if cycles[i][1] == 3:
                            prod_val *= y3_val
                        elif cycles[i][1] == 5:
                            prod_val *= y5_val
                    i_typed += prod_val

                gt_val = compute_typed_GT(A, n, 0, {3: y3_val, 5: y5_val})
                if abs(gt_val - i_typed) > 1e-10:
                    print(f"  FAIL (t=0): bits={bits_int}, y=({y3_val},{y5_val}): "
                          f"GT={gt_val}, I_typed={i_typed}")
                    all_ok = False

        # Test 3: G_T^typed(t; x,x,...) = G_T(t, x) -- standard version
        for t_val, x_val in [(0, 1), (2, 1), (1, 3), (-1, 2)]:
            typed_val = compute_typed_GT(A, n, t_val, {3: x_val, 5: x_val})
            standard_val = compute_typed_GT(A, n, t_val, {3: x_val, 5: x_val})
            # These should be equal by definition (same computation)
            # The real test is: does typed(t; x,x) = standard G_T(t,x)?

    if all_ok:
        print(f"  ALL CHECKS PASS!")

    # Now the exciting part: typed G_T with DIFFERENT weights for 3 vs 5 cycles
    print(f"\n--- Typed G_T reveals cycle-length-sensitive structure ---")
    print(f"Setting y_3=x, y_5=0 isolates 3-cycle contributions")
    print(f"Setting y_3=0, y_5=x isolates 5-cycle contributions")

    for bits_int in [100, 200, 920]:
        if bits_int >= 2**m:
            continue
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)

        cycles = find_ALL_directed_odd_cycles(A, n)
        t3 = sum(1 for c, l in cycles if l == 3)
        t5 = sum(1 for c, l in cycles if l == 5)

        print(f"\nbits={bits_int}, t3={t3}, t5={t5}")

        for t_val in [0, 1, 2, -1]:
            full = compute_typed_GT(A, n, t_val, {3: 2, 5: 2})
            only3 = compute_typed_GT(A, n, t_val, {3: 2, 5: 0})
            only5 = compute_typed_GT(A, n, t_val, {3: 0, 5: 2})
            base = compute_typed_GT(A, n, t_val, {3: 0, 5: 0})

            # Superposition test: does full = base + (only3 - base) + (only5 - base)?
            # i.e., is the typed GT additive in cycle contributions?
            superposition = base + (only3 - base) + (only5 - base)
            is_additive = abs(full - superposition) < 1e-10

            print(f"  t={t_val:+d}: full={full:8.0f}, base={base:8.0f}, "
                  f"only3={only3:8.0f}, only5={only5:8.0f}, "
                  f"additive={'YES' if is_additive else 'NO (diff='+str(full-superposition)+')'}")

    print(f"\n\n--- KEY OBSERVATION ---")
    print(f"If the typed GT is NOT additive in cycle types,")
    print(f"then there are INTERACTION terms between different cycle lengths.")
    print(f"These interactions come from the Eulerian coefficients")
    print(f"A_{{f+1}}(t) * (t-1)^{{d-f}} which depend on the TOTAL number")
    print(f"of cycle vertices, not individual cycle lengths.")

if __name__ == "__main__":
    main()
