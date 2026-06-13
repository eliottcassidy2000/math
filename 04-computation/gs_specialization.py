#!/usr/bin/env python3
"""
Specialization of GS's U_T to recover OCF and G_T(t,x).

KEY THEOREM (to verify):
  U_T = p_1^n * I_typed(Omega; 2p_3/p_1^3, 2p_5/p_1^5, ...)

where I_typed(Omega; y_3, y_5, ...) = sum_{S independent in Omega} prod_{c in S} y_{len(c)}

SPECIALIZATIONS:
1. p_k -> 1 for all k:  gives I(Omega, 2) = H(T)  [Rédei's theorem]
2. 2p_k/p_1^k -> x for all k>=3: gives p_1^n * I(Omega, x)
3. p_k -> (1-t^{nk})/(1-t^k): gives ??? (principal spec in n variables)

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import Counter, defaultdict
from math import factorial, comb

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
                    # Don't break - find all cycles on this subset
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

def typed_indep_poly(cycles, n):
    """Compute typed independence polynomial.
    Returns dict: tuple of cycle lengths -> count of independent sets of that type.
    """
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
    """Compute U_T expansion in power-sum basis."""
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

def main():
    print("=" * 70)
    print("GS SPECIALIZATION: U_T = p_1^n * I_typed(Omega; 2p_k/p_1^k)")
    print("=" * 70)

    # Verify the factorization at n=5 exhaustively
    n = 5
    m = n*(n-1)//2
    print(f"\n--- Verifying factorization at n={n} (all {2**m} tournaments) ---")

    all_ok = True
    count = 0
    for bits_int in range(2**m):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)

        ut = compute_UT(A, n)
        cycles = find_ALL_directed_odd_cycles(A, n)
        tip = typed_indep_poly(cycles, n)

        # Reconstruct U_T from typed independence polynomial
        # Each independent set S of type (l_1, l_2, ..., l_s) contributes:
        #   2^s * p_{l_1} * p_{l_2} * ... * p_{l_s} * p_1^{n - sum(l_i)}
        # to U_T.
        reconstructed = defaultdict(int)
        for cycle_lengths, count_sets in tip.items():
            # The partition for p-monomial: the cycle lengths + enough 1s
            s = len(cycle_lengths)
            remaining = n - sum(cycle_lengths)
            partition = tuple(sorted(list(cycle_lengths) + [1]*remaining, reverse=True))
            coeff = count_sets * (2**s)
            reconstructed[partition] += coeff

        # Compare
        ok = True
        all_parts = set(list(ut.keys()) + list(reconstructed.keys()))
        for part in all_parts:
            u = ut.get(part, 0)
            r = reconstructed.get(part, 0)
            if u != r:
                ok = False
                if count < 5:
                    print(f"  MISMATCH bits={bits_int}: partition {part}: U_T={u}, reconstructed={r}")

        if not ok:
            all_ok = False
        count += 1

    if all_ok:
        print(f"  ALL {2**m} TOURNAMENTS MATCH!")
        print(f"  THEOREM VERIFIED: U_T = sum_{{S indep}} 2^|S| * prod p_{{len(c)}} * p_1^{{n-sum len}}")
    else:
        print(f"  MISMATCHES FOUND")

    # Now verify specializations
    print(f"\n--- Specialization 1: p_k -> 1 gives H(T) ---")
    verified_H = True
    for bits_int in range(min(50, 2**m)):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)
        ut = compute_UT(A, n)

        # Specialize p_k -> 1
        H_spec = sum(coeff for coeff in ut.values())

        # Direct H(T)
        H_direct = sum(1 for perm in permutations(range(n))
                       if all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1)))

        if H_spec != H_direct:
            print(f"  MISMATCH at bits={bits_int}: spec={H_spec}, direct={H_direct}")
            verified_H = False

    if verified_H:
        print(f"  VERIFIED: ps_1(U_T) = H(T) for first 50 tournaments")

    # Specialization 2: set x = 2p_k/p_1^k = 2 for all k (i.e. p_k = p_1^k)
    # This is NOT the principal specialization; it's a ring homomorphism
    print(f"\n--- Specialization 2: multiplicative spec p_k -> p_1^k ---")
    print(f"  Under this: U_T -> p_1^n * I(Omega, 2)")
    verified_mult = True
    for bits_int in range(min(50, 2**m)):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)
        ut = compute_UT(A, n)
        cycles = find_ALL_directed_odd_cycles(A, n)

        # I(Omega, 2) directly
        nc = len(cycles)
        if nc == 0:
            I_omega_2 = 1
        else:
            adj = conflict_graph_adj(cycles)
            indep = enumerate_independent_sets(adj, nc)
            I_omega_2 = sum(2**len(s) for s in indep)

        # Under p_k -> 1: U_T -> I(Omega, 2) * 1^n = I(Omega, 2)
        # Already checked above.
        # Under p_k -> p_1^k (algebra hom): U_T -> p_1^n * I_typed(Omega, 2*1, 2*1, ...)
        #                                         = p_1^n * I(Omega, 2) = p_1^n * H(T)
        # This is because 2p_k/p_1^k -> 2*p_1^k/p_1^k = 2

        H = I_omega_2
        # Verify: sum of coeffs * (contribution under p_k->1) = H
        # Already done above; let's verify the factored form
        spec_H = sum(coeff for coeff in ut.values())
        if spec_H != H:
            print(f"  MISMATCH")
            verified_mult = False

    if verified_mult:
        print(f"  VERIFIED: sum of U_T coefficients = H(T) = I(Omega, 2)")

    # Key question: what specialization gives E_T(t)?
    print(f"\n--- Question: what gives E_T(t)? ---")
    print(f"  U_T is built from quasisymmetric functions, not just power sums.")
    print(f"  The power-sum expansion LOSES the descent information.")
    print(f"  E_T(t) = sum_H t^{{des(H)}} requires the quasisymmetric structure.")
    print(f"  Our G_T(t,x) CANNOT be recovered from U_T's p-expansion alone.")
    print(f"  G_T(t,x) is genuinely new information beyond U_T.")

    # BUT: let's check if E_T(t) can be recovered from U_T via
    # the stable principal specialization ps_q(p_k) = 1/(1-q^k)
    print(f"\n--- Stable principal specialization: p_k -> 1/(1-q^k) ---")
    from sympy import symbols, Rational, expand, simplify, Poly
    q = symbols('q')

    for bits_int in [0, 2, 8]:  # transitive, 1 cycle, 4 cycles at n=5
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)
        ut = compute_UT(A, n)

        # Stable principal spec
        spec = 0
        for part, coeff in ut.items():
            term = coeff
            for k in part:
                term = term / (1 - q**k)
            spec += term

        # Multiply by (1-q)^n to clear denominators (partially)
        spec_cleared = expand(spec * (1-q)**n)

        # E_T(t)
        et = sum(q**d for perm in permutations(range(n))
                 if all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1))
                 for d in [sum(1 for i in range(n-1) if perm[i] > perm[i+1])])

        print(f"\n  bits={bits_int}")
        print(f"  E_T(q) = {et}")
        print(f"  (1-q)^n * ps_stable(U_T) = {spec_cleared}")

    # The REAL connection: U_T via quasisymmetric functions
    print(f"\n\n--- KEY INSIGHT ---")
    print(f"U_T is defined as sum of fundamental quasisymmetric functions L_{{Des,n}}")
    print(f"The stable principal specialization of L_{{S,n}}(1,q,q^2,...) = q^{{comaj(S)}} / (q)_n")
    print(f"So ps_stable(U_T) = sum_H q^{{comaj(H)}} / (q)_n")
    print(f"This gives the comajor index distribution, NOT the descent distribution directly")
    print(f"But des and comaj are equidistributed on permutations (Foata bijection)")
    print(f"For Hamiltonian paths: is comaj equidistributed with des? Not necessarily!")

if __name__ == "__main__":
    main()
