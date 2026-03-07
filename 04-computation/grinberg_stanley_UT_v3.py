#!/usr/bin/env python3
"""
Grinberg-Stanley U_T computation v3 -- exploring normalizations.

Key observations from v2:
1. All phi conventions give the same U_T for tournaments
2. The q-coefficients are NOT always nonneg (T_3 gives q_3 = -1)
3. But C_3 gives q_3 = +1

Hypotheses to test:
A) Maybe Grinberg-Stanley normalize by the Eulerian polynomial A_n(t)?
   Or subtract off the "transitive tournament" contribution?
B) Maybe the right object is U_T - U_{T_n} (difference from transitive)?
C) Maybe we need to look at the ABSOLUTE VALUE and signs encode something else?
D) Maybe the coefficient of q_{(k)} (single part k) counts oriented k-cycles,
   and the sign indicates direction?

Actually, re-reading the task prompt more carefully:
"Their U_T is a symmetric function. When you evaluate U_T at the specialization
 x_i = t^{i-1}, you get the generating function for Hamiltonian paths by descents,
 which should be our E_T(t)."

So let me VERIFY this claim first. Does ps(U_T) = E_T(t)?
From v1, the principal specialization didn't match. But maybe we need
fewer variables, or a different specialization.

Actually, the principal specialization of a symmetric function of degree n
with n variables gives a polynomial in t. For the power-sum symmetric function:
  ps_n(p_lambda) = prod_i (1 + t^{lambda_i} + t^{2*lambda_i} + ... + t^{(n-1)*lambda_i})

The degree-n part of a symmetric function should capture all the information.

Wait -- I think the issue is simpler. The symmetric function U_T has degree n
(it's a sum of p_{type(sigma)} where type(sigma) is a partition of n).
The principal specialization with n variables gives a polynomial of degree n(n-1).
But E_T(t) has degree n-1.

The connection might be through the STABLE principal specialization (infinitely many
variables), which gives p_k -> 1/(1-t^k). Then the degree-n part would be...
no, that diverges.

Let me try a different approach: compute the coefficient of the AUGMENTED monomial
symmetric function. In Grinberg-Stanley, they show that [x_1...x_n] U_T gives
the number of acyclic orientations... no, for tournaments it should be different.

FRESH APPROACH: Let me compute E_T(t) directly from the definition as a
generating function for Hamiltonian paths by descents, and separately compute
what we get from U_T by extracting the coefficient of x_1*x_2*...*x_n
(the "Frobenius character" extraction).

The key identity should be:
  (1/n!) * [x_1 x_2 ... x_n](prod_{i=1}^{n} x_{sigma(i)}) evaluated at x_i -> some specialization

Actually, let me just try the omega involution. In symmetric function theory,
omega sends p_k -> (-1)^{k-1} p_k. Maybe omega(U_T) is the one with nonneg q-coefficients?
"""

from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction


def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def cyclic_tournament(n):
    assert n % 2 == 1
    A = [[0]*n for _ in range(n)]
    half = (n-1) // 2
    for i in range(n):
        for d in range(1, half+1):
            j = (i + d) % n
            A[i][j] = 1
    return A

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def cycle_type(sigma):
    n = len(sigma)
    visited = [False]*n
    lengths = []
    for i in range(n):
        if not visited[i]:
            length = 0; j = i
            while not visited[j]:
                visited[j] = True; j = sigma[j]; length += 1
            lengths.append(length)
    return tuple(sorted(lengths, reverse=True))

def hamiltonian_paths(A, n):
    return [p for p in permutations(range(n))
            if all(A[p[i]][p[i+1]] == 1 for i in range(n-1))]

def descent_count(perm):
    return sum(1 for i in range(len(perm)-1) if perm[i] > perm[i+1])

def find_all_odd_cycles(A, n):
    cycles = []; seen = set()
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(A[perm[i]][perm[(i+1)%length]] == 1 for i in range(length)):
                    min_idx = list(perm).index(min(perm))
                    canonical = perm[min_idx:] + perm[:min_idx]
                    if canonical not in seen:
                        seen.add(canonical); cycles.append(canonical)
    return cycles

def independence_polynomial_coeffs(cycles):
    m = len(cycles)
    cycle_sets = [set(c) for c in cycles]
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = True
    alpha = defaultdict(int)
    for mask in range(2**m):
        verts = [i for i in range(m) if (mask >> i) & 1]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    indep = False; break
            if not indep: break
        if indep:
            alpha[len(verts)] += 1
    max_k = max(alpha.keys()) if alpha else 0
    return [alpha.get(k, 0) for k in range(max_k + 1)]

def compute_UT(A, n):
    UT = defaultdict(int)
    for sigma in permutations(range(n)):
        phi = sum(1 for i in range(n-1) if A[sigma[i]][sigma[i+1]] == 0)
        ct = cycle_type(sigma)
        UT[ct] += (-1)**phi
    return dict(UT)

def format_partition(partition):
    c = Counter(partition)
    parts = []
    for k in sorted(c.keys(), reverse=True):
        parts.append(f"p_{k}" + (f"^{c[k]}" if c[k] > 1 else ""))
    return "*".join(parts)


def apply_omega(UT):
    """Apply the omega involution: p_k -> (-1)^{k-1} p_k."""
    result = {}
    for partition, coeff in UT.items():
        sign = 1
        for k in partition:
            sign *= (-1)**(k-1)
        result[partition] = coeff * sign
    return result


def extract_e_expansion(UT, n):
    """
    Try to extract the e_1^n coefficient of U_T.

    In the power-sum basis, e_1^n = (1/n!) * sum_{sigma in S_n} chi^{(1^n)}(sigma) p_{type(sigma)}
    where chi^{(1^n)} is the sign character: chi^{(1^n)}(sigma) = (-1)^{n - #cycles(sigma)}

    So [e_1^n] U_T = sum_sigma (-1)^{phi(sigma)} * (1/z_lambda) where lambda = type(sigma)
    ... this is getting complicated.

    Actually, the simpler approach: for the monomial symmetric function m_{(1^n)} = e_n,
    we need the coefficient of x_1*x_2*...*x_n in U_T, which picks out the
    contribution from each permutation sigma weighted by (-1)^phi.

    Since p_{type(sigma)} = p_{lambda_1} * ... * p_{lambda_r}, and
    [x_1...x_n] p_{lambda} = n!/z_lambda (where z_lambda is the standard factor),
    we get:
    [x_1...x_n] U_T = sum_sigma (-1)^phi(sigma) * [x_1...x_n] p_{type(sigma)}
    But [x_1...x_n] p_{lambda} is the same for all sigma with the same type lambda,
    so it's:
    = sum_lambda (sum_{sigma: type=lambda} (-1)^phi(sigma)) * [x_1...x_n] p_lambda
    = sum_lambda U_T[lambda] * [x_1...x_n] p_lambda

    [x_1...x_n] p_lambda = n!/z_lambda where z_lambda = prod_i (i^{m_i} * m_i!)
    for a partition with m_i parts of size i.
    """
    total = Fraction(0)
    for partition, coeff in UT.items():
        if coeff == 0:
            continue
        # Compute z_lambda
        c = Counter(partition)
        z = 1
        for k, mult in c.items():
            z *= k**mult
            for j in range(1, mult+1):
                z *= j
        # [x_1...x_n] p_lambda = n!/z_lambda
        import math
        contrib = Fraction(math.factorial(n), z)
        total += Fraction(coeff) * contrib
    return total


def main():
    print("="*70)
    print("GRINBERG-STANLEY U_T: TESTING OMEGA INVOLUTION AND NORMALIZATIONS")
    print("="*70)

    for n, tournaments in [(3, None), (5, None)]:
        print(f"\n{'#'*70}")
        print(f"# n = {n}")
        print(f"{'#'*70}")

        if tournaments is None:
            if n <= 5:
                tournaments = [
                    (transitive_tournament(n), f"T_{n} (transitive)"),
                    (cyclic_tournament(n), f"C_{n} (cyclic)"),
                ]

        for A, name in tournaments:
            print(f"\n{'='*60}")
            print(f"{name}")
            print(f"{'='*60}")

            paths = hamiltonian_paths(A, n)
            ET = defaultdict(int)
            for p in paths:
                ET[descent_count(p)] += 1

            cycles = find_all_odd_cycles(A, n)
            ip = independence_polynomial_coeffs(cycles)

            UT = compute_UT(A, n)
            omega_UT = apply_omega(UT)

            print(f"H(T) = {len(paths)}")
            print(f"E_T(t) = {dict(ET)}")
            print(f"I(Omega,x) = {ip}")

            print(f"\nU_T power-sum expansion:")
            for p in sorted(UT.keys()):
                if UT[p] != 0:
                    print(f"  {UT[p]:+d} * {format_partition(p)}")

            print(f"\nomega(U_T) power-sum expansion:")
            for p in sorted(omega_UT.keys()):
                if omega_UT[p] != 0:
                    print(f"  {omega_UT[p]:+d} * {format_partition(p)}")

            # q-expansion of omega(U_T)
            print(f"\nomega(U_T) q-expansion (odd parts only):")
            for p in sorted(omega_UT.keys()):
                c = omega_UT[p]
                if c == 0 or any(k % 2 == 0 for k in p):
                    continue
                num_big = sum(1 for k in p if k >= 3)
                q = Fraction(c, 2**num_big)
                parts = "*".join(f"q_{k}" for k in p)
                print(f"  {q} * {parts}")

            # [x_1...x_n] extraction
            extract = extract_e_expansion(UT, n)
            extract_omega = extract_e_expansion(omega_UT, n)
            print(f"\n[x_1...x_n] U_T = {extract}")
            print(f"[x_1...x_n] omega(U_T) = {extract_omega}")
            print(f"H(T) = {len(paths)}")

    # Exhaustive n=3: check omega(U_T) q-coefficients
    print(f"\n\n{'#'*70}")
    print(f"# EXHAUSTIVE n=3: omega(U_T) q-coefficients")
    print(f"{'#'*70}")

    for A in all_tournaments(3):
        n = 3
        UT = compute_UT(A, n)
        omega_UT = apply_omega(UT)
        cycles = find_all_odd_cycles(A, n)
        ip = independence_polynomial_coeffs(cycles)
        paths = hamiltonian_paths(A, n)

        q_coeffs = {}
        for p in sorted(omega_UT.keys()):
            c = omega_UT[p]
            if c == 0 or any(k % 2 == 0 for k in p):
                continue
            num_big = sum(1 for k in p if k >= 3)
            q_coeffs[p] = Fraction(c, 2**num_big)

        print(f"\n  H(T)={len(paths)}, I(Omega,x)={ip}")
        print(f"  omega(U_T) q-coeffs: {q_coeffs}")

    # Exhaustive n=3: check if omega makes all q-coefficients nonneg
    print(f"\n\n{'='*60}")
    print(f"Are ALL omega(U_T) q-coefficients nonneg for n=3?")
    print(f"{'='*60}")

    all_nonneg = True
    for A in all_tournaments(3):
        UT = compute_UT(A, 3)
        omega_UT = apply_omega(UT)
        for p, c in omega_UT.items():
            if all(k % 2 == 1 for k in p):
                num_big = sum(1 for k in p if k >= 3)
                q = Fraction(c, 2**num_big)
                if q < 0:
                    all_nonneg = False
                    break
    print(f"  Result: {all_nonneg}")

    # Now test: maybe the right relationship is between U_T and the
    # independence polynomial evaluated as a GENERATING FUNCTION
    #
    # The THM-063 says G_T(0, x) = I(Omega(T), x) and G_T(t, 2) = E_T(t).
    # So maybe the symmetric function analogue is:
    #   U_T evaluated at a specific specialization gives I(Omega(T), x)
    #   or E_T(t) or both.

    # Test: [x_1 x_2 ... x_n] U_T as a function of something
    print(f"\n\n{'='*60}")
    print(f"Testing [x_1...x_n] U_T for all n=3,5 tournaments")
    print(f"{'='*60}")

    for n in [3, 5]:
        print(f"\n  n={n}:")
        for A in all_tournaments(n):
            UT = compute_UT(A, n)
            extract = extract_e_expansion(UT, n)
            paths = hamiltonian_paths(A, n)
            H = len(paths)
            cycles = find_all_odd_cycles(A, n)
            ip = independence_polynomial_coeffs(cycles)
            ip_at_2 = sum(c * 2**k for k, c in enumerate(ip))

            # Check various relationships
            if extract != 0:
                ratio = Fraction(H) / extract
            else:
                ratio = None
            if n == 3:  # Only print for small n
                print(f"    [x1...xn]U_T = {extract}, H(T) = {H}, ratio = {ratio}")

        if n == 5:
            # Just check if there's a universal ratio
            ratios = set()
            for A in all_tournaments(n):
                UT = compute_UT(A, n)
                extract = extract_e_expansion(UT, n)
                paths = hamiltonian_paths(A, n)
                H = len(paths)
                if extract != 0:
                    ratios.add(Fraction(H) / extract)
            print(f"    Distinct ratios H(T)/[x1...xn]U_T: {ratios}")

    # CRUCIAL TEST: Does the sum of coefficients of odd-part partitions in U_T,
    # weighted by z_lambda^{-1} * n!, give H(T)?
    # Or rather: does p_k at x=2 give 2?
    #   p_k(2) for k=1 gives 2, for k>1 gives 2 (if one variable x=2) or 2^k (if x_i=2 for all i)
    #
    # Wait -- I(Omega(T), 2) = H(T). And the OCF says:
    # H(T) = sum_{independent sets S in Omega(T)} 2^{|S|}
    # So at x=2: alpha_0 + 2*alpha_1 + 4*alpha_2 + ...
    #
    # And U_T in the q-basis: c_0 * q_1^n + c_1 * q_3 * q_1^{n-3} + ...
    # Evaluating at q_k = 2 for all k:
    # This gives sum of c_lambda * 2^{#parts(lambda)} = sum c_lambda * 2^r
    # where r is the number of parts.
    #
    # Hmm, that would match I(Omega, 2) if c_lambda is grouped by #parts.

    print(f"\n\n{'='*60}")
    print(f"TEST: sum of q-coefficients * 2^(#parts) vs I(Omega, 2)")
    print(f"{'='*60}")

    for n in [3, 5]:
        print(f"\n  n={n}:")
        all_match = True
        for A in all_tournaments(n):
            UT = compute_UT(A, n)
            paths = hamiltonian_paths(A, n)
            H = len(paths)
            cycles = find_all_odd_cycles(A, n)
            ip = independence_polynomial_coeffs(cycles)

            # Evaluate U_T at p_k = 2 for odd k, p_k = 0 for even k
            # (This is the specialization that selects odd-part partitions)
            val_odd = 0
            for partition, coeff in UT.items():
                if all(k % 2 == 1 for k in partition):
                    val_odd += coeff * (2 ** len(partition))

            # Also try: p_k = 1 for all k
            val_all_1 = 0
            for partition, coeff in UT.items():
                val_all_1 += coeff  # p_k(1) = 1

            # Also try: p_1 = n (principal spec at t=1), p_k = n for all k
            # At t=1, p_k(1,1,...,1) = n for all k
            import math
            val_at_t1 = 0
            for partition, coeff in UT.items():
                val_at_t1 += coeff * (n ** len(partition))

            if n == 3:
                print(f"    H={H}, U_T(p_odd=2,p_even=0)={val_odd}, U_T(p=1)={val_all_1}, U_T(p=n)={val_at_t1}, n!={math.factorial(n)}")
            else:
                pass  # too many

        # Summary for n=5
        if n == 5:
            vals = set()
            for A in all_tournaments(n):
                UT = compute_UT(A, n)
                H = len(hamiltonian_paths(A, n))
                val_odd = sum(c * 2**len(p) for p, c in UT.items() if all(k%2==1 for k in p))
                vals.add((H, val_odd))
            print(f"    Sample (H, U_T|odd=2): {list(vals)[:10]}")

    # Final key test: evaluate U_T at p_k = 2 for all k (not just odd)
    print(f"\n\n{'='*60}")
    print(f"TEST: U_T evaluated at p_k = 2 for ALL k")
    print(f"{'='*60}")

    for n in [3, 5]:
        print(f"\n  n={n}:")
        for A in all_tournaments(n):
            UT = compute_UT(A, n)
            H = len(hamiltonian_paths(A, n))
            val = sum(c * 2**len(p) for p, c in UT.items())
            if n == 3:
                print(f"    H={H}, U_T(p_k=2)={val}")
        if n == 5:
            results = set()
            for A in all_tournaments(n):
                UT = compute_UT(A, n)
                H = len(hamiltonian_paths(A, n))
                val = sum(c * 2**len(p) for p, c in UT.items())
                results.add((H, val))
            print(f"    All (H, U_T|p=2) pairs: {sorted(results)}")


if __name__ == "__main__":
    main()
