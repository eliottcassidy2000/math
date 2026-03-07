#!/usr/bin/env python3
"""
Grinberg-Stanley U_T computation v4 -- CORRECT approach.

After reviewing the paper more carefully:

Grinberg-Stanley Theorem 1.30: For a digraph D on V={1,...,n},
  U_D = sum_{sigma in S_n} (-1)^{asc_D(sigma)} p_{type(sigma)}

where asc_D(sigma) = #{i in {1,...,n-1} : (sigma_i, sigma_{i+1}) is an arc of D}

Key: they use ASCENTS (arcs present in D), not descents/non-arcs!

For a tournament T:
  asc_T(sigma) = #{i : T[sigma_i][sigma_{i+1}} = 1}
               = #{forward arcs in the sequence sigma}

A Hamiltonian path has asc_T = n-1 (all arcs forward), so contributes
(-1)^{n-1} * p_{type(sigma)}.

The permutation sigma = identity for a Hamiltonian path has type (1^n).
BUT sigma in their notation is a permutation, and type = cycle type.

Wait, I need to be more careful about what the "type" means here.
The TYPE of sigma is its CYCLE TYPE as a permutation.
But the Hamiltonian path condition depends on sigma as a SEQUENCE, not as a permutation.

Actually, re-reading: they sum over ALL permutations sigma of V.
sigma = (sigma_1, ..., sigma_n) is a sequence.
type(sigma) = cycle type of sigma VIEWED AS A PERMUTATION (i.e., sigma(i) = sigma_i in one-line notation).

So for the identity permutation sigma = (1,2,...,n), type = (1^n),
and asc_T(id) = #{i : T[i][i+1] = 1}.

OK so my original computation was RIGHT, just with the sign convention flipped.
Let me check: "forward" convention had phi = #{non-arcs of D} = n-1 - asc_D.
So (-1)^phi = (-1)^{n-1-asc_D} = (-1)^{n-1} * (-1)^{-asc_D} = (-1)^{n-1} * (-1)^{asc_D}.

So U_D(forward) = (-1)^{n-1} * U_D(Grinberg-Stanley).

For n=3: (-1)^2 = 1, so they agree.
For n=5: (-1)^4 = 1, so they agree too.

Hmm, so the conventions all agree for these cases. The issue must be elsewhere.

Let me INSTEAD focus on the ACTUAL relationship between U_T and E_T(t).

The Grinberg-Stanley paper, Corollary 1.36:
  When D is acyclic (e.g., transitive tournament), U_D = n! * e_n
  where e_n is the nth elementary symmetric function.

  More generally, for D with a unique acyclic orientation (like a tournament?),
  U_D involves the quasisymmetric function...

Actually, let me just try the CORRECT specialization. The paper says:
for the STABLE principal specialization (infinitely many variables with x_i = t^{i-1}):

  ps(p_k) = 1/(1-t^k)

And for the POLYNOMIAL principal specialization with m variables:
  ps_m(p_k) = (1-t^{mk})/(1-t^k)

The key identity should be something like:
  ps(U_T) / ps(U_{T_n}) = E_T(t) / E_{T_n}(t) = E_T(t) / A_n(t)?

Or maybe: There's a "balanced" version where you divide by (1-t)^n or similar.

Let me just try all reasonable normalizations computationally.
"""

from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction
import math


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
    """U_T = sum_sigma (-1)^{asc_T(sigma)} p_{type(sigma)}"""
    UT = defaultdict(int)
    for sigma in permutations(range(n)):
        asc = sum(1 for i in range(n-1) if A[sigma[i]][sigma[i+1]] == 1)
        ct = cycle_type(sigma)
        UT[ct] += (-1)**asc
    return dict(UT)


def poly_mul(p1, p2):
    """Multiply two polynomials (dicts {exponent: coeff})."""
    result = defaultdict(Fraction)
    for e1, c1 in p1.items():
        for e2, c2 in p2.items():
            result[e1+e2] += c1 * c2
    return {e: c for e, c in result.items() if c != 0}

def poly_add(p1, p2):
    result = defaultdict(Fraction)
    for e, c in p1.items(): result[e] += c
    for e, c in p2.items(): result[e] += c
    return {e: c for e, c in result.items() if c != 0}

def poly_scale(p, s):
    return {e: c*s for e, c in p.items() if c*s != 0}

def pk_poly_stable(k):
    """p_k under stable principal spec: 1/(1-t^k) truncated to degree N."""
    # We'll represent as exact polynomial up to some max degree
    N = 50
    return {i*k: Fraction(1) for i in range(N//k + 1)}

def pk_poly_finite(k, m):
    """p_k(1, t, ..., t^{m-1}) = sum_{i=0}^{m-1} t^{ik}."""
    return {i*k: Fraction(1) for i in range(m)}

def evaluate_UT_poly(UT, pk_func):
    """Evaluate U_T as polynomial in t using given p_k function."""
    total = defaultdict(Fraction)
    for partition, coeff in UT.items():
        if coeff == 0:
            continue
        term = {0: Fraction(1)}
        for k in partition:
            term = poly_mul(term, pk_func(k))
        for e, c in term.items():
            total[e] += Fraction(coeff) * c
    return {e: c for e, c in total.items() if c != 0}

def poly_to_str(p, max_terms=20):
    """Format polynomial."""
    if not p:
        return "0"
    items = sorted(p.items())[:max_terms]
    return " + ".join(f"{c}*t^{e}" for e, c in items if c != 0)

def eulerian_poly(n):
    """Compute A_n(t) = sum_{sigma in S_n} t^{des(sigma)}."""
    result = defaultdict(int)
    for sigma in permutations(range(n)):
        d = descent_count(sigma)
        result[d] += 1
    return dict(result)


def main():
    print("="*70)
    print("GRINBERG-STANLEY U_T: FINDING THE RIGHT SPECIALIZATION")
    print("="*70)

    for n in [3, 5]:
        print(f"\n{'#'*70}")
        print(f"# n = {n}")
        print(f"{'#'*70}")

        An = eulerian_poly(n)
        print(f"\nA_{n}(t) = {An}")

        tournaments = [
            (transitive_tournament(n), f"T_{n}"),
            (cyclic_tournament(n), f"C_{n}"),
        ]

        for A, name in tournaments:
            print(f"\n{'='*50}")
            print(f"{name}")
            print(f"{'='*50}")

            paths = hamiltonian_paths(A, n)
            ET = defaultdict(int)
            for p in paths:
                ET[descent_count(p)] += 1
            ET = dict(ET)

            cycles = find_all_odd_cycles(A, n)
            ip = independence_polynomial_coeffs(cycles)

            UT = compute_UT(A, n)

            print(f"H(T) = {len(paths)}")
            print(f"E_T(t) = {ET}")
            print(f"I(Omega,x) = {ip}")

            # Try stable principal specialization
            ps_stable = evaluate_UT_poly(UT, pk_poly_stable)

            # Try finite principal specialization with n variables
            ps_n = evaluate_UT_poly(UT, lambda k: pk_poly_finite(k, n))

            # Truncate to reasonable degree for display
            max_deg_display = 2*n
            ps_stable_trunc = {e: c for e, c in ps_stable.items() if e <= max_deg_display}
            ps_n_trunc = {e: c for e, c in ps_n.items() if e <= max_deg_display}

            print(f"\nps_stable(U_T) (low degrees): {poly_to_str(ps_stable_trunc)}")
            print(f"ps_{n}(U_T) (low degrees): {poly_to_str(ps_n_trunc)}")

            # Try: ps_stable(U_T) * (1-t)^n
            # Since p_k -> 1/(1-t^k) under stable spec,
            # p_1^n -> 1/(1-t)^n, which diverges.
            # Maybe multiplying by (1-t)^n gives a polynomial?

            one_minus_t = {0: Fraction(1), 1: Fraction(-1)}
            product = {0: Fraction(1)}
            for _ in range(n):
                product = poly_mul(product, one_minus_t)

            ps_times_factor = poly_mul(ps_stable, product)
            ps_times_factor_trunc = {e: c for e, c in ps_times_factor.items()
                                     if e <= max_deg_display and c != 0}
            print(f"\nps_stable(U_T) * (1-t)^{n}: {poly_to_str(ps_times_factor_trunc)}")

            # Check if this is E_T(t) * something
            # E_T(t) for the transitive tournament is just t^0 = 1

            # Try with n-1 powers of (1-t)
            product2 = {0: Fraction(1)}
            for _ in range(n-1):
                product2 = poly_mul(product2, one_minus_t)

            ps_times_factor2 = poly_mul(ps_stable, product2)
            ps_times_factor2_trunc = {e: c for e, c in ps_times_factor2.items()
                                      if e <= max_deg_display and c != 0}
            print(f"ps_stable(U_T) * (1-t)^{n-1}: {poly_to_str(ps_times_factor2_trunc)}")

    # Actually, let me try a completely different approach.
    # Instead of principal specialization, let me look at the
    # MONOMIAL SYMMETRIC FUNCTION expansion and the Frobenius formula.
    #
    # Key insight: U_T in the e-basis or h-basis might be more revealing.
    # For a tournament, E_T(t) = sum_H t^{des(H)}.
    # This is related to the DESCENT polynomial of acyclic orientations.
    #
    # From Stanley's theory of P-partitions:
    # For a poset P on {1,...,n}, the P-Eulerian polynomial is
    #   W_P(t) = (1-t)^n sum_{k>=0} Omega(P, k) t^k
    # where Omega(P,k) = #{order-preserving maps P -> {1,...,k}}.
    #
    # For a tournament T, a Hamiltonian path corresponds to a LINEAR EXTENSION.
    # If T is viewed as a partial order (transitive closure), then
    # H(T) = #{linear extensions}, but T is only a partial order if it's acyclic!
    #
    # For non-acyclic tournaments, we still have Hamiltonian paths, but they
    # don't correspond to linear extensions of a poset.
    #
    # Let me try yet another approach: direct brute-force comparison.

    print(f"\n\n{'#'*70}")
    print("# DIRECT COMPARISON: U_T odd-part coefficients vs I(Omega, x)")
    print(f"{'#'*70}")

    print(f"\nFor each tournament, comparing:")
    print(f"  - q-coefficient of q_{{(k)}} (single part k) in U_T")
    print(f"  - alpha_k = # independent sets of size k in Omega(T)")
    print(f"  - Various groupings of q-coefficients")

    for n in [3, 5]:
        print(f"\n--- n = {n} ---")
        count = 0
        for A in all_tournaments(n):
            count += 1
            UT = compute_UT(A, n)
            cycles = find_all_odd_cycles(A, n)
            ip = independence_polynomial_coeffs(cycles)
            H = len(hamiltonian_paths(A, n))

            # Group q-coefficients by number of parts
            # q_coeff[r] = sum of q-coefficients for partitions with exactly r parts
            q_by_nparts = defaultdict(Fraction)
            q_by_partition = {}
            for partition, coeff in UT.items():
                if coeff == 0 or any(k % 2 == 0 for k in partition):
                    continue
                num_big = sum(1 for k in partition if k >= 3)
                q_coeff = Fraction(coeff, 2**num_big)
                q_by_partition[partition] = q_coeff

                # Number of "big" parts (>= 3) tells us the "weight" in indep poly
                q_by_nparts[num_big] += q_coeff

            if n == 3 or (n == 5 and count <= 5):
                print(f"\n  T={count}, H={H}, I(Omega,x)={ip}")
                print(f"    q-coefficients by partition: {dict(q_by_partition)}")
                print(f"    Sum by #big-parts: {dict(q_by_nparts)}")
                for k in range(len(ip)):
                    alpha_k = ip[k]
                    qk = q_by_nparts.get(k, Fraction(0))
                    print(f"    alpha_{k}={alpha_k}, sum(q-coeffs with {k} big parts)={qk}")

    # THE BIG QUESTION: for all n=3 tournaments, does:
    # sum of q-coefficients with exactly k "big parts" (parts >= 3)
    # equal alpha_k?
    print(f"\n\n{'#'*70}")
    print("# GLOBAL TEST: q-coefficients grouped by #big-parts vs alpha_k")
    print(f"{'#'*70}")

    for n in [3, 5]:
        all_match = True
        mismatch_count = 0
        total = 0
        for A in all_tournaments(n):
            total += 1
            UT = compute_UT(A, n)
            cycles = find_all_odd_cycles(A, n)
            ip = independence_polynomial_coeffs(cycles)

            q_by_nparts = defaultdict(Fraction)
            for partition, coeff in UT.items():
                if coeff == 0 or any(k % 2 == 0 for k in partition):
                    continue
                num_big = sum(1 for k in partition if k >= 3)
                q_by_nparts[num_big] += Fraction(coeff, 2**num_big)

            max_k = max(len(ip)-1, max(q_by_nparts.keys()) if q_by_nparts else 0)
            for k in range(max_k + 1):
                alpha_k = ip[k] if k < len(ip) else 0
                qk = q_by_nparts.get(k, Fraction(0))
                if qk != alpha_k:
                    all_match = False
                    mismatch_count += 1

        print(f"  n={n}: all_match = {all_match}, mismatches = {mismatch_count}/{total}")

    # Let me try a DIFFERENT grouping: by total weight sum(partition) - n
    # or by sum of (part-1)/2 for parts >= 3
    print(f"\n\n{'#'*70}")
    print("# ALTERNATIVE GROUPINGS")
    print(f"{'#'*70}")

    for n in [3]:
        print(f"\n  n={n}:")
        for A in all_tournaments(n):
            UT = compute_UT(A, n)
            cycles = find_all_odd_cycles(A, n)
            ip = independence_polynomial_coeffs(cycles)
            H = len(hamiltonian_paths(A, n))

            # What if we evaluate U_T at p_1 = 1, p_k = x for odd k>=3, p_even = 0?
            # This gives: sum over odd-part partitions lambda of coeff(lambda) * x^{#big-parts}
            # = the q-generating function grouped by big parts
            # With the 2-division: q_k = 2*p_k so p_k = q_k/2
            # So p_1 = 1, p_k = x/2 for odd k >= 3
            # gives: sum coeff(lambda) * (x/2)^{#big-parts} = sum q_coeff * x^{#big-parts}

            # Just evaluate and compare with I(Omega, x) at various x values
            for x_val in [1, 2, 3]:
                ut_val = 0
                for partition, coeff in UT.items():
                    if any(k % 2 == 0 for k in partition):
                        continue
                    num_big = sum(1 for k in partition if k >= 3)
                    ut_val += coeff * (Fraction(x_val, 2))**num_big

                ip_val = sum(c * x_val**k for k, c in enumerate(ip))
                print(f"    H={H}, x={x_val}: U_T(p1=1,p_odd=x/2)={ut_val}, I(Omega,x)={ip_val}, "
                      f"match={ut_val==ip_val}")

    # Let me try: p_1 = 1, p_k = x for odd k >= 3, p_even = 0, no 2-division
    print(f"\n  Without 2-division:")
    for A in all_tournaments(3):
        UT = compute_UT(A, 3)
        cycles = find_all_odd_cycles(A, 3)
        ip = independence_polynomial_coeffs(cycles)
        H = len(hamiltonian_paths(A, 3))

        for x_val in [1, 2]:
            ut_val = 0
            for partition, coeff in UT.items():
                if any(k % 2 == 0 for k in partition):
                    continue
                num_big = sum(1 for k in partition if k >= 3)
                ut_val += coeff * x_val**num_big

            ip_val = sum(c * x_val**k for k, c in enumerate(ip))
            print(f"    H={H}, x={x_val}: U_T(p1=1,p_odd=x)={ut_val}, I(Omega,x)={ip_val}")


if __name__ == "__main__":
    main()
