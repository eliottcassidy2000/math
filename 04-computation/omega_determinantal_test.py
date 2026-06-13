#!/usr/bin/env python3
"""
Test whether I(Omega(T), x) has a determinantal representation.

KEY INSIGHT: A polynomial p(x) = det(I + x*M) has all real roots if M is
a real symmetric positive semidefinite matrix (all eigenvalues of M are
nonneg real, so roots of p are -1/lambda_i, all real negative).

More generally, p(x) = det(I + x*M) for M symmetric has all real roots.

QUESTION: Can I(Omega(T), x) be written as det(I + x*M_T) for some
symmetric matrix M_T derived from the tournament T?

If alpha_k = coefficient of x^k, then:
  det(I + x*M) = sum_k e_k(eigenvalues of M) * x^k

where e_k is the k-th elementary symmetric polynomial.

So the question reduces to: are the alpha_k equal to e_k of some
nonneg real sequence?

This is equivalent to the Newton inequalities being satisfied AND
alpha_k = e_k(lambda_1, ..., lambda_n) for real lambda_i.

The lambda_i can be recovered: they are the roots of the characteristic
polynomial det(M - lambda*I) = 0, which equals
  (-1)^n * (lambda^n - alpha_1 * lambda^{n-1} + alpha_2 * lambda^{n-2} - ...)
if the alpha_k are elementary symmetric polynomials.

Wait — the roots of det(I + x*M) = 0 are x = -1/lambda_i.
So if I(Omega(T), x) = det(I + x*M), the roots of I(Omega,x) are -1/lambda_i
where lambda_i are eigenvalues of M. All real roots of I iff all eigenvalues of M
are real and nonzero, which is automatic for real symmetric M with nonzero eigenvalues.

But alpha_0 = 1, and alpha_k count independent sets, so they're nonneg integers.
The polynomial I(Omega,x) = 1 + alpha_1*x + alpha_2*x^2 + ... + alpha_d*x^d.

For this to equal det(I + x*M) where M is n x n (n = alpha_1 = |V(Omega)|):
  alpha_k = e_k(eigenvalues of M)

This means: the alpha_k ARE the elementary symmetric polynomials of the
eigenvalues of M. This is a strong constraint.

By Newton's identities, this is equivalent to checking that the sequence
(alpha_0, alpha_1, ..., alpha_d, 0, 0, ...) satisfies certain relations.

Actually, the simplest test: if I(Omega(T), x) = det(I + x*M), then the
polynomial has degree exactly rank(M), and all roots are real and negative
(when M is positive semidefinite). So:

1. ALL ROOTS REAL AND NEGATIVE => necessary for determinantal representation
2. alpha_k = e_k(r_1, ..., r_n) where r_i = eigenvalues of M >= 0

Test: given the roots -1/r_1, ..., -1/r_d of I(Omega,x), check that
r_1, ..., r_d > 0 (they should be since roots of I are negative).
Then M would have eigenvalues r_1,...,r_d,0,...,0.

THIS ALWAYS WORKS if all roots of I(Omega,x) are real and negative!
Because: given roots x_1,...,x_d < 0, set r_i = -1/x_i > 0, then
  I(Omega,x) = prod(1 + r_i * x) = det(I + x * diag(r_1,...,r_d))

So a determinantal representation ALWAYS EXISTS if roots are real and negative.
The question is whether M can be chosen to be a NATURAL matrix (e.g., related to
the tournament adjacency matrix or the conflict graph structure).

Let's instead test a STRONGER property: can we find a STRUCTURED determinantal
representation where M is derived from T in a natural way?

CANDIDATE: M = some function of the adjacency matrix of Omega(T).

Actually, for ANY graph G, the MATCHING polynomial mu(G,x) = det(xI - A(G))
where A(G) is the adjacency matrix. Heilmann-Lieb: matching poly has all real roots.

But I(G,x) != matching poly in general.

For LINE GRAPHS, I(G,x) IS related to the matching polynomial of the root graph.
Omega(T) is sometimes a line graph (n<=5) but not always.

NEW IDEA: The independence polynomial of a CLAW-FREE graph can be written as
a "theta function" using the Chudnovsky-Seymour fractional representation.
Maybe there's a matrix M_T such that I(Omega(T), x) = permanent-like formula
involving T's adjacency matrix.

Let me just compute and look for patterns.

Author: opus-2026-03-06-S17
"""

import sys
import os
import random
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, random_tournament,
                             find_odd_cycles, conflict_graph)

random.seed(42)


def indep_poly_coeffs(adj):
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def tournament_adj_matrix(T):
    """Return the antisymmetric matrix A where A[i][j] = T[i][j] - T[j][i]."""
    n = len(T)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i][j] = T[i][j] - T[j][i]
    return A


def test_det_representation(T):
    """
    Test various matrix constructions M_T to see if det(I + x*M_T) = I(Omega(T), x).
    """
    n = len(T)
    cycles = find_odd_cycles(T)
    if not cycles:
        return None

    cg = conflict_graph(cycles)
    ip = indep_poly_coeffs(cg)
    m = len(cg)  # number of odd cycles

    if len(ip) <= 1:
        return None

    # The "eigenvalues" we need: roots of I(Omega,x) are at -1/r_i
    # So det(I + x*diag(r_1,...,r_m)) = I(Omega,x) with r_i > 0
    roots = np.roots(list(reversed(ip)))
    real_roots = all(abs(r.imag) < 1e-8 for r in roots)
    neg_roots = all(r.real < -1e-10 for r in roots)

    if not (real_roots and neg_roots):
        return {'real_roots': False, 'ip': ip}

    r_values = sorted([-1.0/r.real for r in roots])

    # Now test: is there a NATURAL M such that eigenvalues of M are r_values?
    # Candidate 1: M = f(A_Omega) where A_Omega = adjacency matrix of conflict graph
    A_omega = np.array(cg, dtype=float)
    eig_omega = sorted(np.linalg.eigvalsh(A_omega))

    # Candidate 2: M = cI - A_Omega for some c (complement-like)
    # If I(G,x) = det(I + x*(cI - A(G))), eigenvalues are c - lambda_i(G)
    # We need c - lambda_i(G) = r_i, so c = r_i + lambda_i(G) for all i
    # This works only if r_i + lambda_i(G) is constant

    # Check if r_values + eigenvalues of Omega is approximately constant
    if len(r_values) == m:
        sums = [r_values[i] + eig_omega[i] for i in range(m)]
        sum_spread = max(sums) - min(sums)
    else:
        sum_spread = float('inf')

    # Candidate 3: Can the antisymmetric tournament matrix produce the right eigenvalues?
    # A_T is antisymmetric, eigenvalues are purely imaginary.
    # A_T^2 is negative semidefinite, eigenvalues are -(singular values)^2
    # -A_T^2 is PSD, eigenvalues are (singular values)^2
    A_T = tournament_adj_matrix(T)
    neg_A2 = -A_T @ A_T
    eig_A2 = sorted(np.linalg.eigvalsh(neg_A2))

    # Candidate 4: Laplacian of Omega
    deg_omega = np.diag(np.sum(A_omega, axis=0))
    L_omega = deg_omega - A_omega
    eig_L = sorted(np.linalg.eigvalsh(L_omega))

    return {
        'ip': ip,
        'r_values': r_values,
        'eig_omega': eig_omega[:5],  # first 5
        'sum_spread': sum_spread,
        'eig_A2': eig_A2[:min(5, n)],
        'eig_L_omega': eig_L[:5],
        'real_roots': True,
        'm': m,
        'n': n,
    }


# ============================================================
print("=" * 70)
print("DETERMINANTAL REPRESENTATION TEST FOR I(Omega(T), x)")
print("=" * 70)

for n in [5, 6, 7]:
    print(f"\n{'='*60}")
    print(f"  n = {n}")
    print(f"{'='*60}")

    count = min(1 << (n*(n-1)//2), 500) if n <= 6 else 300
    if n <= 5:
        mode = 'exhaustive'
        count = 1 << (n*(n-1)//2)
    else:
        mode = 'random'

    sum_spreads = []
    all_results = []

    for idx in range(count):
        if mode == 'exhaustive':
            T = tournament_from_bits(n, idx)
        else:
            T = random_tournament(n)

        result = test_det_representation(T)
        if result and result.get('real_roots'):
            all_results.append(result)
            if result['sum_spread'] < 100:
                sum_spreads.append(result['sum_spread'])

    print(f"  Analyzed: {len(all_results)} tournaments with real roots")
    if sum_spreads:
        print(f"  Sum spread (r_i + eig_omega_i = const?):")
        print(f"    min={min(sum_spreads):.4f}, max={max(sum_spreads):.4f}, "
              f"avg={sum(sum_spreads)/len(sum_spreads):.4f}")
        print(f"    Fraction with spread < 0.01: "
              f"{sum(1 for s in sum_spreads if s < 0.01)}/{len(sum_spreads)}")
        print(f"    Fraction with spread < 0.1: "
              f"{sum(1 for s in sum_spreads if s < 0.1)}/{len(sum_spreads)}")

    # Show a few examples
    if all_results:
        print(f"\n  Sample results:")
        for res in all_results[:5]:
            print(f"    IP={res['ip']}, r={[f'{r:.3f}' for r in res['r_values'][:4]]}, "
                  f"eig_Omega={[f'{e:.3f}' for e in res['eig_omega'][:4]]}, "
                  f"spread={res['sum_spread']:.4f}")

# ============================================================
# Test a different angle: can I(Omega(T), x) be expressed via
# the SKEW-ADJACENCY MATRIX of T?
# The skew-adjacency matrix S_T has S[i][j] = 1 if i->j, -1 if j->i, 0 if i=j.
# det(I + x*S_T) is a polynomial in x. Is it related to I(Omega(T), x)?
# ============================================================
print(f"\n{'='*60}")
print(f"  TEST: det(I + x*S_T) vs I(Omega(T), x)")
print(f"{'='*60}")

for n in [4, 5, 6]:
    print(f"\n  --- n = {n} ---")
    count = min(1 << (n*(n-1)//2), 200)
    matches = 0
    total = 0

    for bits in range(count):
        T = tournament_from_bits(n, bits)
        cycles = find_odd_cycles(T)
        if not cycles:
            continue
        cg = conflict_graph(cycles)
        ip = indep_poly_coeffs(cg)

        # Skew-adjacency matrix
        S = tournament_adj_matrix(T)
        # det(I + x*S)
        # This is a polynomial in x. Compute at a few points and compare.
        xs = [0, 1, 2, 3, -1]
        ip_vals = [sum(c * x**k for k, c in enumerate(ip)) for x in xs]
        det_vals = [float(np.linalg.det(np.eye(n) + x * S)) for x in xs]

        total += 1
        if all(abs(ip_vals[i] - det_vals[i]) < 1e-6 for i in range(len(xs))):
            matches += 1

    print(f"  det(I+xS_T) = I(Omega,x): {matches}/{total}")

# ============================================================
# What IS det(I + x*S_T)?
# S_T antisymmetric => eigenvalues are purely imaginary: +/- i*sigma_k
# det(I + x*S) = prod(1 + i*sigma_k*x)(1 - i*sigma_k*x) = prod(1 + sigma_k^2 * x^2)
# This is a polynomial in x^2 only! Even-degree terms only.
# So det(I + x*S_T) = 1 + c_2*x^2 + c_4*x^4 + ...
# while I(Omega,x) = 1 + alpha_1*x + alpha_2*x^2 + ...
# They can only agree if alpha_1 = 0, which means no odd cycles. No match.

# Better: try det(I + x * S_T^2) or det(I + x * (-S_T * S_T^T))
# ============================================================
print(f"\n{'='*60}")
print(f"  TEST: det(I + x*(-S_T^2)) vs I(Omega(T), x)")
print(f"{'='*60}")

for n in [4, 5]:
    print(f"\n  --- n = {n} ---")
    m = n * (n - 1) // 2
    count = 1 << m

    for bits in range(min(count, 20)):
        T = tournament_from_bits(n, bits)
        cycles = find_odd_cycles(T)
        if not cycles:
            continue
        cg = conflict_graph(cycles)
        ip = indep_poly_coeffs(cg)

        S = tournament_adj_matrix(T)
        M = -S @ S  # PSD matrix

        # det(I + x*M) as polynomial
        eig_M = sorted(np.linalg.eigvalsh(M), reverse=True)
        # Evaluate det(I + x*M) at several points
        xs = [0, 1, 2, -1, 0.5]
        det_vals = [float(np.linalg.det(np.eye(n) + x * M)) for x in xs]
        ip_vals = [sum(c * x**k for k, c in enumerate(ip)) for x in xs]

        print(f"  bits={bits}: IP={ip}, det_vals={[f'{v:.1f}' for v in det_vals[:3]]}, "
              f"ip_vals={[f'{v:.1f}' for v in ip_vals[:3]]}, "
              f"eig_M={[f'{e:.2f}' for e in eig_M[:4]]}")

# ============================================================
# Let's try: what MATRIX would give I(Omega(T),x) = det(I + x*M)?
# For degree-2 polys: I = 1 + alpha_1*x + alpha_2*x^2
# We need e_1(lambdas) = alpha_1, e_2(lambdas) = alpha_2
# For a 2x2 PSD matrix: eigenvalues lambda_1, lambda_2 >= 0
# lambda_1 + lambda_2 = alpha_1, lambda_1*lambda_2 = alpha_2
# These are the roots of t^2 - alpha_1*t + alpha_2 = 0
# Real positive roots iff alpha_1^2 >= 4*alpha_2 (discriminant) and alpha_1 >= 0
# This is exactly the real-rootedness condition!
# ============================================================
print(f"\n{'='*60}")
print(f"  EIGENVALUE RECOVERY: I(Omega,x) = det(I + x*diag(lambda))")
print(f"{'='*60}")

for n in [5, 6]:
    print(f"\n  --- n = {n} ---")
    count = min(1 << (n*(n-1)//2), 200)

    for bits in range(count):
        T = tournament_from_bits(n, bits)
        cycles = find_odd_cycles(T)
        if not cycles:
            continue
        cg = conflict_graph(cycles)
        ip = indep_poly_coeffs(cg)
        if len(ip) <= 2:
            continue

        # Recover eigenvalues: roots of t^d - alpha_1*t^{d-1} + alpha_2*t^{d-2} - ...
        # where d = degree = alpha(Omega)
        d = len(ip) - 1
        char_coeffs = []
        for k in range(d + 1):
            char_coeffs.append((-1)**k * ip[k] if k < len(ip) else 0)
        # Reverse for np.roots (highest degree first)
        char_poly = [(-1)**d] + [(-1)**(d-k) * ip[k] for k in range(1, d+1)]
        char_poly_rev = list(reversed(char_poly))

        lambdas = np.roots(char_poly_rev)
        all_real = all(abs(l.imag) < 1e-6 for l in lambdas)
        all_pos = all(l.real > -1e-6 for l in lambdas)

        if not all_real or not all_pos:
            print(f"  bits={bits}: IP={ip}, lambdas NOT all real positive: "
                  f"{[f'{l:.3f}' for l in lambdas]}")
            break
    else:
        print(f"  All tested: eigenvalues always real and positive "
              f"(consistent with PSD determinantal representation)")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
