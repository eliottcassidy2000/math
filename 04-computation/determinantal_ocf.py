#!/usr/bin/env python3
"""
CREATIVE EXPLORATION: Determinantal form of OCF.

OCF says H(T) = I(Omega(T), 2) = sum_S 2^|S| where S ranges over
independent sets of vertex-disjoint odd directed cycles.

Björklund's approach: cycle covers via inclusion-exclusion and determinants.
For UNDIRECTED graphs, the independence polynomial I(G, x) = det(I + xA_G)
when G is a LINE GRAPH (Heilmann-Lieb). But Omega(T) is not always a line graph.

However, for ANY graph G, the MATCHING polynomial mu(G, x) = det(xI - A_G).
And for line graphs L(H), I(L(H), x) = mu(H, 1/x) * x^{|V(L(H))|}... not quite.

QUESTION 1: Is there a matrix M_T (derived from tournament T) such that
  H(T) = det(I + 2*M_T)?

If M_T exists, then I(Omega(T), x) = det(I + x*M_T) would give real-rootedness
iff M_T is positive semidefinite. The failure at n=9 means M_T is NOT always PSD,
but the determinant might still give a correct formula.

QUESTION 2: The Grinberg-Stanley proof uses:
  ham(T) = [z^n] exp(sum_{k odd} tr(A^k) z^k / k)
         = [z^n] exp(arctanh_odd(zA))

This is the coefficient of z^n in a formal power series. Can we express this
as a determinant?

QUESTION 3: The Irving-Omar determinantal formula:
  W_D(z) = det(I + zXA_bar) / det(I - zXA)
where X = diag(x_1,...,x_n) and A is the adjacency. Setting all x_i = 1:
  W_D(z) = det(I + zA_bar) / det(I - zA)

For tournaments, A_bar = A^T (transpose/converse), so:
  W_D(z) = det(I + zA^T) / det(I - zA) = det(I + zA^T) / det(I - zA)

The coefficient [z^n] of this rational function might give H(T) / something.

Let's compute all of this for the counterexample and for comparison tournaments.

Author: opus-2026-03-06-S19
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament
import numpy as np
from numpy.linalg import det, inv, eigvals
from itertools import combinations, permutations

# The counterexample
T_ce = [
    [0, 1, 0, 1, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 1, 1, 1, 1, 0],
    [0, 0, 1, 0, 0, 1, 0, 1, 0],
    [1, 1, 0, 1, 0, 0, 0, 1, 0],
    [1, 1, 0, 0, 1, 0, 1, 1, 1],
    [0, 1, 0, 1, 1, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 1, 1, 0],
]
n = 9

def hamiltonian_count(T, n):
    h = 0
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            h += 1
    return h

def irving_omar_W(T, n, z_val):
    """Compute W_D(z) = det(I + z*A^T) / det(I - z*A) at z=z_val."""
    A = np.array(T, dtype=float)
    I_n = np.eye(n)
    num = det(I_n + z_val * A.T)
    denom = det(I_n - z_val * A)
    if abs(denom) < 1e-15:
        return None
    return num / denom

def irving_omar_series(T, n, max_terms=12):
    """Compute W_D(z) as power series in z, return coefficients."""
    A = np.array(T, dtype=float)
    I_n = np.eye(n)
    # W(z) = det(I + zA^T) * det(I - zA)^{-1}
    # = det(I + zA^T) * sum_{k>=0} z^k * something
    # Better: compute numerically at many points, then polynomial interpolation

    # Use Cauchy integral / discrete Fourier
    N = max_terms + 5
    coeffs = np.zeros(N)
    # Evaluate at roots of unity
    for j in range(N):
        z = 0.1 * np.exp(2j * np.pi * 1j * j / N)
        W = irving_omar_W(T, n, z)
        if W is None:
            return None
        coeffs[0] += W
        # This is getting complicated. Use direct power series.

    # Alternative: expand det(I+zA^T)/det(I-zA) via log
    # log W(z) = log det(I+zA^T) - log det(I-zA)
    #          = tr log(I+zA^T) - tr log(I-zA)
    #          = sum_{k>=1} (-1)^{k+1} z^k tr((A^T)^k)/k + sum_{k>=1} z^k tr(A^k)/k
    #          = sum_{k>=1} z^k [tr(A^k)/k + (-1)^{k+1} tr((A^T)^k)/k]
    # For tournaments: tr(A^k) = tr((A^T)^k) (same eigenvalues, A and A^T similar)
    # Wait: tr(A^k) = sum of k-th powers of eigenvalues = sum of k-th powers of eigenvalues of A^T
    # Yes, tr(A^k) = tr((A^T)^k) always.
    # So log W(z) = sum_{k>=1} z^k * tr(A^k)/k * [1 + (-1)^{k+1}]
    #             = sum_{k odd} 2*z^k * tr(A^k)/k
    # Because 1+(-1)^{k+1} = 2 for k odd, 0 for k even!

    # So log W(z) = 2 * sum_{k odd} tr(A^k) z^k / k
    # Therefore W(z) = exp(2 * sum_{k odd} tr(A^k) z^k / k)

    # And the Grinberg-Stanley formula says:
    # H(T) = [z^n] exp(sum_{k odd} tr(A^k) z^k / k) = [z^n] W(z)^{1/2}
    # So W(z)^{1/2} has [z^n] = H(T)!

    # Actually wait. Let me be more careful.
    # Grinberg-Stanley: ham(T) = sum_{sigma, all cycles odd} 2^{psi(sigma)}
    # where psi(sigma) = number of cycles. This equals I(Omega(T), 2).
    # The generating function is:
    # sum_sigma (product over cycles C of sigma) (product over arcs of C) * z^{|support|}
    # weighted by 2^{#cycles} restricted to odd-cycle-only permutations.

    # The exp formula: if we define
    # p_k = tr(A^k) = sum of k-walks returning to start = k * (# directed k-cycles)
    # Then exp(sum_{k>=1} p_k z^k / k) = sum_sigma prod_cycles z^|C| a(C)
    # = product over possible cycle types.
    # Restricting to odd k: exp(sum_{k odd} p_k z^k / k)
    # counts ONLY permutations with all odd cycles.
    # Each cycle of length k contributes z^k and there's a factor of 1 per cycle.
    # The 2^{psi} weighting... hmm, that's from the independence polynomial.

    # Let me just compute the power series numerically.
    return None

def compute_series_via_traces(T, n, max_k=None):
    """Compute exp(sum_{k odd} tr(A^k) z^k / k) as power series."""
    if max_k is None:
        max_k = n + 2
    A = np.array(T, dtype=float)

    # Compute tr(A^k) for k=1,...,max_k
    traces = []
    Ak = np.eye(n)
    for k in range(1, max_k + 1):
        Ak = Ak @ A
        traces.append(np.trace(Ak))

    # Build the odd-cycle generating function coefficients
    # log f(z) = sum_{k odd, k>=1} tr(A^k) * z^k / k
    # f(z) = exp(log f(z))
    # Compute coefficients of f(z) up to z^n using the exp formula

    # log_coeffs[k] = coefficient of z^k in log f(z)
    log_coeffs = np.zeros(max_k + 1)
    for k in range(1, max_k + 1):
        if k % 2 == 1:  # odd k
            log_coeffs[k] = traces[k-1] / k

    # Exponentiate: if log f = sum c_k z^k, then f = exp(sum c_k z^k)
    # Use the recursive formula: f[0]=1, f[m] = (1/m) sum_{k=1}^{m} k*c_k*f[m-k]
    f = np.zeros(max_k + 1)
    f[0] = 1.0
    for m in range(1, max_k + 1):
        s = 0.0
        for k in range(1, m + 1):
            s += k * log_coeffs[k] * f[m - k]
        f[m] = s / m

    return f, traces, log_coeffs

print("=" * 70)
print("DETERMINANTAL FORM OF OCF — DEEP EXPLORATION")
print("=" * 70)

# Part 1: The Irving-Omar W(z) and its square root
print("\n--- Part 1: Irving-Omar W(z) decomposition ---")

A = np.array(T_ce, dtype=float)
I_n = np.eye(n)

# log W(z) = 2 * sum_{k odd} tr(A^k) z^k / k
# So W(z)^{1/2} = exp(sum_{k odd} tr(A^k) z^k / k)
# And [z^n] W(z)^{1/2} should give something related to H(T).

f, traces, log_coeffs = compute_series_via_traces(T_ce, n)

print("  tr(A^k) for k=1,...,9:")
for k in range(1, n+1):
    print(f"    k={k}: tr(A^{k}) = {traces[k-1]:.0f}")

print(f"\n  Odd-trace GF: exp(sum_{{k odd}} tr(A^k) z^k / k)")
print(f"  Coefficients [z^0]...[z^{n}]:")
for k in range(n + 1):
    print(f"    [z^{k}] = {f[k]:.6f}")

print(f"\n  [z^{n}] = {f[n]:.6f}")

H_ce = hamiltonian_count(T_ce, n)
print(f"  H(T) = {H_ce}")
print(f"  Match: {abs(f[n] - H_ce) < 0.5}")

# Part 2: The Grinberg-Stanley formula with 2^psi weighting
print("\n--- Part 2: Grinberg-Stanley with 2^psi ---")
# The actual OCF formula is:
# H(T) = sum_{sigma: all odd cycles} 2^{#cycles(sigma)}
# = sum_{sigma: all odd cycles} prod_{cycle C} 2
# = [z^n] exp(sum_{k odd} 2 * tr(A^k) z^k / k)  ← factor of 2 per cycle
# = [z^n] exp(2 * sum_{k odd} tr(A^k) z^k / k)
# = [z^n] W(z)

# So [z^n] W(z) = H(T)!

# Compute W(z) coefficients: W = exp(2 * sum_{k odd} tr(A^k) z^k / k)
log_W = np.zeros(n + 3)
for k in range(1, n + 2):
    if k % 2 == 1:
        log_W[k] = 2 * traces[k-1] / k

W = np.zeros(n + 3)
W[0] = 1.0
for m in range(1, n + 2):
    s = 0.0
    for k in range(1, m + 1):
        s += k * log_W[k] * W[m - k]
    W[m] = s / m

print(f"  W(z) = exp(2 * sum_{{k odd}} tr(A^k) z^k / k)")
print(f"  Coefficients [z^0]...[z^{n}]:")
for k in range(n + 1):
    print(f"    [z^{k}] = {W[k]:.6f}")

print(f"\n  [z^{n}] W(z) = {W[n]:.6f}")
print(f"  H(T) = {H_ce}")
print(f"  Match: {abs(W[n] - H_ce) < 0.5}")

# Part 3: det(I + 2A^T) / det(I - 2A) at z=... wait, W(z) = det(I+zA^T)/det(I-zA)
print("\n--- Part 3: Determinantal verification ---")
# W(z) = det(I+zA^T)/det(I-zA)
# Let's verify: [z^n] det(I+zA^T)/det(I-zA) should be H(T)

# Compute via partial fractions / Laurent series
# det(I-zA) = prod (1-z*lambda_i) where lambda_i are eigenvalues of A
eigs_A = eigvals(A)
print(f"  Eigenvalues of A: {np.sort_complex(eigs_A)}")

# det(I-zA) = prod(1 - z*lambda_i)
# det(I+zA^T) = prod(1 + z*lambda_i^*) ... no, A^T has eigenvalues conj(lambda)
# For real matrix A: eigenvalues of A^T = eigenvalues of A (same char poly)
# So det(I+zA^T) = prod(1 + z*lambda_i)

# Therefore W(z) = prod((1+z*lambda_i)/(1-z*lambda_i))
# = prod((1+z*lambda_i)/(1-z*lambda_i))

# log W(z) = sum_i log((1+z*lambda_i)/(1-z*lambda_i))
#           = sum_i 2*arctanh(z*lambda_i)
#           = sum_i 2*sum_{k odd} (z*lambda_i)^k / k
#           = 2*sum_{k odd} z^k/k * sum_i lambda_i^k
#           = 2*sum_{k odd} z^k * tr(A^k) / k  ✓

print(f"\n  Verified: log W(z) = 2*sum_{{k odd}} tr(A^k) z^k / k")
print(f"  This is 2*arctanh(zA) (operator-level)")

# Part 4: Can we write H(T) = det(something)?
print("\n--- Part 4: Determinantal expression for H(T) ---")
# W(z) = det(I+zA^T)/det(I-zA)
# = det((I+zA^T)(I-zA)^{-1})  (if I-zA invertible)
# = det(Cayley_z)
# where Cayley_z = (I+zA^T)(I-zA)^{-1}

# At z=1: Cayley_1 = (I+A^T)(I-A)^{-1}
# But det(I-A) may be 0!
det_I_minus_A = det(I_n - A)
print(f"  det(I - A) = {det_I_minus_A:.6f}")

if abs(det_I_minus_A) > 1e-10:
    Cayley = (I_n + A.T) @ inv(I_n - A)
    eigs_Cayley = eigvals(Cayley)
    print(f"  Cayley eigenvalues: {np.sort_complex(eigs_Cayley)}")
    det_Cayley = det(Cayley)
    print(f"  det(Cayley) = {det_Cayley:.6f}")
    print(f"  W(1) should be infinity or H-related...")
else:
    print(f"  I-A is singular! Can't form Cayley at z=1.")

# W(z) is a rational function. The poles are at z = 1/lambda_i.
# [z^n] of W(z) is computed by residues.
print(f"\n  Poles of W(z) at z = 1/lambda_i:")
for lam in eigs_A:
    if abs(lam) > 1e-10:
        print(f"    1/({lam:.4f}) = {1/lam:.4f}")
    else:
        print(f"    1/({lam:.4f}) = infinity")

# Part 5: W(z) = det(Cayley(z)). Can we find a DIFFERENT matrix M such that
# H(T) = det(I + 2M)?
print("\n--- Part 5: Search for H(T) = det(I + 2M) ---")
# If I(Omega,x) = det(I+xM), then expanding:
# det(I+xM) = 1 + x*tr(M) + x^2*(tr(M)^2-tr(M^2))/2 + ...
# Matching with I(Omega,x) = 1 + 94x + 10x^2 + x^3:
# tr(M) = 94
# (tr(M)^2 - tr(M^2))/2 = 10 => tr(M^2) = 94^2 - 20 = 8816 - 20 = 8816
# Wait: (94^2 - tr(M^2))/2 = 10 => tr(M^2) = 94^2 - 20 = 8816 - 20 = 8816
# Hmm: 94^2 = 8836. 8836 - 20 = 8816. tr(M^2) = 8816.
# det(I+xM) for 3x3 M: 1 + s1*x + s2*x^2 + s3*x^3
# where s1 = tr(M), s2 = (tr(M)^2-tr(M^2))/2, s3 = det(M)
# s1 = 94, s2 = 10, s3 = 1

# eigenvalues of M: sum = 94, sum of products of pairs = 10, product = 1
# These are the ROOTS of t^3 - 94t^2 + 10t - 1 = 0
# Wait: char poly of M is t^3 - s1*t^2 + s2*t - s3 = t^3 - 94t^2 + 10t - 1
M_eigs = np.roots([1, -94, 10, -1])
print(f"  If I(Omega,x) = det(I+xM) for 3x3 matrix M:")
print(f"    M eigenvalues would be roots of t^3 - 94t^2 + 10t - 1 = 0:")
print(f"    {M_eigs}")
print(f"    All positive? {all(e.real > 0 and abs(e.imag) < 1e-8 for e in M_eigs)}")

# For I(Omega,x) to have real roots, M must have real eigenvalues.
# The eigenvalues of M are t^3-94t^2+10t-1=0.
# Check: this polynomial has discriminant...
# Its roots are the RECIPROCALS of the roots of I(Omega,x) (with sign flip)
# I(Omega,x) = 1+94x+10x^2+x^3 has roots r1,r2,r3
# det(I+xM) = (1+x*m1)(1+x*m2)(1+x*m3)
# So -1/m_i = r_i => m_i = -1/r_i
# Since I has complex roots, M has complex eigenvalues => M is NOT real symmetric.
# This is consistent with the counterexample failing real-rootedness.

print(f"\n  Since I(Omega,x) has complex roots, no REAL symmetric M exists.")
print(f"  The 'determinantal' form det(I+xM) = I(Omega,x) requires complex M.")

# Part 6: For comparison, check a real-rooted case
print("\n--- Part 6: Comparison with real-rooted tournament ---")
np.random.seed(42)
for trial in range(20):
    T2 = random_tournament(n)
    f2, traces2, _ = compute_series_via_traces(T2, n)
    H2 = int(round(f2[n]))

    # Quick check: H from OCF
    if H2 > 100:
        print(f"\n  Tournament #{trial}: H = {H2}")
        print(f"  tr(A^k) odd: ", end="")
        for k in [1,3,5,7,9]:
            print(f"tr(A^{k})={traces2[k-1]:.0f}", end="  ")
        print()

        A2 = np.array(T2, dtype=float)
        eigs2 = eigvals(A2)
        print(f"  |eigenvalues|: {sorted([abs(e) for e in eigs2], reverse=True)[:5]}")

        # W(z) decomposition
        log_W2 = np.zeros(n + 3)
        for k in range(1, n + 2):
            if k % 2 == 1:
                log_W2[k] = 2 * traces2[k-1] / k
        W2 = np.zeros(n + 3)
        W2[0] = 1.0
        for m in range(1, n + 2):
            s = 0.0
            for k in range(1, m + 1):
                s += k * log_W2[k] * W2[m - k]
            W2[m] = s / m

        print(f"  W coeffs: {[round(W2[k]) for k in range(n+1)]}")
        break

# Part 7: The KEY insight — connection between W(z) and independence polynomial
print("\n--- Part 7: W(z) and I(Omega, x) connection ---")
print("  Established:")
print("  W(z) = exp(2 * sum_{k odd} tr(A^k) z^k / k)")
print("  = det(I + zA^T) / det(I - zA)")
print()
print("  [z^k] of W(z) = sum over odd-cycle-only permutations sigma on k vertices")
print("                   of 2^{#cycles(sigma)} * product of arc weights")
print()
print("  In particular, [z^n] W(z) = H(T) = I(Omega(T), 2)")
print()
print("  The LOWER coefficients [z^k] for k < n count 'partial' cycle covers:")
print("  [z^k] = sum over vertex sets S of size k,")
print("          sum over odd-cycle-only permutations of S,")
print("          2^{#cycles} * product of arcs")
print()
print("  This is I(Omega(T[S]), 2) summed over all k-element subsets S!")
print("  (where T[S] is the induced subtournament on S)")

# Verify this at small k
print("\n  Verifying [z^k] = sum_{|S|=k} I(Omega(T[S]), 2):")
for k in [3, 5]:
    total = 0
    count_subsets = 0
    for subset in combinations(range(n), k):
        sub_T = [[T_ce[i][j] for j in subset] for i in subset]
        fk, _, _ = compute_series_via_traces(sub_T, k, max_k=k+2)
        total += fk[k]
        count_subsets += 1
    print(f"  k={k}: [z^{k}] = {W[k]:.1f}, sum of I(Omega(T[S]),2) = {total:.1f}, match = {abs(W[k]-total)<0.5}")

# Part 8: The Z-function approach
print("\n--- Part 8: Formal power series Z(T,x,z) ---")
print("  Define Z(T, x, z) = exp(x * sum_{k odd} tr(A^k) z^k / k)")
print("  Then Z(T, 2, z) = W(z), and [z^n] Z(T, 2, z) = H(T)")
print()
print("  Z(T, x, z) = det((I + zA^T)^{x/2}) / det((I - zA)^{x/2})")
print("  (formal fractional power of determinant)")
print()
print("  I(Omega(T), x) should be related to [z^n] Z(T, x, z)!")
print("  But Z(T, x, z) involves ALL subtournaments, not just cycle structure.")

# Compute [z^n] Z(T, x, z) for general x
# Z(T,x,z) = exp(x * sum_{k odd} t_k z^k / k) where t_k = tr(A^k)
# [z^n] Z = sum over compositions of n into odd parts (k1,...,km)
#           of x^m * prod(t_{ki}/ki) / m!  ... no, more complicated.

# Actually: exp(x * g(z)) = sum_{m>=0} x^m g(z)^m / m!
# So [z^n] exp(x*g(z)) = sum_{m>=0} x^m/m! * [z^n] g(z)^m
# where g(z) = sum_{k odd} t_k z^k / k

# g(z)^m: coefficient of z^n is sum over m-tuples of odd integers (k1,...,km)
# summing to n, of prod(t_{ki}/ki).
# This counts collections of m directed odd cycles covering n vertices (with multiplicity).

# For x=2: [z^n] Z(T,2,z) = sum_{m>=0} 2^m/m! * [z^n] g(z)^m = H(T) ✓

# For general x: [z^n] Z(T,x,z) is a polynomial in x!
# The coefficient of x^m is [z^n] g(z)^m / m!
# = (1/m!) * (number of ordered m-tuples of directed odd cycles with total length n)
# divided by cycle lengths
# = (1/m!) * (number of odd-cycle-only permutations of [n] with exactly m cycles)
# / product of cycle lengths...

# Actually, [z^n] g(z)^m / m! = sum over UNORDERED m-tuples of odd cycle lengths
# summing to n, weighted by prod(t_k/k). This counts m-cycle permutations
# of [n] with all odd cycles, divided by m!... which is the number of
# such permutations with cycle-type partition, divided by the number of
# such permutations with that partition type.

# Let's just compute it numerically.
print("\n  Computing [z^n] Z(T, x, z) as polynomial in x:")
# We need g(z)^m for m = 0, 1, 2, ..., up to n/3 (since smallest odd cycle = 3)

g_coeffs = np.zeros(n + 3)
for k in range(1, n + 2):
    if k % 2 == 1 and k <= n:
        g_coeffs[k] = traces[k-1] / k

# g_coeffs: coeff of z^k in g(z)

# g(z)^m: use convolution
max_m = n // 3 + 1
g_power_n = np.zeros(max_m + 1)  # g_power_n[m] = [z^n] g(z)^m
g_power = np.zeros(n + 3)  # current g(z)^m
g_power[0] = 1.0  # g^0 = 1

for m in range(max_m + 1):
    g_power_n[m] = g_power[n] if n < len(g_power) else 0
    # Multiply by g(z) for next iteration
    new_g = np.zeros(n + 3)
    for i in range(n + 1):
        for j in range(n + 1 - i):
            if i + j < n + 3:
                new_g[i + j] += g_power[i] * g_coeffs[j]
    g_power = new_g

# [z^n] Z(T,x,z) = sum_m x^m * g_power_n[m] / m!
import math
poly_x = np.zeros(max_m + 1)
for m in range(max_m + 1):
    poly_x[m] = g_power_n[m] / math.factorial(m)

print(f"  [z^n] Z(T, x, z) = ", end="")
terms = []
for m in range(max_m + 1):
    if abs(poly_x[m]) > 0.001:
        terms.append(f"{poly_x[m]:.4f}*x^{m}")
print(" + ".join(terms))

# Evaluate at x=2
val_at_2 = sum(poly_x[m] * 2**m for m in range(max_m + 1))
print(f"  Value at x=2: {val_at_2:.4f}")
print(f"  H(T) = {H_ce}")

# THIS polynomial in x IS I(Omega(T), x)!?
# Let's check: I(Omega, x) = 1 + 94x + 10x^2 + x^3
print(f"\n  Compare with I(Omega, x) = 1 + 94x + 10x^2 + x^3:")
for m in range(max_m + 1):
    ip_coeff = [1, 94, 10, 1][m] if m <= 3 else 0
    print(f"    x^{m}: Z-coeff={poly_x[m]:.4f}, I(Omega)-coeff={ip_coeff}")

print("\n" + "=" * 70)
print("KEY DISCOVERY CHECK:")
print("Is [z^n] Z(T,x,z) = I(Omega(T), x)?")
print("=" * 70)
