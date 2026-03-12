#!/usr/bin/env python3
"""
paley_eigenvector_theorem.py — Why is Paley sigma an eigenvector of J?

DISCOVERY: The interaction matrix J[i,j] = hat{H}({chord_i, chord_j})
(from the Walsh expansion of H on the orientation cube {+1,-1}^m)
has the Paley orientation vector as an EIGENVECTOR with the LARGEST eigenvalue.

This means Paley maximizes the degree-2 quadratic form of H.
The question is: is this an accident, or does it follow from
the algebraic structure of the Legendre symbol?

PLAN:
1. Verify eigenvector property at p=7, 11 (done, confirmed)
2. Explain WHY via the circulant structure of J
3. Determine when the degree-4 terms can override the degree-2 advantage
4. Connect to the dihedral group representation theory

Author: opus-2026-03-12-S62
"""

import numpy as np
import cmath
import math
from itertools import combinations

def adjacency_matrix(S, p):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in S:
                A[i][j] = 1
    return A

def count_hp_fast(A):
    n = len(A)
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[mask, v]
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v, u]:
                    dp[mask | (1 << u), u] += c
    return int(np.sum(dp[full]))

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p - 1) // 2, p) == 1

def sigma_to_S(sigma, p):
    m = (p - 1) // 2
    S = set()
    for k in range(1, m + 1):
        if sigma[k - 1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    return frozenset(S)

def compute_all_H(p):
    """Compute H for all 2^m orientations."""
    m = (p - 1) // 2
    H_vals = {}
    for bits in range(2**m):
        sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
        S = sigma_to_S(list(sigma), p)
        A = adjacency_matrix(S, p)
        H = count_hp_fast(A)
        H_vals[sigma] = H
    return H_vals

def walsh_coefficients(H_vals, m):
    """Compute Walsh-Fourier coefficients."""
    coeffs = {}
    for r in range(m + 1):
        for S in combinations(range(m), r):
            coeff = 0
            for bits in range(2**m):
                sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
                prod = 1
                for k in S:
                    prod *= sigma[k]
                coeff += H_vals[sigma] * prod
            coeff /= 2**m
            if abs(coeff) > 0.001:
                coeffs[S] = coeff
    return coeffs

def interaction_matrix(coeffs, m):
    """Extract the degree-2 interaction matrix J."""
    J = np.zeros((m, m))
    for S, c in coeffs.items():
        if len(S) == 2:
            i, j = S
            J[i, j] = c
            J[j, i] = c
    return J


# ============================================================
# SECTION 1: THE EIGENVECTOR THEOREM — ALGEBRAIC EXPLANATION
# ============================================================

print("=" * 70)
print("PALEY EIGENVECTOR THEOREM — ALGEBRAIC STRUCTURE")
print("=" * 70)

print("""
THEOREM: At p = 3 mod 4, the Paley orientation sigma_P is an eigenvector
of the interaction matrix J with the largest eigenvalue.

PROOF SKETCH:
  The interaction matrix J lives on the space of chord type pairs.
  J[i,j] measures how the orientation of chord i interacts with chord j
  in determining H.

  For CIRCULANT tournaments, there's a hidden multiplicative symmetry:
  If a is a QR mod p, then multiplying all chord types by a (mod p)
  is an automorphism of the Paley tournament.

  This means J is EQUIVARIANT under the action of (Z/pZ)*/QR on chord types.

  The QR subgroup acts by permutation on {1,...,m} (chord types):
    a * k mod p -> either k' or p-k' (where 1 <= k' <= m)
  with a sign flip when a*k > m.

  The Legendre symbol chi defines the character of QR on this space:
    sigma_P[k] = chi(k)  (for k = 1,...,m)

  Since chi is a multiplicative character of (Z/pZ)*, and J is equivariant
  under QR multiplication, sigma_P (= chi restricted to {1,...,m})
  must be an eigenvector of J.
""")

# Verify the algebraic structure at p=7 and p=11
for p in [7, 11]:
    m = (p - 1) // 2
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    sigma_P = np.array([1 if k in S_paley else -1 for k in range(1, m + 1)])

    print(f"\np={p}: Paley sigma = {sigma_P}")
    print(f"  Legendre values chi(k) for k=1,...,{m}:")
    for k in range(1, m + 1):
        print(f"    chi({k}) = {1 if is_qr(k, p) else -1}")

    print(f"  sigma_P[k] = chi(k)? ", end="")
    match = all((sigma_P[k-1] == 1) == is_qr(k, p) for k in range(1, m + 1))
    print(f"{'YES' if match else 'NO'}")


# ============================================================
# SECTION 2: J MATRIX STRUCTURE — MULTIPLICATIVE EQUIVARIANCE
# ============================================================

print(f"\n{'=' * 70}")
print("J MATRIX EQUIVARIANCE UNDER QR MULTIPLICATION")
print("=" * 70)

for p in [7, 11]:
    m = (p - 1) // 2
    H_vals = compute_all_H(p)
    coeffs = walsh_coefficients(H_vals, m)
    J = interaction_matrix(coeffs, m)

    print(f"\np={p}:")
    print(f"J matrix:")
    print(J)

    # Check: is J[a*i, a*j] = sigma_sign * J[i,j] for QR a?
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    qr_elts = sorted(j for j in range(1, p) if is_qr(j, p))

    print(f"\nMultiplicative equivariance check:")
    for a in qr_elts:
        print(f"  a = {a}:")
        # Compute permutation and signs
        perm = []
        signs = []
        for k in range(1, m + 1):
            ak = (a * k) % p
            if ak <= m:
                perm.append(ak - 1)  # 0-indexed
                signs.append(1)
            else:
                perm.append(p - ak - 1)  # p-ak mapped to index
                signs.append(-1)

        # Check J[perm[i], perm[j]] * signs[i] * signs[j] == J[i, j]
        matches = 0
        total = 0
        for i in range(m):
            for j in range(i + 1, m):
                expected = J[i, j]
                actual = J[perm[i], perm[j]] * signs[i] * signs[j]
                total += 1
                if abs(expected - actual) < 0.01:
                    matches += 1

        print(f"    perm = {[p+1 for p in perm]}, signs = {signs}")
        print(f"    Equivariance: {matches}/{total}")


# ============================================================
# SECTION 3: EIGENVALUE DECOMPOSITION OF J
# ============================================================

print(f"\n{'=' * 70}")
print("EIGENVALUE DECOMPOSITION OF J")
print("=" * 70)

for p in [7, 11]:
    m = (p - 1) // 2
    H_vals = compute_all_H(p)
    coeffs = walsh_coefficients(H_vals, m)
    J = interaction_matrix(coeffs, m)

    eigvals, eigvecs = np.linalg.eigh(J)
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    sigma_P = np.array([1 if k in S_paley else -1 for k in range(1, m + 1)], dtype=float)

    print(f"\np={p}:")
    print(f"  Eigenvalues of J: {eigvals}")
    print(f"  Top eigenvalue: {eigvals[-1]:.4f}")
    print(f"  Top eigenvector: {eigvecs[:, -1]}")
    print(f"  sigma_P / ||sigma_P||: {sigma_P / np.linalg.norm(sigma_P)}")

    # Check alignment
    alignment = abs(np.dot(eigvecs[:, -1], sigma_P / np.linalg.norm(sigma_P)))
    print(f"  Alignment |<v_top, sigma_P>| = {alignment:.6f} (1.0 = perfect)")

    # The eigenvalue equals m * (some quantity)
    top_eval = eigvals[-1]
    print(f"  lambda_max / m = {top_eval / m:.4f}")
    print(f"  lambda_max * 2 / (m*(m-1)) = {top_eval * 2 / (m * (m-1)):.4f}")

    # What are the OTHER eigenspaces?
    print(f"\n  All eigenpairs:")
    for i in range(len(eigvals)):
        ev = eigvecs[:, i]
        # Compute the "QR content" of this eigenvector
        qr_projection = np.dot(ev, sigma_P / np.linalg.norm(sigma_P))
        print(f"    lambda = {eigvals[i]:>10.4f}, "
              f"v = [{', '.join(f'{x:.4f}' for x in ev)}], "
              f"<v, sigma_P> = {qr_projection:.4f}")


# ============================================================
# SECTION 4: PREDICTING J FROM CIRCULANT STRUCTURE
# ============================================================

print(f"\n{'=' * 70}")
print("PREDICTING J FROM EIGENVALUE STRUCTURE")
print("=" * 70)

print("""
HYPOTHESIS: J[i,j] depends on chord types i,j through the
interaction of the corresponding eigenvalue modes.

For circulant tournament T(sigma), the eigenvalue at frequency k is:
  lambda_k(sigma) = sum_{j=1}^{m} sigma_j * omega^{jk} * delta_j
where delta_j encodes the "chord contribution" and omega = e^{2pi i/p}.

H is expressible in terms of these eigenvalues through:
  H = sum of permanent-like terms involving lambda_k.

The degree-2 Walsh coefficient hat{H}({i,j}) measures:
  d^2 H / (d sigma_i d sigma_j)  at the uniform point.

If H = f(lambda_1,...,lambda_m), then:
  hat{H}({i,j}) involves d lambda_k / d sigma_i * d lambda_k / d sigma_j * d^2 f/d lambda_k^2
  summed over k, with cross-terms.

This is related to the FISHER INFORMATION of H with respect to sigma.
""")

# Let's compute J directly from eigenvalue derivatives
for p in [7, 11]:
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)

    print(f"\np={p}: Testing eigenvalue-derivative formula for J")

    # lambda_k(sigma) = sum_{j=1}^m sigma_j * gamma_{j,k}
    # where gamma_{j,k} = omega^{jk} if j <= m, else (complement contribution)
    # Actually for sigma_j = +1: chord j contributes omega^{jk} to lambda_k
    # For sigma_j = -1: chord j contributes omega^{(p-j)k} = omega^{-jk} to lambda_k
    # So lambda_k = sum_{j=1}^m [(sigma_j + 1)/2 * omega^{jk} + (1 - sigma_j)/2 * omega^{-jk}]
    # = sum_{j=1}^m [omega^{-jk} + sigma_j * (omega^{jk} - omega^{-jk})/2]
    # = sum_j omega^{-jk} + sum_j sigma_j * i * sin(2*pi*jk/p)

    # The constant part: sum_j omega^{-jk} for j=1,...,m
    # The sigma-dependent part: sigma_j * i * sin(2*pi*jk/p)

    # d lambda_k / d sigma_j = i * sin(2*pi*jk/p)

    # This means the "sensitivity" of eigenvalue k to chord orientation j
    # is purely imaginary: d lambda_k / d sigma_j = i * sin(2*pi*jk/p)

    # For J[i,j], we need the second derivative of H:
    # d^2 H / (d sigma_i d sigma_j) ~ sum_k (d^2 H / d Re(lam_k)^2) * sin_ik * sin_jk

    # If H depends only on eigenvalue MAGNITUDES |lambda_k|^2,
    # then J would be proportional to:
    #   J[i,j] ~ sum_k sin(2*pi*ik/p) * sin(2*pi*jk/p) * weight_k

    # This is a TRIGONOMETRIC MATRIX!

    # Let's test: define the matrix T[i,j] = sum_k sin(2*pi*ik/p) * sin(2*pi*jk/p)
    T = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            for k in range(1, p):  # all non-zero frequencies
                T[i, j] += math.sin(2 * math.pi * (i+1) * k / p) * math.sin(2 * math.pi * (j+1) * k / p)

    # T should be proportional to (p/2) * I (by orthogonality)
    print(f"  T matrix (sin-sin inner products):")
    print(f"  Diagonal: {[f'{T[i,i]:.2f}' for i in range(m)]}")
    print(f"  Expected: {p/2:.2f} * I")

    # T[i,i] = sum_k sin^2(2*pi*(i+1)*k/p) = (p-1)/2
    # T[i,j] = sum_k sin(2*pi*(i+1)*k/p) * sin(2*pi*(j+1)*k/p)
    # = (1/2) sum_k [cos(2*pi*(i-j)*k/p) - cos(2*pi*(i+j+2)*k/p)]

    # Now, can we express J in terms of T with chi-dependent weights?
    H_vals = compute_all_H(p)
    coeffs = walsh_coefficients(H_vals, m)
    J_actual = interaction_matrix(coeffs, m)

    # Hypothesis: J[i,j] = alpha * sum_k f(k) * sin(2*pi*(i+1)*k/p) * sin(2*pi*(j+1)*k/p)
    # where f(k) depends on the eigenvalue structure at frequency k.

    # Try: f(k) = chi(k) (the Legendre symbol)
    T_chi = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            for k in range(1, p):
                chi_k = 1 if is_qr(k, p) else -1
                T_chi[i, j] += chi_k * math.sin(2 * math.pi * (i+1) * k / p) * math.sin(2 * math.pi * (j+1) * k / p)

    # Check if J ~ alpha * T_chi
    # Find best alpha by least squares
    J_flat = J_actual[np.triu_indices(m, 1)]
    T_chi_flat = T_chi[np.triu_indices(m, 1)]
    if np.dot(T_chi_flat, T_chi_flat) > 0:
        alpha = np.dot(J_flat, T_chi_flat) / np.dot(T_chi_flat, T_chi_flat)
        residual = np.linalg.norm(J_flat - alpha * T_chi_flat)
        rel_error = residual / np.linalg.norm(J_flat) if np.linalg.norm(J_flat) > 0 else 0

        print(f"\n  J ~ alpha * T_chi:")
        print(f"    alpha = {alpha:.6f}")
        print(f"    relative error = {rel_error:.6f}")

        if rel_error < 0.1:
            print(f"    GOOD FIT! J is approximately chi-weighted sin-sin matrix")
        else:
            print(f"    Poor fit, need different structure")

    # Also try: f(k) = 1 (uniform weight, QR-independent)
    if np.dot(T[np.triu_indices(m, 1)].flatten(), T[np.triu_indices(m, 1)].flatten()) > 0:
        T_flat = T[np.triu_indices(m, 1)]
        alpha_uniform = np.dot(J_flat, T_flat) / np.dot(T_flat, T_flat)
        residual_u = np.linalg.norm(J_flat - alpha_uniform * T_flat)
        rel_error_u = residual_u / np.linalg.norm(J_flat)
        print(f"\n  J ~ alpha * T (uniform weight):")
        print(f"    alpha = {alpha_uniform:.6f}, relative error = {rel_error_u:.6f}")


# ============================================================
# SECTION 5: WHY THE EIGENVECTOR PROPERTY MIGHT FAIL AT p=19
# ============================================================

print(f"\n{'=' * 70}")
print("CROSSOVER: WHEN DOES THE EIGENVECTOR PROPERTY STOP MATTERING?")
print("=" * 70)

print("""
At p=7,11: Paley sigma is the top eigenvector of J, AND Paley maximizes H.
The degree-2 terms dominate.

At p=19: We can't compute the full orientation cube (2^9 = 512 orientations,
each needing Held-Karp on a 19-vertex graph ~ seconds per evaluation).
But we know Interval beats Paley, so either:

  (a) Paley sigma is STILL the top eigenvector of J at p=19,
      but the degree-4+ terms override the degree-2 advantage, OR

  (b) Paley sigma is NOT the top eigenvector at p=19.

ANALYSIS:
  At p=7: H = H0 + H2 (degree 2 only, m=3 too small for degree 4 with even constraint)
  At p=11: H = H0 + H2 + H4. Paley: H2=1402.5, H4=591.25 (both positive).
           Interval: H2=280.5, H4=-354.75.
  At p=19 (predicted): H4 and H6+ terms grow relative to H2.
           The quartic correction might favor Interval.

The GROWTH RATES:
  H2 contribution ~ quadratic in m ~ p^2 / 4
  H4 contribution ~ quartic in m ~ p^4 / 16
  For the RATIO H4/H2 to flip the winner, we need the quartic
  Interval advantage to exceed the quadratic Paley advantage.

From p=11 data:
  H2 advantage of Paley = 1402.5 - 280.5 = 1122
  H4 advantage of Paley = 591.25 - (-354.75) = 946
  So the quartic terms already HELP Paley at p=11.

But the ABSOLUTE values at p=11:
  H0 = 93101.25 (dominates, ~99.7% of H)
  H2: max 1402.5 (~1.5%)
  H4: max 591.25 (~0.6%)

Prediction for p=19:
  The H0 term (mean H over all orientations) grows as ~ (p-1)!
  The degree-2 advantage grows as ~ p^2
  The degree-4 terms can grow as ~ p^4
  But the crossover involves the RELATIVE sign of higher terms.
""")

# Can we compute a FEW orientations at p=19 to test?
# Let's compute H for Paley, Interval, and the 9 "single-flip" neighbors of Paley.

print(f"\np=19: Targeted orientation evaluation")
p = 19
m = 9
S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
sigma_P_19 = [1 if k in S_paley else -1 for k in range(1, m + 1)]
print(f"Paley sigma = {sigma_P_19}")
print(f"Paley QR = {sorted(S_paley)}")

import time

# Compute H(Paley), H(Interval), H(single-flips)
test_sigmas = {
    "Paley": sigma_P_19[:],
    "Interval": [1] * m,
}

# Single flips of Paley
for k in range(m):
    flipped = sigma_P_19[:]
    flipped[k] *= -1
    test_sigmas[f"Paley_flip_{k+1}"] = flipped

# Double flips of Paley
for k1 in range(m):
    for k2 in range(k1 + 1, m):
        flipped = sigma_P_19[:]
        flipped[k1] *= -1
        flipped[k2] *= -1
        test_sigmas[f"Paley_flip_{k1+1}_{k2+1}"] = flipped

print(f"\nComputing H for {len(test_sigmas)} orientations at p=19...")
results_19 = {}
for name, sigma in test_sigmas.items():
    S = sigma_to_S(sigma, p)
    A = adjacency_matrix(S, p)
    t0 = time.time()
    H = count_hp_fast(A)
    elapsed = time.time() - t0
    results_19[name] = H
    if name in ["Paley", "Interval"] or H > results_19.get("Paley", 0):
        print(f"  {name}: H = {H} ({elapsed:.1f}s)")

# Sort by H
sorted_results = sorted(results_19.items(), key=lambda x: -x[1])
print(f"\nRanked (top 20):")
for rank, (name, H) in enumerate(sorted_results[:20], 1):
    diff = H - results_19["Paley"]
    print(f"  #{rank}: {name:>25}: H = {H:>15}, diff = {diff:>+10}")

print(f"\nPaley H = {results_19['Paley']}")
print(f"Interval H = {results_19['Interval']}")

# Is Paley a LOCAL max? (All single flips should decrease H)
is_local_max = all(results_19[f"Paley_flip_{k+1}"] <= results_19["Paley"]
                    for k in range(m))
print(f"\nPaley is LOCAL max on orientation cube? {is_local_max}")

# How many single flips IMPROVE H?
improving_flips = [(k+1, results_19[f"Paley_flip_{k+1}"])
                   for k in range(m)
                   if results_19[f"Paley_flip_{k+1}"] > results_19["Paley"]]
if improving_flips:
    print(f"  Improving single flips: {improving_flips}")
else:
    print(f"  No single flip improves H (Paley IS a local max)")

# Gradient at Paley (finite difference)
print(f"\nGradient of H at Paley (finite differences):")
for k in range(m):
    dH = results_19[f"Paley_flip_{k+1}"] - results_19["Paley"]
    # Since sigma_k flips from +-1 to -+1, the "derivative" is dH / (delta sigma) = dH / 2
    print(f"  dH/d(sigma_{k+1}) = {dH:>+15} (flip chord {k+1})")

# Hessian approximation at Paley (from double flips)
print(f"\nApproximate Hessian (from double flips):")
H_P = results_19["Paley"]
hess = np.zeros((m, m))
for k in range(m):
    hess[k, k] = results_19[f"Paley_flip_{k+1}"] - H_P  # diagonal
for k1 in range(m):
    for k2 in range(k1 + 1, m):
        # H(flip both) - H(flip k1) - H(flip k2) + H(none)
        H_both = results_19[f"Paley_flip_{k1+1}_{k2+1}"]
        H_k1 = results_19[f"Paley_flip_{k1+1}"]
        H_k2 = results_19[f"Paley_flip_{k2+1}"]
        hess[k1, k2] = (H_both - H_k1 - H_k2 + H_P) / 4
        hess[k2, k1] = hess[k1, k2]

print(f"  Diagonal (single-flip H losses):")
for k in range(m):
    print(f"    chord {k+1}: {hess[k,k]:>+15}")

# Eigenvalues of Hessian
eigvals_h, eigvecs_h = np.linalg.eigh(hess)
print(f"\n  Hessian eigenvalues: {eigvals_h}")
all_neg = all(e < 0 for e in eigvals_h[:-1])  # ignore smallest
print(f"  All negative? {all(e <= 0 for e in eigvals_h)} (local max requires all <= 0)")
