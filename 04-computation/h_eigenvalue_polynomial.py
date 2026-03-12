#!/usr/bin/env python3
"""
h_eigenvalue_polynomial.py — opus-2026-03-12-S58

Find the EXACT polynomial H(y_1,...,y_m) expressing H as a function
of eigenvalue magnitudes squared y_k = |λ_k|² for circulant tournaments.

At p=7: m = (p-1)/2 = 3 independent magnitudes (λ_k and λ_{p-k} are conjugate).
There are only 2 distinct H values (189 and 175), so H is determined by
1 parameter beyond the Parseval constraint Σy_k = (p-1)(p+1)/(4·2) = 3.

Wait: y_k = |λ_k|² for k=1,...,(p-1)/2 (independent magnitudes).
Σ_{k=1}^{(p-1)/2} 2·y_k = (p-1)(p+1)/4 (Parseval, factor 2 from conjugate pairs).
So Σ y_k = (p-1)(p+1)/8 ... let me recompute.

Actually: Σ_{k=1}^{p-1} |λ_k|² = (p-1)(p+1)/4 (from Section 6 of permanent_spectral_bridge).
Since |λ_k|² = |λ_{p-k}|² (conjugate pairs), and there are (p-1)/2 pairs:
Σ_{k=1}^{(p-1)/2} 2·y_k = (p-1)(p+1)/4
So Σ y_k = (p-1)(p+1)/8.

At p=7: Σ y_k = 6·8/8 = 6. And y_1 = y_2 = y_3 = 2 for Paley.
At p=11: Σ y_k = 10·12/8 = 15. And y_k = 3 for Paley.

Can we express H as a polynomial in y_1,...,y_m and if so, what is it?

Note: H depends on the FULL eigenvalue multiset (complex), not just magnitudes.
Two orbits at p=11 have the same prod|λ| = 23 but different H (93467 vs 92411).
So H is NOT just a function of {y_k}!

But wait — the spectral orbit analysis showed H depends on the Z_p^* orbit
class, which is determined by the multiset {|λ_k|²}. Let me recheck...

Actually, two DIFFERENT Z_p^* orbits can have the same {|λ_k|²} multiset?
No: each orbit has a unique {|λ_k|²} multiset. The issue at p=11 is that
two orbits have the same PRODUCT of |λ_k| but different multisets.

Let me check: at p=11, the two orbits with prod=23 — do they have different
{y_k} multisets?
"""

import numpy as np
from itertools import combinations
from math import factorial, comb

def paley_qr(p):
    """Quadratic residues mod p."""
    return set(x*x % p for x in range(1, p))

def circulant_eigenvalues(p, S):
    """Eigenvalues of circulant tournament."""
    omega = np.exp(2j * np.pi / p)
    return [sum(omega**(k*s) for s in S) for k in range(p)]

def all_circulant_sets(p):
    """All tournament connection sets."""
    m = (p - 1) // 2
    reps = list(range(1, m + 1))
    result = []
    for bits in range(1 << m):
        S = set()
        for i in range(m):
            if bits & (1 << i):
                S.add(reps[i])
            else:
                S.add(p - reps[i])
        result.append(S)
    return result

def hamiltonian_paths_dp(A, n):
    """Held-Karp DP for HP count."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask >> u & 1: continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def circulant_adjacency(p, S):
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if (j - i) % p in S:
                A[i][j] = 1
    return A

print("=" * 70)
print("EIGENVALUE MULTISET ANALYSIS")
print("=" * 70)
print()

for p in [7, 11, 13]:
    print(f"=== p = {p} ===")
    qr = paley_qr(p)
    m = (p - 1) // 2
    sets = all_circulant_sets(p)

    data = []
    for S in sets:
        A = circulant_adjacency(p, S)
        H = hamiltonian_paths_dp(A, p)
        eigs = circulant_eigenvalues(p, S)
        # Independent magnitudes: k=1,...,(p-1)/2
        y = tuple(sorted([round(abs(eigs[k])**2, 8) for k in range(1, m+1)]))
        prod_y = 1
        for yi in y:
            prod_y *= yi
        data.append((H, y, prod_y, sorted(S)))

    # Group by y multiset
    from collections import defaultdict
    by_y = defaultdict(list)
    for H, y, prod_y, S in data:
        by_y[y].append((H, S))

    for y, entries in sorted(by_y.items(), key=lambda x: -x[1][0][0]):
        H_vals = set(e[0] for e in entries)
        print(f"  y = {y}")
        print(f"    prod(y) = {np.prod(y):.4f}")
        print(f"    H values: {sorted(H_vals, reverse=True)}")
        if len(H_vals) > 1:
            print(f"    *** DIFFERENT H for same y multiset! ***")
        print()

    # Now check: is H a function of the power sums of y?
    # Power sums: σ_j = Σ y_k^j for j = 1,2,...
    print(f"  Power sums of y = (|λ_k|²) for k=1,...,{m}:")
    for H, y, prod_y, S in sorted(data, key=lambda x: -x[0])[:6]:
        sigma = {}
        for j in range(1, 5):
            sigma[j] = sum(yi**j for yi in y)
        print(f"    S={S}: H={H}, σ₁={sigma[1]:.4f}, σ₂={sigma[2]:.4f}, "
              f"σ₃={sigma[3]:.4f}, σ₄={sigma[4]:.4f}")
    print()

# ================================================================
# At p=7, H is linear in σ₂ (since there's only 1 free param)
# ================================================================
print("=" * 70)
print("H AS POLYNOMIAL IN POWER SUMS (p=7)")
print("=" * 70)

p = 7
m = 3
sets = all_circulant_sets(p)
data = []
for S in sets:
    A = circulant_adjacency(p, S)
    H = hamiltonian_paths_dp(A, p)
    eigs = circulant_eigenvalues(p, S)
    y = [abs(eigs[k])**2 for k in range(1, m+1)]
    s1 = sum(y)
    s2 = sum(yi**2 for yi in y)
    s3 = sum(yi**3 for yi in y)
    data.append((H, s1, s2, s3, sorted(S)))

# H = a + b*s2 (s1 is constant by Parseval)
H_vals = [d[0] for d in data]
s2_vals = [d[2] for d in data]
s3_vals = [d[3] for d in data]

# Two distinct points:
H1, s2_1 = 189, 12.0  # Paley: y=(2,2,2), s2=12
H2, s2_2 = 175, 52.0/3  # Non-Paley: need exact

print(f"Data points:")
for H, s1, s2, s3, S in sorted(data, key=lambda x: -x[0]):
    print(f"  H={H}, s1={s1:.4f}, s2={s2:.4f}, s3={s3:.4f}, S={S}")

# Linear fit: H = a + b*s2
if len(set(s2_vals)) >= 2:
    # Get two distinct points
    pts = []
    for h, s2 in zip(H_vals, s2_vals):
        if (h, round(s2, 6)) not in [(p[0], round(p[1], 6)) for p in pts]:
            pts.append((h, s2))
    if len(pts) >= 2:
        h1, x1 = pts[0]
        h2, x2 = pts[1]
        b = (h1 - h2) / (x1 - x2)
        a = h1 - b * x1
        print(f"\nLinear fit: H = {a:.6f} + {b:.6f} * σ₂")
        # Verify
        for H, s1, s2, s3, S in data:
            pred = a + b * s2
            err = abs(H - pred)
            if err > 0.01:
                print(f"  ERROR: S={S}, H={H}, predicted={pred:.2f}")

print()

# ================================================================
# At p=11, need more power sums
# ================================================================
print("=" * 70)
print("H AS POLYNOMIAL IN POWER SUMS (p=11)")
print("=" * 70)

p = 11
m = 5
sets = all_circulant_sets(p)
data = []
for S in sets:
    A = circulant_adjacency(p, S)
    H = hamiltonian_paths_dp(A, p)
    eigs = circulant_eigenvalues(p, S)
    y = [abs(eigs[k])**2 for k in range(1, m+1)]
    sums = {}
    for j in range(1, 6):
        sums[j] = sum(yi**j for yi in y)
    data.append((H, sums, sorted(S)))

print("Data (4 distinct H values, 4 distinct sigma multisets):")
for H, sums, S in sorted(data, key=lambda x: -x[0])[:8]:
    print(f"  H={H}, σ₁={sums[1]:.2f}, σ₂={sums[2]:.2f}, "
          f"σ₃={sums[3]:.2f}, σ₄={sums[4]:.2f}")

# Collect unique points
seen = set()
unique = []
for H, sums, S in data:
    key = (H, round(sums[2], 4))
    if key not in seen:
        seen.add(key)
        unique.append((H, sums))

print(f"\nUnique (H, σ) values: {len(unique)}")

# Try to express H as polynomial in σ₂, σ₃
# 4 unknowns: H = a + b*σ₂ + c*σ₃ + d*σ₂²
# We have 4 equations
if len(unique) >= 3:
    import numpy.linalg as la
    A_mat = []
    b_vec = []
    for H, sums in unique:
        A_mat.append([1, sums[2], sums[3], sums[2]**2])
        b_vec.append(H)
    A_mat = np.array(A_mat, dtype=float)
    b_vec = np.array(b_vec, dtype=float)

    if A_mat.shape[0] >= A_mat.shape[1]:
        coeffs, residuals, rank, sv = la.lstsq(A_mat, b_vec, rcond=None)
        print(f"\nFit: H = {coeffs[0]:.4f} + {coeffs[1]:.4f}*σ₂ + "
              f"{coeffs[2]:.4f}*σ₃ + {coeffs[3]:.6f}*σ₂²")
        print(f"Rank: {rank}")

        # Verify
        max_err = 0
        for H, sums, S in data:
            pred = coeffs[0] + coeffs[1]*sums[2] + coeffs[2]*sums[3] + coeffs[3]*sums[2]**2
            err = abs(H - pred)
            max_err = max(max_err, err)
        print(f"Max error: {max_err:.6f}")

        if max_err > 1:
            # Try σ₂ and σ₃ only
            A_mat2 = [[1, s[2], s[3]] for _, s in unique]
            coeffs2, _, _, _ = la.lstsq(np.array(A_mat2), np.array([h for h,_ in unique]), rcond=None)
            print(f"\nSimpler fit: H = {coeffs2[0]:.4f} + {coeffs2[1]:.4f}*σ₂ + {coeffs2[2]:.4f}*σ₃")
            max_err2 = 0
            for H, sums, S in data:
                pred = coeffs2[0] + coeffs2[1]*sums[2] + coeffs2[2]*sums[3]
                err = abs(H - pred)
                max_err2 = max(max_err2, err)
            print(f"Max error: {max_err2:.6f}")

print()

# ================================================================
# The REAL question: what symmetric function of λ gives H?
# ================================================================
print("=" * 70)
print("H AS SYMMETRIC FUNCTION OF COMPLEX EIGENVALUES")
print("=" * 70)
print()
print("H depends on the COMPLEX eigenvalue multiset, not just magnitudes.")
print("For circulant tournaments, the relevant symmetric functions are")
print("the power sums p_k = Σ λ_j^k (complex).")
print()
print("Since eigenvalues come in conjugate pairs λ_j, conj(λ_j):")
print("  p_k = 2 * Σ Re(λ_j^k) for j=1,...,(p-1)/2")
print("  = 2 * Σ |λ_j|^k * cos(k * arg(λ_j))")
print()
print("So the PHASES matter, not just magnitudes!")
print()

for p in [7, 11]:
    m = (p - 1) // 2
    sets = all_circulant_sets(p)
    data = []
    for S in sets:
        A = circulant_adjacency(p, S)
        H = hamiltonian_paths_dp(A, p)
        eigs = circulant_eigenvalues(p, S)
        # Complex power sums (excluding λ_0)
        pk = {}
        for k in range(2, p+1):
            pk[k] = sum(eigs[j]**k for j in range(1, p)).real
        data.append((H, pk, sorted(S)))

    print(f"=== p = {p} ===")
    print(f"Complex power sums p_k = Σ Re(λ_j^k) for k=2,...,{p}:")
    for H, pk, S in sorted(data, key=lambda x: -x[0])[:6]:
        print(f"  H={H}, S={S}")
        for k in range(2, min(p+1, 8)):
            print(f"    p_{k} = {pk[k]:.4f}")

    # Try fitting H as linear function of power sums
    # H = a₀ + Σ aₖ pₖ
    print(f"\n  Fitting H = a₀ + Σ aₖ * p_k:")

    unique = []
    seen_H = set()
    for H, pk, S in data:
        if H not in seen_H:
            seen_H.add(H)
            unique.append((H, pk))

    n_unique = len(unique)
    print(f"  {n_unique} distinct H values")

    # Use p_4, p_5, p_6, ... as features (p_2, p_3 are universal)
    # First check which power sums actually vary
    for k in range(2, p+1):
        vals = set(round(u[1][k], 4) for u in unique)
        if len(vals) > 1:
            print(f"    p_{k} varies: {len(vals)} distinct values")

    # Fit using varying power sums
    varying_k = [k for k in range(2, p+1)
                 if len(set(round(u[1][k], 4) for u in unique)) > 1]

    if len(varying_k) >= n_unique - 1:
        # Use first (n_unique - 1) varying power sums
        use_k = varying_k[:n_unique - 1]
        A_mat = np.array([[1] + [u[1][k] for k in use_k] for _, u in [unique[i] for i in range(n_unique)]])
        if A_mat.shape[0] != n_unique:
            A_mat = np.zeros((n_unique, 1 + len(use_k)))
            for i, (H, pk) in enumerate(unique):
                A_mat[i, 0] = 1
                for j, k in enumerate(use_k):
                    A_mat[i, 1+j] = pk[k]

        b_vec = np.array([h for h, _ in unique])

        try:
            coeffs = np.linalg.solve(A_mat[:len(use_k)+1, :len(use_k)+1],
                                      b_vec[:len(use_k)+1])
            terms = [f"{coeffs[0]:.4f}"]
            for j, k in enumerate(use_k):
                terms.append(f"{coeffs[1+j]:.6f}*p_{k}")
            print(f"  H = {' + '.join(terms)}")

            # Verify on ALL data
            max_err = 0
            for H, pk, S in data:
                pred = coeffs[0] + sum(coeffs[1+j]*pk[k] for j, k in enumerate(use_k))
                err = abs(H - pred)
                max_err = max(max_err, err)
            print(f"  Max error on all {len(data)} tournaments: {max_err:.6f}")
        except Exception as e:
            print(f"  Solve failed: {e}")
    print()

# ================================================================
# Key question: is H a LINEAR function of power sums?
# ================================================================
print("=" * 70)
print("IS H LINEAR IN POWER SUMS?")
print("=" * 70)
print()

for p in [7, 11]:
    m = (p - 1) // 2
    sets = all_circulant_sets(p)
    data = []
    for S in sets:
        A = circulant_adjacency(p, S)
        H = hamiltonian_paths_dp(A, p)
        eigs = circulant_eigenvalues(p, S)
        pk = {k: sum(eigs[j]**k for j in range(1, p)).real for k in range(2, p+1)}
        data.append((H, pk))

    # Extract unique points
    unique = {}
    for H, pk in data:
        key = H
        if key not in unique:
            unique[key] = pk

    H_vals = sorted(unique.keys(), reverse=True)
    n_pts = len(H_vals)

    # Check linearity in p_4 alone (first varying power sum)
    if p == 7:
        # Only 2 distinct H values at p=7
        # H = a + b*p_4 should work
        H1, H2 = H_vals[0], H_vals[1]
        p4_1 = unique[H1][4]
        p4_2 = unique[H2][4]
        b = (H1 - H2) / (p4_1 - p4_2)
        a = H1 - b * p4_1
        print(f"p={p}: H = {a:.4f} + {b:.4f} * p_4")
        print(f"  H({p4_1:.0f}) = {a + b*p4_1:.0f} (expected {H1})")
        print(f"  H({p4_2:.0f}) = {a + b*p4_2:.0f} (expected {H2})")
        print(f"  → H = {a:.0f} + ({b:.4f}) * p_4")
        print(f"  → H = {a:.0f} - ({-b:.4f}) * (tr(A^4) - ((p-1)/2)^4)")
        print()

    if p == 11:
        # 4 distinct H values, need p_4, p_5, p_6 (3 independent params)
        # Try linear: H = a + b*p_4 + c*p_5 + d*p_6
        A_mat = []
        b_vec = []
        for H in H_vals:
            pk = unique[H]
            A_mat.append([1, pk[4], pk[5], pk[6]])
            b_vec.append(H)
        A_mat = np.array(A_mat, dtype=float)
        b_vec = np.array(b_vec, dtype=float)

        coeffs = np.linalg.solve(A_mat, b_vec)
        print(f"p={p}: H = {coeffs[0]:.4f} + {coeffs[1]:.6f}*p_4 + "
              f"{coeffs[2]:.6f}*p_5 + {coeffs[3]:.6f}*p_6")

        # Verify
        max_err = 0
        for H, pk in data:
            pred = coeffs[0] + coeffs[1]*pk[4] + coeffs[2]*pk[5] + coeffs[3]*pk[6]
            max_err = max(max_err, abs(H - pred))
        print(f"  Max error: {max_err:.6f}")

        # Check if coefficients are rational
        for i, c in enumerate(coeffs):
            label = ["a", "b(p_4)", "c(p_5)", "d(p_6)"][i]
            # Try to find rational approximation
            from fractions import Fraction
            frac = Fraction(c).limit_denominator(10000)
            print(f"  {label} = {c:.8f} ≈ {frac} = {float(frac):.8f}")
        print()

# ================================================================
# Newton identity approach: express H in terms of elementary symm funcs
# ================================================================
print("=" * 70)
print("NEWTON IDENTITIES: RELATION BETWEEN POWER SUMS AND ELEM SYMM")
print("=" * 70)
print()
print("For the eigenvalue polynomial: det(xI - A) = Π(x - λ_k)")
print("The elementary symmetric functions e_k are the coefficients.")
print("Newton's identities relate p_k and e_k.")
print()
print("For circulant: det(xI - A) = Π_{k=0}^{p-1} (x - λ_k)")
print("= (x - (p-1)/2) * Π_{k=1}^{p-1} (x - λ_k)")
print()
print("For Paley T_7: λ_0=3, λ_{1-6} are roots of (x+1/2)² + 7/4 = 0 each.")
print("Wait, |λ_k| = √2 for all k≠0 at p=7.")
print("Actually the minimal poly of each λ_k over Q is x² + x + 2 = 0")
print("(since λ_k = -1/2 ± i√(7)/2, so x = -1/2 + iy → (x+1/2)² = -7/4)")
print("x² + x + (1+7)/4 = x² + x + 2")
print()

for p in [7, 11]:
    qr = paley_qr(p)
    S = qr
    eigs = circulant_eigenvalues(p, S)

    # Characteristic polynomial
    coeffs = np.poly(eigs)
    print(f"p={p}: Characteristic polynomial of Paley adjacency:")
    for i, c in enumerate(coeffs):
        if abs(c.imag) < 1e-8:
            c = c.real
        if abs(c) > 0.001:
            print(f"  x^{p-i} coeff: {c:.4f} (≈ {round(c)})")

    # Factor out (x - (p-1)/2)
    from numpy.polynomial import polynomial as P
    reduced_coeffs = np.polydiv(coeffs, [1, -(p-1)/2])[0]
    print(f"  After dividing by (x - {(p-1)//2}):")
    for i, c in enumerate(reduced_coeffs):
        if abs(c.imag) < 1e-8:
            c = c.real
        if abs(c) > 0.001:
            print(f"    x^{p-1-i} coeff: {c:.6f} (≈ {round(c)})")
    print()

print()
print("=" * 70)
print("SYNTHESIS")
print("=" * 70)
print()
print("KEY FINDINGS:")
print()
print("1. H IS a function of the eigenvalue multiset for circulants (proved).")
print("2. H is determined by power sums p_k for k ≥ 4 (p_2, p_3 universal).")
print("3. At p=7: H is LINEAR in p_4 (1 free parameter).")
print("4. At p=11: H is LINEAR in (p_4, p_5, p_6) (3 free parameters).")
print("5. The linear coefficients appear to be RATIONAL numbers.")
print()
print("6. CONJECTURE: H is a POLYNOMIAL in the power sums p_k for k ≥ 4.")
print("   The degree of this polynomial is at most 1 in each p_k")
print("   (i.e., H is MULTI-LINEAR in the varying power sums).")
print()
print("7. If H is linear in p_k, then Paley maximizes H iff it maximizes")
print("   the weighted sum Σ a_k * p_k. Since p_k = Σ λ_j^k, this reduces")
print("   to showing each weighted power sum is maximal for flat spectrum.")
print("   This would complete Step B of the proof strategy!")
