#!/usr/bin/env python3
"""
pi_gauss_spectral_89c.py — Gauss sums, spectral structure, and π in Paley tournaments
opus-2026-03-14-S89c

The skew-adjacency matrix S of the Paley tournament P_p has eigenvalues
related to Gauss sums. We investigate:

1. The spectral decomposition S = Σ λ_j P_j
2. How eigenvalue phases relate to π
3. The determinant formula det(S) and its connection to class numbers
4. Trace formulas: Tr(S^k) counts closed directed walks of length k
"""

import cmath
import math
from fractions import Fraction

def legendre(a, p):
    """Legendre symbol (a/p)"""
    if a % p == 0:
        return 0
    return pow(a, (p-1)//2, p) if pow(a, (p-1)//2, p) != p-1 else -1

def gauss_sum(p):
    """Compute the quadratic Gauss sum g = Σ_{t=0}^{p-1} (t/p) ζ^t"""
    zeta = cmath.exp(2j * cmath.pi / p)
    g = sum(legendre(t, p) * zeta**t for t in range(p))
    return g

def paley_adjacency(p):
    """Skew-adjacency matrix: S[i][j] = (j-i)/p"""
    S = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j:
                S[i][j] = legendre(j - i, p)
    return S

def paley_eigenvalues(p):
    """Eigenvalues of skew-adjacency matrix of P_p via Gauss sum formula"""
    # The circulant matrix has eigenvalues λ_k = Σ_{j=1}^{p-1} (j/p) ω^{jk}
    # where ω = e^{2πi/p}
    omega = cmath.exp(2j * cmath.pi / p)
    eigs = []
    for k in range(p):
        lam = sum(legendre(j, p) * omega**(j*k) for j in range(1, p))
        eigs.append(lam)
    return eigs

print("=" * 70)
print("GAUSS SUMS AND π IN PALEY TOURNAMENT SPECTRA")
print("=" * 70)

for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    if not all(p % i != 0 for i in range(2, int(p**0.5)+1)):
        continue
    if p % 4 != 3:  # Paley tournaments need p ≡ 3 mod 4
        continue

    print(f"\n--- p = {p} ---")

    # Gauss sum
    g = gauss_sum(p)
    print(f"  Gauss sum g = {g.real:.6f} + {g.imag:.6f}i")
    print(f"  |g| = {abs(g):.6f}, √p = {math.sqrt(p):.6f}")
    print(f"  g² = {(g**2).real:.6f} + {(g**2).imag:.6f}i")
    # For p ≡ 3 mod 4: g² = -p (imaginary Gauss sum)
    print(f"  -p = {-p}")
    print(f"  arg(g)/π = {cmath.phase(g)/math.pi:.8f}")

    # Eigenvalues
    eigs = paley_eigenvalues(p)

    # λ_0 = 0 always (sum of Legendre symbols)
    print(f"\n  Eigenvalue λ_0 = {eigs[0].real:.6f} + {eigs[0].imag:.6f}i (should be 0)")

    # Non-zero eigenvalues
    nonzero = [(k, eigs[k]) for k in range(1, p)]
    # Group by absolute value
    print(f"\n  Non-zero eigenvalues (grouped by |λ|):")

    # For Paley, all |λ_k| = √p for k ≠ 0
    for k, lam in nonzero[:min(6, len(nonzero))]:
        phase = cmath.phase(lam) / math.pi
        print(f"    k={k}: λ = {lam.real:+.6f} {lam.imag:+.6f}i, |λ|={abs(lam):.4f}, arg/π={phase:+.8f}")

    # The phases arg(λ_k)/π for k=1,...,p-1
    phases = sorted([cmath.phase(eigs[k])/math.pi for k in range(1, p)])
    print(f"\n  Sorted phases (arg/π):")
    for ph in phases:
        print(f"    {ph:+.8f}")

    # The eigenvalues come in conjugate pairs λ_k, λ_{p-k}
    # For Paley with p ≡ 3 mod 4:
    # λ_k = (k/p) × g × ω^{...} but the exact formula involves characters
    # Actually λ_k = χ(k) × g where χ is the Legendre symbol? No.
    # The circulant eigenvalues are λ_k = Σ_{j=1}^{p-1} (j/p) ω^{jk}
    # = (k/p) × g when k ≠ 0 (this follows from the Gauss sum multiplication formula)

    # Verify: λ_k = (k/p) × g
    print(f"\n  Verification: λ_k = (k/p) × g?")
    all_match = True
    for k in range(1, min(p, 8)):
        predicted = legendre(k, p) * g
        actual = eigs[k]
        diff = abs(predicted - actual)
        match = diff < 1e-8
        if not match:
            all_match = False
        if k <= 5:
            print(f"    k={k}: (k/p)={legendre(k,p):+d}, predicted={predicted.real:+.4f}{predicted.imag:+.4f}i, actual={actual.real:+.4f}{actual.imag:+.4f}i, {'✓' if match else '✗'}")
    print(f"    All match: {'✓' if all_match else '✗'}")

    # So all eigenvalues are ±g, with multiplicities (p-1)/2 each!
    # This means the spectrum has exactly 2 distinct non-zero eigenvalues: +g and -g.
    # det(S) = 0^1 × g^{(p-1)/2} × (-g)^{(p-1)/2}
    # = g^{(p-1)/2} × (-1)^{(p-1)/2} × g^{(p-1)/2}
    # = (-1)^{(p-1)/2} × g^{p-1}

    count_plus = sum(1 for k in range(1,p) if legendre(k,p) == 1)
    count_minus = sum(1 for k in range(1,p) if legendre(k,p) == -1)
    print(f"\n  Multiplicities: +g appears {count_plus} times, -g appears {count_minus} times")

    # det(S) = g^{(p-1)/2} × (-g)^{(p-1)/2}  [ignoring 0 eigenvalue]
    # Wait, λ_0 = 0, so det(S) = 0 always! That's right.
    # But the product of nonzero eigenvalues:
    prod_nonzero = 1
    for k in range(1, p):
        prod_nonzero *= eigs[k]
    print(f"\n  Product of nonzero eigenvalues: {prod_nonzero.real:.4f} + {prod_nonzero.imag:.4f}i")
    # Should be g^{(p-1)/2} × (-g)^{(p-1)/2} = (-1)^{(p-1)/2} × g^{p-1}
    # g^{p-1} = (g²)^{(p-1)/2} = (-p)^{(p-1)/2}
    # For p ≡ 3 mod 4: (p-1)/2 is odd, so (-1)^{(p-1)/2} = -1
    # and (-p)^{(p-1)/2} = (-1)^{(p-1)/2} × p^{(p-1)/2} = -p^{(p-1)/2}
    # So product = (-1) × (-p^{(p-1)/2}) = p^{(p-1)/2}
    predicted_prod = p**((p-1)//2)
    print(f"  Predicted: p^((p-1)/2) = {predicted_prod}")
    print(f"  Match: {'✓' if abs(prod_nonzero.real - predicted_prod) < 1 else '✗'}")

    # Connection to π: each eigenvalue has argument ±arg(g) where
    # arg(g) = π/2 + ε for small ε.
    # More precisely, for p ≡ 3 mod 4:
    # g = i√p (up to sign), so arg(g) = ±π/2.
    # The actual value: g = i^{...} × √p with a specific phase.
    arg_g = cmath.phase(g)
    print(f"\n  arg(g) = {arg_g:.8f} = {arg_g/math.pi:.8f} × π")
    print(f"  arg(g)/π - 1/2 = {arg_g/math.pi - 0.5:.8f}")

print("\n" + "=" * 70)
print("THE π CONNECTION: GAUSS SUM PHASE")
print("=" * 70)

print("""
  For p ≡ 3 (mod 4), the quadratic Gauss sum satisfies:
    g² = (-1)^{(p-1)/2} × p = -p

  So g = ±i√p, meaning arg(g) = ±π/2.

  ALL non-zero eigenvalues of P_p have |λ| = √p and arg = ±π/2.

  The skew-adjacency matrix is pure imaginary (up to the zero eigenvalue).
  This is because P_p is a "self-complementary" tournament:
    (a/p) = -(−a/p) when p ≡ 3 mod 4.

  The factor π appears because:
  1. Gauss sums g = i√p have phase π/2 (EXACT, not approximate)
  2. The eigenvalues sit at ±i√p on the imaginary axis
  3. The spectral radius is √p (Ramanujan-like bound)
  4. det(S_{00}) = p^{(p-3)/2} (Pfaffian involves p^{(p-3)/4})
""")

# Trace formula: connection to cycle counts
print("=" * 70)
print("TRACE FORMULA: CLOSED WALKS AND CYCLE COUNTS")
print("=" * 70)

for p in [3, 7, 11]:
    if p % 4 != 3:
        continue

    S = paley_adjacency(p)
    print(f"\n  p = {p}:")

    # Tr(S^k) = sum of k-th powers of eigenvalues
    # = 0^k + (p-1)/2 × (i√p)^k + (p-1)/2 × (-i√p)^k
    # When k is odd: (i√p)^k + (-i√p)^k = 0 (since i^k + (-i)^k = 0 for odd k? No.)
    # i^k + (-i)^k = i^k(1 + (-1)^k) = 0 for odd k. ✓
    # When k is even: (i√p)^k + (-i√p)^k = 2×(i√p)^k for k≡0 mod 4, or -2×... hmm

    # Let me just compute directly
    for k in range(1, 8):
        # Matrix power trace
        Sk = [[0]*p for _ in range(p)]
        # Initialize to identity
        for i in range(p):
            Sk[i][i] = 1
        # Multiply k times
        for _ in range(k):
            Sk_new = [[0]*p for _ in range(p)]
            for i in range(p):
                for j in range(p):
                    for l in range(p):
                        Sk_new[i][j] += Sk[i][l] * S[l][j]
            Sk = Sk_new
        tr = sum(Sk[i][i] for i in range(p))

        # Spectral prediction
        eigs = paley_eigenvalues(p)
        tr_spectral = sum(e**k for e in eigs)

        print(f"    Tr(S^{k}) = {tr}, spectral = {tr_spectral.real:.1f}+{tr_spectral.imag:.1f}i")

    # Tr(S^3) counts directed 3-cycles: each 3-cycle ABC is counted 3 times
    # (starting at A, B, or C), times 2 orientations.
    # Actually Tr(S^3) = Σ_{i,j,k} S[i][j]S[j][k]S[k][i]
    # Each directed 3-cycle i→j→k→i contributes S[i][j]S[j][k]S[k][i] = 1×1×1 = 1
    # Each COUNTER-directed 3-cycle contributes (-1)³ = -1
    # But for TOURNAMENTS, every triple forms either a 3-cycle or a transitive triple.
    # Directed 3-cycles: if (i→j, j→k, k→i), trace term = S[i][j]·S[j][k]·S[k][i]
    # But S[i][j] = (j-i)/p for Paley, which is ±1.
    # Hmm this gets complex. The key point is:
    # Tr(S^k) = 0 for odd k (because eigenvalues come in ±pairs for Paley)
    # This is a SYMMETRY property!

print("\n" + "=" * 70)
print("SPECTRAL DETERMINANT OF MINOR S_00")
print("=" * 70)

for p in [3, 7, 11, 19, 23]:
    if p % 4 != 3:
        continue

    S = paley_adjacency(p)

    # Minor S_00: delete row 0, col 0
    S00 = [[S[i][j] for j in range(1, p)] for i in range(1, p)]

    # Compute det via numpy-free method (expansion for small p)
    # For efficiency, use the formula we proved: det(S_00) = p^{(p-3)/2}
    predicted = p**((p-3)//2)

    # Direct computation via row reduction for small p
    import copy
    M = [row[:] for row in S00]
    n = len(M)
    det_val = 1
    for col in range(n):
        # Find pivot
        pivot = None
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            det_val = 0
            break
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
            det_val *= -1
        det_val *= M[col][col]
        for row in range(col+1, n):
            if M[row][col] != 0:
                factor = Fraction(M[row][col], M[col][col])
                for j in range(col, n):
                    M[row][j] = Fraction(M[row][j]) - factor * Fraction(M[col][j])

    print(f"  p={p}: det(S_00) = {det_val}, predicted p^((p-3)/2) = {predicted}, match: {'✓' if int(det_val) == predicted else '✗'}")
    if p <= 11:
        # Pfaffian: since S_00 is skew-symmetric of even size (p-1 is even),
        # det = Pf²
        pf_sq = int(det_val)
        pf = int(round(abs(pf_sq)**0.5))
        if pf * pf == abs(pf_sq):
            print(f"    Pf(S_00)² = {pf_sq}, |Pf| = {pf} = p^((p-3)/4) = {p**((p-3)//4)}")

print("\n" + "=" * 70)
print("SUMMARY: WHERE π LIVES IN THE SPECTRUM")
print("=" * 70)

print("""
  The Paley tournament P_p (p ≡ 3 mod 4) has spectrum:

  • One zero eigenvalue (from all-ones vector)
  • (p-1)/2 eigenvalues at +i√p  (QR characters)
  • (p-1)/2 eigenvalues at -i√p  (NQR characters)

  The eigenvalues lie on the IMAGINARY axis at ±i√p.

  arg(eigenvalue) = ±π/2    (EXACT, not approximate)

  This π/2 phase comes from the Gauss sum identity:
    g² = -p  ⟹  g = ±i√p  ⟹  arg(g) = ±π/2

  This is why Gauss sums are called "imaginary" for p ≡ 3 mod 4.
  The factor π enters through the complex exponential:
    ζ_p = e^{2πi/p}
  which is the ONLY place π appears in the definition.

  DEEPER: The Gauss sum g = Σ (t/p) e^{2πit/p} is a discrete
  Fourier transform of the Legendre symbol. The resulting phase
  π/2 is a consequence of quadratic reciprocity + Euler's formula.
""")

print("Done!")
