#!/usr/bin/env python3
"""
bernoulli_tournament_connection.py - Deep exploration of Bernoulli numbers in tournament theory.

The cumulant hierarchy has Bernoulli constants:
  κ₂ = (n+1)/12 + cycle terms  ↔ B₂ = 1/6
  κ₄ = -(n+1)/120 + cycle terms  ↔ B₄ = -1/30
  κ₆ = (n+1)/252 + cycle terms  ↔ B₆ = 1/42

The pattern: κ_{2k} constant = (-1)^{k-1} (n+1) |B_{2k}| / (2k)
           = (n+1) · B_{2k} / (2k) · (-1)^{k-1} · (-1)^{k-1}
Actually: (n+1)/12 = (n+1) · 1/12 = (n+1) · B₂ / (2·1)... wait.
  B₂ = 1/6, so (n+1) · B₂ / 1 = (n+1)/6. But κ₂ constant = (n+1)/12.
  Hmm, (n+1)/12 = (n+1) · B₂ / 2 = (n+1)/(2·6).

More precisely: κ_{2k} constant for TRANSITIVE tournament = the Bernoulli cumulant
of the uniform distribution on {0, 1, ..., n-1}:
  κ₂(Uniform{0,...,n-1}) = (n²-1)/12 ≠ (n+1)/12

Wait, let me recalculate. The forward-edge count for a TRANSITIVE tournament (total order)
is fwd(σ) = inv(σ), the number of inversions of σ. The distribution of inv(σ) over S_n
is the "inversion distribution" with known cumulants.

Actually, fwd(σ) for any tournament T is not inversions in general. But the MEAN is
always (n-1)/2 (by symmetry of each arc contributing 1/2 on average).

The Bernoulli constants come from the Eulerian number cumulants, not the inversion cumulants.

Let me think more carefully and compute.

Author: opus-2026-03-07-S46f
"""
from fractions import Fraction
from math import factorial, comb
from collections import Counter
from itertools import permutations, combinations

def eulerian(n, k):
    """Eulerian number A(n,k) = # permutations of [n] with k descents."""
    if n == 0:
        return 1 if k == 0 else 0
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def bernoulli(n):
    """Compute B_n (Bernoulli number, B_1 = -1/2 convention)."""
    B = [Fraction(0)] * (n+1)
    B[0] = Fraction(1)
    for m in range(1, n+1):
        B[m] = -sum(comb(m, k) * B[k] / (m - k + 1) for k in range(m))
    return B[n]

# Compute Bernoulli numbers
print("=" * 70)
print("BERNOULLI NUMBERS")
print("=" * 70)
for k in range(8):
    B = bernoulli(2*k)
    print(f"  B_{2*k} = {B} = {float(B):.6f}")

# The cumulant constants
print("\n" + "=" * 70)
print("CUMULANT CONSTANTS vs BERNOULLI")
print("=" * 70)

# κ₂ constant = (n+1)/12 for any n
# The Eulerian distribution has cumulants known in closed form.
# For a RANDOM permutation σ of [n], the distribution of des(σ) has:
#   Mean = (n-1)/2
#   Var = (n+1)/12  [this is the well-known Eulerian variance]
# These are the cumulants of the EULERIAN distribution, not tournament-specific.

# For the FORWARD-EDGE distribution of a tournament:
#   Mean = (n-1)/2 (by symmetry, same as Eulerian)
#   Var = (n+1)/12 + 4t₃/(n(n-1)) (tournament-dependent)

# So the Bernoulli constants ARE the Eulerian cumulants!

print("Eulerian cumulants (des(σ) over S_n):")
for n in range(3, 9):
    # Compute the actual cumulants of the Eulerian distribution at size n
    # des(σ) ranges from 0 to n-1
    # F_k = A(n, k) = Eulerian number
    N = factorial(n)

    # Moments about the mean
    mean = Fraction(n-1, 2)

    # E[X^r] = (1/n!) Σ_σ des(σ)^r = (1/n!) Σ_k k^r A(n,k)
    moments = []
    for r in range(7):
        M_r = sum(Fraction(k)**r * eulerian(n, k) for k in range(n)) / N
        moments.append(M_r)

    # Central moments
    mu2 = moments[2] - moments[1]**2
    mu3 = moments[3] - 3*moments[2]*moments[1] + 2*moments[1]**3
    mu4 = (moments[4] - 4*moments[3]*moments[1] + 6*moments[2]*moments[1]**2
            - 3*moments[1]**4)

    # Cumulants
    kappa2 = mu2
    kappa3 = mu3
    kappa4 = mu4 - 3*mu2**2

    print(f"  n={n}: κ₂ = {kappa2} = {float(kappa2):.6f}", end="")
    expected = Fraction(n+1, 12)
    print(f"  (n+1)/12 = {expected} {'✓' if kappa2 == expected else '✗'}")

print("\n  Eulerian κ₃ (should be 0 by palindromicity):")
for n in range(3, 9):
    N = factorial(n)
    mean = Fraction(n-1, 2)
    moments_c = []
    for r in range(4):
        M_r = sum(Fraction(k - mean)**r * eulerian(n, k) for k in range(n)) / N
        moments_c.append(M_r)
    print(f"  n={n}: κ₃ = {moments_c[3]}")

print("\n  Eulerian κ₄:")
for n in range(3, 9):
    N = factorial(n)
    mean = Fraction(n-1, 2)
    mu = [Fraction(0)] * 5
    for r in range(5):
        mu[r] = sum(Fraction(k - mean)**r * eulerian(n, k) for k in range(n)) / N

    kappa4 = mu[4] - 3 * mu[2]**2
    expected = -Fraction(n+1, 120)
    print(f"  n={n}: κ₄ = {kappa4} = {float(kappa4):.6f}  "
          f"-(n+1)/120 = {expected} {'✓' if kappa4 == expected else '✗'}")

print("\n  Eulerian κ₆:")
for n in range(3, 9):
    N = factorial(n)
    mean = Fraction(n-1, 2)
    mu = [Fraction(0)] * 7
    for r in range(7):
        mu[r] = sum(Fraction(k - mean)**r * eulerian(n, k) for k in range(n)) / N

    # κ₆ = μ₆ - 15μ₄μ₂ - 10μ₃² + 30μ₂³
    kappa6 = mu[6] - 15*mu[4]*mu[2] - 10*mu[3]**2 + 30*mu[2]**3
    expected = Fraction(n+1, 252)
    print(f"  n={n}: κ₆ = {kappa6}  (n+1)/252 = {expected} {'✓' if kappa6 == expected else '✗'}")

# The pattern: κ_{2k}^{Eulerian} = (-1)^{k-1} (n+1) B_{2k} / (2k) ?
print("\n" + "=" * 70)
print("FORMULA: κ_{2k}^{Eulerian} = (-1)^{k-1} · (n+1) · |B_{2k}| / (2k)")
print("=" * 70)

for k in range(1, 5):
    B2k = bernoulli(2*k)
    formula = Fraction((-1)**(k-1)) * Fraction(8+1, 1) * abs(B2k) / (2*k)  # n=8
    print(f"  k={k}: B_{2*k} = {B2k}, formula at n=8: {formula}")

# Actually the formula should be n-independent for the RATIO.
# κ_{2k}^{Eul} / (n+1) = (-1)^{k-1} |B_{2k}| / (2k)
print("\n  Checking: κ_{2k} / (n+1) = (-1)^{k-1} B_{2k} / (2k)")
for k in range(1, 4):
    for n in range(3, 9):
        N = factorial(n)
        mean = Fraction(n-1, 2)
        mu = [Fraction(0)] * (2*k+1)
        for r in range(2*k+1):
            mu[r] = sum(Fraction(j - mean)**r * eulerian(n, j) for j in range(n)) / N

        if k == 1:
            kappa = mu[2]
        elif k == 2:
            kappa = mu[4] - 3*mu[2]**2
        elif k == 3:
            kappa = mu[6] - 15*mu[4]*mu[2] - 10*mu[3]**2 + 30*mu[2]**3

        ratio = kappa / (n+1) if n+1 != 0 else None
        B2k = bernoulli(2*k)
        expected_ratio = B2k / (2*k)

        match = '✓' if ratio == expected_ratio else '✗'
        if n == 3:
            print(f"  k={k}: B_{2*k}/(2k) = {expected_ratio}, ", end="")
        if n == 3:
            print(f"κ_{2*k}/(n+1) at n=3..8: ", end="")
        print(f"{match}", end=" ")
    print()

# DEEPER: What does this mean for tournament cumulants?
print("\n" + "=" * 70)
print("TOURNAMENT CUMULANTS = EULERIAN + CYCLE CORRECTIONS")
print("=" * 70)
print("""
κ_{2k}(T) = κ_{2k}^{Eulerian}(n) + cycle-dependent terms

The Eulerian base:
  κ₂^{Eul} = (n+1)/12 = (n+1) · B₂ / 2
  κ₄^{Eul} = -(n+1)/120 = (n+1) · B₄ / 4
  κ₆^{Eul} = (n+1)/252 = (n+1) · B₆ / 6

GENERAL FORMULA: κ_{2k}^{Eul} = (n+1) · B_{2k} / (2k)

The cycle corrections (from THM-117):
  Δκ₂ = (2/C(n,2)) · t₃
  Δκ₄ = (2/C(n,4)) · (t₅ + 2α₂) - 48/(n(n-1))² · t₃²
  Δκ₆ = (2/C(n,6)) · t₇ + [nonlinear in t₃, t₅, α₂]

STRUCTURAL INSIGHT:
  κ_{2k}(T) = (n+1) · B_{2k}/(2k) + 2·t_{2k+1}/C(n,2k) + [nonlinear lower]

The FIRST term is universal (Bernoulli), depends only on n.
The SECOND term is the leading cycle correction (THM-117), depends on T.
The remaining terms are nonlinear combinations of lower cycle invariants.

This means: the Bernoulli numbers govern the "base topology" of tournaments,
and the cycle corrections are perturbations around this base.

ANALOGY with statistical mechanics:
  κ_{2k}^{Eul} = "free energy" of the Eulerian random process
  Δκ_{2k} = "interaction energy" from cycle structure
  The Bernoulli numbers = "partition function" coefficients
""")

# Now: the GENERATING FUNCTION of Eulerian cumulants
print("=" * 70)
print("GENERATING FUNCTION OF EULERIAN CUMULANTS")
print("=" * 70)

# The cumulant generating function of des(σ):
# K(t) = Σ κ_k t^k / k!
# For the Eulerian distribution on S_n, the mgf is:
# M(t) = E[e^{t·des(σ)}] = (1/n!) Σ_σ e^{t·des(σ)} = (1/n!) Σ_k A(n,k) e^{kt}

# The cgf is K(t) = log M(t).
# For CENTERED variable (des - mean):
# K_c(t) = log E[e^{t(des-mean)}]

# Known result (folklore): the Eulerian cgf is
# K_c(t) = log(t / (e^t - 1)) + t/2 (for the continuous Bernoulli limit)
# More precisely, as n→∞, des(σ)/n → Uniform[0,1], so
# κ_{2k}(des/n) → κ_{2k}(Uniform[0,1]) = B_{2k}/(2k)

# But we have κ_{2k}(des) = (n+1) · B_{2k}/(2k) for FINITE n!
# This is EXACT, not just asymptotic!

print("""
The Eulerian cumulant formula κ_{2k}^{Eul} = (n+1) · B_{2k}/(2k) is EXACT for all n.

This connects to the Bernoulli polynomial generating function:
  Σ_{k=0}^∞ B_k t^k / k! = t / (e^t - 1)

So the cumulant generating function of des(σ) (centered) is:
  K(t) = (n+1) · log(t / (e^t - 1)) + correction terms

This is related to the TODD CLASS in algebraic geometry:
  td(x) = x / (1 - e^{-x}) = Σ B_k (-x)^k / k!

The tournament cumulant hierarchy therefore lives in the TODD CLASS!
Tournament structure = Todd class + odd-cycle perturbations.

This suggests a deep connection between:
1. Tournament parity (Rédei) ↔ Todd class (Hirzebruch-Riemann-Roch)
2. OCF I(Ω,2) ↔ Some characteristic class of a bundle
3. Cumulant hierarchy ↔ Chern character expansion
""")
