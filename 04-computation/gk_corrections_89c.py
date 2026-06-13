#!/usr/bin/env python3
"""
gk_corrections_89c.py — Analyze corrections to the naive formula
opus-2026-03-15-S89c

Naive: g_k(m) = Σ_{r=1}^k C(k-1,r-1)·C(m,r)·2^{r-1}
True g_k is degree 3 for k≥3.
Correction E_k(m) = true - naive.

For k=4: E_4(m) = -(m-1)(m-2)(m-3)(m-4)/3 = -8C(m-1,4)
What about k=5,6,7?
"""

from math import comb, factorial
from fractions import Fraction
from sympy import symbols, interpolate, Rational, factor, expand, Poly

x = symbols('x')

def g_naive(k, m):
    return sum(comb(k-1, r-1) * comb(m, r) * 2**(r-1) for r in range(1, k+1))

# True g_k polynomials: 3·g_k(m) = a·m³+b·m²+c·m+d
gk_coeffs = {
    3: (2, 0, 1, 0),
    4: (10, -33, 50, -24),
    5: (388, -2040, 3431, -1776),
    6: (69660, -380445, 653748, -342960),
    7: (19826270, -109486152, 189674605, -100014720),
    8: (7309726742, -40641958545, 70757788486, -37425556680),
    9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
}

def g_true(k, m):
    if k == 1: return m
    if k == 2: return m*m
    if k in gk_coeffs:
        a, b, c, d = gk_coeffs[k]
        return (a*m**3 + b*m**2 + c*m + d) // 3
    return None

# Compute corrections in binomial basis: C(m,r) for r=0,...,k
print("="*70)
print("BINOMIAL BASIS OF CORRECTIONS E_k(m) = true - naive")
print("="*70)

for k in range(3, 10):
    # Compute E_k(m) for m=0,...,k+3
    E_vals = []
    for m in range(0, max(k+4, 12)):
        t = g_true(k, m)
        n = g_naive(k, m)
        if t is not None:
            E_vals.append(t - n)
        else:
            break

    print(f"\nk={k}: E_{k}(m) = {E_vals[:10]}")

    # Forward differences to get binomial coefficients
    # E_k(m) = Σ Δ^r E_k(0) · C(m,r)
    curr = list(E_vals)
    deltas = []
    for r in range(len(curr)):
        deltas.append(curr[0])
        curr = [curr[i+1] - curr[i] for i in range(len(curr)-1)]
        if not curr:
            break

    # Trim trailing zeros
    while deltas and deltas[-1] == 0:
        deltas.pop()

    print(f"  Δ^r E_{k}(0): {deltas}")

    # Express as polynomial
    if len(E_vals) >= 2:
        pts = [(Rational(m), Rational(E_vals[m])) for m in range(min(len(E_vals), k+2))]
        poly = interpolate(pts, x)
        poly = expand(poly)
        print(f"  E_{k}(m) = {poly}")
        print(f"  E_{k}(m) = {factor(poly)}")

# Analyze the naive formula itself
print("\n" + "="*70)
print("NAIVE FORMULA: g_k^naive(m) = Σ C(k-1,r-1)·C(m,r)·2^{r-1}")
print("="*70)

# The naive formula = (1+2)^{k-1} evaluated in the "binomial-C(m,r)" basis
# By the binomial theorem: Σ_{r=1}^k C(k-1,r-1)·C(m,r)·2^{r-1}
# = Σ_{r=0}^{k-1} C(k-1,r)·C(m,r+1)·2^r
# Using Vandermonde: Σ C(k-1,r)·C(m,r+1)·2^r = ... hmm, not standard Vandermonde.
# But Σ_r C(k-1,r)·C(m,r+1)·2^r = C(m,1)·Σ_r C(k-1,r)·C(m-1,r)·2^r / ... no.
# Actually this is the coefficient extraction from (1+x)^{k-1}·(1+2x)^... no.

# The DEGREE of the naive formula: since C(m,r) is degree r in m,
# and r goes up to k, the naive formula has degree k.
# The true g_k has degree 3 for k≥3.
# So the correction E_k must cancel degrees 4,...,k.

# Let's look at the LEADING term of E_k in the monomial basis
for k in range(3, 10):
    # Leading coefficient of naive in m^k
    # C(m,k) has leading coeff 1/k!, and its coefficient in naive is C(k-1,k-1)·2^{k-1} = 2^{k-1}
    # So leading term of naive = 2^{k-1}/k! × m^k
    lc_naive = Fraction(2**(k-1), factorial(k))

    # Leading coeff of true: a_k/3 × m³ for k≥3
    if k in gk_coeffs:
        lc_true = Fraction(gk_coeffs[k][0], 3)
    elif k <= 2:
        lc_true = Fraction(1)
    else:
        lc_true = "?"

    # Leading correction = lc_true - lc_naive (different degrees!)
    print(f"  k={k}: naive leading = {lc_naive}/m^{k}, true leading = {lc_true}/m³")
    if k >= 4:
        print(f"    correction must kill the m^{k} term: -2^{{{k-1}}}/{k}! = {-lc_naive}")

# EGF connection: the naive formula has EGF interpretation
print("\n" + "="*70)
print("EGF OF NAIVE FORMULA")
print("="*70)
print("g_k^naive(m) = Σ C(k-1,r-1)·C(m,r)·2^{r-1}")
print("= coefficient of t^{k-1} in (1+t)^{k-1} summed with C(m,r)·2^{r-1}")
print()
print("Using the identity Σ_r C(m,r)·z^r = (1+z)^m:")
print("  Σ_{r≥1} C(m,r)·2^{r-1} = ((1+2)^m - 1)/2 = (3^m - 1)/2")
print("  But we also sum over compositions, so it's not just 3^m.")
print()
print("Actually: g_k^naive(m) = Σ_r C(k-1,r-1)·C(m,r)·2^{r-1}")
print("= (1/2) Σ_r C(k-1,r-1)·C(m,r)·2^r")
print("Let u^r replace 2^r, we get:")
print("Σ_r C(k-1,r-1)·C(m,r)·u^r = u·Σ_r C(k-1,r)·C(m,r+1)·u^r")
print()

# By Vandermonde-Chu generalization:
# Σ_r C(a,r)·C(b,r+c)·z^r = sum involving hypergeometric
# For c=1: Σ C(k-1,r)·C(m,r+1)·z^r
# This can be evaluated using the identity:
# Σ_r C(n,r)·C(m,r+s)·z^r = C(m,s)·₂F₁(-n, s+1; s-m+1; -z)

# Let's just verify numerically
print("Numerical check of naive = binomial sum:")
for k in [4, 5]:
    for m in [3, 5, 7]:
        # Direct computation
        val = g_naive(k, m)
        # Try (1+2)^something
        print(f"  k={k}, m={m}: naive={val}, 3^m={(3**m-1)//2}")

# Interestingly, for k=1: g_1(m) = m. And C(0,0)·C(m,1)·1 = m. ✓
# For k→∞: the naive formula Σ C(k-1,r-1)·C(m,r)·2^{r-1}
# ≈ Σ_r C(m,r)·2^{r-1}·(k-1)^{r-1}/(r-1)!
# ≈ (1/2) Σ_r C(m,r)·(2(k-1))^{r-1}·/(r-1)!
# This grows exponentially in k.

# KEY INSIGHT: look at the correction polynomial factored form
print("\n" + "="*70)
print("CORRECTION POLYNOMIALS (FACTORED)")
print("="*70)

for k in range(4, 10):
    E_vals = [g_true(k, m) - g_naive(k, m) for m in range(k+4)]
    pts = [(Rational(m), Rational(E_vals[m])) for m in range(min(len(E_vals), k+1))]
    poly = interpolate(pts, x)
    poly = expand(poly)
    f = factor(poly)
    print(f"k={k}: E_{k}(m) = {f}")

    # Check if it has factor (m-1)(m-2)
    val_at_1 = poly.subs(x, 1)
    val_at_2 = poly.subs(x, 2)
    print(f"  E_{k}(1) = {val_at_1}, E_{k}(2) = {val_at_2}")

# Also: the correction should be expressible in terms of falling factorials
# or as a combination of C(m,r) for r ≥ 4
print("\n" + "="*70)
print("CORRECTION IN C(m,r) BASIS")
print("="*70)

for k in range(4, 10):
    E_vals = [g_true(k, m) - g_naive(k, m) for m in range(k+4)]
    # Forward differences
    curr = list(E_vals)
    deltas = []
    for r in range(len(curr)):
        deltas.append(curr[0])
        curr = [curr[i+1] - curr[i] for i in range(len(curr)-1)]
        if not curr:
            break

    # Remove trailing zeros
    while deltas and deltas[-1] == 0:
        deltas.pop()

    terms = []
    for r, d in enumerate(deltas):
        if d != 0:
            terms.append(f"{d}·C(m,{r})")
    print(f"k={k}: E_{k}(m) = {' + '.join(terms)}")

    # Compare with naive C(m,r) coefficients
    naive_deltas = [0]  # r=0
    for r in range(1, k+1):
        naive_deltas.append(comb(k-1, r-1) * 2**(r-1))

    # Total = naive + correction
    total_deltas = [0] * max(len(deltas), len(naive_deltas))
    for i, d in enumerate(naive_deltas):
        total_deltas[i] += d
    for i, d in enumerate(deltas):
        total_deltas[i] += d

    # Remove trailing zeros
    while total_deltas and total_deltas[-1] == 0:
        total_deltas.pop()

    print(f"  true C(m,r) coeffs: {total_deltas}")
    print(f"  naive C(m,r) coeffs: {naive_deltas}")
    print(f"  correction C(m,r) coeffs: {deltas}")

print("\nDone!")
