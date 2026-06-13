#!/usr/bin/env python3
"""
q1_limit_exact.py — opus-2026-03-13-S67h

Determine the EXACT limit of Q_1/ΣQ as p → ∞ for Interval tournaments.

Q_k = sin²(mπk/p) / sin²(πk/p), m = (p-1)/2, ΣQ = m(m+1)/2.

We compute Q_1/ΣQ analytically:
  Q_1 = sin²(mπ/p) / sin²(π/p)
  ΣQ = m(m+1)/2

As p → ∞ with m = (p-1)/2:
  sin(mπ/p) = sin((p-1)π/(2p)) = sin(π/2 - π/(2p)) = cos(π/(2p))
  sin(π/p) = π/p - (π/p)³/6 + ...

So Q_1 = cos²(π/(2p)) / sin²(π/p)
       = (1 - π²/(8p²) + ...) / (π²/p² - 2π⁴/(3p⁴) + ...)
       = p²/π² · (1 - π²/(8p²)) · (1 + 2π²/(3p²) + ...) + O(1)
       = p²/π² - 1/8 + 2/(3) + O(1/p²)
Wait, this needs more care.

Let me just compute to very high precision and fit.
"""

import numpy as np
from math import pi, sin, cos, log

def is_prime(n):
    if n < 2: return False
    if n < 4: return n > 1
    if n % 2 == 0 or n % 3 == 0: return False
    d = 5
    while d*d <= n:
        if n % d == 0 or (n+2) % d == 0: return False
        d += 6
    return True

# EXACT ANALYSIS
# Q_1 = sin²(mπ/p) / sin²(π/p) where m = (p-1)/2
#
# Let x = π/p. Then:
#   sin(mπ/p) = sin((p-1)π/(2p)) = sin(π/2 - x/2) = cos(x/2)
#   sin(π/p) = sin(x)
#
# So Q_1 = cos²(x/2) / sin²(x) where x = π/p → 0
#
# Using sin(x) = 2sin(x/2)cos(x/2):
#   Q_1 = cos²(x/2) / (4sin²(x/2)cos²(x/2)) = 1/(4sin²(x/2))
#
# WAIT! This simplifies beautifully!
# Q_1 = 1/(4sin²(π/(2p)))
#
# And ΣQ = m(m+1)/2 = (p-1)(p+1)/8 = (p²-1)/8
#
# So Q_1/ΣQ = 8 / (4sin²(π/(2p)) · (p²-1))
#            = 2 / (sin²(π/(2p)) · (p²-1))
#
# As p → ∞: sin(π/(2p)) ≈ π/(2p), so sin² ≈ π²/(4p²)
# Q_1/ΣQ ≈ 2 / (π²/(4p²) · p²) = 2 · 4 / π² = 8/π²   ✓
#
# But we need the EXACT subleading correction!
# sin²(π/(2p)) = (π/(2p))² - (π/(2p))⁴/3 + ... = π²/(4p²)(1 - π²/(12p²) + ...)
#
# Q_1/ΣQ = 2 / ((π²/(4p²))(1 - π²/(12p²) + ...) · (p²-1))
#         = 8/(π²) · 1/((1 - π²/(12p²))(1 - 1/p²))
#         = 8/π² · 1/(1 - 1/p² - π²/(12p²) + O(p⁻⁴))
#         = 8/π² · (1 + (1 + π²/12)/p² + O(p⁻⁴))
#
# So Q_1/ΣQ = 8/π² + 8/π² · (1 + π²/12)/p² + O(p⁻⁴)
#            = 8/π² + 8(12 + π²)/(12π²p²) + O(p⁻⁴)

print("EXACT DERIVATION OF Q_1/ΣQ LIMIT")
print("=" * 60)
print()
print("THEOREM: For Interval tournament C_p^{1..m} with m=(p-1)/2:")
print()
print("  Q_1 = 1/(4·sin²(π/(2p)))")
print()
print("Proof: Let x = π/p. Then")
print("  sin(mπ/p) = sin((p-1)π/(2p)) = sin(π/2 - x/2) = cos(x/2)")
print("  sin(π/p)  = sin(x) = 2·sin(x/2)·cos(x/2)")
print()
print("  Q_1 = cos²(x/2) / sin²(x)")
print("       = cos²(x/2) / [4·sin²(x/2)·cos²(x/2)]")
print("       = 1 / [4·sin²(x/2)]")
print("       = 1 / [4·sin²(π/(2p))]   □")
print()

# Verify
print("Verification:")
for p in [7, 11, 13, 17, 23, 53, 97, 197, 499, 997]:
    if not is_prime(p): continue
    m = (p-1)//2

    # Direct computation
    Q1_direct = sin(m*pi/p)**2 / sin(pi/p)**2

    # Formula
    Q1_formula = 1 / (4 * sin(pi/(2*p))**2)

    # Ratio to check
    print(f"  p={p:4d}: Q1_direct = {Q1_direct:.8f}, formula = {Q1_formula:.8f}, match = {abs(Q1_direct - Q1_formula) < 1e-10}")

print()
print("COROLLARY: Q_1/ΣQ = 2 / [(p²-1)·sin²(π/(2p))]")
print()
print("  Since ΣQ = m(m+1)/2 = (p²-1)/8,")
print("  Q_1/ΣQ = [1/(4sin²(π/(2p)))] / [(p²-1)/8]")
print("         = 2 / [(p²-1)·sin²(π/(2p))]")
print()

print("LIMIT (PROVED):")
print(f"  Q_1/ΣQ → 8/π² = {8/pi**2:.12f} as p → ∞")
print()
print("Proof: sin(π/(2p)) = π/(2p) - O(p⁻³)")
print("  sin²(π/(2p)) = π²/(4p²) · (1 - π²/(12p²) + O(p⁻⁴))")
print("  (p²-1)·sin²(π/(2p)) = π²/4 · (1-1/p²)(1-π²/(12p²)+...) → π²/4")
print("  Q_1/ΣQ = 2/(π²/4) = 8/π²  □")
print()

# Full asymptotic expansion
print("ASYMPTOTIC EXPANSION:")
print("  Q_1/ΣQ = 8/π² + 8(12+π²)/(12π²) · p⁻² + O(p⁻⁴)")
c2 = 8*(12+pi**2)/(12*pi**2)
print(f"  = {8/pi**2:.10f} + {c2:.10f}/p² + O(p⁻⁴)")
print()

# Verify expansion
print("Verification of asymptotic expansion (should be O(p⁻⁴) residual):")
for p in [7, 11, 23, 53, 97, 197, 499, 997]:
    if not is_prime(p): continue
    m = (p-1)//2
    exact = 2 / ((p**2 - 1) * sin(pi/(2*p))**2)
    approx_0 = 8/pi**2
    approx_2 = 8/pi**2 + c2/p**2

    resid_0 = (exact - approx_0) * p**2  # Should → c2
    resid_2 = (exact - approx_2) * p**4  # Should → constant

    print(f"  p={p:4d}: exact-8/π² = {exact-approx_0:.2e}, (exact-8/π²)·p² = {resid_0:.6f} (expect {c2:.6f}), "
          f"O(p⁻⁴) residual·p⁴ = {resid_2:.4f}")

print()
print()
print("=" * 60)
print("NOW: Q_1/m² EXACT LIMIT")
print("=" * 60)
print()
print("Q_1/m² = 1/(4m²·sin²(π/(2p)))")
print("       = 1/(4·((p-1)/2)²·sin²(π/(2p)))")
print("       = 1/((p-1)²·sin²(π/(2p)))")
print()
print("As p→∞: (p-1)²·sin²(π/(2p)) → p²·(π/(2p))² = π²/4")
print(f"So Q_1/m² → 4/π² = {4/pi**2:.10f}")
print()

# But the convergence is slow — verify
print("Convergence of Q_1/m²:")
for p in [7, 11, 23, 53, 97, 197, 499, 997, 1999, 4999, 9973]:
    if not is_prime(p): continue
    m = (p-1)//2
    Q1 = 1 / (4 * sin(pi/(2*p))**2)
    ratio = Q1 / m**2
    diff = ratio - 4/pi**2
    print(f"  p={p:5d}: Q_1/m² = {ratio:.10f}, 4/π² = {4/pi**2:.10f}, diff = {diff:.2e}, diff·p = {diff*p:.4f}")

# The diff·p should converge (O(1/p) correction)
print()
print("Observation: diff ∝ 1/p. The next-order correction:")
print("  Q_1/m² = 4/π² + c₁/p + c₂/p² + ...")
print()

# More precise: Q_1/m² = 1/((p-1)²·sin²(π/(2p)))
# = 1/(p²(1-1/p)²·(π/(2p))²·(1-π²/(12p²)+...))
# = 4/(π²(1-1/p)²·(1-π²/(12p²)+...))
# = 4/π²·(1+2/p+O(1/p²))·(1+π²/(12p²)+...)
# = 4/π²·(1+2/p+(1+π²/12)/p²+...)
# So c₁ = 8/π² and diff·p → 8/π²

c1_pred = 8/pi**2
print(f"Predicted c₁ = 8/π² = {c1_pred:.6f}")
for p in [997, 1999, 4999, 9973]:
    if not is_prime(p): continue
    m = (p-1)//2
    Q1 = 1 / (4 * sin(pi/(2*p))**2)
    ratio = Q1 / m**2
    diff_p = (ratio - 4/pi**2) * p
    print(f"  p={p:5d}: diff·p = {diff_p:.6f} (expect {c1_pred:.6f})")

print()
print()
print("=" * 60)
print("TROPICAL CONSTANT κ REVISITED")
print("=" * 60)
print()
print("κ = log(max subset prod) / log(F_p)")
print("max subset prod = ∏{Q_k : Q_k > 1}")
print()
print("How many Q_k > 1?")
print("Q_k > 1 ⟺ sin²(mπk/p) > sin²(πk/p)")
print("       ⟺ |sin(mπk/p)| > |sin(πk/p)|")
print()

for p in [23, 53, 97, 197, 499, 997]:
    if not is_prime(p): continue
    m = (p-1)//2

    count_above = 0
    log_prod_above = 0
    for k in range(1, m+1):
        Q = sin(m*pi*k/p)**2 / sin(pi*k/p)**2
        if Q > 1:
            count_above += 1
            log_prod_above += log(Q)

    frac = count_above / m
    Fp = 0  # Can't compute for large p, use log approximation
    log_Fp = (p-1) * log((1+5**0.5)/2) - 0.5*log(5)  # Binet's formula approx
    kappa = log_prod_above / log_Fp

    print(f"  p={p:4d}: {count_above}/{m} modes > 1 ({frac:.4f}), "
          f"κ = {kappa:.6f}")

print()
print("Fraction of Q_k > 1 → 1/3 (confirmed from previous analysis)")
print("κ → ??? Let's check convergence more carefully:")

# Actually compute κ for many primes
kappas = []
for p in range(5, 5000):
    if not is_prime(p): continue
    m = (p-1)//2

    log_prod_above = 0
    for k in range(1, m+1):
        Q = sin(m*pi*k/p)**2 / sin(pi*k/p)**2
        if Q > 1:
            log_prod_above += log(Q)

    log_Fp = (p-1) * log((1+5**0.5)/2) - 0.5*log(5)
    kappa = log_prod_above / log_Fp
    kappas.append((p, kappa))

# Show last few
print("\nKappa values (last 10):")
for p, k in kappas[-10:]:
    print(f"  p={p:4d}: κ = {k:.8f}")

print(f"\nFinal value: κ = {kappas[-1][1]:.8f}")
print(f"2/3 = {2/3:.8f}")
print(f"Diff from 2/3: {kappas[-1][1] - 2/3:.2e}")

# Check if κ → 2/3 exactly using asymptotic analysis
# log(∏ Q_k>1 Q_k) = Σ_{Q_k>1} log Q_k
# = Σ_{Q_k>1} [log sin²(mπk/p) - log sin²(πk/p)]
# In the limit, this becomes an integral:
# m · ∫_{θ : sin²(mθ)>sin²(θ)} [log sin²(mθ) - log sin²(θ)] dθ/(2π/p)
# Hmm this is getting complicated. Let me try a different approach.

print()
print()
print("=" * 60)
print("THE FRACTION Q_k > 1: EXACT ANALYSIS")
print("=" * 60)
print()
print("Q_k > 1 ⟺ |sin(mπk/p)| > |sin(πk/p)|")
print("Let θ = πk/p ∈ (0, π/2) for k = 1,...,m.")
print("Then Q_k > 1 ⟺ |sin(mθ)| > |sin(θ)|")
print()
print("For large m, sin(mθ) oscillates rapidly. The fraction of θ ∈ (0,π/2)")
print("where |sin(mθ)| > |sin(θ)| converges to:")
print("  ∫_0^{π/2} P(|sin(mθ)| > sin(θ)) dθ · (2/π)")
print()
print("For large m, sin(mθ) is approximately equidistributed on [-1,1]")
print("with density 1/(π√(1-x²)). So:")
print("  P(|sin(mθ)| > sin(θ)) = P(|X| > sin(θ)) where X ~ arcsin distr.")
print("  = 2·(1/π)·arccos(sin(θ)) = 2·(1/π)·(π/2 - θ) = 1 - 2θ/π")
print()
print("Average over θ ~ Uniform(0, π/2):")
print("  ∫_0^{π/2} (1 - 2θ/π) · (2/π) dθ = (2/π)[θ - θ²/π]_0^{π/2}")
print("  = (2/π)(π/2 - π/4) = (2/π)(π/4) = 1/2")
print()
print("But we observed 1/3, not 1/2! The discrepancy is because the")
print("equidistribution assumption breaks down: for small θ, mθ is NOT")
print("large enough for equidistribution. The θ ≈ 0 region contributes")
print("disproportionately because Q_k > 1 there always (Q_1 is huge).")
print()
print("Let me re-examine numerically, separating the contribution of")
print("Q_1 (always > 1) from the bulk.")

# Count fraction excluding k=1
fracs_bulk = []
for p in range(5, 2000):
    if not is_prime(p): continue
    m = (p-1)//2
    if m < 3: continue

    count_above = 0
    for k in range(2, m+1):  # Skip k=1
        Q = sin(m*pi*k/p)**2 / sin(pi*k/p)**2
        if Q > 1:
            count_above += 1

    frac = count_above / (m-1)  # Out of m-1 bulk modes
    fracs_bulk.append((p, m, frac))

print("\nBulk (k≥2) fraction with Q_k > 1:")
for p, m, f in fracs_bulk[-10:]:
    print(f"  p={p:4d}: {f:.6f}")

print(f"\nBulk limit: {fracs_bulk[-1][2]:.6f}")
print(f"Is it 1/3 too? No, with k=1 included:")

# Total fraction
fracs_total = []
for p in range(5, 2000):
    if not is_prime(p): continue
    m = (p-1)//2
    if m < 3: continue

    count = sum(1 for k in range(1, m+1)
                if sin(m*pi*k/p)**2 / sin(pi*k/p)**2 > 1)
    fracs_total.append((p, count/m))

print("\nTotal fraction (all k):")
for p, f in fracs_total[-10:]:
    print(f"  p={p}: {f:.6f}")

# These should still converge to 1/3
# The actual condition is more subtle due to the Fejér kernel structure
