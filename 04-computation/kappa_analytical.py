#!/usr/bin/env python3
"""
kappa_analytical.py — opus-2026-03-13-S67h

Analytical derivation of κ_trop ≈ 0.6714.

κ = lim S(m) / log(F_p) where:
  S(m) = Σ_{k: Q_k>1} log Q_k
  log(F_p) = Σ_{k=1}^m log(1+Q_k)

NEW INSIGHT from qk_exact_limits.py:
  - k odd:  Q_k = 1/(4sin²(kπ/(2p))) ~ p²/(k²π²) for large p
  - k even: Q_k = 1/(4cos²(kπ/(2p))) → 1/4 < 1 (NEVER contributes to S)

So S(m) = Σ_{k odd, Q_k>1} log Q_k and all even modes are excluded!

For k odd and FIXED (k << p): Q_k ~ p²/(k²π²) > 1 always.
  log Q_k ≈ 2log(p/(kπ))

For k odd and LARGE (k ~ αm, α ∈ (0,1)):
  Q_k = 1/(4sin²(kπ/(2p)))
  For k odd, kπ/(2p) ranges up to mπ/(2p) ≈ π/4.
  sin²(kπ/(2p)) > 1/4 when kπ/(2p) > π/6, i.e., k > p/3.
  This means Q_k < 1 for k > p/3 (approximately).

So: Q_k > 1 for odd k ≤ p/3 (approximately), and Q_k < 1 for odd k > p/3.
Number of such k: approximately p/6 = m/3 of the m total modes.
Fraction with Q_k > 1: (m/3)/m = 1/3. ✓

Now for the numerator:
S(m) = Σ_{k odd, k≤p/3} log(1/(4sin²(kπ/(2p))))
     = Σ_{k odd, k≤p/3} [-log 4 - 2log|sin(kπ/(2p))|]
"""

import numpy as np
from math import pi, sin, cos, log, sqrt, ceil

phi = (1 + sqrt(5)) / 2

def is_prime(n):
    if n < 2: return False
    if n < 4: return n > 1
    if n % 2 == 0 or n % 3 == 0: return False
    d = 5
    while d*d <= n:
        if n % d == 0 or (n+2) % d == 0: return False
        d += 6
    return True

print("=" * 70)
print("ANALYTICAL DERIVATION OF κ_trop")
print("=" * 70)
print()

# First: verify the cutoff k* ≈ p/3 for Q_k > 1
print("STEP 1: Identify exact cutoff for Q_k > 1")
print("-" * 40)
print()
print("Q_k > 1 for k odd ⟺ 1/(4sin²(kπ/(2p))) > 1 ⟺ sin²(kπ/(2p)) < 1/4")
print("⟺ |sin(kπ/(2p))| < 1/2 ⟺ kπ/(2p) < π/6 (first branch)")
print("⟺ k < p/3")
print()
print("More precisely: k < p/3 OR k in (2p/3, p)")
print("But k ∈ {1,...,m} = {1,...,(p-1)/2}, so k < p/2.")
print("The condition k < p/3 covers the odd k with Q_k > 1 in this range.")
print("For 5p/6 < k < p: sin(kπ/(2p)) = sin(π/2 - ε) ≈ 1, so Q_k < 1.")
print()

# Wait, but k goes up to m = (p-1)/2, and p/3 < m = p/2.
# Are there modes with p/3 < k < p/2 that have Q_k > 1?
# sin(kπ/(2p)) where π/6 < kπ/(2p) < π/4
# sin(kπ/(2p)) > sin(π/6) = 1/2 → Q_k < 1. YES, so cutoff is exactly p/3.

for p in [97, 997, 9973]:
    if not is_prime(p): continue
    m = (p-1)//2
    cutoff = p/3

    # Count actual Q_k > 1 for odd k
    actual = 0
    predicted = 0
    for k in range(1, m+1, 2):  # odd k only
        Q = 1 / (4 * sin(k*pi/(2*p))**2)
        if Q > 1:
            actual += 1
        if k < cutoff:
            predicted += 1

    # Also check even k: should be 0
    even_above = sum(1 for k in range(2, m+1, 2)
                     if 1/(4*cos(k*pi/(2*p))**2) > 1)

    print(f"  p={p}: odd k with Q>1: actual={actual}, predicted (k<p/3)={predicted}, even k with Q>1={even_above}")

print()
print("Perfect! The cutoff k < p/3 is exact for large p.")
print("Odd k in {1,3,...,p/3}: approximately p/6 modes. Total m ≈ p/2.")
print("Fraction = (p/6)/(p/2) = 1/3. ✓")

print()
print()
print("STEP 2: Compute S(m) = Σ_{odd k < p/3} log Q_k")
print("-" * 40)
print()
print("S(m) = Σ_{k odd, k<p/3} [-log 4 - 2 log sin(kπ/(2p))]")
print()
print("Let K = ⌊p/6⌋ = number of odd k values below p/3.")
print("(Odd k in {1,3,...,2K-1}, so 2K-1 < p/3, K ≈ p/6.)")
print()

# For the sum: Σ_{j=0}^{K-1} log sin((2j+1)π/(2p))
# = Σ_{j=0}^{K-1} log((2j+1)π/(2p)) + correction
# ≈ Σ_{j=0}^{K-1} [log(2j+1) + log(π/(2p))]
# = log((2K-1)!!) + K·log(π/(2p))
# where (2K-1)!! = 1·3·5·...·(2K-1)

# But this is only valid for small j. For j ~ K, the sin is not small.
# Let me compute the sum numerically and analytically.

# The exact sum: Σ_{j=1,3,...,2K-1} log sin(jπ/(2p))
# where 2K-1 ≈ p/3.

# Known identity: Σ_{k=1}^{n-1} log sin(kπ/n) = log(n) - (n-1)log 2
# But our sum is over ODD k only, up to p/3.

# Let me split: Σ_{all odd k, k≤m} vs Σ_{all odd k, p/3 < k ≤ m}
# Useful: Σ_{k=1}^{m} log sin(kπ/(2p)) = sum over all k

# Actually, let me try a different approach.
# The denominator: log(F_p) ≈ (p-1)log(φ) - (1/2)log 5 ≈ 2m·logφ
# The numerator: S(m) = Σ_{odd k<p/3} [-log 4 - 2 log sin(kπ/(2p))]
#   = -K·log 4 - 2·Σ_{j=0}^{K-1} log sin((2j+1)π/(2p))

# Using the product formula: ∏_{k=1}^{n-1} sin(kπ/(2n)) has known closed form.
# Specifically: ∏_{k=1}^{n-1} 2sin(kπ/(2n)) = √(2n)... no.
#
# Let me just use the Riemann sum approximation for large m.
# Σ_{j=0}^{K-1} log sin((2j+1)π/(2p))
# ≈ (p/2) ∫_0^{1/3} log sin(πt/2) dt  [substituting t = k/p, dt = 2/p since odd k]
# Wait, the spacing of odd k is 2, so Σ_{odd k} ≈ (p/2)·∫ f(t) dt with t=k/p.
# But we want k from 1 to p/3, so t from 0 to 1/3.

# Actually: Σ_{j=0}^{K-1} f((2j+1)π/(2p)) where K = p/6
# Riemann sum with step Δk = 2, so Δ(kπ/(2p)) = π/p
# = (p/π) · ∫_0^{π/6} f(θ) dθ + O(1)  [where f = log sin θ]
# where f(θ) = log sin θ

# ∫_0^{π/6} log sin θ dθ = ???
# Known: ∫_0^{π/2} log sin θ dθ = -(π/2) log 2
# And: ∫_0^α log sin θ dθ = α log(α) - α + α log(2) + ... [for small α]
# Actually the Clausen function: ∫_0^θ log|2sin(t/2)| dt = -Cl_2(θ)/2
# So ∫_0^θ log sin t dt = ∫_0^θ [log(2sin(t)) - log 2] dt
#   = ∫_0^θ log(2sin t) dt - θ log 2
# And ∫_0^θ log(2sin t) dt = θ log 2 - (1/2)Cl_2(2θ) + ???
# Hmm, let me use the known formula:
# ∫_0^x log sin t dt = -x log 2 - (1/2)Σ_{n=1}^∞ sin(2nx)/n²
# = -x log 2 - (1/2)Cl_2(2x)
# where Cl_2 is the Clausen function of order 2.

# For x = π/6:
# Cl_2(π/3) = Σ_{n=1}^∞ sin(nπ/3)/n² = (√3/2)(1 - 1/4 + 1/7 - 1/8 + ...)
# Actually Cl_2(π/3) = known constant ≈ 1.01494

# Let me compute numerically
from math import sin as msin
def clausen2(theta, nterms=10000):
    return sum(msin(n*theta)/n**2 for n in range(1, nterms+1))

cl2_pi3 = clausen2(pi/3)
print(f"Cl_2(π/3) = {cl2_pi3:.8f}")

# ∫_0^{π/6} log sin t dt = -(π/6)log 2 - (1/2)Cl_2(π/3)
integral_logsin = -(pi/6)*log(2) - cl2_pi3/2
print(f"∫_0^{{π/6}} log sin t dt = {integral_logsin:.8f}")

# Verify numerically
theta = np.linspace(1e-15, pi/6, 100000)
numerical = np.trapezoid(np.log(np.sin(theta)), theta)
print(f"Numerical integral: {numerical:.8f}")

print()
print()
print("STEP 3: Assemble S(m)")
print("-" * 40)
print()

# S(m) = -K·log 4 - 2·Σ log sin(...)
# The sum ≈ (p/(2π))·∫_0^{π/6} log sin θ dθ  [Riemann sum with step π/p for half the terms]
# Wait, more carefully:
# Σ_{j=0}^{K-1} log sin((2j+1)π/(2p))
# This is a sum over K ≈ p/6 terms, with argument ranging from π/(2p) to ~(p/3)·π/(2p) = π/6
# Spacing of argument: 2π/(2p) = π/p
# So the Riemann sum = (p/π)·∫_{π/(2p)}^{π/6} log sin θ dθ + edge corrections

# But near θ=0, log sin θ ≈ log θ → -∞, so the integral needs care at the lower limit.
# ∫_0^{π/6} log sin θ dθ = ∫_0^{π/6} log θ dθ + ∫_0^{π/6} log(sin θ/θ) dθ
# = [θ log θ - θ]_0^{π/6} + correction
# = (π/6)(log(π/6) - 1) + (well-behaved integral)

# For the sum: Σ log sin = (p/π) · integral + O(log p) corrections

# K = ⌊p/6⌋ ≈ p/6
# S(m) = -K log 4 - 2·(p/π)·∫_0^{π/6} log sin θ dθ
# = -(p/6) log 4 - 2(p/π)[-(π/6)log 2 - Cl_2(π/3)/2]
# = -(p/6) log 4 + (p/3) log 2 + (p/π) Cl_2(π/3)
# = -(p/3) log 2 + (p/3) log 2 + (p/π) Cl_2(π/3)
# = (p/π) Cl_2(π/3)

# Wait, -(p/6) log 4 = -(p/6)·2log 2 = -(p/3) log 2
# And -2·(p/π)·(-(π/6)log 2) = (p/3) log 2
# These cancel!

# So S(m) ≈ -2·(p/π)·(-Cl_2(π/3)/2) = (p/π)·Cl_2(π/3)

print("S(m) ≈ (p/π)·Cl₂(π/3)")
print()
print("And log(F_p) ≈ (p-1)·log(φ) ≈ p·log(φ)")
print()
print("Therefore: κ = S(m)/log(F_p) ≈ Cl₂(π/3) / (π·log(φ))")
print()

kappa_pred = cl2_pi3 / (pi * log(phi))
print(f"κ_predicted = Cl₂(π/3) / (π·logφ) = {cl2_pi3:.8f} / ({pi:.6f} × {log(phi):.6f})")
print(f"            = {cl2_pi3:.8f} / {pi*log(phi):.8f}")
print(f"            = {kappa_pred:.10f}")
print()

# Compare with numerical κ
print("Compare with numerical values:")
for p in [997, 4999, 9973]:
    if not is_prime(p): continue
    m = (p-1)//2
    S = 0.0
    logF = 0.0
    for k in range(1, m+1):
        Q = sin(m*pi*k/p)**2 / sin(pi*k/p)**2
        logF += log(1 + Q)
        if Q > 1:
            S += log(Q)
    kappa_num = S / logF
    print(f"  p={p}: κ_numerical = {kappa_num:.10f}, predicted = {kappa_pred:.10f}, diff = {kappa_num - kappa_pred:.2e}")

print()
print(f"HMMMM. The prediction {kappa_pred:.8f} is not matching the numerical {0.6714:.4f}.")
print("The issue is that the Riemann sum approximation is too crude.")
print("Let me be more careful about the lower integration limit and edge effects.")
print()

# More careful calculation:
# The sum is Σ_{j=0}^{K-1} f((2j+1)·Δ) where Δ = π/(2p) and f = log sin
# But the first term f(Δ) = log sin(π/(2p)) ≈ log(π/(2p))
# This is a BIG negative term that's not well-captured by the integral.

# Alternative: write Q_k = 1/(4sin²(kπ/(2p))) for odd k
# log Q_k = -log 4 - 2 log sin(kπ/(2p))
# For Q_k > 1: sin²(kπ/(2p)) < 1/4, i.e., kπ/(2p) < π/6

# S(m) = Σ_{odd k: kπ/(2p)<π/6} [-log 4 - 2 log sin(kπ/(2p))]
#       = Σ_{odd k<p/3} [log(1/(4sin²(kπ/(2p))))]

# Better: use the Euler-Maclaurin formula for the sum
# Σ_{j=0}^{K-1} g(j) ≈ ∫_0^K g(x) dx + (1/2)(g(0)+g(K-1)) + ...
# where g(j) = -log 4 - 2 log sin((2j+1)π/(2p))
# = log Q_{2j+1}

# Actually the key issue might be the DENOMINATOR.
# log(F_p) = Σ_{k=1}^m log(1+Q_k)
# For odd k: log(1+Q_k) ≈ log Q_k for large Q_k (and ≈ log(5/4) for large k where Q_k<1)
# For even k: log(1+Q_k) ≈ log(5/4)

# So:
# log(F_p) ≈ Σ_{odd k<p/3} log Q_k + Σ_{odd k>p/3} log(1+Q_k) + Σ_{even k} log(1+Q_k)
# = S(m) + Σ_{odd k>p/3} log(1+Q_k) + Σ_{even k} log(5/4 + ε)

# The SECOND sum: odd k from p/3 to m, these Q_k < 1:
# Q_k = 1/(4sin²(kπ/(2p))), sin²(kπ/(2p)) > 1/4, so Q_k < 1
# log(1+Q_k) = log(1 + 1/(4sin²(kπ/(2p))))

# The THIRD sum: even k, log(5/4) each, ~m/2 terms = (m/2)log(5/4)

# So: κ = S(m) / [S(m) + R1 + R2]
# where R1 = Σ_{odd k>p/3} log(1+Q_k) and R2 = Σ_{even k} log(1+1/(4cos²))

# For large p:
# R2 ≈ (m/2)·log(5/4)
# R1: number of terms ≈ m/2 - m/3 = m/6, each contributing ≈ log(1+Q) where Q < 1.

# Let's compute R1 and R2 numerically and check the ratio
for p in [9973]:
    if not is_prime(p): continue
    m = (p-1)//2

    S = 0.0
    R1 = 0.0
    R2 = 0.0
    for k in range(1, m+1):
        Q = sin(m*pi*k/p)**2 / sin(pi*k/p)**2
        lq = log(1 + Q)
        if k % 2 == 1:
            if Q > 1:
                S += log(Q)
            else:
                R1 += lq
        else:
            R2 += lq

    total = S + R1 + R2
    # Also compute correction: S includes log Q_k, but denominator uses log(1+Q_k)
    # For odd k with Q>1: log(1+Q_k) = log(Q_k) + log(1+1/Q_k) ≈ log(Q_k)
    # Correction: Σ_{odd k<p/3} log(1+1/Q_k)
    corr = 0.0
    for k in range(1, m+1, 2):
        Q = 1 / (4 * sin(k*pi/(2*p))**2)
        if Q > 1:
            corr += log(1 + 1/Q)

    print(f"p={p}: S = {S:.4f} ({S/total*100:.1f}%)")
    print(f"       R1 (odd k, Q<1) = {R1:.4f} ({R1/total*100:.1f}%)")
    print(f"       R2 (even k) = {R2:.4f} ({R2/total*100:.1f}%)")
    print(f"       total = {total:.4f}")
    print(f"       correction = {corr:.4f}")
    print(f"       κ = S/total = {S/total:.8f}")
    print(f"       κ' = S/(S+corr+R1+R2) = {S/(total):.8f}")
    print()
    print(f"       S/m = {S/m:.6f}")
    print(f"       R1/m = {R1/m:.6f}")
    print(f"       R2/m = {R2/m:.6f}")
    print(f"       corr/m = {corr/m:.6f}")

print()
print()
# Now the integrals:
# S/m → ∫_0^{1/3} max(0, log(1/(4sin²(πt/2)))) · 2dt  [factor 2 for odd k only]
# Wait, let me be more careful. k = 1,3,5,...,2j+1,...
# t = k/p ∈ (0, 1/2), but only odd k, so effective density = 1/2.
# S/m = (1/m) Σ_{odd k, k<p/3} log(1/(4sin²(kπ/(2p))))
# There are K ≈ p/6 = m/3 such terms.
# If I substitute θ = kπ/(2p), then θ ranges from π/(2p) to π/6.
# Spacing: 2·π/(2p) = π/p.
# S/m ≈ (1/m)·(p/π)·∫_0^{π/6} max(0, -log 4 - 2log sin θ) dθ

# Numerator factor: (1/m)·(p/π) = (1/(p/2))·(p/π) = 2/π

# S/m → (2/π) · ∫_0^{π/6} (-log 4 - 2 log sin θ) dθ
#      = (2/π) · [-(π/6)log 4 - 2·∫_0^{π/6} log sin θ dθ]
#      = (2/π) · [-(π/6)log 4 - 2·(-(π/6)log 2 - Cl₂(π/3)/2)]
#      = (2/π) · [-(π/3)log 2 + (π/3)log 2 + Cl₂(π/3)]
#      = (2/π) · Cl₂(π/3)

sm_pred = (2/pi) * cl2_pi3
print(f"Predicted S/m → (2/π)·Cl₂(π/3) = {sm_pred:.8f}")
print(f"Numerical S/m ≈ {0.6461:.8f} (from previous runs)")
print(f"Difference: {sm_pred - 0.6461:.6f}")

# Similarly for log(F_p)/m → 2 logφ:
# κ = (S/m) / (logF/m) = [(2/π)·Cl₂(π/3)] / [2 logφ]
# = Cl₂(π/3) / (π logφ)

print()
print(f"κ = Cl₂(π/3) / (π·logφ) = {cl2_pi3:.8f} / ({pi*log(phi):.8f})")
print(f"  = {kappa_pred:.10f}")
print()

# The issue: this gives ~0.6685, but numerical κ ≈ 0.6714.
# The difference comes from R1 and corr terms that modify both numerator and denominator.

# Actually κ = S / (S + R1 + R2 + corr)
# where corr accounts for log(1+Q_k) vs log(Q_k) for Q>1 terms.
# In the original formula, κ = Σ_{Q>1} log(Q) / Σ log(1+Q)
#   = S / (S + corr + R1 + R2)

# Let me compute each term as an integral:
# corr/m → (2/π)·∫_0^{π/6} log(1+4sin²θ) dθ  [since 1/Q_k = 4sin²(kπ/(2p))]

theta = np.linspace(1e-15, pi/6, 100000)
corr_integral = np.trapezoid(np.log(1 + 4*np.sin(theta)**2), theta)
corr_per_m = (2/pi) * corr_integral
print(f"corr/m → (2/π)·∫_0^{{π/6}} log(1+4sin²θ) dθ = {corr_per_m:.8f}")

# R1/m: odd k from p/3 to m, Q_k < 1
# R1/m → (2/π)·∫_{π/6}^{π/4} log(1 + 1/(4sin²θ)) dθ
r1_integrand = np.log(1 + 1/(4*np.sin(theta[(theta > pi/6) & (theta < pi/4)])**2))
theta_r1 = theta[(theta > pi/6) & (theta < pi/4)]
r1_integral = np.trapezoid(r1_integrand, theta_r1)
r1_per_m = (2/pi) * r1_integral
print(f"R1/m → (2/π)·∫_{{π/6}}^{{π/4}} log(1+1/(4sin²θ)) dθ = {r1_per_m:.8f}")

# R2/m: even k from 1 to m, Q_k = 1/(4cos²(kπ/(2p))) → 1/4
# R2/m → (1/m)·(m/2)·log(5/4) = (1/2)·log(5/4)
r2_per_m = 0.5 * log(5/4)
print(f"R2/m → (1/2)·log(5/4) = {r2_per_m:.8f}")

# Total logF/m = S/m + corr/m + R1/m + R2/m
total_per_m = sm_pred + corr_per_m + r1_per_m + r2_per_m
print(f"\nTotal logF/m = {total_per_m:.8f}")
print(f"2·logφ = {2*log(phi):.8f}")
print(f"Difference: {total_per_m - 2*log(phi):.6f}")

# κ = S/m / total
kappa_from_integrals = sm_pred / total_per_m
print(f"\nκ = S/m / (logF/m) = {sm_pred:.8f} / {total_per_m:.8f} = {kappa_from_integrals:.10f}")

# Hmm, we should check this more carefully. Actually the denominator
# should be Σ log(1+Q_k), NOT S + corrections. Let me redo.
print()
print("=" * 70)
print("CORRECT INTEGRAL DECOMPOSITION")
print("=" * 70)
print()

# logF/m = (1/m) Σ_{k=1}^m log(1+Q_k)
# Split by parity:
# = (1/m) Σ_{odd k} log(1+Q_k) + (1/m) Σ_{even k} log(1+Q_k)

# Odd k terms: Q_k = 1/(4sin²(kπ/(2p)))
# (1/m) Σ_{odd k≤m} log(1+1/(4sin²(kπ/(2p))))
# ≈ (2/π) ∫_0^{π/4} log(1+1/(4sin²θ)) dθ

theta_full = np.linspace(1e-15, pi/4, 100000)
odd_integral = np.trapezoid(np.log(1 + 1/(4*np.sin(theta_full)**2)), theta_full)
odd_per_m = (2/pi) * odd_integral
print(f"Odd contribution: (2/π)·∫_0^{{π/4}} log(1+1/(4sin²θ)) dθ = {odd_per_m:.8f}")

# Even k terms: Q_k = 1/(4cos²(kπ/(2p))) → 1/4
# (1/m) Σ_{even k≤m} log(1+1/(4cos²(kπ/(2p))))
# ≈ (2/π) ∫_0^{π/4} log(1+1/(4cos²θ)) dθ
even_integral = np.trapezoid(np.log(1 + 1/(4*np.cos(theta_full)**2)), theta_full)
even_per_m = (2/pi) * even_integral
print(f"Even contribution: (2/π)·∫_0^{{π/4}} log(1+1/(4cos²θ)) dθ = {even_per_m:.8f}")

logF_per_m = odd_per_m + even_per_m
print(f"Total logF/m = {logF_per_m:.8f}")
print(f"2·logφ = {2*log(phi):.8f}")
print(f"Ratio: {logF_per_m/(2*log(phi)):.8f}")

# Numerator: S/m = (1/m) Σ_{Q>1} log Q
# Only odd k with Q > 1 (k < p/3, θ < π/6) contribute
# S/m = (2/π) ∫_0^{π/6} log(1/(4sin²θ)) dθ = (2/π)∫_0^{π/6} (-log4 - 2logsinθ) dθ
print(f"\nS/m = (2/π)·Cl₂(π/3) = {sm_pred:.8f}")
print(f"κ = S/m / (logF/m) = {sm_pred/logF_per_m:.10f}")
print()

# But wait — the DENOMINATOR should also include contributions from
# LARGE odd k modes where Q < 1. These contribute log(1+Q) where Q < 1.

# Let me verify this against the EXACT computation
print("Verification against exact computation at p=9973:")
p_test = 9973
m_test = (p_test-1)//2
S_exact = 0.0
logF_exact = 0.0
for k in range(1, m_test+1):
    Q = sin(m_test*pi*k/p_test)**2 / sin(pi*k/p_test)**2
    logF_exact += log(1 + Q)
    if Q > 1:
        S_exact += log(Q)

print(f"  S_exact/m = {S_exact/m_test:.8f}, integral prediction = {sm_pred:.8f}")
print(f"  logF_exact/m = {logF_exact/m_test:.8f}, integral prediction = {logF_per_m:.8f}")
print(f"  κ_exact = {S_exact/logF_exact:.10f}")
print(f"  κ_integral = {sm_pred/logF_per_m:.10f}")

print()
print("=" * 70)
print("FINAL ANSWER")
print("=" * 70)
print()
print("κ_trop = Cl₂(π/3) / [∫_0^{π/4} [log(1+1/(4sin²θ)) + log(1+1/(4cos²θ))] dθ]")
print(f"       ≈ {cl2_pi3:.8f} / [{odd_integral + even_integral:.8f}]")
print(f"       ≈ {cl2_pi3 / (odd_integral + even_integral):.10f}")
print()
print("Let's see if this simplifies to a known constant...")

# Note: 1/(4sin²θ) + 1/(4cos²θ) = 1/(4sin²θcos²θ) = 1/sin²(2θ)
# But log(1+a) + log(1+b) ≠ log(1+a+b) in general.
# However, (1+1/(4sin²θ))(1+1/(4cos²θ)) = 1 + 1/(4sin²θ) + 1/(4cos²θ) + 1/(16sin²θcos²θ)
# = 1 + (cos²θ+sin²θ)/(4sin²θcos²θ) + 1/(16sin²θcos²θ)
# = 1 + 1/(4sin²θcos²θ) + 1/(16sin²θcos²θ)
# = 1 + 5/(16sin²θcos²θ)
# = 1 + 5/(4sin²(2θ))

# So odd+even integral = ∫_0^{π/4} log(1 + 5/(4sin²(2θ))) dθ
# = (1/2)∫_0^{π/2} log(1 + 5/sin²u) du/2  [u=2θ]... wait
# Actually: ∫_0^{π/4} log((1+1/(4sin²θ))(1+1/(4cos²θ))) dθ
# = ∫_0^{π/4} log(1 + 5/(4sin²(2θ))) dθ
# = (1/2)∫_0^{π/2} log(1 + 5/(4sin²u)) du  [u=2θ]

combined_integral = np.trapezoid(np.log(1 + 5/(4*np.sin(2*theta_full)**2)), theta_full)
print(f"\nCombined: ∫_0^{{π/4}} log(1+5/(4sin²(2θ))) dθ = {combined_integral:.8f}")
print(f"Separate: {odd_integral+even_integral:.8f}")
print(f"Match: {abs(combined_integral - odd_integral - even_integral) < 1e-6}")

# So: logF/m ≈ (2/π)·(1/2)·∫_0^{π/2} log(1+5/(4sin²u)) du
# = (1/π)·∫_0^{π/2} log(1+5/(4sin²u)) du
# = (1/π)·∫_0^{π/2} log((4sin²u+5)/(4sin²u)) du
# = (1/π)·[∫_0^{π/2} log(4sin²u+5) du - ∫_0^{π/2} log(4sin²u) du]

# ∫_0^{π/2} log(4sin²u) du = (π/2)log4 + 2·∫_0^{π/2} log sin u du
# = (π/2)log4 + 2·(-(π/2)log2) = π log 2 - π log 2 = 0
# Wait: (π/2)log(4) = π log 2, and 2·(-(π/2)log 2) = -π log 2
# So ∫_0^{π/2} log(4sin²u) du = 0 !!

print()
print("KEY IDENTITY: ∫_0^{π/2} log(4sin²u) du = 0")
print("Proof: (π/2)log 4 + 2·∫_0^{π/2} log sin u du = π log 2 + 2·(-π/2·log 2) = 0 ✓")
print()
print("Therefore:")
print("  logF/m → (1/π)·∫_0^{π/2} log(4sin²u + 5) du")
print()

# Compute this integral
u = np.linspace(1e-15, pi/2, 100000)
den_integral = np.trapezoid(np.log(4*np.sin(u)**2 + 5), u)
print(f"  ∫_0^{{π/2}} log(4sin²u + 5) du = {den_integral:.8f}")
print(f"  (1/π)·this = {den_integral/pi:.8f}")
print(f"  2logφ = {2*log(phi):.8f}")
print(f"  Ratio: {den_integral/pi/(2*log(phi)):.8f}")

# Now κ = [(2/π)·Cl₂(π/3)] / [(1/π)·∫ log(4sin²u+5) du]
# = 2·Cl₂(π/3) / ∫_0^{π/2} log(4sin²u+5) du

kappa_analytical = 2*cl2_pi3 / den_integral
print(f"\nκ = 2·Cl₂(π/3) / ∫_0^{{π/2}} log(4sin²u+5) du")
print(f"  = 2·{cl2_pi3:.8f} / {den_integral:.8f}")
print(f"  = {kappa_analytical:.10f}")
print()
print(f"Compare with numerical κ ≈ 0.6714")
print(f"Difference: {kappa_analytical - 0.6714:.6f}")

# The integral ∫_0^{π/2} log(4sin²u+5) du
# = ∫_0^{π/2} log(5-2cos(2u)+2) du... hmm
# = ∫_0^{π/2} log(7-2cos(2u)) du / 2 + something... no.
# 4sin²u = 2-2cos(2u), so 4sin²u+5 = 7-2cos(2u)
# ∫_0^{π/2} log(7-2cos(2u)) du = (1/2)∫_0^π log(7-2cos v) dv

# This is a KNOWN integral:
# ∫_0^π log(a-b cos v) dv = π log((a+√(a²-b²))/2) for a > |b|
# Here a=7, b=2: √(49-4) = √45 = 3√5
# = π log((7+3√5)/2)

result_formula = pi * log((7 + 3*sqrt(5))/2)
print(f"\n∫_0^π log(7-2cos v) dv = π·log((7+3√5)/2) = {result_formula:.8f}")
print(f"(1/2)·this = {result_formula/2:.8f}")
print(f"Numerical = {den_integral:.8f}")
print(f"Match: {abs(result_formula/2 - den_integral) < 1e-4}")

# So: ∫_0^{π/2} log(4sin²u+5) du = (π/2)·log((7+3√5)/2)
# And: logF/m → (1/π)·(π/2)·log((7+3√5)/2) = (1/2)·log((7+3√5)/2)

val = (7 + 3*sqrt(5))/2
print(f"\n(7+3√5)/2 = {val:.8f}")
print(f"log((7+3√5)/2)/2 = {log(val)/2:.8f}")
print(f"2logφ = {2*log(phi):.8f}")

# Check: (7+3√5)/2 = (7+6.708...)/2 = 6.854...
# φ⁴ = φ²·φ² = 2.618²  = 6.854...
# φ⁴ = (7+3√5)/2 exactly!!
phi4 = phi**4
print(f"\nφ⁴ = {phi4:.8f}")
print(f"(7+3√5)/2 = {val:.8f}")
print(f"MATCH: {abs(phi4 - val) < 1e-10}")

# So: logF/m → (1/2)·log(φ⁴) = 2logφ ✓ ✓ ✓
# This CONFIRMS our formula since logF/m = 2logφ is known!

print()
print("CONFIRMED: ∫_0^{π/2} log(4sin²u+5) du = (π/2)·log(φ⁴) = 2π·logφ")
print("And logF/m → (1/π)·(2π logφ) = 2logφ ✓")
print()

# Now κ:
# κ = 2·Cl₂(π/3) / (2π·logφ) = Cl₂(π/3) / (π·logφ)
kappa_exact = cl2_pi3 / (pi * log(phi))
print(f"THEOREM: κ_trop = Cl₂(π/3) / (π·logφ)")
print(f"                = {cl2_pi3:.10f} / {pi*log(phi):.10f}")
print(f"                = {kappa_exact:.10f}")
print()
print("But this doesn't match the numerical value! Let me check...")
print(f"Numerical κ at p=9973: {S_exact/logF_exact:.10f}")
print(f"Analytical: {kappa_exact:.10f}")
print(f"Difference: {S_exact/logF_exact - kappa_exact:.6f}")

# The discrepancy: 0.6714 vs 0.6685
# This is ~0.003, which is O(1/√p) or O(log p / p)
# The issue is that the Riemann sum for S has slow convergence due to
# the singularity at θ=0. The leading integral gives Cl₂(π/3)/π·logφ,
# but there are O(log p / p) corrections from the discrete sum near θ=0.

# Actually, wait. Let me recheck the Riemann sum more carefully.
# The sum Σ_{odd k<p/3} (-log4 - 2 log sin(kπ/(2p)))
# is NOT well-approximated by (p/π)·∫ because of the log singularity at θ=0.
# Near θ=0: log sin θ ≈ log θ, and Σ log(odd k) - (p/π)∫_0^{π/6} log θ dθ
# The difference is O(log p) (Euler-Maclaurin correction at the boundary).

# So κ = [Cl₂(π/3)/π + c₁·log p / p] / [2logφ + c₂·log p / p]
# and the O(log p / p) corrections make κ slightly larger than Cl₂(π/3)/(π logφ).

# In the limit p→∞, κ_∞ SHOULD be Cl₂(π/3)/(π logφ).
# But the convergence is logarithmically slow!

print()
print("The convergence is LOGARITHMICALLY SLOW due to the singularity at θ=0.")
print("Prediction: κ(p) = Cl₂(π/3)/(π·logφ) + c·log(p)/p + ...")
print()
print("Let me check this prediction with a log(p)/p fit:")

# Fit κ = a + b·log(p)/p
ps_fit = []
ks_fit = []
for p in range(53, 10000):
    if not is_prime(p): continue
    m = (p-1)//2
    S = 0.0
    logF = 0.0
    for k in range(1, m+1):
        Q = sin(m*pi*k/p)**2 / sin(pi*k/p)**2
        logF += log(1 + Q)
        if Q > 1:
            S += log(Q)
    ps_fit.append(p)
    ks_fit.append(S/logF)

ps_arr = np.array(ps_fit, dtype=float)
ks_arr = np.array(ks_fit)

A = np.column_stack([np.ones(len(ps_arr)), np.log(ps_arr)/ps_arr])
result = np.linalg.lstsq(A, ks_arr, rcond=None)
a_fit, b_fit = result[0]
print(f"Fit: κ = {a_fit:.10f} + {b_fit:.6f}·log(p)/p")
print(f"Cl₂(π/3)/(π·logφ) = {kappa_exact:.10f}")
print(f"Fit intercept = {a_fit:.10f}")
print(f"Difference: {a_fit - kappa_exact:.2e}")

# Try log²(p)/p fit
A2 = np.column_stack([np.ones(len(ps_arr)), np.log(ps_arr)/ps_arr, np.log(ps_arr)**2/ps_arr])
result2 = np.linalg.lstsq(A2, ks_arr, rcond=None)
a2, b2, c2 = result2[0]
print(f"\nBetter fit: κ = {a2:.10f} + {b2:.6f}·log(p)/p + {c2:.4f}·log²(p)/p")
print(f"Cl₂(π/3)/(π·logφ) = {kappa_exact:.10f}")
