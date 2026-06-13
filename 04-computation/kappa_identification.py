#!/usr/bin/env python3
"""
kappa_identification.py — opus-2026-03-13-S67h

Identify the exact value of the tropical constant κ ≈ 0.67143.

κ = log(∏{Q_k>1} Q_k) / log(F_p) for Interval tournament as p → ∞

Strategy: compute κ to high precision, check against known constants,
and try to derive analytically from the Fejér kernel structure.

We know:
  Q_k = sin²(mπk/p) / sin²(πk/p)  where m = (p-1)/2
  F_p = ∏_{k=1}^m (1 + Q_k) = Fibonacci(p)
  κ = Σ_{k: Q_k>1} log Q_k / Σ_{k=1}^m log(1+Q_k)

In the continuum limit (p → ∞, θ = 2πk/p):
  Q(θ) = sin²(mθ/2) / sin²(θ/2) = Fejér kernel F_m(θ)

  numerator: ∫_0^π 1_{Q(θ)>1} · log Q(θ) dθ/(2π/p)
  denominator: ∫_0^π log(1 + Q(θ)) dθ/(2π/p)

Both scale as m times an integral, so κ = I_num / I_den where:
  I_num = ∫_0^π 1_{F_m(θ)>1} · log F_m(θ) dθ · m/π
  I_den = ∫_0^π log(1+F_m(θ)) dθ · m/π

Hmm, but these integrals depend on m. Let me think differently.

For large m, log F_p ≈ 2m·log(φ) and the log-product of Q>1 also
scales linearly in m. So κ = lim Σ_{Q_k>1} log Q_k / (2m log φ).
"""

import numpy as np
from math import pi, sin, cos, log, sqrt, factorial
# from scipy import integrate  # not available

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

# High-precision computation of κ for large p
print("HIGH-PRECISION κ COMPUTATION")
print("=" * 60)

kappas = []
for p in range(5, 10000):
    if not is_prime(p): continue
    m = (p-1)//2

    log_prod_above = 0.0
    log_Fp = 0.0
    for k in range(1, m+1):
        smk = sin(m*pi*k/p)
        sk = sin(pi*k/p)
        Q = smk**2 / sk**2
        log_Fp += log(1 + Q)
        if Q > 1:
            log_prod_above += log(Q)

    kappa = log_prod_above / log_Fp
    kappas.append((p, m, kappa))

# Show convergence
for p, m, k in kappas[-20:]:
    print(f"  p={p:5d}, m={m:4d}: κ = {k:.10f}")

kappa_final = kappas[-1][2]
print(f"\nBest estimate: κ = {kappa_final:.10f}")

# Check known constants
print(f"\n  2/3 = {2/3:.10f}  (diff = {kappa_final - 2/3:.2e})")
print(f"  log(2) = {log(2):.10f}  (diff = {kappa_final - log(2):.2e})")
print(f"  1 - 1/(2e) = {1 - 1/(2*np.e):.10f}  (diff = {kappa_final - (1-1/(2*np.e)):.2e})")
print(f"  2log(φ)/π = {2*log(phi)/pi:.10f}  (diff = {kappa_final - 2*log(phi)/pi:.2e})")
print(f"  3/π + 1/6 = {3/pi + 1/6:.10f}")
print(f"  1 - 1/e = {1 - 1/np.e:.10f}")
print(f"  1 - 1/3 = {2/3:.10f}")
print(f"  π/4 - 1/12 = {pi/4 - 1/12:.10f}")
print(f"  1/log(4) = {1/log(4):.10f}")

# Try OEIS-style: look at 0.67143
print(f"\n  0.67143... Let's get more digits and try to identify")

# Richardson extrapolation for better estimate
# κ(p) = κ_∞ + c/p + d/p² + ...
# Use last several values
ps = [p for p, m, k in kappas[-100:]]
ks = [k for p, m, k in kappas[-100:]]

# Fit κ = a + b/p + c/p²
import numpy as np
A = np.column_stack([np.ones(len(ps)), 1/np.array(ps, dtype=float), 1/np.array(ps, dtype=float)**2])
result = np.linalg.lstsq(A, ks, rcond=None)
a, b, c = result[0]
print(f"\nRichardson extrapolation: κ_∞ = {a:.12f}")
print(f"  Correction: {b:.6f}/p + {c:.4f}/p²")

# Try inverse-square-root fit
A2 = np.column_stack([np.ones(len(ps)), 1/np.sqrt(np.array(ps, dtype=float))])
result2 = np.linalg.lstsq(A2, ks, rcond=None)
print(f"  Alt fit (1/√p): κ_∞ = {result2[0][0]:.12f}")

# Now try to derive analytically
print()
print("=" * 60)
print("ANALYTICAL DERIVATION ATTEMPT")
print("=" * 60)
print()

# Key formula: κ = lim_{m→∞} S(m) / log(F_{2m+1})
# where S(m) = Σ_{k: Q_k(m)>1} log Q_k(m)
# and Q_k(m) = sin²(mπk/(2m+1)) / sin²(πk/(2m+1))
# (using p = 2m+1)
#
# In the large-m limit, θ_k = πk/(2m+1) ≈ πk/(2m) and the sum becomes
# a Riemann sum of an integral.
#
# S(m) = Σ_{k=1}^m 1_{Q_k>1} · [log sin²(mθ_k) - log sin²(θ_k)]
# where θ_k = πk/p.
#
# Substituting t = k/m, so t ranges over {1/m, 2/m, ..., 1}:
#   θ_k = πk/p ≈ πt/2
#   mθ_k ≈ mπt/2
#   sin(mθ_k) = sin(mπt/2)
#
# S(m) ≈ Σ_{t} [log sin²(mπt/2) - log sin²(πt/2)] · 1_{...}
#
# For the denominator:
# log F_p = Σ_{k=1}^m log(1 + Q_k) ≈ 2m log φ
#
# So we need S(m) / (2m log φ)
# Hmm, S(m) grows linearly in m (checked below)

print("Growth of S(m) = Σ_{Q>1} log Q:")
for p, m, kap in kappas:
    if m in [10, 20, 50, 100, 200, 500, 1000, 2000, 5000]:
        S_m = kap * (p-1) * log(phi)  # approx log F_p ≈ (p-1)logφ
        print(f"  m={m:4d}: S(m) = {S_m:.4f}, S(m)/m = {S_m/m:.6f}, S(m)/(m*logφ) = {S_m/(m*log(phi)):.6f}")

print()

# More precisely: let's compute S(m)/m directly
print("S(m)/m convergence:")
s_over_m = []
for p, m, kap in kappas:
    if m < 10: continue
    # S(m) = kappa * log(F_p)
    # log(F_p) = Σ log(1+Q_k), already computed implicitly
    # Instead recompute
    S_m = 0.0
    for k in range(1, m+1):
        Q = sin(m*pi*k/p)**2 / sin(pi*k/p)**2
        if Q > 1:
            S_m += log(Q)
    s_over_m.append((m, S_m/m))

for m, ratio in s_over_m[-10:]:
    print(f"  m={m:4d}: S(m)/m = {ratio:.8f}")

# Does S(m)/m converge?
print(f"\nS(m)/m at m={s_over_m[-1][0]}: {s_over_m[-1][1]:.8f}")
print(f"2*κ*log(φ) = {2*a*log(phi):.8f}")  # Should match

# Actually, κ = S(m) / log F_p, and log F_p ≈ 2m log φ
# So S(m)/m ≈ 2κ log φ
print(f"Predicted S(m)/m → 2κ_∞·logφ = {2*a*log(phi):.8f}")

print()
print("=" * 60)
print("THE INTEGRAL REPRESENTATION")
print("=" * 60)
print()

# In the continuum limit:
# S(m)/m → ∫_0^1 1_{|sin(πmt)|>|sin(πt/2)|} · [2log|sin(πmt)| - 2log|sin(πt/2)|] dt
#
# But sin(πmt) oscillates wildly. The standard approach:
# Use equidistribution of {mt mod 1}: for irrational multiplier m,
# sin(πmt) is equidistributed.
#
# But m here IS the variable going to ∞ and we sum over k, so
# the equidistribution is in k (not m). For each fixed t = k/m,
# as m → ∞, we need sin(mπt/2) which is sin(πk/2) — FIXED!
#
# Wait, that's not right. Let me redo more carefully.
# θ_k = πk/p = πk/(2m+1), so
# mθ_k = mπk/(2m+1) = πk(1 - 1/(2m+1))/2 → πk/2 as m → ∞
# But k ranges from 1 to m, so k/m ranges from 1/m to 1.
#
# For k << m: mθ_k ≈ πk/2 which grows.
# For k = αm: mθ_k ≈ παm/2, which is O(m) — oscillates rapidly.
#
# The condition Q_k > 1 becomes |sin(mπk/(2m+1))| > |sin(πk/(2m+1))|
# For k = O(1): Q_k ≈ sin²(πk/2) / sin²(πk/(2m)) ≈ sin²(πk/2) · (2m/(πk))²
# → ∞ for k odd, → 0 for k even!
#
# For k = αm, α ∈ (0,1): Q_k = sin²(mπαm/(2m+1)) / sin²(παm/(2m+1))
#   ≈ sin²(πα m²/(2m)) / sin²(πα/2)
#   = sin²(παm/2) / sin²(πα/2)
#   This oscillates between 0 and 1/sin²(πα/2) ≤ 1 for α > 1/3 (roughly).
#
# Hmm, this is quite complex. Let me try numerical integration instead.

# For very large m, compute the "integral" directly
m_test = 10000
p_test = 2*m_test + 1  # Not prime, but OK for the limit analysis
S = 0.0
N_above = 0
for k in range(1, m_test+1):
    Q = sin(m_test*pi*k/p_test)**2 / sin(pi*k/p_test)**2
    if Q > 1:
        S += log(Q)
        N_above += 1

logF = sum(log(1 + sin(m_test*pi*k/p_test)**2/sin(pi*k/p_test)**2)
           for k in range(1, m_test+1))

kappa_integral = S / logF
print(f"m = {m_test}: κ = {kappa_integral:.10f}")
print(f"  S/m = {S/m_test:.8f}")
print(f"  logF/m = {logF/m_test:.8f}")
print(f"  2·log(φ) = {2*log(phi):.8f}")
print(f"  fraction above: {N_above/m_test:.6f}")

# Also try with the EXACT Fejér kernel interpretation
# Q_k = F_m(2πk/p) where F_m is the Fejér kernel.
# log F_p = Σ log(1 + F_m(2πk/p))
# This is a Riemann sum for (p/(2π)) ∫_0^{2π} log(1+F_m(θ)) dθ
# But p/(2π) = (2m+1)/(2π) ≈ m/π
# So logF ≈ (m/π) ∫_0^{2π} log(1+F_m(θ)) dθ

# The integral ∫_0^{2π} log(1+F_m(θ)) dθ / (2π)
# For the Fejér kernel F_m(θ) = (sin(mθ/2)/sin(θ/2))²:

# This integral is known! It equals 2·log something.
# Actually, ∫_0^{2π} log(1 + |D_m(θ)|²) dθ where D_m is Dirichlet kernel...

# Let me compute it numerically for comparison
def integral_log_fejer(m, npts=100000):
    theta = np.linspace(1e-10, pi, npts)
    Q = (np.sin(m*theta) / np.sin(theta))**2  # Note: this is Dirichlet²
    # Actually Fejér is sin²(mθ/2)/sin²(θ/2)
    Q_fejer = (np.sin(m*theta/2) / np.sin(theta/2))**2
    integrand = np.log(1 + Q_fejer)
    return np.trapz(integrand, theta) / pi  # ∫_0^π / π = average over half-period

for m in [10, 50, 100, 500, 1000]:
    val = integral_log_fejer(m)
    print(f"  m={m:4d}: (1/π)∫_0^π log(1+F_m(θ))dθ = {val:.6f}, val/log(m) = {val/log(m):.6f}, val/(2logφ) = {val/(2*log(phi)):.6f}")

print()
print("Observation: the integral grows as 2·log(m) + const, not 2·log(φ)!")
print("This is because logF_p ≈ 2m·logφ but the RIEMANN SUM has m terms,")
print("each contributing ~2logφ/1 in the average, while the INTEGRAL is")
print("dominated by the peak at θ=0 which contributes ~log(m²).")
print()
print("The resolution: the Riemann sum ≠ integral for sharply peaked functions.")
print("The peak of F_m at θ=0 has width O(1/m), capturing ~O(1) sample points")
print("near θ=0, each contributing O(log m²). But there are ~m points total,")
print("each contributing O(1) on average.")

print()
print("=" * 60)
print("ALTERNATIVE: κ AS AN INTEGRAL")
print("=" * 60)

# κ = Σ_{Q_k>1} log Q_k / Σ log(1+Q_k)
# Both sums have m terms. In the limit, this ratio is NOT an integral/integral.
# It's a discrete ratio. But we can write:
# κ = <log Q · 1_{Q>1}> / <log(1+Q)>
# where <·> is average over k.
#
# The denominator: <log(1+Q_k)> = (1/m)Σ log(1+Q_k) = logF_p/m ≈ 2logφ
#
# The numerator: (1/m)Σ_{Q_k>1} log Q_k → some limit L
# κ_∞ = L / (2logφ)
#
# From the data: S(m)/m → 2κ_∞·logφ ≈ 2·0.6714·0.4812 ≈ 0.646

S_over_m_values = []
logF_over_m_values = []
for p_test in [100003, 50021, 20011, 10007, 5003, 2003]:
    # Use closest prime
    m = (p_test-1)//2
    p_act = p_test

    S = 0.0
    logF = 0.0
    for k in range(1, m+1):
        Q = sin(m*pi*k/p_act)**2 / sin(pi*k/p_act)**2
        logF += log(1 + Q)
        if Q > 1:
            S += log(Q)

    S_over_m_values.append((m, S/m))
    logF_over_m_values.append((m, logF/m))

print("\nConvergence of averages:")
for (m, sm), (_, lm) in zip(S_over_m_values, logF_over_m_values):
    print(f"  m={m:5d}: S/m = {sm:.8f}, logF/m = {lm:.8f}, κ = {sm/lm:.8f}")

print(f"\n  2·logφ = {2*log(phi):.8f}")

# The key: logF/m → 2logφ (this is known from F_p = φ^p/√5)
# So κ = lim S(m)/m / (2logφ) = lim S(m)/(2m·logφ)

print()
print("=" * 60)
print("TRYING TO IDENTIFY κ EXACTLY")
print("=" * 60)
print()
print(f"κ_∞ = {a:.12f}")

# Try: is κ related to a known integral involving sin functions?
#
# Specifically:
# S(m)/m = (1/m) Σ_{k: Q_k>1} log Q_k
# = (1/m) Σ_k max(0, log Q_k)  [since log Q > 0 iff Q > 1]
# + some correction from modes with Q exactly 1 (measure zero)
#
# Actually NOT. Some Q_k < 1 have log Q_k < 0, we're summing only the positive ones.
# But it's NOT the same as max(0, log Q_k) because there may be Q_k in (0,1) that we skip.
# Well, 1_{Q>1}·log Q = max(0, log Q) for Q > 0, so yes it IS max(0, log Q_k).

# So S(m)/m = (1/m) Σ_{k=1}^m max(0, log Q_k)
# = (1/m) Σ max(0, 2log|sin(mπk/p)| - 2log|sin(πk/p)|)

# Using the substitution θ = πk/p, this is a Riemann sum for
# (p/π) · (1/m) · ∫_0^{π/2} max(0, 2log|sin(mθ)| - 2log|sin(θ)|) dθ/(p/(πm))
# ≈ ∫_0^{π/2} max(0, 2log|sin(mθ)| - 2log|sin(θ)|) dθ / (π/(2m))

# This is getting hairy. Let me try ISC (Inverse Symbolic Calculator) approach.
# κ = 0.671429...

# Let me try more combinations
import itertools

target = a  # Our best κ estimate
print(f"Target: {target:.10f}")
print()

# Combinations of π, φ, e, log, sqrt
candidates = {
    "2/3": 2/3,
    "log(2)": log(2),
    "1-1/(2e)": 1-1/(2*np.e),
    "2logφ/π": 2*log(phi)/pi,
    "π/(2e)": pi/(2*np.e),
    "1-1/π": 1-1/pi,
    "2/(1+√2)": 2/(1+sqrt(2)),
    "3-e": 3-np.e,
    "φ-1+1/π": phi-1+1/pi,
    "2φ/π-1/e": 2*phi/pi - 1/np.e,
    "ln(2)+1/42": log(2)+1/42,
    "2/π+1/π²": 2/pi+1/pi**2,
    "8/(3π)": 8/(3*pi),
    "π/4-1/12": pi/4-1/12,
    "2log(φ²)/π": 2*log(phi**2)/pi,
    "4logφ/π": 4*log(phi)/pi,
    "log(φ)/log(√5)": log(phi)/log(sqrt(5)),
    "(1+8/π²)/3": (1+8/pi**2)/3,
    "8/(π²+4)": 8/(pi**2+4),
    "2·(8/π²-1/2)": 2*(8/pi**2-0.5),
    "(4/π)·logφ": 4*log(phi)/pi,
    "2/(3·log(1+√2))": 2/(3*log(1+sqrt(2))),
}

sorted_cands = sorted(candidates.items(), key=lambda x: abs(x[1]-target))
print("Closest matches:")
for name, val in sorted_cands[:10]:
    print(f"  {name:30s} = {val:.10f}  diff = {abs(val-target):.2e}")

# Let me try rational multiples of irrational constants
print("\nTrying p/q · C for small p,q and C ∈ {π, φ, e, √2, √3, √5, log2}:")
constants = {"π": pi, "φ": phi, "e": np.e, "√2": sqrt(2), "√3": sqrt(3),
             "√5": sqrt(5), "log2": log(2), "logφ": log(phi), "1/π": 1/pi}

best = []
for cname, cval in constants.items():
    for p in range(1, 20):
        for q in range(1, 20):
            val = p * cval / q
            if abs(val - target) < 0.001:
                best.append((f"{p}·{cname}/{q}", val, abs(val-target)))

best.sort(key=lambda x: x[2])
for name, val, diff in best[:15]:
    print(f"  {name:30s} = {val:.10f}  diff = {diff:.2e}")

# Also check sums
print("\nTrying a·C₁ + b·C₂ for small integers a,b:")
best2 = []
for c1name, c1 in constants.items():
    for c2name, c2 in constants.items():
        if c1name >= c2name: continue
        for a in range(-5, 6):
            for b in range(-5, 6):
                val = a*c1 + b*c2
                if abs(val - target) < 0.0001:
                    best2.append((f"{a}·{c1name}+{b}·{c2name}", val, abs(val-target)))

best2.sort(key=lambda x: x[2])
for name, val, diff in best2[:10]:
    print(f"  {name:30s} = {val:.10f}  diff = {diff:.2e}")
