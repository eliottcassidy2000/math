#!/usr/bin/env python3
"""
beyond_natural_base_s20bd.py -- kind-pasteur-2026-03-22-S20bd

BEYOND NATURAL BASES: What happens when b is not a natural number?

The fiber fraction at base b is f_{1/b}(k) = (1/b)_k / k!
with generating function (1-x)^{-1/b}.

The "generalized pi" for base b is Gamma(1/b)^b.
  b=2: Gamma(1/2)^2 = pi
  b=3: Gamma(1/3)^3 = 19.23
  b=infinity: Gamma(0)^inf = ???

What about non-integer b? Fractional? Irrational? Complex?

THE FUNKIEST QUESTION: What is the "tournament" at b = golden ratio?
At b = e? At b = pi? At b = 2+i?

Author: kind-pasteur-2026-03-22-S20bd
"""
import sys
import numpy as np
from math import gamma, lgamma, pi, sqrt, e, log, factorial
import cmath
sys.stdout.reconfigure(line_buffering=True)

def generalized_pi(b):
    """Gamma(1/b)^b -- the generalized pi for base b."""
    if b <= 0: return float('inf')
    try:
        return gamma(1/b)**b
    except:
        return float('inf')

def fiber_fraction(b, k):
    """(1/b)_k / k! -- the fiber fraction at base b, index k."""
    result = 1.0
    for j in range(k):
        result *= (1/b + j)
    result /= factorial(k)
    return result

print("=" * 70)
print("  BEYOND NATURAL BASES: THE GAMMA(1/b)^b LANDSCAPE")
print("=" * 70)

# ================================================================
# 1. THE LANDSCAPE OF GAMMA(1/b)^b
# ================================================================
print(f"\n{'='*70}")
print(f"  1. GAMMA(1/b)^b FOR ALL POSITIVE REAL b")
print(f"{'='*70}\n")

print(f"  {'b':>8s} {'Gamma(1/b)':>12s} {'Gamma(1/b)^b':>14s} {'log':>8s} {'interpretation':>25s}")

special_b = [
    (0.5, "sub-binary"),
    (1.0, "unary (trivial)"),
    ((1+sqrt(5))/2, "golden ratio"),
    (2.0, "BINARY = TOURNAMENTS"),
    (e, "Euler's number"),
    (3.0, "TERNARY"),
    (pi, "pi itself"),
    (4.0, "quaternary"),
    (5.0, "quinary"),
    (7.0, "septenary"),
    (10.0, "decimal"),
    (12.0, "duodecimal"),
    (100.0, "centenary"),
    (1000.0, "millenary"),
]

for b, name in special_b:
    g = gamma(1/b)
    lg = b * log(g)
    gb_str = f"{g**b:.4f}" if lg < 50 else f"~10^{lg/log(10):.1f}"
    print(f"  {b:>8.4f} {g:>12.6f} {gb_str:>14s} {lg:>8.4f} {name:>25s}")

# ================================================================
# 2. THE MINIMUM OF GAMMA(1/b)^b
# ================================================================
print(f"\n{'='*70}")
print(f"  2. WHERE IS THE MINIMUM OF Gamma(1/b)^b?")
print(f"{'='*70}\n")

# Compute Gamma(1/b)^b for fine grid of b
bs = np.linspace(0.1, 20, 1000)
gbs = np.array([gamma(1/b)**b for b in bs])

min_idx = np.argmin(gbs)
b_min = bs[min_idx]
gb_min = gbs[min_idx]

print(f"  Minimum of Gamma(1/b)^b:")
print(f"    b_min = {b_min:.4f}")
print(f"    Gamma(1/b_min)^b_min = {gb_min:.6f}")
print(f"    For comparison: pi = {pi:.6f}")

# Refine with finer grid
bs_fine = np.linspace(b_min - 0.5, b_min + 0.5, 10000)
gbs_fine = np.array([gamma(1/b)**b for b in bs_fine])
min_idx_fine = np.argmin(gbs_fine)
b_min_fine = bs_fine[min_idx_fine]
gb_min_fine = gbs_fine[min_idx_fine]

print(f"    Refined: b_min = {b_min_fine:.6f}, value = {gb_min_fine:.8f}")

# Is this e? Or something else?
print(f"    e = {e:.6f}")
print(f"    b_min / e = {b_min_fine / e:.6f}")
print(f"    b_min is {'close to e' if abs(b_min_fine - e) < 0.1 else 'NOT close to e'}")

# ================================================================
# 3. THE LIMITS
# ================================================================
print(f"\n{'='*70}")
print(f"  3. THE LIMITS: b -> 0, b -> 1, b -> infinity")
print(f"{'='*70}\n")

# b -> 0: Gamma(1/b) -> Gamma(inf) -> inf, and inf^0 = 1? No.
# Actually Gamma(1/b) for b->0 means Gamma(inf) which diverges.
# But Gamma(1/b)^b with b->0: let x = 1/b -> inf.
# Gamma(x)^{1/x} for x -> inf: by Stirling, Gamma(x) ~ sqrt(2pi/x)(x/e)^x
# So Gamma(x)^{1/x} ~ (sqrt(2pi/x))^{1/x} * x/e -> x/e -> inf.
# So Gamma(1/b)^b -> inf as b -> 0.

# b -> 1: Gamma(1) = 1, so Gamma(1)^1 = 1.
print(f"  b -> 1: Gamma(1)^1 = {gamma(1)**1:.6f}")

# b -> inf: Gamma(1/b) -> Gamma(0+) -> +inf, but Gamma(1/b)^b ...
# 1/b -> 0, Gamma(1/b) ~ b (since Gamma(x) ~ 1/x for x -> 0+)
# So Gamma(1/b)^b ~ b^b which -> inf.
# But more precisely: Gamma(x) = (1/x - gamma_euler + O(x)) for x -> 0
# So Gamma(1/b) ~ b for large b.
# Gamma(1/b)^b ~ b^b -> inf.

for b in [1, 2, 5, 10, 50, 100]:
    g = gamma(1/b)
    ratio = g / b
    print(f"  b={b:>6d}: Gamma(1/b) = {g:.4f}, Gamma(1/b)/b = {ratio:.4f}")

print(f"\n  Gamma(1/b)/b -> 1 as b -> inf (since Gamma(x) ~ 1/x near 0)")
print(f"  So Gamma(1/b)^b ~ b^b -> infinity")

# ================================================================
# 4. NON-INTEGER b: WHAT DOES IT MEAN?
# ================================================================
print(f"\n{'='*70}")
print(f"  4. FRACTIONAL b: WHAT IS A phi-ARY TOURNAMENT?")
print(f"{'='*70}\n")

phi = (1 + sqrt(5)) / 2

print(f"""  At integer b, a b-ary tournament assigns one of b outcomes to each pair.
  At fractional b, there is no literal "b-ary tournament."
  But the FIBER FRACTION f_{{1/b}}(k) = (1/b)_k / k! is well-defined for all b > 0.

  The GF (1-x)^{{-1/b}} is an ANALYTIC FUNCTION of b.
  It interpolates smoothly between integer values.

  PHYSICAL INTERPRETATION of fractional b:
  A system with b = 1.5 "outcomes per comparison" could be:
  - A MIXED system: some pairs binary, some ternary (average 1.5)
  - A PROBABILISTIC system: each comparison has ~1.5 effective outcomes
    (e.g., 50% chance of binary, 50% chance of ternary)
  - A FRACTIONAL ENTROPY system: log2(b) bits per comparison

  For b = phi = {phi:.6f} (golden ratio):
  - Fiber fraction at k=1: f = 1/phi = {1/phi:.6f} = phi - 1 = 1/phi
  - The golden ratio fiber fraction IS 1/phi = the golden ratio complement!
  - Gamma(1/phi)^phi = {generalized_pi(phi):.6f}
  - This is the "golden pi" -- the generalized circle constant for
    a comparison system with phi effective outcomes.
""")

# Fiber fractions at b = phi
print(f"  GOLDEN RATIO FIBER FRACTIONS (b = phi = {phi:.4f}):")
for k in range(8):
    f = fiber_fraction(phi, k)
    f_binary = fiber_fraction(2, k)
    print(f"    k={k}: f_phi = {f:.6f}, f_binary = {f_binary:.6f}, ratio = {f/f_binary:.4f}")

# ================================================================
# 5. COMPLEX b: THE ANALYTIC CONTINUATION
# ================================================================
print(f"\n{'='*70}")
print(f"  5. COMPLEX b: THE FULL ANALYTIC CONTINUATION")
print(f"{'='*70}\n")

# For complex b, Gamma(1/b)^b involves complex Gamma function.
# The complex Gamma function is meromorphic on all of C.
# Gamma(1/b)^b = exp(b * log(Gamma(1/b)))

def complex_generalized_pi(b):
    """Gamma(1/b)^b for complex b, using log-Gamma."""
    z = 1/b
    # Use Stirling for large |z| or numerical for moderate
    # lgamma only works for real, use cmath for complex
    # For now: only real part matters for magnitude
    if isinstance(b, complex):
        # Approximate: |Gamma(z)|^2 = pi*z / (sin(pi*z) * Gamma(1-z)) ... too complex
        # Just compute magnitude via lgamma on real part
        try:
            # Gamma(z) for complex z: use scipy if available, else approximate
            from scipy.special import gamma as cgamma
            val = cgamma(z)**b
            return val
        except ImportError:
            return None
    return gamma(1/b)**b

# Special complex values
print(f"  Gamma(1/b)^b at special complex b:")
print(f"    b = 2:      {generalized_pi(2):.6f} (= pi)")
print(f"    b = 2+0i:   {generalized_pi(2):.6f}")

# What about b = 2i (purely imaginary)?
# 1/b = 1/(2i) = -i/2
# Gamma(-i/2) is complex.
# |Gamma(-i/2)|^2 = pi / ((-i/2) * sin(pi*(-i/2)) * ... ) -- complicated.

# For b -> infinity along the real axis: Gamma(1/b)^b -> b^b -> inf
# For b -> infinity along imaginary axis: |Gamma(1/b)^b| -> ???

# Let's just compute the generalized pi for the first few primes
print(f"\n  THE PRIME GENERALIZED PIs:")
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
    gp = generalized_pi(p)
    print(f"    b={p:>2d}: Gamma(1/{p})^{p} = {gp:>14.4f}, log = {log(gp):>8.4f}")

# ================================================================
# 6. THE SPECIAL VALUES
# ================================================================
print(f"\n{'='*70}")
print(f"  6. SPECIAL VALUES AND IDENTITIES")
print(f"{'='*70}\n")

# Gamma(1/2)^2 = pi (well-known)
# Gamma(1/3)^3 = ? Is this a known constant?
# Gamma(1/4)^4 = ? Related to the lemniscate constant?

# The lemniscate constant: omega = 2*integral_0^1 dt/sqrt(1-t^4)
# = 2*pi*sqrt(2) / Gamma(1/4)^2 / (something)
# Actually: Gamma(1/4) = sqrt(2*pi^3 / AGM(1, sqrt(2)))^{1/2}... complicated.
# But Gamma(1/4) = 3.6256... and Gamma(1/4)^4 = 172.79.

# The Chowla-Selberg formula connects Gamma at rational arguments
# to periods of elliptic curves.

g14 = gamma(0.25)
g13 = gamma(1/3)
g12 = gamma(0.5)

print(f"  Gamma(1/2) = sqrt(pi) = {g12:.6f}")
print(f"  Gamma(1/3) = {g13:.6f}")
print(f"  Gamma(1/4) = {g14:.6f}")
print()
print(f"  GENERALIZED PIs:")
print(f"    Gamma(1/2)^2 = pi = {g12**2:.6f}")
print(f"    Gamma(1/3)^3 = {g13**3:.6f}")
print(f"    Gamma(1/4)^4 = {g14**4:.6f}")
print()

# Ratio pattern
print(f"  RATIOS Gamma(1/b)^b / Gamma(1/(b-1))^(b-1):")
prev = pi
for b in range(3, 12):
    curr = gamma(1/b)**b
    ratio = curr / prev
    print(f"    b={b}: ratio = {ratio:.6f}")
    prev = curr

# ================================================================
# 7. THE b -> infinity LIMIT
# ================================================================
print(f"\n{'='*70}")
print(f"  7. THE ASYMPTOTICS: b -> infinity")
print(f"{'='*70}\n")

# For large b: Gamma(1/b) ~ b (since Gamma(x) ~ 1/x near 0)
# So log(Gamma(1/b)^b) = b * log(Gamma(1/b)) ~ b * log(b)
# And Gamma(1/b)^b ~ b^b

# More precisely: Gamma(x) = 1/x - gamma_euler + O(x) for x -> 0
# So Gamma(1/b) = b - gamma_euler/b + O(1/b^2)? No.
# Gamma(x) for x -> 0: Gamma(x) = 1/x - gamma + (gamma^2/2 + pi^2/12)*x + ...
# So Gamma(1/b) = b - gamma + O(1/b)
# Gamma(1/b)^b = (b - gamma)^b ~ b^b * (1 - gamma/b)^b ~ b^b * e^{-gamma}

euler_gamma = 0.5772156649

for b in [10, 20, 50, 100, 1000]:
    g = gamma(1/b)
    actual = g**b
    approx = b**b * np.exp(-euler_gamma)
    ratio = actual / approx
    print(f"  b={b:>5d}: Gamma(1/b)^b = {actual:.4e}, b^b*e^{{-gamma}} = {approx:.4e}, ratio = {ratio:.6f}")

print(f"""
  ASYMPTOTIC: Gamma(1/b)^b ~ b^b * e^{{-gamma}} as b -> infinity

  where gamma = {euler_gamma:.6f} is the EULER-MASCHERONI CONSTANT.

  So the generalized pi grows as b^b (super-exponential).
  The "correction factor" is e^{{-gamma}} = {np.exp(-euler_gamma):.6f}.

  THIS CONNECTS TOURNAMENT THEORY TO THE EULER-MASCHERONI CONSTANT!

  The hierarchy:
    b=2: Gamma(1/2)^2 = pi                (the circle constant)
    b=3: Gamma(1/3)^3 = 19.23             (the ternary constant)
    b->inf: Gamma(1/b)^b ~ b^b * e^(-gamma)  (Euler-Mascheroni!)

  The three great constants pi, e, gamma ALL appear in the
  fiber fraction hierarchy at different positions:
    pi at b=2 (binary = tournaments)
    e in the base of the asymptotic growth
    gamma as the correction factor at large b
""")

# ================================================================
# 8. THE REFLECTION FORMULA
# ================================================================
print(f"{'='*70}")
print(f"  8. THE REFLECTION: Gamma(1/b) * Gamma(1-1/b) = pi/sin(pi/b)")
print(f"{'='*70}\n")

# The Gamma reflection formula: Gamma(x)*Gamma(1-x) = pi/sin(pi*x)
# Set x = 1/b: Gamma(1/b)*Gamma(1-1/b) = pi/sin(pi/b)

for b_val in [2, 3, 4, 5, 6, phi, e, pi]:
    x = 1/b_val
    lhs = gamma(x) * gamma(1-x)
    rhs = pi / np.sin(pi*x)
    label = {2:'binary', 3:'ternary', 4:'quat', 5:'quin', 6:'sen'}.get(b_val, f'{b_val:.3f}')
    print(f"  b={b_val:>7.3f} ({label:>7s}): Gamma(1/b)*Gamma(1-1/b) = {lhs:.6f}, pi/sin(pi/b) = {rhs:.6f}")

print(f"""
  THE REFLECTION FORMULA says:
  Gamma(1/b) = pi / (sin(pi/b) * Gamma(1-1/b))

  For b=2: Gamma(1/2) = pi / (sin(pi/2) * Gamma(1/2)) = pi / Gamma(1/2)
  => Gamma(1/2)^2 = pi. (The classical result.)

  For general b: Gamma(1/b)^b = (pi / sin(pi/b))^b / Gamma(1-1/b)^b

  This connects the "generalized pi" Gamma(1/b)^b to sin(pi/b),
  which is a CHORD LENGTH on the unit circle at angle pi/b.

  THE GEOMETRIC MEANING:
  sin(pi/b) is the half-chord subtended by angle 2*pi/b on the unit circle.
  This is the side length of a regular b-gon inscribed in a unit circle,
  divided by 2.

  So: Gamma(1/b)^b ~ (pi/sin(pi/b))^b ~ (b/pi * pi)^b = b^b for large b
  (using sin(pi/b) ~ pi/b).

  THE BEAUTIFUL CHAIN:
  REGULAR b-GON -> chord length sin(pi/b) -> reflection formula ->
  Gamma(1/b) -> fiber fraction -> tournament theory at base b

  THE TOURNAMENT AT BASE b IS CONTROLLED BY THE REGULAR b-GON.
""")
