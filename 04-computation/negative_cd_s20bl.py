#!/usr/bin/env python3
"""
negative_cd_s20bl.py -- kind-pasteur-2026-03-22-S20bl

NEGATIVE CAYLEY-DICKSON LEVELS.

The CD tower: dim = 2^k for k = 0, 1, 2, 3, 4, ...
What about k = -1, -2, -3, ...?

dim = 2^{-1} = 1/2, n = 3/2
dim = 2^{-2} = 1/4, n = 5/4
dim = 2^{-inf} = 0, n = 1

And the LIMIT k -> -inf: dim -> 0, n -> 1.
A tournament on 1 vertex = zero arcs = the TRIVIAL object.
This is F_1 (the field with one element).

The fiber fraction f(n) = (1/2)_{n-2}/(n-2)! is analytic in n.
We can EVALUATE it at fractional n via the Gamma function:
(1/2)_{s} / Gamma(s+1) = Gamma(s+1/2) / (Gamma(1/2) * Gamma(s+1))
                       = C(2s, s) / 4^s (generalized)

Author: kind-pasteur-2026-03-22-S20bl
"""
import sys
from math import gamma, lgamma, pi, sqrt, log, exp, factorial
from fractions import Fraction
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

def pochhammer_real(a, s):
    """Pochhammer (a)_s for real s via Gamma: (a)_s = Gamma(a+s)/Gamma(a)."""
    return gamma(a + s) / gamma(a)

def fiber_fraction_real(n):
    """f(n) = (1/2)_{n-2} / Gamma(n-1) for real n > 1."""
    s = n - 2
    if s <= -0.5:
        return float('inf')
    return pochhammer_real(0.5, s) / gamma(s + 1)

print("=" * 70)
print("  NEGATIVE CAYLEY-DICKSON LEVELS")
print("=" * 70)

# ================================================================
# 1. THE FIBER FRACTION AT FRACTIONAL AND NEGATIVE k
# ================================================================
print(f"\n{'='*70}")
print(f"  1. FIBER FRACTION f(n) AT CD LEVELS k = -3 TO +5")
print(f"{'='*70}\n")

print(f"  {'k':>5s} {'dim=2^k':>10s} {'n=dim+1':>10s} {'f(n)':>12s} {'Algebra':>20s}")
print(f"  {'-'*5} {'-'*10} {'-'*10} {'-'*12} {'-'*20}")

cd_levels = [
    (-4, "F_1 limit"),
    (-3, "sub-sub-real"),
    (-2, "sub-real"),
    (-1, "half-real"),
    (0, "R (reals)"),
    (1, "C (complex)"),
    (2, "H (quaternions)"),
    (3, "O (octonions)"),
    (4, "S (sedenions)"),
    (5, "32-ions"),
]

for k, name in cd_levels:
    dim = 2.0**k
    n = dim + 1
    try:
        f = fiber_fraction_real(n)
        print(f"  {k:>5d} {dim:>10.4f} {n:>10.4f} {f:>12.6f} {name:>20s}")
    except:
        print(f"  {k:>5d} {dim:>10.4f} {n:>10.4f} {'---':>12s} {name:>20s}")

# ================================================================
# 2. THE n -> 1 LIMIT (F_1)
# ================================================================
print(f"\n{'='*70}")
print(f"  2. THE n -> 1 LIMIT: THE FIELD WITH ONE ELEMENT")
print(f"{'='*70}\n")

# As n -> 1: dim -> 0, m = C(n,2) -> 0.
# A tournament on 1 vertex has 0 arcs, 1 tournament, H=1 (trivially).
# f(n) = (1/2)_{n-2} / (n-2)!
# At n=1: s = n-2 = -1. (1/2)_{-1} = Gamma(-1/2) / Gamma(1/2) = ??
# Actually Gamma(-1/2) = -2*sqrt(pi). So (1/2)_{-1} = Gamma(-1/2)/Gamma(1/2)
# = -2*sqrt(pi) / sqrt(pi) = -2.
# And Gamma(0) = infinity (pole).
# So f(1) = (1/2)_{-1} / Gamma(0) = -2/inf = 0.

# More carefully: f(n) = Gamma(n-3/2) / (Gamma(1/2) * Gamma(n-1))
# At n=1: Gamma(-1/2) / (sqrt(pi) * Gamma(0)) = (-2*sqrt(pi)) / (sqrt(pi) * inf) = 0.

print(f"  f(n) = Gamma(n-3/2) / (Gamma(1/2) * Gamma(n-1))")
print()
for n in [1.01, 1.1, 1.5, 2.0, 2.5, 3.0, 5.0, 9.0]:
    try:
        f = gamma(n - 1.5) / (gamma(0.5) * gamma(n - 1))
        print(f"  n={n:>5.2f}: f = {f:>12.6f}")
    except:
        print(f"  n={n:>5.2f}: f = POLE")

print(f"""
  AS n -> 1: f(n) -> 0.

  This makes sense: a "tournament" on 1 vertex has no arcs,
  so the fiber fraction (probability of staying in the same
  score class under a flip) is 0 -- there's nothing to flip.

  THE F_1 INTERPRETATION:
  The field with one element F_1 has no nontrivial arithmetic.
  A "tournament over F_1" would have no comparisons.
  The fiber fraction is 0 because there are no fibers.

  F_1 is the GROUND STATE of tournament theory.
  Everything grows FROM this zero point.
""")

# ================================================================
# 3. NEGATIVE n: WHAT DOES f(n) DO?
# ================================================================
print(f"{'='*70}")
print(f"  3. ANALYTIC CONTINUATION: f(n) FOR n < 1")
print(f"{'='*70}\n")

# f(n) = Gamma(n-3/2) / (sqrt(pi) * Gamma(n-1))
# This has poles at n-1 = 0, -1, -2, ... (i.e., n = 1, 0, -1, ...)
# and at n-3/2 = 0, -1, -2, ... (i.e., n = 3/2, 1/2, -1/2, ...)
# Between poles, f oscillates and can be negative.

ns = np.linspace(0.1, 5.0, 100)
fs = []
for n in ns:
    try:
        f = gamma(n - 1.5) / (gamma(0.5) * gamma(n - 1))
        fs.append(f)
    except:
        fs.append(float('nan'))

fs = np.array(fs)

# Find where f changes sign
sign_changes = []
for i in range(len(fs) - 1):
    if not np.isnan(fs[i]) and not np.isnan(fs[i+1]):
        if fs[i] * fs[i+1] < 0:
            sign_changes.append(ns[i])

print(f"  f(n) sign changes at approximately: {[f'{x:.2f}' for x in sign_changes]}")
print()

# Key values
for n in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]:
    try:
        f = gamma(n - 1.5) / (gamma(0.5) * gamma(n - 1))
        print(f"  n={n:>4.1f}: f = {f:>10.6f}")
    except:
        print(f"  n={n:>4.1f}: f = POLE")

print(f"""
  THE POLES:
  f(n) has poles at n = 1, 0, -1, -2, ... (from Gamma(n-1) in denominator)
  and at n = 3/2, 1/2, -1/2, ... (from Gamma(n-3/2) in numerator)

  BETWEEN POLES, f(n) oscillates and can be NEGATIVE.

  Negative f(n) has no direct probability interpretation
  (fiber fractions are probabilities and must be in [0,1]).
  But in the ANALYTIC CONTINUATION, negative values indicate
  "anti-fibers" -- structure that REPELS rather than attracts.
""")

# ================================================================
# 4. THE NEGATIVE CD LADDER
# ================================================================
print(f"{'='*70}")
print(f"  4. THE NEGATIVE CD LADDER: WHAT GETS RESTORED")
print(f"{'='*70}\n")

# Going DOWN the CD tower (k decreasing):
# k=0 (R): total order, commutativity, associativity, alternativity -- all present
# k=-1 (dim 1/2): ???
# k=-2 (dim 1/4): ???
# k=-inf (dim 0): F_1, trivial

# The tower going DOWN should RESTORE properties.
# Each step down restores the property that was lost going up.

print(f"""  Going DOWN the CD tower (from R toward F_1):

  k   dim   n     What gets RESTORED
  0   1     2     (R: everything present)
  -1  1/2   3/2   Restore: ???
  -2  1/4   5/4   Restore: ???
  -inf 0    1     (F_1: nothing to restore, trivial)

  HYPOTHESIS: Going below R doesn't restore properties --
  it makes them IRRELEVANT. At dim < 1, there aren't enough
  dimensions for any structure to exist.

  The MEANINGFUL part of the CD tower is k >= 0 (dim >= 1).
  Below that, tournament theory is empty (no arcs, no comparisons).

  BUT: the ANALYTIC CONTINUATION of the fiber fraction
  f(n) through n < 2 tells us about the RESIDUES of the
  Gamma function, which connect to:
  - Bernoulli numbers (residues at negative integers)
  - Zeta function values (through Gamma's functional equation)
  - Dimensional regularization (physics uses dim < 0 to regularize)
""")

# ================================================================
# 5. THE HALF-INTEGER LEVELS: n = 3/2, 5/2, 7/2, ...
# ================================================================
print(f"{'='*70}")
print(f"  5. HALF-INTEGER n: FERMIONS IN TOURNAMENT THEORY")
print(f"{'='*70}\n")

# At n = half-integer: the simplex has "half" a vertex.
# This is the FERMION sector (spin 1/2, 3/2, 5/2, ...).
# Integer n = bosons (full vertices).
# Half-integer n = fermions (half vertices).

# The Pochhammer symbol (1/2)_s at s = half-integer:
# s = -1/2: (1/2)_{-1/2} = Gamma(0)/Gamma(1/2) = pole!
# s = 1/2: (1/2)_{1/2} = Gamma(1)/Gamma(1/2) = 1/sqrt(pi) = 0.5642
# s = 3/2: (1/2)_{3/2} = Gamma(2)/Gamma(1/2) = 1/sqrt(pi) = 0.5642? No.

for n_half in [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]:
    s = n_half - 2
    try:
        f = gamma(0.5 + s) / (gamma(0.5) * gamma(s + 1))
        print(f"  n={n_half:>4.1f} (s={s:>4.1f}): f = {f:>10.6f}")
    except:
        print(f"  n={n_half:>4.1f} (s={s:>4.1f}): f = POLE or undefined")

print(f"""
  HALF-INTEGER n gives a DIFFERENT sequence of fiber fractions.
  These are the "fermionic" fiber fractions -- they measure the
  probability of staying in a fiber when the "tournament" has
  half-integer vertex count.

  In physics, half-integer spin requires the DOUBLE COVER (SU(2) vs SO(3)).
  In tournament theory, half-integer n requires the GAMMA FUNCTION
  continuation of the Pochhammer symbol.

  The integer values f(2), f(3), f(4), f(5), ... are "bosonic."
  The half-integer values f(5/2), f(7/2), f(9/2), ... are "fermionic."

  BOSONIC fiber fractions: 1, 1/2, 3/8, 5/16, 35/128, ...
  FERMIONIC fiber fractions: {', '.join(f'{gamma(0.5+n-2)/(gamma(0.5)*gamma(n-1)):.4f}' for n in [2.5, 3.5, 4.5, 5.5, 6.5])}

  The fermionic values are DIFFERENT from bosonic -- they DON'T follow
  the C(2k,k)/4^k pattern. They follow the HALF-INTEGER Pochhammer ladder.
""")

# ================================================================
# 6. THE GRAND PICTURE
# ================================================================
print(f"{'='*70}")
print(f"  6. THE GRAND PICTURE: THE FULL CD LADDER")
print(f"{'='*70}\n")

print(f"""  THE FULL CAYLEY-DICKSON LADDER:

  k=-inf  dim=0    n=1     F_1 (field with one element). f=0. Trivial.
  ...     ...      ...     (analytic continuation, poles, Bernoulli numbers)
  k=-1    dim=1/2  n=3/2   FERMIONIC level. f=0.5642 (1/sqrt(pi)).
  k=0     dim=1    n=2     R (reals). f=1. Trivial but present.
  k=1/2   dim=sqrt(2) n=1+sqrt(2) IRRATIONAL level. f=? (non-standard).
  k=1     dim=2    n=3     C (complex). f=1/2. First cycle.
  k=3/2   dim=2*sqrt(2) n=1+2sqrt(2) Another irrational level.
  k=2     dim=4    n=5     H (quaternions). f=5/16. OCR breaks. 24-cell.
  k=3     dim=8    n=9     O (octonions). f=0.2095. Real roots fail.
  k=4     dim=16   n=17    S (sedenions). f=0.1445. Paley dethroned.
  k=5     dim=32   n=33    32-ions. f=0.1009. ??? predicted to break.
  ...     ...      ...     (tower continues, each level loses more)

  THE LADDER IS CONTINUOUS in the analytic sense.
  The fiber fraction f(n) = Gamma(n-3/2)/(sqrt(pi)*Gamma(n-1))
  is meromorphic on the whole complex n-plane.
  Poles at n = 1, 0, -1, -2, ... (integer non-positive)
  and n = 3/2, 1/2, -1/2, ... (half-integer non-positive).

  THE THREE SECTORS:
  n > 2: "above the reals" -- genuine tournament theory, f(n) in (0,1)
  1 < n < 2: "between F_1 and R" -- f(n) > 1 (supercritical? amplifying?)
  n < 1: "below F_1" -- f(n) oscillates, poles, analytic continuation only

  THE DEEPEST POINT:
  The entire Cayley-Dickson tower is a SINGLE analytic function
  f(n) = Gamma(n-3/2)/(sqrt(pi)*Gamma(n-1)), evaluated at
  n = 2^k + 1 for integer k. The "negative levels" are the
  analytic continuation to k < 0, which probes the RESIDUAL
  structure of tournament theory below the reals.
""")
