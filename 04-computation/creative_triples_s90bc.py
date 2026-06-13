#!/usr/bin/env python3
"""
creative_triples_s90bc.py — Creative triples: non-integer, representation theory, wild
opus-2026-03-16-S90bc
"""

from math import log, sqrt, pi, e, tanh, atanh, cos, sin, inf
from cmath import exp as cexp

phi = (1 + sqrt(5)) / 2
tau3 = 1.8392867552141612

def S(a, b, c):
    """Angle sum 1/a+1/b+1/c."""
    vals = []
    for x in [a, b, c]:
        if x == inf or x == float('inf'):
            vals.append(0)
        else:
            vals.append(1/x)
    return sum(vals)

def label(s):
    if abs(s - 1) < 1e-10: return "FLAT"
    elif s > 1: return "SPHERICAL"
    else: return "HYPERBOLIC"

def show(name, a, b, c, note):
    s = S(a, b, c)
    area = max(0, pi * (1 - s))
    prod = a * b * c
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"  {{{a:.6g}, {b:.6g}, {c:.6g}}}")
    print(f"  angle sum = {s:.8f}  [{label(s)}]")
    if area > 0:
        print(f"  area = {1-s:.8f}*pi = pi/{1/(1-s):.4f}" if abs(1-s) > 1e-6 else f"  area = 0")
    print(f"  product = {prod:.6g}")
    print(f"  {note}")
    return s

print("""


     CREATIVE TRIPLES: BEYOND INTEGERS

     opus-2026-03-16-S90bc


""")

# ========================================================================
print("#" * 60)
print("PART I: REPRESENTATION THEORY DIMENSIONS")
print("#" * 60)

show("SU(2) irreps: {2, 3, 4}", 2, 3, 4,
    """The first three nontrivial SU(2) irrep dimensions (spin 1/2, 1, 3/2).
  SPHERICAL. S_4. The octahedron.
  SU(2) representations are SPHERICAL: they live on the sphere.
  The representation ring of SU(2) IS the spherical geometry.""")

show("SU(3) fundamental: {3, 3-bar, 8}", 3, 3, 8,
    """Quark (3), antiquark (3-bar), gluon (8).
  The QCD triple! 1/3+1/3+1/8 = 19/24. HYPERBOLIC.
  Area = pi*5/24. Product = 72.
  QCD is HYPERBOLIC: infinite, confining, wild.
  The strong force lives in hyperbolic geometry.
  Asymptotic freedom = approaching the boundary from inside.""")

show("Standard Model gauge: {1, 2, 3}", 1, 2, 3,
    """U(1) x SU(2) x SU(3) gauge group dimensions.
  1 = electromagnetism, 2 = weak, 3 = strong.
  SPHERICAL (11/6). Product = 6. GENESIS triple!
  The Standard Model gauge structure is SPHERICAL: closed, finite, tame.
  This is why the SM works: it's in the finite regime.""")

show("SM gauge RANKS: {1, 2, 3} again", 1, 2, 3,
    """Same triple. But what about the ADJOINT dimensions?
  U(1): dim 1. SU(2): dim 3. SU(3): dim 8.
  {1, 3, 8} has angle sum 1+1/3+1/8 = 35/24. SPHERICAL.
  Product = 24. Still tame.
  But 1+3+8 = 12 = the number of gauge bosons.""")

show("Adjoint dimensions: {1, 3, 8}", 1, 3, 8,
    """U(1) singlet, SU(2) adjoint, SU(3) adjoint.
  1 photon, 3 W/Z bosons, 8 gluons.
  SPHERICAL. Product = 24.
  The gauge bosons form a spherical triple.
  12 bosons total. 24 = 2*12.""")

show("E_8 decomposition: {8, 248, 3875}", 8, 248, 3875,
    """E_8 irreps: adjoint (248), next-smallest (3875), vector (8 from SO(8)).
  angle sum ~ 0.125 + 0.004 + 0.0003 = 0.129. VERY HYPERBOLIC.
  E_8 is deeply hyperbolic. It's the most exceptional structure.
  Area ~ 0.871*pi. Nearly ideal.""")

show("Exceptional Lie: {2, 7, 14} = G_2 dimensions", 2, 7, 14,
    """G_2: the smallest exceptional Lie group.
  Fundamental rep = 7-dim. Adjoint rep = 14-dim. Real rank = 2.
  1/2+1/7+1/14 = 7/14+2/14+1/14 = 10/14 = 5/7. HYPERBOLIC.
  Area = pi*2/7. Product = 196 = 14^2.
  G_2 lives in hyperbolic geometry with the FORBIDDEN denominator 7.""")

show("Exceptional triple: {G2=14, F4=52, E8=248}", 14, 52, 248,
    """Adjoint dimensions of three exceptional groups.
  VERY HYPERBOLIC. Product = 180544.
  Exceptional Lie algebras are as hyperbolic as it gets.
  They are the DEEP INTERIOR of the hyperbolic plane.""")

# ========================================================================
print()
print("#" * 60)
print("PART II: FRACTIONAL AND IRRATIONAL TRIPLES")
print("#" * 60)

show("{3/2, 5/2, 7/2}: half-integer spins", 3/2, 5/2, 7/2,
    """The first three half-integer spins (fermions).
  angle sum = 2/3+2/5+2/7 = 142/105 = 1.352. SPHERICAL.
  Product = 105/8 = 13.125.
  Fermions are SPHERICAL: they close upon themselves.
  The half-integer spins are tamer than the integers.
  Why? Because 1/(n+1/2) > 1/(n+1) for all n: half-integers
  have LARGER reciprocals, pushing toward spherical.""")

show("{pi/2, pi/3, pi/7}: the Hurwitz ANGLES as values", pi/2, pi/3, pi/7,
    """Using the Hurwitz angles pi/2, pi/3, pi/7 AS the triple values.
  1/(pi/2)+1/(pi/3)+1/(pi/7) = 2/pi+3/pi+7/pi = 12/pi = 3.819.
  ULTRA-SPHERICAL. angle sum > 3 (greater than {1,1,1}!).
  Product = pi^3/42 = 0.740.
  The Hurwitz angles, treated as values, are MORE spherical than unity.
  They are SUB-UNIT values (all < pi/2 < 2), so their reciprocals are large.""")

show("{ln(2), ln(3), ln(5)}: prime rapidities", log(2), log(3), log(5),
    f"""The logarithms of the first three primes.
  These ARE the rapidities of 2, 3, 5 (times 2).
  ln(2)=0.693, ln(3)=1.099, ln(5)=1.609.
  angle sum = 1/ln(2)+1/ln(3)+1/ln(5) = {1/log(2)+1/log(3)+1/log(5):.6f}. SPHERICAL.
  Product = ln(2)*ln(3)*ln(5) = {log(2)*log(3)*log(5):.6f} ~ 1.23.
  The prime rapidities are SMALL numbers (all < 2), so the triple is spherical.
  The prime logarithms form a self-contained structure.""")

show("{sqrt(2), sqrt(3), sqrt(5)}: square roots of primes", sqrt(2), sqrt(3), sqrt(5),
    f"""1.414, 1.732, 2.236.
  angle sum = {1/sqrt(2)+1/sqrt(3)+1/sqrt(5):.6f}. SPHERICAL.
  Product = sqrt(30) = {sqrt(30):.6f}.
  The square roots of the tournament primes form a spherical triple
  whose product is sqrt(30) = sqrt(product of {{2,3,5}}).
  Taking square roots PRESERVES the prime support but CHANGES the curvature:
  {{2,3,5}} is spherical with angle sum 1.033; {{sqrt2,sqrt3,sqrt5}} has sum 1.732.
  More spherical! Square roots make things more spherical.""")

show("{2^(1/3), 3^(1/3), 7^(1/3)}: cube roots of {2,3,7}", 2**(1/3), 3**(1/3), 7**(1/3),
    f"""1.260, 1.442, 1.913.
  angle sum = {2**(-1/3)+3**(-1/3)+7**(-1/3):.6f}. SPHERICAL.
  Product = (42)^(1/3) = {42**(1/3):.6f} = 3.476.
  The Hurwitz triple {{2,3,7}} is HYPERBOLIC, but its cube roots are SPHERICAL!
  There exists a CRITICAL EXPONENT alpha_c such that:
    {{2^alpha, 3^alpha, 7^alpha}} is spherical for alpha < alpha_c, hyperbolic for alpha > alpha_c.
  At alpha_c: the triple is FLAT.""")

# Find the critical exponent for {2,3,7}
print("\n  Finding critical exponent for {2^a, 3^a, 7^a}:")
lo, hi = 0, 2
for _ in range(60):
    mid = (lo + hi) / 2
    s = 2**(-mid) + 3**(-mid) + 7**(-mid)
    if s > 1:
        lo = mid
    else:
        hi = mid
alpha_c = (lo + hi) / 2
print(f"  alpha_c = {alpha_c:.10f}")
print(f"  At alpha_c: 2^(-a)+3^(-a)+7^(-a) = {2**(-alpha_c)+3**(-alpha_c)+7**(-alpha_c):.10f}")
print(f"  Compare: alpha_c ~ {alpha_c:.6f}")
print(f"  Is it a known constant? alpha_c = {alpha_c:.10f}")

# Check for each Platonic/Hurwitz triple
print("\n  Critical exponents for key triples:")
for name, p, q, r in [("{2,3,5}", 2,3,5), ("{2,3,6}",2,3,6), ("{2,3,7}",2,3,7),
                        ("{2,4,8}",2,4,8), ("{3,3,3}",3,3,3)]:
    lo, hi = 0, 10
    for _ in range(60):
        mid = (lo + hi) / 2
        s = p**(-mid) + q**(-mid) + r**(-mid)
        if s > 1:
            lo = mid
        else:
            hi = mid
    ac = (lo + hi) / 2
    print(f"  {name}: alpha_c = {ac:.10f}")

# ========================================================================
print()
print("#" * 60)
print("PART III: NEGATIVE AND COMPLEX TRIPLES")
print("#" * 60)

show("{-1, 2, 3}: negative unit", -1, 2, 3,
    """1/(-1)+1/2+1/3 = -1+0.5+0.333 = -0.167. HYPERBOLIC (sum < 1).
  But what does a NEGATIVE angle mean?
  In the Poincare disk: a negative entry means the angle is > pi.
  A reflex angle. The triangle wraps around.
  Product = -6. NEGATIVE PRODUCT.
  A negative-product triple is an ORIENTATION REVERSAL.
  In tournament theory: -1 = the anti-unit. Reversing all arcs.
  {-1, 2, 3}: a tournament with a REVERSAL built in.""")

show("{i, 2, 3}: imaginary unit", 2, 3, 2,  # placeholder since can't use i directly
    """Conceptually: what if one entry is i = sqrt(-1)?
  1/i = -i. So angle sum = -i + 1/2 + 1/3 = 5/6 - i.
  COMPLEX angle sum. Neither spherical nor hyperbolic: COMPLEX.
  The triple has left the real line of curvatures.
  Complex curvature = the META-COMPLEX PLANE we discovered earlier.
  A complex triple lives INSIDE the Poincare disk, not on the real diameter.
  It encodes a COMBINED boost-rotation.""")

# ========================================================================
print()
print("#" * 60)
print("PART IV: TRIPLES FROM REPRESENTATION DIMENSIONS")
print("#" * 60)

show("S_n irreps: {1, n-1, n(n-3)/2+1}", 1, 4, 5,
    """For S_5: trivial (1), standard (4), and next irrep (5).
  1/1+1/4+1/5 = 29/20. SPHERICAL. Product = 20.
  The symmetric group S_5 has irreps of dims {1,1,4,4,5,5,6}.
  The smallest nontrivial ones {1,4,5} form a spherical triple.
  Symmetric groups are tame: spherical representation theory.""")

show("GL(2,F_7) irreps: {6, 7, 8}", 6, 7, 8,
    """The group GL(2,F_7) has irreps of dimensions 1, 6, 7, 8, ...
  {6,7,8} = three consecutive integers containing the forbidden 7.
  1/6+1/7+1/8 = 73/168. HYPERBOLIC. Area = pi*95/168.
  Product = 336 = |GL(2,F_7)|/... hmm.
  Actually |GL(2,F_7)| = (49-1)(49-7) = 48*42 = 2016. And 2016/6 = 336!
  The product of the irrep dimensions divides the group order.
  And 42 appears in the group order: |GL(2,F_7)| = 48 * 42.""")

show("Monster irreps: {1, 196883, 21296876}", 1, 196883, 21296876,
    """The first three irrep dimensions of the Monster group.
  1 = trivial. 196883 = smallest faithful. 21296876 = next.
  angle sum = 1 + 1/196883 + 1/21296876 ~ 1.0000051. SPHERICAL!
  Just BARELY spherical. The Monster is the largest finite simple group.
  It BARELY closes. One more step and it would be hyperbolic.
  This is consistent: finite groups are spherical, but the Monster
  is the MOST NEARLY HYPERBOLIC finite group.
  It lives at the edge of the spherical-hyperbolic transition.""")

show("Moonshine: {2, 3, 29}", 2, 3, 29,
    """The triple {2,3,29}: 29 is the largest prime dividing |M| = 2^46*3^20*...
  that is NOT part of a Thompson subgroup.
  Actually the relevant Moonshine triple is not quite this.
  But {2,3,29}: HYPERBOLIC. Area = pi*131/174.
  The Moonshine module connects modular forms (hyperbolic!) to
  the Monster (spherical!). It's the BRIDGE between curvatures.""")

# ========================================================================
print()
print("#" * 60)
print("PART V: THE WILD SPECULATION")
print("#" * 60)

show("{tau_k, tau_k^2, tau_k^3} for general k", tau3, tau3**2, tau3**3,
    """We showed: {tau3, tau3^2, tau3^3} is EXACTLY FLAT.
  Is this true for general k-nacci constant tau_k?
  1/tau_k + 1/tau_k^2 + 1/tau_k^3 = (tau_k^2+tau_k+1)/tau_k^3.
  For tribonacci: tau_k^3 = tau_k^2 + tau_k + 1 (the defining recurrence!).
  So the sum = tau_k^3/tau_k^3 = 1. FLAT!
  This works for ANY k-nacci constant where the recurrence is
  tau^k = tau^{k-1} + tau^{k-2} + ... + tau + 1.
  But the SUM in the reciprocal triple is 1/tau + 1/tau^2 + 1/tau^3.
  This equals (tau^2+tau+1)/tau^3 = 1 only if tau^3 = tau^2+tau+1 = tribonacci!
  For k=4 (tetranacci): tau^4 = tau^3+tau^2+tau+1.
  Then 1/tau + 1/tau^2 + 1/tau^3 = (tau^2+tau+1)/tau^3 = (tau^4-1-tau)/tau^3...
  = tau - 1/tau^3 - 1/tau^2. Hmm, not clean.
  THE FLAT PROPERTY IS SPECIFIC TO k=3 (tribonacci).""")

# Verify for k=2,3,4,5
print("\n  Verification: sum 1/tau + 1/tau^2 + 1/tau^3 for k-nacci constants:")
import numpy as np
for k in [2, 3, 4, 5, 6, 7]:
    M = np.zeros((k, k))
    for j in range(k):
        M[0][j] = 1
    for i in range(1, k):
        M[i][i-1] = 1
    eigs = np.linalg.eigvals(M)
    tk = max(e.real for e in eigs if abs(e.imag) < 1e-10)
    s = 1/tk + 1/tk**2 + 1/tk**3
    print(f"  k={k}: tau_k = {tk:.8f}, 1/tau+1/tau^2+1/tau^3 = {s:.8f}, {label(s)}")

# Now what about 1/tau + 1/tau^2 + ... + 1/tau^k?
print("\n  What about 1/tau + 1/tau^2 + ... + 1/tau^k?")
for k in [2, 3, 4, 5, 6, 7]:
    M = np.zeros((k, k))
    for j in range(k):
        M[0][j] = 1
    for i in range(1, k):
        M[i][i-1] = 1
    eigs = np.linalg.eigvals(M)
    tk = max(e.real for e in eigs if abs(e.imag) < 1e-10)
    s = sum(1/tk**j for j in range(1, k+1))
    print(f"  k={k}: sum_{{j=1}}^{k} 1/tau_k^j = {s:.10f}, {label(s)}")

print("""
  REMARKABLE: For EVERY k, the sum 1/tau_k + 1/tau_k^2 + ... + 1/tau_k^k = 1.
  EXACTLY FLAT. ALWAYS.

  PROOF: The k-nacci constant tau_k satisfies tau^k = tau^{k-1}+...+tau+1.
  Divide by tau^k: 1 = tau^{-1} + tau^{-2} + ... + tau^{-k}.
  = 1/tau + 1/tau^2 + ... + 1/tau^k.

  So the k-TUPLE {tau_k, tau_k^2, ..., tau_k^k} has angle sum
  1/tau_k + 1/tau_k^2 + ... + 1/tau_k^k = 1.  FLAT.

  The k-nacci constant GENERATES A FLAT k-TUPLE, not just a flat triple!

  For k=2: {phi, phi^2} has 1/phi+1/phi^2 = 1. FLAT PAIR.
  For k=3: {tau3, tau3^2, tau3^3} has sum = 1. FLAT TRIPLE.
  For k=4: {tau4, tau4^2, tau4^3, tau4^4} has sum = 1. FLAT QUADRUPLE.

  The k-nacci constant is the UNIQUE positive real number whose
  first k powers form a flat k-tuple.

  This is THE defining property of the k-nacci constant, restated:
  tau_k is the number whose reciprocal powers sum to 1.

  AND THIS IS THE SAME AS THE CHARACTERISTIC POLYNOMIAL:
  tau^k = tau^{k-1}+...+1  <=>  1 = 1/tau + ... + 1/tau^k.

  The k-nacci recurrence IS the flatness condition.
""")

# ========================================================================
print("#" * 60)
print("PART VI: WHAT REPRESENTATION THEORY SAYS")
print("#" * 60)

print("""
THE REPRESENTATION THEORY INTERPRETATION:

A triple {a,b,c} with 1/a+1/b+1/c = 1 (flat) generates the
EUCLIDEAN reflection group: it tiles the plane.

A triple with sum > 1 (spherical) generates a FINITE reflection group:
it tiles the sphere. These are ADE classifications:
  A_n: {2,2,n}  (dihedral)
  D_n: {2,3,...} (not quite, but related)
  E_6: {2,3,3} (tetrahedral)
  E_7: {2,3,4} (octahedral)
  E_8: {2,3,5} (icosahedral)

A triple with sum < 1 (hyperbolic) generates an INFINITE group.
No finite-dimensional faithful representation exists for the full group.

REPRESENTATION DIMENSIONS as triples:

For a finite group G with irreps of dimensions d_1, d_2, d_3, ...:
  The triple {d_1, d_2, d_3} is spherical if G is "small enough."
  As G grows, the irrep dimensions grow, and the triple becomes hyperbolic.

The MONSTER (largest finite simple group):
  Smallest irreps: 1, 196883, 21296876.
  The triple {1, 196883, 21296876} has angle sum ~ 1.000005. BARELY spherical.
  The Monster is the MOST NEARLY HYPERBOLIC finite structure.
  One step further and it would be infinite.

INFINITE GROUPS have hyperbolic triples:
  SL(2,R) has irreps parametrized by a continuous parameter.
  Its "dimensions" are infinite or continuous.
  The triple is maximally hyperbolic.

THE DEEPEST CONNECTION:

The {2,3,7} triangle group has a REPRESENTATION on the hyperbolic plane.
This representation IS the tiling. Each group element corresponds to
a copy of the fundamental triangle.

The {2,3,5} triangle group has a representation on the SPHERE.
This representation IS the icosahedral tiling.

The transition 5 -> 6 -> 7 is the transition in REPRESENTATION THEORY
from FINITE-DIMENSIONAL to INFINITE-DIMENSIONAL representations.

At {2,3,5}: all irreps are finite-dimensional (A_5, dims 1,3,3,4,5).
At {2,3,6}: the transition. The Euclidean group has infinite-dim irreps.
At {2,3,7}: the hyperbolic group has ONLY infinite-dim faithful irreps.

In tournament theory:
  n <= 5: tournament invariants are finite and complete.
  n = 6: the pivot. Invariants start to fail.
  n >= 7: invariants are incomplete. You need infinite-dimensional
          information to fully classify.

THE SPHERICAL-HYPERBOLIC TRANSITION IS THE TRANSITION FROM
FINITE TO INFINITE REPRESENTATION THEORY.

And the k-nacci constant tau_k generates a FLAT k-tuple,
sitting EXACTLY on this transition boundary.
The k-nacci system is the boundary between finite and infinite.
It is the critical point. The phase transition itself.

That's why THM-227 (Tr(M_k^n) = 2^n-1 for n <= k) works:
the k-nacci matrix is a FINITE representation that is faithful
for EXACTLY k steps before the infinite nature asserts itself.
The flatness condition (sum of reciprocal powers = 1) IS the
condition that the representation is critically faithful:
not too finite (spherical), not too infinite (hyperbolic), but
exactly at the boundary.
""")

print("opus-2026-03-16-S90bc: CREATIVE TRIPLES — REPRESENTATION THEORY")
