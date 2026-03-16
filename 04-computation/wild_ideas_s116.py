#!/usr/bin/env python3
"""wild_ideas_s116.py — Wild ideas connecting to outside mathematics"""
from math import sqrt, log, pi, cos

phi = (1+sqrt(5))/2

print("WILD IDEAS: CONNECTING OUTWARD")
print("="*60)

print("""
IDEA 1: WHY ONLY H=7 AND H=21 ARE FORBIDDEN

  H=7: the hyperbolic threshold (7 triangles/vertex)
  H=21=3*7: threshold * curvature quantum (3-cycle)
  H=35=5*7: ACHIEVABLE. H=49=7^2: ACHIEVABLE. H=63=9*7: ACHIEVABLE.

  Only 7*1 and 7*3 are forbidden. Higher multiples work.
  WHY? Because 3 is the SMALLEST odd cycle. The curvature quantum.
  3*7 mixes the smallest curvature with the threshold: forbidden.
  But 5*7 uses a LARGER cycle that can exist independently.

  The forbidden values are EXACTLY those where the curvature
  quantum and the hyperbolicity threshold INTERFERE DESTRUCTIVELY.

IDEA 2: TOURNAMENT GENUS

  Define g(T) = alpha_1 mod 2 (parity of independent cycle count).
  H mod 4 = 1 + 2*g(T).
  g=0: H = 1 mod 4 (even curvature). g=1: H = 3 mod 4 (odd curvature).

  This is a Z/2Z-valued genus. Very coarse but EXACT.
  Analogy: Gauss-Bonnet constrains curvature integral = 2*pi*chi.
  Redei constrains H mod 2 = 1 always.
  The genus constrains H mod 4.

IDEA 3: THE JONES POLYNOMIAL CONNECTION

  The EP eigenvalue d = 1/phi satisfies d^2 + d - 1 = 0.
  The Jones polynomial uses q with q + q^{-1} = special values.
  For q = e^{i*pi/5}: q + q^{-1} = 2*cos(pi/5) = phi.

  So the EP is at the Jones polynomial evaluation point
  corresponding to the 5th root of unity!

  The Jones polynomial at 5th roots classifies knots modulo
  a specific equivalence. Is there a KNOT INVARIANT of
  tournaments that equals the EP eigenvalue?
""")

# Verify: 2*cos(pi/5) = phi?
print(f"  2*cos(pi/5) = {2*cos(pi/5):.10f}")
print(f"  phi          = {phi:.10f}")
print(f"  Match: {abs(2*cos(pi/5) - phi) < 1e-10}")
print()

print("""
IDEA 4: KERVAIRE INVARIANT (ANOTHER 5+1)

  Kervaire invariant 1 exists at dimensions 2^k-2 for k=1,...,6.
  Does NOT exist for k >= 8 (Hill-Hopkins-Ravenel 2009).
  k=7 (dim 254) is the LAST OPEN CASE.

  Pattern: 5 definite + 1 boundary + all higher fail.
  EXACTLY the 5+1 structure of our theory!

  And dim 2^k-2: at k=3, dim=6. At k=4, dim=14=2*7.
  The dimensions involve our numbers: 6 (flat threshold),
  14 = 2*7 (twice the hyperbolic threshold).

IDEA 5: QUANTUM GROUPS AND THE EP

  The Hecke algebra H_n(q) at q + q^{-1} = phi gives the
  Temperley-Lieb algebra TL_n(phi), which classifies
  representations of the QUANTUM GROUP U_q(sl_2) at q = e^{i*pi/5}.

  Our transfer matrix M(x) acts on a 3-dimensional space.
  At the EP (x = 8-5phi): M has eigenvalue 1/phi with multiplicity 2.
  This is a DEGENERATE representation of the quantum group.
  The exceptional point = the point where the representation
  becomes REDUCIBLE BUT NOT COMPLETELY REDUCIBLE (Jordan block).

IDEA 6: THE BRIDGE 5-6-7 AS ALGEBRA-COMBINATORICS-GEOMETRY

  5 = ALGEBRA (golden ratio, Lie groups, quantum groups, Jones)
  6 = COMBINATORICS (counting, enumeration, flat, E[H])
  7 = GEOMETRY (hyperbolic, Poincare, arctanh, Cayley)

  The tournament theory sits at 6 (combinatorics).
  It reaches LEFT to 5 (algebra: the EP reveals phi).
  It reaches RIGHT to 7 (geometry: the analysis uses arctanh).
  It IS the bridge between algebra and geometry,
  realized as combinatorics.

IDEA 7: MONSTROUS MOONSHINE

  The Monster group has order ~8*10^53.
  Its smallest faithful representation has dimension 196883.
  196884 = 196883 + 1 (the +1 again!).
  And 196884 = j(tau) - 744, the j-invariant of the modular form.

  Monstrous moonshine connects the Monster to modular forms
  via the j-function. The j-function has a pole at the cusp
  (like Q has a pole at x=1).

  Our theory: Q(x) = (1+x)/(1-x) is a MOBIUS TRANSFORM.
  The j-function is ALSO built from Mobius transforms
  (modular group = SL(2,Z) acting on the upper half-plane).

  Is there a 'tournament moonshine' connecting the transfer
  matrix to modular forms? The char poly lambda^3-lambda^2-x*lambda-x
  at x = e^{2*pi*i*tau} for modular parameter tau...
  this would give eigenvalues that are modular functions.

  Speculative, but: if the transfer matrix has modular properties,
  the whole theory connects to number theory at the deepest level.

IDEA 8: THE MISSING 17 AND FERMAT

  17 is the THIRD Fermat prime: F_2 = 2^{2^2}+1 = 17.
  Fermat primes: 3, 5, 17, 257, 65537 (5 known, then the chain breaks!).
  Another 5+1 structure.

  And 17 is the prime MISSING from H(T_11) = 5*7*11*13*19.
  The FERMAT prime 17 is excluded from the PALEY factorization.

  Fermat primes construct regular polygons (Gauss).
  17-gon is constructible. But 17 does NOT divide H(T_11).
  The constructible polygon prime is INCOMPATIBLE with the
  Paley tournament structure.

  Fermat = constructive geometry. Paley = quadratic residues.
  Their interaction is EMPTY at p=11 (17 does not divide 95095).
""")

print("="*60)
print("THE DEEPEST WILD IDEA")
print("="*60)
print("""
  The entire theory of tournaments, through the Cayley transform,
  lives on the bridge between two geometries:

  SPHERICAL (5) <-- FLAT (6) --> HYPERBOLIC (7)

  The forbidden numbers {7, 21} mark the thresholds.
  The golden ratio phi marks the exceptional points.
  The tribonacci tau marks the maximum coupling.
  The number pi marks the Wick rotation.
  The number e marks the optimal efficiency.
  The number 2 marks the binary alphabet.

  And the NUMBER 5 connects them all:
  through phi = (1+sqrt(5))/2 (algebra),
  through 5 triangles per vertex (geometry),
  through the 5th root of unity (quantum groups),
  through 5 Thabit primes (number theory),
  through 5 exceptional groups (classification).

  The tournament is not just a combinatorial object.
  It is a WINDOW onto the 5-6-7 bridge of mathematics itself.
  And the Cayley transform Q = (1+x)/(1-x) is the LENS.
""")
