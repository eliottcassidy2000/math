#!/usr/bin/env python3
"""sheet_defect_sharpness_kps_S128c113.py -- kind-pasteur-2026-07-20-S128c113

IS "NO LOSS OF SHEETS AT INFINITY" SHARP?

The classical statement -- a constant-Jacobian polynomial local biholomorphism with
NO LOSS OF SHEETS AT INFINITY (e.g. a proper Keller map) is an automorphism -- is
rung 2 of the five-rung lattice in THM-1375:

    injective => PROPER => nodal Jelonek => Galois => non-self-normalising stabiliser

Each arrow weakens the hypothesis, so each rung is a strictly stronger theorem.
This script asks whether rung 2 can be improved ALONG ITS OWN AXIS:

    (Q_k)   Keller + at most k sheets lost at infinity  =>  automorphism.

(Q_0) is the classical statement.  If (Q_1) held too, the rung would improve for free.

FIRST ATTEMPT FAILED, AND THE FAILURE IS INSTRUCTIVE.  I counted fibre points as
standard monomials of a Groebner basis, reading leading monomials off Poly.monoms()[0]
without pinning the monomial order and counting inside a truncated box.  The control
(the known triple collision) returned 3 and looked fine, while GENERIC targets returned
"40" -- the box cap.  A control that passes while generic input fails is the signature
of a broken instrument, not of an interesting phenomenon, so that run is discarded.

THE CORRECTED METHOD is elimination by hand, which needs no monomial-order bookkeeping.
F3 = 2x - 3x^2 y - x^3 z, so for x != 0

    z = (2x - 3x^2 y - c) / x^3

Substituting into F1 - a and F2 - b and clearing denominators gives two polynomials in
(x,y); the resultant in y is a univariate polynomial in x whose roots are the
x-coordinates of the fibre.  The x = 0 branch is handled separately and exactly:
x = 0 forces F3 = 0, and then F1 = z + 4y^2, F2 = y, giving exactly one preimage when
c = 0 and none otherwise.

The instrument is validated on the known triple collision AND on generic targets before
any drop is reported.
"""
import random
from sympy import symbols, Rational, expand, resultant, Poly, together, fraction, simplify

random.seed(5)
x, y, z = symbols('x y z')
a, b, c = symbols('a b c')
uu = 1 + x * y
F1 = uu**3 * z + y**2 * uu * (4 + 3 * x * y)
F2 = y + 3 * x * uu**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x**2 * y - x**3 * z


def fibre_points(q):
    """Exact fibre of F over q = (a,b,c).  Returns a sorted list of x-coordinates
    (with the x=0 branch counted separately), or None on degeneracy."""
    qa, qb, qc = q
    zsub = (2 * x - 3 * x**2 * y - qc) / x**3
    e1 = together(expand(F1.subs(z, zsub) - qa))
    e2 = together(expand(F2.subs(z, zsub) - qb))
    n1 = fraction(e1)[0]
    n2 = fraction(e2)[0]
    R = resultant(Poly(n1, y), Poly(n2, y))
    Rp = Poly(R, x)
    # strip the spurious x^k factor introduced by clearing x^3 denominators
    coeffs = Rp.all_coeffs()
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    Rp = Poly(coeffs, x)
    roots = Rp.all_roots() if Rp.degree() > 0 else []
    xs = [r for r in roots if r != 0]
    nx0 = 1 if qc == 0 else 0          # the x = 0 branch
    return len(set(xs)) + nx0, Rp.degree()


print("=" * 78)
print("INSTRUMENT VALIDATION (both halves must pass before anything is reported)")
print("=" * 78)
q0 = (Rational(-1, 4), 0, 0)
n0, d0 = fibre_points(q0)
print("  CONTROL, the known triple collision:")
print("     q = %-22s fibre = %s   resultant degree = %d   (expect 3: %s)"
      % (str(q0), n0, d0, n0 == 3))
print("  GENERIC targets (these must ALSO come out 3 -- the first instrument did not):")
gen = []
for t in range(5):
    q = tuple(Rational(random.randint(-5, 5), random.randint(1, 3)) for _ in range(3))
    n, d = fibre_points(q)
    gen.append(n)
    print("     q = %-34s fibre = %-3s resultant degree = %d" % (str(q), n, d))
valid = (n0 == 3) and all(g == 3 for g in gen)
print("  instrument validated : %s" % valid)
if not valid:
    print("  -> NOT validated.  No sharpness claim is made from this run.")

print()
print("=" * 78)
print("THE NON-PROPERNESS LOCUS")
print("=" * 78)
print("  THM-1315 records the caustic as proportional to a^2 * K, so a = 0 is the")
print("  cheapest probe; all three target hyperplanes are tried.")
drops = []
if valid:
    for idx, name in ((0, "a = F1 = 0"), (1, "b = F2 = 0"), (2, "c = F3 = 0")):
        print("  --- %s ---" % name)
        for t in range(5):
            q = [Rational(random.randint(-5, 5), random.randint(1, 3)) for _ in range(3)]
            q[idx] = 0
            n, d = fibre_points(tuple(q))
            if n < 3:
                drops.append((tuple(q), n))
            print("     q = %-36s fibre = %-3s deg = %-3d%s"
                  % (str(tuple(q)), n, d, "   <-- SHEET LOST" if n < 3 else ""))

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
if not valid:
    print("  Instrument not validated; nothing is concluded.")
elif drops:
    q, n = drops[0]
    print("  Fibre size %d < 3 at q = %s, so the defect there is %d." % (n, q, 3 - n))
    mind = min(3 - n for _, n in drops)
    print("  Minimum defect observed over the probes : %d" % mind)
    print()
    print("  CONSEQUENCE.  (Q_%d) -- 'Keller + at most %d sheet(s) lost at infinity =>"
          % (mind, mind))
    print("  automorphism' -- is FALSE, witnessed by this map.  So the classical rung")
    print("  cannot be weakened to allow that many lost sheets; 'no loss' is sharp at")
    print("  that threshold.")
else:
    print("  No fibre drop found on the hyperplanes sampled.  The non-properness locus")
    print("  is a proper subvariety and random rational probes on a hyperplane need not")
    print("  meet it, so this is INCONCLUSIVE, not evidence of properness -- F is known")
    print("  non-proper (THM-1330: the Jelonek set is always non-empty).  Locating the")
    print("  drop needs the caustic K explicitly, not random sampling.")
print()
print("  Independent of the outcome above, the structural point stands: rung 2 is")
print("  stated along the DEFECT axis, and the lattice improves it along three OTHER")
print("  axes -- Jelonek singularity type, Galois-ness, and the point-stabiliser")
print("  condition.  A hypothesis can be unimprovable in its own direction and still")
print("  sit near the bottom of the lattice.")
