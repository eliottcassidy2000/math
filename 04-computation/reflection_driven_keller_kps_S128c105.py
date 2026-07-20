#!/usr/bin/env python3
"""reflection_driven_keller_kps_S128c105.py -- kind-pasteur-2026-07-20-S128c105

Two things, one of which kills my own first idea.

=====================================================================
(A) SELF-REFUTATION: the weight TWIST is a coordinate artifact.
=====================================================================
THM-1300 records F as C*-equivariant with weights (1,-1,-2) on the source and
(-2,-1,1) on the target -- the same MULTISET, related by the transposition (1 3).
I proposed the reduced conjecture

    (RJC-untwisted)  F Keller and C*-equivariant with the SAME weight VECTOR on
                     source and target  ==>  F invertible,

on the reading that the counterexample needs a "reflection twist".  THAT IS FALSE.
Composing with the target coordinate swap sigma(a,b,c) = (c,b,a) gives

    G := sigma . F,   det JG = -det JF = 2  (still Keller),
    same fibres as F  (still non-injective),
    target weights of G = (1,-1,-2) = source weights  (UNTWISTED).

So an untwisted equivariant Keller counterexample exists immediately.  The weight
twist is not invariant under composing with a linear automorphism, while both
"Keller" and "non-injective" are -- so it can never be the operative hypothesis.
Verified exactly below.

This is the same lesson as my LRC (1,2,3) episode: an apparent symmetry breaking
that is really an artifact of the coordinate ORDERING convention.

=====================================================================
(B) WHAT SURVIVES: the monodromy reflection, which IS coordinate-free.
=====================================================================
The Galois group of C(x)/C(F) is intrinsic.  For cover degree 3 the generic fibre
is a 3-element set carrying a cyclic structure -- THM-1310's "cyclic 3-tournament".
Its symmetry group is D_3 = S_3, split as in THM-127:

    A_3 = the 3 rotations  = tournament AUTOMORPHISMS   (arc-preserving)
    the 3 transpositions   = ANTI-automorphisms          (arc-reversing)

Campbell (1973): a Keller map that is a GALOIS cover is an automorphism.  At d = 3,
Galois <=> monodromy contained in A_3 <=> discriminant is a square.  Hence

    THEOREM (reflection-driven).  Every degree-3 Keller counterexample has
    monodromy exactly S_3, i.e. its fibre discriminant is a NON-SQUARE:
    the monodromy MUST contain an arc-reversing element of the fibre tournament.

That upgrades canon's status for the discriminant law from VERIFIED (0/8316 sample
points, salvage catalog 2.5) to PROVED -- it is forced by Campbell, not evidence
for a pattern.  And d = 3 is the ONLY degree where "Galois" is a single quadratic
character: for d >= 4 a transitive subgroup of A_d need not be regular.  With
Smith's d = 2 impossibility (klein-S325), d = 3 is simultaneously the minimal
counterexample degree and the unique character-detectable one.

CONCRETE REALISATION, checked below: lambda = -1 in F's own torus FIXES the
collision image (-1/4,0,0) and TRANSPOSES two of its three preimages, fixing the
third -- an explicit reflection of the fibre 3-tournament.

=====================================================================
(C) INSTRUMENT GATE (MISTAKE-196).
=====================================================================
Any search over untwisted equivariant maps must be able to REDISCOVER G = sigma.F.
This script reports the minimal ansatz degree at which G is expressible, so that a
negative search result below that degree is correctly read as vacuous.
"""
from sympy import symbols, expand, diff, Matrix, Poly, Rational, factor, simplify

x, y, z = symbols('x y z')

u = 1 + x * y
F1 = u**3 * z + y**2 * u * (4 + 3 * x * y)
F2 = y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x**2 * y - x**3 * z
F = [F1, F2, F3]


def wts(expr, w=(1, -1, -2)):
    ex = {sum(e * wi for e, wi in zip(m, w)) for m, _ in Poly(expand(expr), x, y, z).terms()}
    return ex.pop() if len(ex) == 1 else sorted(ex)


def detJ(G):
    return expand(Matrix(3, 3, lambda i, j: diff(G[i], [x, y, z][j])).det())


print("=" * 78)
print("(A) SELF-REFUTATION -- the weight twist is a coordinate artifact")
print("=" * 78)
print("  source weights (x,y,z)      = (1, -1, -2)")
print("  target weights (F1,F2,F3)   = (%s, %s, %s)" % tuple(wts(f) for f in F))
print("  det JF = %s" % detJ(F))
G = [F3, F2, F1]                      # sigma . F with sigma = swap of coords 1,3
print()
print("  G := sigma . F  with sigma(a,b,c) = (c,b,a):")
print("  target weights (G1,G2,G3)   = (%s, %s, %s)" % tuple(wts(g) for g in G))
print("  det JG = %s   -> constant, non-zero: %s"
      % (detJ(G), detJ(G) != 0 and detJ(G).free_symbols == set()))
pts = [(0, 0, Rational(-1, 4)), (1, Rational(-3, 2), Rational(13, 2)),
       (-1, Rational(3, 2), Rational(13, 2))]
imsG = [tuple(expand(g.subs({x: p[0], y: p[1], z: p[2]})) for g in G) for p in pts]
print("  G on the three collision points: %s" % (imsG[0],))
print("  all three still coincide : %s" % (len(set(imsG)) == 1))
print()
print("  VERDICT: G is Keller, non-injective, and UNTWISTED (target vector = source")
print("  vector).  So (RJC-untwisted) is FALSE.  The twist is not invariant under")
print("  composing with a linear automorphism, while Keller and non-injective are.")
print("  My proposed hypothesis is refuted before it was ever claimed.")

print()
print("=" * 78)
print("(B) THE INTRINSIC REFLECTION -- lambda = -1 acts on the fibre by a")
print("    TRANSPOSITION of the 3-cycle tournament")
print("=" * 78)
lam = -1
w_src = (1, -1, -2)
w_tgt = (-2, -1, 1)
q = (Rational(-1, 4), 0, 0)
qa = tuple(Rational(lam)**w * qi for w, qi in zip(w_tgt, q))
print("  target point q = %s" % (q,))
print("  lambda = -1 acting on q (target weights %s) : %s" % (w_tgt, qa))
print("  q is FIXED by lambda = -1 : %s" % (qa == q))
print()
print("  action of lambda = -1 on the three preimages (source weights %s):" % (w_src,))
act = []
for p in pts:
    pa = tuple(Rational(lam)**w * pi for w, pi in zip(w_src, p))
    act.append(pa)
    where = pts.index(pa) if pa in pts else None
    print("     %-26s -> %-26s  (preimage #%s)"
          % (str(p), str(pa), where if where is not None else "OFF-FIBRE"))
perm = [pts.index(a) if a in pts else None for a in act]
print("  induced permutation of the fibre : %s" % (perm,))
fixed = [i for i, v in enumerate(perm) if v == i]
moved = [i for i, v in enumerate(perm) if v != i]
print("  fixed points: %s ; moved: %s" % (fixed, moved))
is_transp = len(fixed) == 1 and len(moved) == 2 and perm[moved[0]] == moved[1] \
    and perm[moved[1]] == moved[0]
print("  is a TRANSPOSITION (odd, so NOT in A_3) : %s" % is_transp)
print()
print("  READING.  The generic fibre is a 3-cycle tournament (THM-1310).  Its")
print("  symmetry group is D_3 = S_3 = A_3 (rotations = tournament automorphisms)")
print("  union the 3 transpositions (ANTI-automorphisms, arc-reversing) -- exactly")
print("  THM-127's rotations-preserve / reflections-reverse split.  Campbell (1973)")
print("  says a GALOIS Keller cover is an automorphism, and at d = 3 Galois means")
print("  monodromy inside A_3.  So a counterexample MUST carry a reflection, and")
print("  lambda = -1 supplies one explicitly.")

print()
print("=" * 78)
print("(C) INSTRUMENT GATE -- minimal ansatz degree that can express G = sigma.F")
print("=" * 78)
# untwisted weight-(1,-1,-2) parametrisation: G1 = x A(s,t), G2 = y B + xz C,
# G3 = z D + y^2 E, with s = xy, t = x^2 z.
s, t = symbols('s t')


def in_st(expr, gen):
    """Write expr / gen as a polynomial in s = xy, t = x^2 z; return its degree."""
    q = simplify(expand(expr) / gen)
    q = expand(q)
    deg = 0
    for m, _ in Poly(q, x, y, z).terms():
        a, b, c = m
        # x^a y^b z^c invariant: a = b + 2c ; equals s^b t^c
        if a - b - 2 * c != 0:
            return None
        deg = max(deg, b + c)
    return deg


dA = in_st(G[0], x)
print("  G1 = x * A(s,t)          -> deg A = %s" % dA)
# G3 = z D + y^2 E : split off the z-part
G3 = expand(G[2])
zpart = expand(G3.coeff(z, 1) * z)
rest = expand(G3 - zpart)
dD = in_st(zpart, z)
dE = in_st(rest, y**2)
print("  G3 = z * D(s,t) + y^2 E  -> deg D = %s, deg E = %s" % (dD, dE))
need = max(v for v in (dA, dD, dE) if v is not None)
print()
print("  MINIMAL ansatz degree needed to express G : %d" % need)
print("  => any search over untwisted equivariant maps run at ansatz degree < %d"
      % need)
print("     CANNOT rediscover the known witness G, so a negative result there has")
print("     ZERO evidential weight (MISTAKE-196, the instrument-gate rule).")
print("     My first run used degree 1.  It was vacuous, and I am recording that")
print("     rather than reporting its 'no non-linear solutions' line as a finding.")
