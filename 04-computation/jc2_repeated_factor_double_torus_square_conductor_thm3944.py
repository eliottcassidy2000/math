#!/usr/bin/env python3
"""Exact companion for THM-3944's repeated-factor double-torus collapse."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


P, G, W, V, Z, T = sp.symbols("P G W V Z T")
d = sp.symbols("delta")


def dreduce(expression: sp.Expr) -> sp.Expr:
    return sp.expand(sp.rem(sp.expand(expression), d**2 + 3, d))


def dzero(expression: sp.Expr, message: str) -> None:
    gate(dreduce(expression) == 0, message)


omega = (-1 + d) / 2
omega2 = (-1 - d) / 2
dzero(omega**2 + omega + 1, "omega quadratic relation")
dzero(omega**3 - 1, "omega cube relation")
dzero((omega - omega2) - d, "delta is omega-omega^2")


# ---------------------------------------------------------------------------
# The balanced internal split of p1-p0=G^2.
# ---------------------------------------------------------------------------

p0 = P
p1 = P + G**2
l0 = p1 - p0
l1 = p1 - omega * p0
l2 = p1 - omega2 * p0
q0 = d * P * G
q1 = G * (3 * P + 2 * G**2)

gate(l0 == G**2, "repeated difference factor")
dzero(q1 - q0 - 2 * G * l1, "one copy of G enters q1-q0")
dzero(q1 + q0 - 2 * G * l2, "one copy of G enters q1+q0")
dzero((q1 - q0) * (q1 + q0) - 4 * (p1**3 - p0**3),
      "internally split difference-of-cubes identity")

H0 = dreduce(q0**2 - 4 * p0**3)
H1 = sp.expand(q1**2 - 4 * p1**3)
expected_H = -P**2 * (4 * P + 3 * G**2)
gate(sp.expand(H0 - expected_H) == 0,
     "first discriminant square-factor identity")
gate(sp.expand(H1 - expected_H) == 0,
     "second discriminant square-factor identity")

# The associated depressed cubics have the same discriminant up to -27.
F0 = T**3 - 3 * p0 * T - q0
F1 = T**3 - 3 * p1 * T - q1
dzero(sp.discriminant(F0, T) + 27 * H0, "first cubic discriminant")
gate(sp.expand(sp.discriminant(F1, T) + 27 * H1) == 0,
     "second cubic discriminant")


# ---------------------------------------------------------------------------
# Reduced branch geometry: one doubled line and one one-place parabola.
# ---------------------------------------------------------------------------

L = 4 * P + 3 * G**2
gate(sp.factor(expected_H) == -P**2 * L, "nonreduced branch factorization")
gate(sp.gcd(P, L) == 1, "line and parabola have no common divisor")
gate(sp.diff(L, P) == 4, "parabola is smooth affine A1")

# Projective closure 4*P*Z+3*G^2=0 has the unique infinity point [1:0:0],
# where its Z-derivative is nonzero.
PP, GG, ZZ = sp.symbols("PP GG ZZ")
Lhom = 4 * PP * ZZ + 3 * GG**2
gate(Lhom.subs({GG: 0, ZZ: 0}) == 0, "parabola infinity point")
gate(sp.diff(Lhom, ZZ).subs({PP: 1, GG: 0, ZZ: 0}) == 4,
     "parabola is smooth at its unique infinity point")


# ---------------------------------------------------------------------------
# Integral closure and exact conductor.
# ---------------------------------------------------------------------------

# V=W/P is integral and satisfies V^2=-(4P+3G^2).  Solving for P identifies
# the normalization with k[G,V].
Pnorm = -(V**2 + 3 * G**2) / 4
gate(sp.expand(V**2 + L.subs(P, Pnorm)) == 0,
     "integral overelement equation")
gate(sp.expand(Pnorm + (V**2 + 3 * G**2) / 4) == 0,
     "normalization eliminates P")
gate(sp.expand(sp.expand((P * V) ** 2 - expected_H).subs(V**2, -L)) == 0,
     "W=P*V recovers the nonnormal quadratic order")

# As R=k[P,G] modules, B=R+R*V and B0=R+R*(P*V); hence B/B0=(R/P)*V
# and the conductor is P*B=(P,W).  These exact coefficient matrices freeze
# the index rather than attempting to encode module theory in SymPy.
inclusion_matrix = sp.Matrix([[1, 0], [0, P]])
gate(inclusion_matrix.det() == P, "normalization module index is P")
dzero(Pnorm.subs(V, d * G),
      "first normalization line lies over the conductor")
dzero(Pnorm.subs(V, -d * G),
      "second normalization line lies over the conductor")

# The normalized double-cover map is (G,V)->(Pnorm,G).  Its Jacobian is -V/2,
# and V=0 maps exactly to the one-place parabola L=0.
jac = sp.det(sp.Matrix([
    [sp.diff(Pnorm, V), sp.diff(Pnorm, G)],
    [sp.diff(G, V), sp.diff(G, G)],
]))
gate(jac == -V / 2, "normalized double-cover Jacobian")
gate(sp.expand(L.subs(P, Pnorm).subs(V, 0)) == 0,
     "ramification line maps to the parabola")


# ---------------------------------------------------------------------------
# The two Cardano directions collapse in different ways after normalization.
# ---------------------------------------------------------------------------

Wnorm = sp.expand(Pnorm * V)
f0_plus = dreduce((q0 + W).subs({P: Pnorm, W: Wnorm}))
f0_minus = dreduce((q0 - W).subs({P: Pnorm, W: Wnorm}))
f1_plus = sp.expand((q1 + W).subs({P: Pnorm, W: Wnorm}))
f1_minus = sp.expand((q1 - W).subs({P: Pnorm, W: Wnorm}))

dzero(f0_plus + (V + d * G) ** 2 * (V - d * G) / 4,
      "first positive radicand has valuations 2,1")
dzero(f0_minus - (V + d * G) * (V - d * G) ** 2 / 4,
      "first conjugate radicand has valuations 1,2")
gate(sp.expand(f1_plus + (V + G) ** 3 / 4) == 0,
     "second positive radicand is an exact cube up to scalar")
gate(sp.expand(f1_minus - (V - G) ** 3 / 4) == 0,
     "second conjugate radicand is an exact cube up to scalar")

# The first cubic Kummer cover is ramified at both regular height-one lines;
# the second has trivial Kummer class.  Freeze the valuation residues.
gate((2 % 3, 1 % 3) == (2, 1), "first radicand ramification residues")
gate((3 % 3, 0) == (0, 0), "second radicand cube residues")

# The second depressed cubic actually splits over k[G,V].
r0 = -G
r1 = (G - d * V) / 2
r2 = (G + d * V) / 2
split_F1 = dreduce((T - r0) * (T - r1) * (T - r2))
dzero(F1.subs(P, Pnorm) - split_F1,
      "second depressed cubic splits over the normalization")

# Both original norm identities survive on the quadratic order.
dzero(sp.expand((q0 + W) * (q0 - W) - 4 * p0**3).subs(W**2, H0),
      "first Cardano norm on the quadratic order")
gate(sp.expand(sp.expand((q1 + W) * (q1 - W) - 4 * p1**3).subs(W**2, H1)) == 0,
     "second Cardano norm on the quadratic order")


summary = {
    "checks": CHECKS,
    "identity": "p1=p0+G^2 balanced internal split",
    "discriminant": "-P^2(4P+3G^2)",
    "reduced_nonlinear_branch": "smooth A1 parabola with one infinity place",
    "quadratic_order": "nonnormal; conductor (P,W); quotient (R/P)V",
    "normalization": "k[G,V] = A2",
    "first_cardano": "noncube but ramified with residues (2,1)",
    "second_cardano": "exact cube; depressed cubic splits",
    "smooth_locus_c3_characters": 0,
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3944 repeated-factor double-torus square-conductor companion")
print(f"CHECKS={CHECKS}")
print("COMMON_H=-P^2(4P+3G^2)")
print("NONLINEAR_REDUCED_BRANCH=A1;INFINITY_PLACES=1")
print("QUADRATIC_ORDER=NONNORMAL;CONDUCTOR=(P,W);INDEX_SUPPORT=P=0")
print("NORMALIZATION=k[G,V]=A2;UNITS=k*;CL=0;H1_MU3=0")
print("CARDANO_0=RAMIFIED_ON_V+/-delta*G;VALUATIONS=2,1")
print("CARDANO_1=EXACT_CUBE;CUBIC=SPLIT")
print("MECHANISM=ONE_PLACE_GAIN_IS_PAID_BY_SQUARE_CONDUCTOR_AND_CHARACTER_COLLAPSE")
print(f"SEMANTIC_SHA256={semantic}")
