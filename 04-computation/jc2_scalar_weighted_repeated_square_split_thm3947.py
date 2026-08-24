#!/usr/bin/env python3
"""Exact companion for THM-3947's repeated-square split trichotomy.

Reproduction:
  python3 04-computation/jc2_scalar_weighted_repeated_square_split_thm3947.py
  python3 -O 04-computation/jc2_scalar_weighted_repeated_square_split_thm3947.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


P, G, alpha, r = sp.symbols("P G alpha r", nonzero=True)
z = sp.symbols("z")
omega = (-1 + sp.sqrt(-3))/2
delta = sp.expand(omega-omega**2)

gate(sp.simplify(omega**2+omega+1) == 0, "cube-root relation")
gate(sp.simplify(delta**2+3) == 0, "delta square")

p0 = P
p1 = P+alpha*G**2
L0 = sp.expand(p1-p0)
L1 = sp.expand(p1-omega*p0)
L2 = sp.expand(p1-omega**2*p0)

# Split the two copies of G in L0=alpha G^2 between the complementary
# difference-of-cubes packets.  The nonzero scalar r is the only weight.
D = sp.expand(r*G*L1)
S = sp.expand(alpha*G*L2/r)
q0 = sp.expand(S-D)
q1 = sp.expand(S+D)
H = sp.expand(q0**2-4*p0**3)

gate(sp.expand(L0*L1*L2-(p1**3-p0**3)) == 0,
     "difference-of-cubes packet")
gate(sp.simplify(D*S-(p1**3-p0**3)) == 0,
     "scalar-weighted internal split")
gate(sp.simplify(H-(q1**2-4*p1**3)) == 0,
     "common double-torus discriminant")

A = sp.simplify(alpha*(1-omega**2)/r-r*(1-omega))
B = sp.simplify(alpha*(alpha/r-r))
gate(sp.simplify(q0-G*(A*P+B*G**2)) == 0,
     "two scalar moduli A and B control q0")
gate(sp.simplify(H-(G**2*(A*P+B*G**2)**2-4*P**3)) == 0,
     "weighted-homogeneous discriminant")

F = sp.expand((A*z+B)**2-4*z**3)
gate(sp.simplify(H-G**6*F.subs(z, P/G**2)) == 0,
     "three-parabola dehomogenization")
disc_F = sp.factor(sp.discriminant(F, z))
gate(sp.simplify(disc_F+16*B**3*(A**3+27*B)) == 0,
     "cubic root discriminant")
gate(sp.simplify(r*B-alpha*(alpha-r**2)) == 0,
     "first collision seam is r^2=alpha")
gate(sp.simplify(A**3+27*B
                 - 3*delta*(r**2+omega*alpha)**3/r**3) == 0,
     "second collision seam is r^2=-omega alpha")
gate(sp.simplify(1+omega) != 0,
     "the two collision seams are disjoint")

# First collision: the doubled component is the line p0=P=0, and the
# remaining component is the THM-3944 parabola.
H_first = sp.factor(H.subs(alpha, r**2))
first_expected = -P**2*(4*P+3*r**2*G**2)
gate(sp.expand(H_first-first_expected) == 0,
     "first seam is double p0 plus one parabola")
gate(sp.simplify(A.subs(alpha, r**2)**2+3*r**2) == 0,
     "first seam nonzero simple root")

# Second collision: the doubled component is p1=0.  It is the index-swapped
# copy of the first conductor geometry.
alpha_second = -omega**2*r**2
p1_second = sp.expand(P+alpha_second*G**2)
H_second = sp.factor(H.subs(alpha, alpha_second))
second_expected = sp.expand(
    -4*p1_second**2*(P+alpha_second*G**2/4)
)
gate(sp.simplify(H_second-second_expected) == 0,
     "second seam is double p1 plus one parabola")
gate(sp.simplify(B.subs(alpha, alpha_second)) != 0,
     "second seam has no zero cubic root")

# On the second discriminant seam B=-A^3/27, the double and simple roots are
# A^2/9 and A^2/36.  This also proves there is no hidden triple-root case.
u = sp.symbols("u")
universal_second = sp.factor((u-sp.Rational(1, 27))**2-4*u**3)
gate(sp.expand(
    universal_second
    + (9*u-1)**2*(36*u-1)/sp.Integer(729)
) == 0, "universal double-root factorization")
gate(sp.Rational(1, 9) != sp.Rational(1, 36),
     "second seam has two distinct components")
gate(sp.simplify(A.subs(alpha, r**2)) != 0,
     "B=0 never produces a triple component")

# Every nonzero root zeta gives a smooth one-place parabola P=zeta G^2.
# Its projective conic P Z-zeta G^2 has determinant zeta/4 and exactly one
# point at infinity, where d/dZ=P is nonzero.  The zero root is the affine
# line P=0 and is also one-place.
croot = sp.symbols("croot", nonzero=True)
conic_matrix = sp.Matrix(((0, 0, sp.Rational(1, 2)),
                          (0, -croot, 0),
                          (sp.Rational(1, 2), 0, 0)))
gate(sp.simplify(conic_matrix.det()-croot/4) == 0,
     "nonzero parabola closure is a smooth conic")
Z = sp.symbols("Z")
projective_parabola = P*Z-croot*G**2
gate(projective_parabola.subs({Z: 0, G: 0, P: 1}) == 0,
     "parabola has its unique infinity point")
gate(sp.diff(projective_parabola, Z).subs({Z: 0, G: 0, P: 1}) == 1,
     "the unique infinity point is smooth")

summary = {
    "checks": CHECKS,
    "family": "p1-p0=alpha G^2;one G assigned to each complementary factor",
    "generic": "three distinct smooth one-place parabolas;full H reduced reducible",
    "seam_r2_eq_alpha": "double p0 line plus THM3944 parabola",
    "seam_r2_eq_minus_omega_alpha": "double p1 parabola plus simple parabola",
    "triple_component": "impossible for alpha*r nonzero",
    "scope": "factorization geometry only;generic normal order and Cardano lattice open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3947 scalar-weighted repeated-square split exact companion")
print(f"CHECKS={CHECKS}")
print("GENERIC=THREE_DISTINCT_SMOOTH_ONE_PLACE_PARABOLAS;FULL_H=REDUCED_REDUCIBLE")
print("SEAM_R2_EQ_ALPHA=DOUBLE_P0_LINE_PLUS_PARABOLA")
print("SEAM_R2_EQ_MINUS_OMEGA_ALPHA=DOUBLE_P1_PARABOLA_PLUS_PARABOLA")
print("TRIPLE_COMPONENT=NONE_FOR_NONZERO_PARAMETERS")
print("CONDUCTOR_SEAMS=INDEX_SWAPPED_THM3944_GEOMETRY")
print("OPEN=GENERIC_NORMAL_ORDER_AND_CARDANO_LATTICE")
print(f"SEMANTIC_SHA256={semantic}")
