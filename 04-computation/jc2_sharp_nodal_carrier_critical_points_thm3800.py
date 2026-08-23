#!/usr/bin/env python3
"""Exact companion for THM-3800's sharp-carrier critical-point correction."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


r, z, e = sp.symbols("r z e")
c = sp.symbols("c", nonzero=True)
variables = (r, z, e)
surface = r**2 * e - z**3 + c**3 * r
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 * c**3 + 6 * r * e],
        [-9 * z**2, -3 * c**3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


A = e**2 - z / (3 * c**3) + 4 * z**3 * e**3 / c**6
K = c**3 + 2 * r * e
A_z = sp.factor(sp.diff(A, z))
A_e = sp.factor(sp.diff(A, e))
critical = [sp.factor(bracket(A, q)) for q in variables]
zero(critical[0] - (-3 * r**2 * A_z - 9 * z**2 * A_e),
     "r Hamiltonian")
zero(critical[1] + 3 * K * A_e, "z Hamiltonian")
zero(critical[2] - 3 * K * A_z, "e Hamiltonian")
zero(bracket(A, surface), "surface Casimir")


# The K=0 torus was designed to be escaped.  Its restriction is linear, so
# the remaining Hamiltonian is nonzero there.
e_torus = -c**6 / (4 * z**3)
r_torus = 2 * z**3 / c**3
A_torus = sp.factor(A.subs(e, e_torus))
zero(A_torus + z / (3 * c**3), "sharp carrier torus restriction")
zero(sp.diff(A_torus, z) + 1 / (3 * c**3),
     "sharp carrier torus derivative")
zero(surface.subs({r: r_torus, e: e_torus}), "torus lies on surface")
zero(K.subs({r: r_torus, e: e_torus}), "torus lies on K=0")
zero(critical[0].subs({r: r_torus, e: e_torus}) - r_torus**2 / c**3,
     "K=0 remaining Hamiltonian nonzero")


# On K!=0, criticality forces A_z=A_e=0.  The two equations are exactly the
# triangular z/e laws below.
z_law = -6 * c**3 * e**2
e_law = 1296 * c**3 * e**7 - 1
zero(
    A_z.subs(z, z_law) - e_law / (3 * c**3),
    "A_z equals e-law multiple",
)
zero(
    A_e.subs(z, z_law) + 2 * e * e_law,
    "A_e equals e-law multiple",
)


# Explicit two branches over each of the seven e-roots.  Work with an
# algebraic symbol w satisfying w^2=3 and reduce every numerator modulo it.
w = sp.symbols("w", nonzero=True)
eta = sp.symbols("eta", nonzero=True)


def reduce_branch(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    numerator = sp.rem(sp.Poly(numerator, w), sp.Poly(w**2 - 3, w)).as_expr()
    numerator = sp.rem(
        sp.Poly(sp.expand(numerator), eta),
        sp.Poly(1296 * c**3 * eta**7 - 1, eta),
    ).as_expr()
    return sp.factor(numerator)


for sign in (-1, 1):
    r0 = c**3 * (-1 + sign / w) / (2 * eta)
    z0 = -6 * c**3 * eta**2
    substitutions = {r: r0, z: z0, e: eta}
    gate(reduce_branch(surface.subs(substitutions)) == 0,
         f"branch {sign} lies on surface")
    for index, component in enumerate(critical):
        gate(reduce_branch(component.subs(substitutions)) == 0,
             f"branch {sign} Hamiltonian {index}")
    K0 = sp.factor(K.subs(substitutions))
    zero(K0 - sign * c**3 / w, f"branch {sign} K value")


# General-c Groebner basis over Q(c).  Its monic leading coefficients make
# specialization safe for every c!=0; c=1 is replayed separately as the
# hostile control that exposed the failed mate-support premise.
domain = sp.QQ.frac_field(c)
groebner_general = sp.groebner(
    [surface] + critical, r, z, e, order="lex", domain=domain
)
expected_general = [
    r**2 + 1296 * c**6 * e**6 * r + 216 * c**9 * e**5,
    z + 6 * c**3 * e**2,
    e**7 - 1 / (1296 * c**3),
]
gate(len(groebner_general.polys) == 3, "general Groebner length")
for index, expected in enumerate(expected_general):
    zero(groebner_general.polys[index].as_expr() - expected,
         f"general Groebner row {index}")

groebner_c1 = sp.groebner(
    [surface.subs(c, 1)] + [component.subs(c, 1) for component in critical],
    r, z, e, order="lex", domain=sp.QQ,
)
expected_c1 = [
    r**2 + 1296 * e**6 * r + 216 * e**5,
    z + 6 * e**2,
    e**7 - sp.Rational(1, 1296),
]
gate(len(groebner_c1.polys) == 3, "c=1 Groebner length")
for index, expected in enumerate(expected_c1):
    zero(groebner_c1.polys[index].as_expr() - expected,
         f"c=1 Groebner row {index}")


# Reducedness and exact count.  The e-polynomial has seven simple nonzero
# roots; z is linear; the r-quadratic has nonzero discriminant on every root.
e_polynomial = 1296 * c**3 * e**7 - 1
gate(sp.gcd(sp.Poly(e_polynomial, e), sp.Poly(sp.diff(e_polynomial, e), e)) == 1,
     "seven e-roots are simple")
r_quadratic = expected_general[0]
r_discriminant = sp.factor(sp.discriminant(r_quadratic, r))
field_c = sp.QQ.frac_field(c)
discriminant_remainder = sp.rem(
    sp.Poly(r_discriminant, e, domain=field_c),
    sp.Poly(e_law, e, domain=field_c),
).as_expr()
zero(discriminant_remainder - 432 * c**9 * e**5,
     "r discriminant on e-law")
gate(
    sp.gcd(
        sp.Poly(discriminant_remainder, e, domain=field_c),
        sp.Poly(e_law, e, domain=field_c),
    ) == 1,
    "r discriminant is nonzero on all e-roots",
)
gate(432 != 0, "r branches remain distinct in characteristic zero")

derivative_minor = sp.factor(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(A_z, z), sp.diff(A_z, e)],
                [sp.diff(A_e, z), sp.diff(A_e, e)],
            ]
        )
    ).subs(z, z_law)
)
reduced_minor = sp.rem(
    sp.Poly(sp.together(derivative_minor).as_numer_denom()[0], e,
            domain=field_c),
    sp.Poly(e_law, e, domain=field_c),
).as_expr()
zero(reduced_minor + 1008 * e**5, "z/e transverse minor on e-law")
gate(
    sp.gcd(
        sp.Poly(reduced_minor, e, domain=field_c),
        sp.Poly(e_law, e, domain=field_c),
    ) == 1,
    "z/e intersection is transverse at all seven roots",
)


semantic = {
    "correction": "reserved four-weight mate lane fails because carrier is critical",
    "carrier": "A=e^2-z/(3c^3)+(4/c^6)z^3e^3",
    "critical_points": "1296c^3e^7=1; z=-6c^3e^2; r=c^3(-1+/-1/sqrt(3))/(2e)",
    "count": "14 distinct reduced points for every c!=0",
    "K_boundary": "K=0 escaped exactly; critical points have K=+/-c^3/sqrt(3)",
    "consequence": "no regular Darboux mate; polynomial support enumeration superseded",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3800-sharp-torus-escaping-nodal-carrier-has-fourteen-critical-points")
print("surface=r^2e-z^3+c^3r;field=algebraically_closed_characteristic_zero;c!=0")
print("carrier=A=e^2-z/(3c^3)+(4/c^6)z^3e^3")
print("failed_assumption=torus_escape_does_not_imply_smooth_carrier")
print("K0_control=A|K=-z/(3c^3);Hamiltonian_r=r^2/c^3")
print("critical_law=1296c^3e^7=1;z=-6c^3e^2;r=c^3(-1+/-1/sqrt(3))/(2e)")
print("critical_count=14;scheme=reduced")
print("K_values=+/-c^3/sqrt(3)")
print("c1_groebner=r^2+1296e^6r+216e^5;z+6e^2;e^7-1/1296")
print("consequence=no_regular_Darboux_mate;reserved_four_weight_support_lane_superseded")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
