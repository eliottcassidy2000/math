#!/usr/bin/env python3
"""Exact companion for THM-3790's cubic pseudo-plane arm gate."""

from __future__ import annotations

import ast
from pathlib import Path

import sympy as sp


r, z, e, c, t, s, zeta = sp.symbols("r z e c t s zeta", nonzero=True)
variables = (r, z, e)
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 * c**3 + 6 * r * e],
        [-9 * z**2, -3 * c**3 - 6 * r * e, 0],
    ]
)

CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


def zero(expr: sp.Expr, message: str) -> None:
    gate(sp.cancel(expr) == 0, message)


# The first infinitesimal arm is exactly C[e,z]/(z^2): over Q(c), the
# surface relation plus (r,z)^2 has Groebner basis (r,z^2).
domain = sp.QQ.frac_field(c)
surface = r**2 * e - z**3 + c**3 * r
arm_square = sp.groebner(
    [r**2, r * z, z**2, surface], r, z, e, order="lex", domain=domain
)
target_arm_square = sp.groebner([r, z**2], r, z, e, order="lex", domain=domain)
for generator in (r, z**2):
    zero(arm_square.reduce(generator)[1], f"arm quotient contains {generator}")
for generator in (r**2, r * z, z**2, surface):
    zero(target_arm_square.reduce(generator)[1],
         f"arm quotient reverse containment {generator}")

# Generic first-normal jets reduce the Poisson bracket to the Bezout row.
aa0 = sp.Function("aa0")(e)
aa1 = sp.Function("aa1")(e)
bb0 = sp.Function("bb0")(e)
bb1 = sp.Function("bb1")(e)
generic_jet = sp.expand(bracket(aa0 + z * aa1, bb0 + z * bb1).subs({r: 0, z: 0}))
expected_jet = 3 * c**3 * (
    aa1 * sp.diff(bb0, e) - sp.diff(aa0, e) * bb1
)
zero(generic_jet - expected_jet, "generic conormal Poisson reduction")


# The minimal nodal immersion and its first-normal Bezout completion.
a0 = t**2
b0 = t**3 - t
a1 = -1 / (3 * c**3)
b1 = -t / (2 * c**3)
bezout = sp.factor(3 * c**3 * (a1 * sp.diff(b0, t) - sp.diff(a0, t) * b1))
gate(bezout == 1, "nodal Bezout completion")
gate(sp.gcd(sp.diff(a0, t), sp.diff(b0, t)) == 1, "nodal immersion")
zero(a0.subs(t, 1) - a0.subs(t, -1), "nodal A collision")
zero(b0.subs(t, 1) - b0.subs(t, -1), "nodal C collision")

# Two quadratics collide away from the diagonal only when their centers
# agree; that is exactly when their derivatives have a common root.
qa, qb, qc, qd, qe, qf = sp.symbols("qa qb qc qd qe qf", nonzero=True)
quad_a = qa * t**2 + qb * t + qc
quad_b = qd * t**2 + qe * t + qf
collision_centres = sp.resultant(
    qa * (t + s) + qb, qd * (t + s) + qe, t
)
derivative_centres = sp.resultant(sp.diff(quad_a, t), sp.diff(quad_b, t), t)
zero(collision_centres - (qa * qe - qd * qb), "quadratic collision centers")
zero(derivative_centres - 2 * (qa * qe - qd * qb),
     "quadratic common derivative center")

# Universal first-normal coefficient of the pulled-back nodal equation.
delta_normal = sp.factor(
    2 * b0 * b1 - (a0 - 1) * (3 * a0 - 1) * a1
)
zero(delta_normal + (t**2 - 1) / (3 * c**3),
     "chosen nodal normal coefficient")

# The coefficient is independent of the chosen Bezout completion.
u1, v1 = sp.symbols("u1 v1")
generic_bezout = 3 * c**3 * (u1 * (3 * t**2 - 1) - 2 * t * v1)
generic_delta_normal = (t**2 - 1) * (2 * t * v1 - (3 * t**2 - 1) * u1)
zero(generic_delta_normal + (t**2 - 1) * generic_bezout / (3 * c**3),
     "universal nodal normal coefficient")

# The smallest global lift and its exact failed bracket.
A = e**2 - z / (3 * c**3)
C = e**3 - e - e * z / (2 * c**3)
lifted = sp.factor(bracket(A, C))
expected = (c**3 + 2 * e * r) * (2 * c**3 + z) / (2 * c**6)
zero(lifted - expected, "lifted nodal bracket")

# Critical equations and the seven-point family.
critical = [sp.factor(bracket(A, q)) for q in variables]
critical_factor = c**3 + 2 * r * e
zero(critical[0] - (r**2 / c**3 - 18 * e * z**2),
     "critical r equation")
zero(critical[1] + 6 * e * critical_factor, "critical z equation")
zero(critical[2] + critical_factor / c**3, "critical e equation")
crit_poly = 8 * zeta**7 + 9 * c**15
r0 = 2 * zeta**3 / c**3
e0 = -c**6 / (4 * zeta**3)


def reduce_at_critical(expr: sp.Expr) -> sp.Expr:
    substituted = sp.cancel(expr.subs({r: r0, z: zeta, e: e0}))
    numerator = sp.together(substituted).as_numer_denom()[0]
    return sp.factor(sp.rem(numerator, crit_poly, zeta))


gate(reduce_at_critical(surface) == 0, "displayed critical point on surface")
gate(reduce_at_critical(c**3 + 2 * r * e) == 0, "displayed K factor")
for index, expression in enumerate(critical):
    gate(reduce_at_critical(expression) == 0,
         f"displayed critical bracket {index}")
gate(sp.gcd(crit_poly, sp.diff(crit_poly, zeta)) == 1,
     "critical degree-seven polynomial squarefree")

# Exhaustiveness, not merely seven examples.  The complete critical ideal on
# the surface equals the triangular ideal displayed below over Q(c).
critical_ideal = sp.groebner(
    [surface, c**3 + 2 * r * e, sp.together(critical[0] * c**3)],
    e, r, z, order="lex", domain=domain
)
critical_triangular = sp.groebner(
    [9 * c**9 * e - 2 * z**4,
     c**3 * r - 2 * z**3,
     8 * z**7 + 9 * c**15],
    e, r, z, order="lex", domain=domain
)
for generator in (9 * c**9 * e - 2 * z**4,
                  c**3 * r - 2 * z**3,
                  8 * z**7 + 9 * c**15):
    zero(critical_ideal.reduce(generator)[1],
         f"critical ideal forward {generator}")
for generator in (surface, c**3 + 2 * r * e,
                  sp.together(critical[0] * c**3)):
    zero(critical_triangular.reduce(generator)[1],
         f"critical ideal reverse {generator}")

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("THM3790_ARM_QUOTIENT B/I^2=C[e,z]/(z^2); r=0 because c^3*r=z^3-r^2*e")
print(f"THM3790_BEZOUT {bezout}")
print("THM3790_NODAL gcd(2*t,3*t^2-1)=1; gamma(1)=gamma(-1)=(1,0)")
print(f"THM3790_NODAL_NORMAL_COEFFICIENT {delta_normal}")
print(f"THM3790_LIFTED_BRACKET {lifted}")
print("THM3790_CRITICAL_POLY 8*zeta^7+9*c^15; squarefree=YES; count=7")
print("THM3790_CRITICAL_POINT r=2*zeta^3/c^3; e=-c^6/(4*zeta^3)")
print("THM3790_CRITICAL_IDEAL exact_triangular_exhaustion=YES")
print("THM3790_CRITICAL_REMAINDERS surface=0 K=0 Dr=0 Dz=0 De=0")
print(f"THM3790_CHECKS {CHECKS}")
print("RESULT=PASS")
