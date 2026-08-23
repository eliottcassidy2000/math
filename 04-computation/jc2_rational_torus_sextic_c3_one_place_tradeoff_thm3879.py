#!/usr/bin/env python3
"""Exact companion for THM-3879's rational torus-sextic near miss.

The script reconstructs the trinodal quartic and its dual sextic, verifies
the torus decomposition, resolves the complete singular packet, checks the
connected Cardano/Kummer layer, and proves the projective one-place line
system empty.
"""

from __future__ import annotations

import hashlib

import sympy as sp


S, T = sp.symbols("S T")
A, B, C, U = sp.symbols("A B C U")
CHECKS = 0


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def equal(label: str, left: object, right: object) -> None:
    global CHECKS
    CHECKS += 1
    if left != right:
        raise AssertionError(f"{label}: {left!r} != {right!r}")


def gate(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def monic(poly: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    p = sp.Poly(poly, variable)
    return sp.expand(p.as_expr() / p.LC())


# A rational trinodal quartic: the two roots of each q_i map to one vertex.
qx = S**2 - T**2
qy = S**2 - 2 * T**2
qz = S * T
X = sp.expand(qy * qz)
Y = sp.expand(qx * qz)
Z = sp.expand(qx * qy)
zero("primal_inverse", S * X * Y - T * (2 * Y - X) * Z)
gate(sp.resultant(qx.subs(S, 1), qy.subs(S, 1), T) != 0,
     "qx and qy have disjoint roots")
gate(qx.subs({S: 1, T: 0}) * qy.subs({S: 1, T: 0}) != 0,
     "finite qz root is disjoint")
gate(qx.subs({S: 0, T: 1}) * qy.subs({S: 0, T: 1}) != 0,
     "infinite qz root is disjoint")

primal_s = sp.Matrix([sp.diff(v, S) for v in (X, Y, Z)])
primal_t = sp.Matrix([sp.diff(v, T) for v in (X, Y, Z)])
dual_raw = primal_s.cross(primal_t)
dual_A = sp.expand(-(S**2 - T**2) ** 2 * (S**2 + 2 * T**2))
dual_B = sp.expand((S**2 - 2 * T**2) ** 2 * (S**2 + T**2))
dual_C = 2 * S**3 * T**3
for label, raw, expected in zip(
    ("dual_A", "dual_B", "dual_C"), dual_raw, (dual_A, dual_B, dual_C)
):
    zero(label, raw - 4 * expected)

# The exact torus decomposition.
Q2 = 2 * A**2 + 3 * A * B + B**2 + 11 * C**2
Q3 = C * (4 * A**2 + 5 * A * B + 2 * B**2 + 14 * C**2)
F = sp.expand(4 * Q2**3 - 27 * Q3**2)
H = S * T * (S**4 + S**2 * T**2 + 2 * T**4)
pull = {A: dual_A, B: dual_B, C: dual_C}
zero("contact_conic_pullback", Q2.subs(pull) - 3 * H**2)
zero("contact_cubic_pullback", Q3.subs(pull) - 2 * H**3)
zero("sextic_equation", F.subs(pull))

# Bidual recovery proves that the degree-six parametrization is birational.
common = -108 * S**3 * T**3 * (S**4 + S**2 * T**2 + 2 * T**4) ** 3
common *= S**8 + 6 * S**6 * T**2 - 19 * S**4 * T**4 + 12 * S**2 * T**6 + 4 * T**8
for label, variable, primal in zip(
    ("bidual_X", "bidual_Y", "bidual_Z"), (A, B, C), (X, Y, Z)
):
    zero(label, sp.diff(F, variable).subs(pull) - common * primal)

# Complete singular locus.  In C=1 the singular Groebner basis has one
# linear x-row and the exact inner^2 times outer eliminant.
x, y = sp.symbols("x y")
f_aff = sp.expand(F.subs({A: x, B: y, C: 1}))
singular_gb = sp.groebner(
    [f_aff, sp.diff(f_aff, x), sp.diff(f_aff, y)], x, y, order="lex"
)
equal("finite_singular_gb_length", len(singular_gb.polys), 2)
inner = y**4 - 13 * y**2 + 128
outer = y**4 + 272 * y**2 + 64
equal(
    "finite_singular_eliminant",
    monic(singular_gb.polys[-1].as_expr(), y),
    sp.expand(inner**2 * outer),
)
gate(sp.Poly(singular_gb.polys[0].as_expr(), x).degree() == 1,
     "finite singular x coordinate is unique")
gate(sp.gcd(sp.Poly(inner, y), sp.Poly(outer, y)).degree() == 0,
     "inner and outer packets are disjoint")
gate(sp.gcd(sp.Poly(inner, y), sp.diff(sp.Poly(inner, y))).degree() == 0,
     "four finite inner points are reduced")
gate(sp.gcd(sp.Poly(outer, y), sp.diff(sp.Poly(outer, y))).degree() == 0,
     "four outer points are reduced")

def numerator_mod(expr: sp.Expr, modulus: sp.Expr) -> sp.Expr:
    numerator = sp.cancel(expr).as_numer_denom()[0]
    return sp.rem(sp.Poly(numerator, y), sp.Poly(modulus, y)).as_expr()


x_inner = (y**3 - 13 * y) / 16
x_outer = y * (y**2 + 136) / 192
for name, expression in (
    ("F", f_aff),
    ("F_x", sp.diff(f_aff, x)),
    ("F_y", sp.diff(f_aff, y)),
):
    zero(f"inner_{name}", numerator_mod(expression.subs(x, x_inner), inner))
    zero(f"outer_{name}", numerator_mod(expression.subs(x, x_outer), outer))

q_aff = Q2.subs({A: x, B: y, C: 1})
g_aff = Q3.subs({A: x, B: y, C: 1})
inner_gb = sp.groebner([q_aff, g_aff], x, y, order="lex")
equal("inner_coordinate_row", inner_gb.polys[0].as_expr(), 16 * x - y**3 + 13 * y)
equal("inner_eliminant", inner_gb.polys[1].as_expr(), inner)
inner_jac = sp.det(sp.Matrix([
    [sp.diff(q_aff, x), sp.diff(q_aff, y)],
    [sp.diff(g_aff, x), sp.diff(g_aff, y)],
]))
inner_jac_rem = numerator_mod(inner_jac.subs(x, x_inner), inner)
gate(sp.gcd(sp.Poly(inner_jac_rem, y), sp.Poly(inner, y)).degree() == 0,
     "all four finite inner points are transverse A2 cusps")

outer_hessian = sp.det(sp.hessian(f_aff, (x, y)))
outer_hessian_rem = numerator_mod(outer_hessian.subs(x, x_outer), outer)
gate(sp.gcd(sp.Poly(outer_hessian_rem, y), sp.Poly(outer, y)).degree() == 0,
     "all four outer points have nondegenerate A1 quadratic part")

# At C=0 there are exactly the two remaining transverse torus intersections.
c = sp.symbols("c")
F_inf = sp.expand(F.subs({A: x, B: 1, C: c}))
zero("infinity_section", F_inf.subs(c, 0) - 4 * (x + 1) ** 3 * (2 * x + 1) ** 3)
for root in (-1, sp.Rational(-1, 2)):
    for name, expression in (
        ("F", F_inf),
        ("F_x", sp.diff(F_inf, x)),
        ("F_c", sp.diff(F_inf, c)),
    ):
        zero(f"infinity_{root}_{name}", expression.subs({x: root, c: 0}))
q_inf = Q2.subs({A: x, B: 1, C: c})
g_inf = Q3.subs({A: x, B: 1, C: c})
inf_jac = sp.det(sp.Matrix([
    [sp.diff(q_inf, x), sp.diff(q_inf, c)],
    [sp.diff(g_inf, x), sp.diff(g_inf, c)],
]))
equal("infinity_cusp_one_transverse", inf_jac.subs({x: -1, c: 0}), -1)
equal("infinity_cusp_two_transverse", inf_jac.subs({x: -sp.Rational(1, 2), c: 0}), sp.Rational(1, 2))

# The natural cubic is irreducible S3, and its quadratic resolvent has the
# explicit codimension-one-unramified C3 Kummer radicand.
cubic = U**3 - Q2 * U - Q3
zero("cubic_discriminant", sp.discriminant(cubic, U) - F)
aa, bb, cc = sp.symbols("aa bb cc")
linear_root = aa * A + bb * B + cc * C
root_coefficients = sp.Poly(sp.expand(linear_root**3 - Q2 * linear_root - Q3), A, B, C).coeffs()
equal("no_linear_polynomial_root", [p.as_expr() for p in sp.groebner(
    root_coefficients, aa, bb, cc, order="lex"
).polys], [sp.Integer(1)])

# The monogenic cubic surface is already normal: its singular locus consists
# of four closed points over the outer nodes, hence has codimension two in a
# hypersurface surface.
cubic_aff = sp.expand(cubic.subs({A: x, B: y, C: 1}))
cubic_singular_gb = sp.groebner(
    [
        cubic_aff,
        sp.diff(cubic_aff, U),
        sp.diff(cubic_aff, x),
        sp.diff(cubic_aff, y),
    ],
    U, x, y,
    order="lex",
)
equal("cubic_singular_u_row", cubic_singular_gb.polys[0].as_expr(), 48 * U + y**2 + 88)
equal("cubic_singular_x_row", cubic_singular_gb.polys[1].as_expr(), 192 * x - y**3 - 136 * y)
equal("cubic_singular_outer_row", cubic_singular_gb.polys[2].as_expr(), outer)

w, rho = sp.symbols("w rho")
kummer_norm = sp.expand((rho * Q3 + w) * (rho * Q3 - w))
kummer_norm = kummer_norm.subs(w**2, 27 * Q3**2 - 4 * Q2**3).subs(rho**2, 27)
zero("kummer_norm_cube", kummer_norm - 4 * Q2**3)

# No projective line has a pullback supported at one normalization address.
la, lb, lc, kap, r, inv = sp.symbols("la lb lc kap r inv")
line_pullback = sp.expand(la * dual_A + lb * dual_B + lc * dual_C)
finite_sixth = sp.Poly(sp.expand(line_pullback - kap * (S - r * T) ** 6), S, T)
finite_equations = finite_sixth.coeffs() + [inv * kap - 1]
equal("finite_one_place_system", [p.as_expr() for p in sp.groebner(
    finite_equations, la, lb, lc, kap, r, inv, order="lex"
).polys], [sp.Integer(1)])
pure_t_equations = sp.Poly(sp.expand(line_pullback - kap * T**6), S, T).coeffs()
equal("infinite_one_place_system", sp.linsolve(
    pure_t_equations, (la, lb, lc, kap)
), sp.FiniteSet((0, 0, 0, 0)))
equal("two_place_boundary", dual_C, 2 * S**3 * T**3)

print("THM3879_PARAM", "[-(S2-T2)^2(S2+2T2):(S2-2T2)^2(S2+T2):2S3T3]")
print("THM3879_TORUS", "F=4Q2^3-27Q3^2;Q2_pull=3H^2;Q3_pull=2H^3")
print("THM3879_SINGULARITIES", "6A2(inner torus intersections)+4A1(outer nodes);genus=0")
print("THM3879_CUBIC", "U^3-Q2*U-Q3 irreducible;disc=F;Galois=S3")
print("THM3879_COMPLETION", "normal monogenic cubic;singular locus=four isolated points")
print("THM3879_KUMMER", "Norm(rho*Q3+w)=4Q2^3;rho^2=27;w^2=27Q3^2-4Q2^3")
print("THM3879_ONE_PLACE", "no line pulls back to a sixth power;line C=0 gives exactly two places")
print("THM3879_SCOPE", "explicit embedded torus sextic near miss;JC2 remains open")
semantic = (
    "dual trinodal rational sextic",
    "six transverse torus cusps and four outer nodes",
    "explicit connected S3 cubic and unramified C3 Kummer layer",
    "normal monogenic completion with isolated singularities",
    "all projective one-place charts excluded",
    "C3 gluing versus one-place infinity tradeoff",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic).encode()).hexdigest())
print("CHECKS", CHECKS)
