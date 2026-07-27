#!/usr/bin/env python3
"""Exact companion for THM-2468."""

from __future__ import annotations

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_unit_basis(groebner_basis: s.GroebnerBasis) -> bool:
    return (
        len(groebner_basis.polys) == 1
        and groebner_basis.polys[0].as_expr() == 1
    )


p, v, zeta, lam = s.symbols("p v zeta lam")
a, b, c, d, e, h = s.symbols("a b c d e h")

f1 = (
    -2981440 * p * v
    + 819896 * p * zeta
    + 24640 * p
    + 3689532 * v**2
    - 1449459 * v * zeta
    - 101640 * v
    + 83853 * zeta
    + 252
)
f2 = (
    -1319329792 * lam * p**3
    + 1443016960 * p * v**2
    - 71554560 * p * v
    + 65591680 * p * zeta
    + 98560 * p
    - 1190488992 * v**3
    + 147581280 * v**2
    - 162339408 * v * zeta
    - 1219680 * v
    + 15944049 * zeta**2
    + 2236080 * zeta
    + 672
)

raw_resultant = s.resultant(f1, f2, zeta)
R = s.Poly(raw_resultant, p, v, lam).primitive()[1].as_expr()
require(s.cancel(raw_resultant / R) == 2**4 * 11**6, "resultant content changed")
require(len(s.Poly(R, p, v, lam).terms()) == 24, "resultant support changed")
require(s.degree(R, p) == 5 and s.degree(R, v) == 5, "resultant bidegree changed")

A = 616 * p - 1089 * v + 63
require(
    s.factor(s.diff(R, lam) / (p**3 * A**2)) == -82458112,
    "B-axis pencil law changed",
)

L5 = (
    155624547606 * v**5
    + 3215383215 * v**4
    - 1700698560 * v**3
    + 58124770 * v**2
    - 855470 * v
    + 2583
)
constant_in_p = s.Poly(R, p).coeff_monomial(1)
linear_in_p = s.Poly(R, p).coeff_monomial(p)
require(s.expand(constant_in_p + 567 * L5) == 0, "fixed quintic section changed")
require(
    s.gcd(s.Poly(L5, v), s.Poly(s.diff(L5, v), v)).degree() == 0,
    "L5 lost squarefreeness",
)
require(
    s.gcd(s.Poly(constant_in_p, v), s.Poly(linear_in_p, v)).degree() == 0,
    "p-zero section lost first-order separation",
)
require(
    s.expand(s.diff(R, v).subs(p, 0) + 567 * s.diff(L5, v)) == 0,
    "p-zero smooth-section derivative changed",
)

top = s.Integer(0)
for (ip, iv, il), coefficient in s.Poly(R, p, v, lam).terms():
    if ip + iv == 5:
        top += coefficient * p**ip * v**iv * lam**il
expected_top = -38974342 * (56 * p - 99 * v) ** 2 * (
    256 * lam * p**3 - 280 * p * v**2 + 231 * v**3
)
require(s.expand(top - expected_top) == 0, "top homogeneous factorization changed")
require(
    s.factor(top.subs(v, 0)) == -31289225347072 * lam * p**5,
    "top form unexpectedly acquired a vertical factor",
)
require(
    s.factor(top.subs(p, 0)) == -88239118492602 * v**5,
    "p-zero section acquired a point at infinity",
)

# Any reducible total-degree-five polynomial has a line or conic factor.
# Since v does not divide the top form when lambda is nonzero, each such
# factor has nonzero leading p coefficient and the following monic shapes.
line_remainder = s.Poly(s.expand(R.subs(p, a * v + b)), v)
gb_line = s.groebner(
    [coefficient for _, coefficient in line_remainder.terms()],
    a,
    b,
    lam,
    order="grevlex",
)
require(is_unit_basis(gb_line), "line-factor coefficient ideal stopped being unit")

quadratic = p**2 + (a * v + b) * p + c * v**2 + d * v + e
quadratic_remainder = s.Poly(
    s.rem(s.Poly(R, p), s.Poly(quadratic, p)).as_expr(),
    p,
    v,
)
quadratic_equations = [coefficient for _, coefficient in quadratic_remainder.terms()]

# Type I: both copies of the repeated wall line 56p-99v.
type_i_substitution = {
    a: -s.Rational(99, 28),
    c: s.Rational(9801, 3136),
}
type_i = [s.factor(eq.subs(type_i_substitution)) for eq in quadratic_equations]
gb_type_i = s.groebner(type_i, b, d, e, lam, order="grevlex")
require(is_unit_basis(gb_type_i), "quadratic top type I stopped being impossible")

# Write h for a root of 256*lambda*h^3-280*h+231.  For lambda!=0,
# h and 280h-231 are both nonzero.
lambda_h = (280 * h - 231) / (256 * h**3)
Q = 256 * lam * p**3 - 280 * p * v**2 + 231 * v**3
require(s.factor(Q.subs({lam: lambda_h, p: h * v})) == 0, "cubic-root chart changed")


def cleared_primitive_equations(substitution: dict[s.Symbol, s.Expr]) -> list[s.Expr]:
    equations: list[s.Expr] = []
    for equation in quadratic_equations:
        numerator = s.cancel(equation.subs(substitution)).as_numer_denom()[0]
        primitive = s.Poly(numerator, b, d, e, h).primitive()[1].as_expr()
        if primitive != 0:
            equations.append(primitive)
    return equations


# Type II: one wall line and one cubic-root line.
type_ii_substitution = {
    lam: lambda_h,
    a: -(s.Rational(99, 56) + h),
    c: s.Rational(99, 56) * h,
}
type_ii_top = p**2 + type_ii_substitution[a] * p * v + type_ii_substitution[c] * v**2
require(
    s.factor(
        s.rem(
            s.Poly(top.subs(lam, lambda_h), p),
            s.Poly(type_ii_top, p),
        ).as_expr()
    )
    == 0,
    "quadratic top type II stopped dividing the top form",
)
type_ii = cleared_primitive_equations(type_ii_substitution)
gb_type_ii = s.groebner(type_ii, b, d, e, h, order="grevlex")
require(
    gb_type_ii.reduce(h)[1] == 0,
    "quadratic top type II no longer forces the forbidden h=0 boundary",
)

# Type III: the product of the other two cubic-root lines, Q/(p-hv).
type_iii_substitution = {
    lam: lambda_h,
    a: h,
    c: -231 * h**2 / (280 * h - 231),
}
type_iii_top = p**2 + type_iii_substitution[a] * p * v + type_iii_substitution[c] * v**2
require(
    s.factor(
        s.rem(
            s.Poly(Q.subs(lam, lambda_h), p),
            s.Poly(type_iii_top, p),
        ).as_expr()
    )
    == 0,
    "quadratic top type III stopped dividing the cubic top factor",
)
type_iii = cleared_primitive_equations(type_iii_substitution)
gb_type_iii = s.groebner(type_iii, b, d, e, h, order="grevlex")
require(
    gb_type_iii.reduce(h**4)[1] == 0,
    "quadratic top type III no longer forces the forbidden h=0 boundary",
)

# The y=0 boundary of the original first flux on C=D=E=0.
B, y, u, Z = s.symbols("B y u Z")
first_flux_A = 616 * B - 1089 * u + 63 * y**2
first_flux_K = (
    -745360 * B * u * y
    + 6160 * B * y**3
    + 922383 * u**2 * y
    - 25410 * u * y**3
    + 63 * y**5
)
N1 = 1331 * first_flux_A * Z + 4 * first_flux_K
require(
    s.expand(N1.subs(y, 0) - 1331 * (616 * B - 1089 * u) * Z) == 0,
    "y=0 first-flux boundary changed",
)

# The connected quadratic lift Y^2=1/p branches at the five simple
# p-zero places.  The branch divisor has even degree, so it has at least six
# places.  Riemann-Hurwitz then gives genus at least two without knowing the
# genus of the plane-quintic normalization.
finite_simple_branch_points = s.degree(L5, v)
minimum_even_branch_count = finite_simple_branch_points + (
    finite_simple_branch_points % 2
)
minimum_lift_genus = -1 + minimum_even_branch_count // 2
require(finite_simple_branch_points == 5, "square-lift simple-zero count changed")
require(minimum_even_branch_count == 6, "square-lift parity floor changed")
require(minimum_lift_genus == 2, "square-lift genus floor changed")

print("THM-2468 exact companion")
print("resultant_content=2^4*11^6")
print("resultant_terms=24")
print("resultant_bidegree=(5,5)")
print("pencil=R_0-82458112*lambda*p^3*A^2")
print("top_factor=(56p-99v)^2*(256lambda*p^3-280p*v^2+231v^3)")
print("line_factor_unit_ideal=YES")
print("quadratic_top_types=3")
print("quadratic_type_I_unit_ideal=YES")
print("quadratic_type_II_forces_h=0")
print("quadratic_type_III_forces_h^4=0")
print("uniform_absolute_irreducibility=lambda_nonzero")
print("p_zero_section=-567*L5")
print("L5_degree=5")
print("L5_squarefree=YES")
print("square_lift_minimum_branch_points=6")
print("square_lift_minimum_genus=2")
print("y_zero_boundary=Z_zero_on_open_first_flux")
print("ALL CHECKS PASSED")
