#!/usr/bin/env python3
"""Exact symbolic audit of the 3+1+3 Keller-map depth-two fiber."""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x, y, z = sp.symbols("x y z")
a, b, c = sp.symbols("a b c")
u = 1 + x * y
F = (
    u**3 * z + y**2 * u * (4 + 3 * x * y),
    y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
)
L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2

v_star = (-sp.Rational(1, 4), 0, 0)
P_minus = (-1, sp.Rational(3, 2), sp.Rational(13, 2))
P_zero = (0, 0, -sp.Rational(1, 4))
P_plus = (1, -sp.Rational(3, 2), sp.Rational(13, 2))


def groebner_at(target):
    equations = [F[i] - target[i] for i in range(3)]
    return sp.groebner(equations, z, y, x, order="lex")


def univariate_x_polynomial(gbasis):
    candidates = [
        sp.Poly(poly.as_expr(), x, domain=sp.QQ)
        for poly in gbasis.polys
        if not poly.as_expr().has(y, z)
    ]
    require(len(candidates) == 1, "expected one x-only Groebner polynomial")
    return candidates[0].monic()


def main():
    jacobian = sp.factor(sp.det(sp.Matrix(F).jacobian((x, y, z))))
    require(jacobian == -2, "Jacobian determinant is not -2")

    g_v = groebner_at(v_star)
    p_v = univariate_x_polynomial(g_v)
    require(p_v.as_expr() == x**3 - x,
            "v* does not have x-roots 0,+1,-1")

    g_zero = groebner_at(P_zero)
    p_zero = univariate_x_polynomial(g_zero)
    require(p_zero.as_expr() == x + sp.Rational(1, 8),
            "P0 does not have the unique x-root -1/8")
    require([poly.as_expr() for poly in g_zero.polys] ==
            [z, y, x + sp.Rational(1, 8)],
            "P0 Groebner basis is not the unique point (-1/8,0,0)")

    g_minus = groebner_at(P_minus)
    g_plus = groebner_at(P_plus)
    p_minus = univariate_x_polynomial(g_minus)
    p_plus = univariate_x_polynomial(g_plus)
    expected_minus = x**3 - sp.Rational(404, 21119) * x \
        - sp.Rational(208, 21119)
    expected_plus = x**3 + sp.Rational(532, 20929) * x \
        - sp.Rational(208, 20929)
    require(p_minus.as_expr() == expected_minus, "wrong P- fiber cubic")
    require(p_plus.as_expr() == expected_plus, "wrong P+ fiber cubic")
    require(sp.discriminant(p_minus.as_expr(), x) != 0,
            "P- cubic is inseparable")
    require(sp.discriminant(p_plus.as_expr(), x) != 0,
            "P+ cubic is inseparable")

    L_zero = sp.expand(L.subs({a: P_zero[0], b: P_zero[1],
                               c: P_zero[2]}))
    require(L_zero == 0, "P0 is not on the Jelonek surface L=0")
    core_at_zero = sp.expand(
        (L * x**3 + (4 - 3 * b * c) * x - 2 * c).subs(
            {a: P_zero[0], b: P_zero[1], c: P_zero[2]}))
    require(core_at_zero == 4 * x + sp.Rational(1, 2),
            "core cubic does not drop directly to a linear polynomial")

    print("det J_F =", jacobian)
    print("F^-1(v*) x-polynomial =", p_v.as_expr())
    print("P- fiber cubic =", sp.expand(21119 * p_minus.as_expr()))
    print("P0 fiber Groebner basis =",
          [poly.as_expr() for poly in g_zero.polys])
    print("P+ fiber cubic =", sp.expand(20929 * p_plus.as_expr()))
    print("L(P0) =", L_zero)
    print("core E(x;P0) =", core_at_zero)
    print("depth-two finite count = 3 + 1 + 3 = 7")
    print("interpretation: two sheets escape at infinity; no finite collision")


if __name__ == "__main__":
    main()
