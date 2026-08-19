#!/usr/bin/env python3
"""Exact audit of the cusp-square packet as an inverse-cubic owner.

For arbitrary U(v,y), the packet (L,T,U,S) satisfies
S^2=T^3+27*L*U^2.  The associated cubic L*X^3+T*X+2*U factors into one
marked rational root and a quadratic Kummer pair.  The script verifies this
factorization, all six natural projection obstructions, and the existing
everywhere-immersive explicit four-coordinate packet.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    """Keep every truth-bearing check active under ``python -O``."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def jac(first: sp.Expr, second: sp.Expr, v: sp.Symbol,
        y: sp.Symbol) -> sp.Expr:
    return sp.factor(
        sp.diff(first, v) * sp.diff(second, y)
        - sp.diff(first, y) * sp.diff(second, v)
    )


def main() -> None:
    v, y, X = sp.symbols("v y X")
    U = sp.Function("U")(v, y)

    T = y**2 - 6 * v * U
    S = y**3 - 9 * v * U * y
    L = v**2 * (8 * v * U - y**2)
    E = sp.expand(L * X**3 + T * X + 2 * U)

    identity = sp.factor(S**2 - T**3 - 27 * L * U**2)
    require(identity == 0, "cusp-square identity")

    qa = v * (8 * v * U - y**2)
    qb = y**2 - 8 * v * U
    qc = 2 * U
    quadratic = qa * X**2 + qb * X + qc
    factorization = sp.expand((v * X + 1) * quadratic)
    require(sp.expand(E - factorization) == 0, "marked-root factorization")

    quadratic_disc = sp.factor(qb**2 - 4 * qa * qc)
    require(sp.expand(quadratic_disc - y**2 * (y**2 - 8 * v * U)) == 0,
            "quadratic Kummer discriminant")
    resultant = sp.factor(sp.resultant(v * X + 1, quadratic, X))
    require(sp.expand(resultant + 2 * v * (y**2 - 9 * v * U)) == 0,
            "marked-root/quadratic resultant")

    cubic_disc = sp.factor(sp.discriminant(E, X))
    require(sp.expand(cubic_disc + 4 * L * S**2) == 0,
            "inverse-cubic discriminant square")
    require(sp.expand(cubic_disc - quadratic_disc * resultant**2) == 0,
            "product discriminant factorization")

    # All six natural two-coordinate projections fail the Keller test.
    jac_TS = jac(T, S, v, y)
    require(sp.factor(jac_TS - 54 * v * U * (U + v * sp.diff(U, v))) == 0,
            "Jac(T,S)")
    divisibility = {}
    for first, second, label in [(L, T, "LT"), (L, U, "LU"),
                                 (L, S, "LS")]:
        value = jac(first, second, v, y)
        require(sp.factor(value.subs(v, 0)) == 0,
                f"v divides Jac({label})")
        divisibility[label] = sp.factor(value / v)

    jac_TU = jac(T, U, v, y)
    require(sp.factor(jac_TU
                      + 2 * (y * sp.diff(U, v)
                             + 3 * U * sp.diff(U, y))) == 0,
            "Jac(T,U) transport equation")
    jac_US_axis = sp.factor(jac(U, S, v, y).subs(v, 0))
    expected_US_axis = 3 * y * (
        y * sp.Subs(sp.diff(U, v), v, 0)
        + 3 * U.subs(v, 0) * sp.diff(U.subs(v, 0), y)
    )
    require(sp.simplify(jac_US_axis - expected_US_axis) == 0,
            "Jac(U,S) axis factor")

    # Positive four-coordinate control from the incoming exact scan.
    Up = (1 + y - sp.Rational(1, 2) * y**2
          - sp.Rational(3, 2) * v * y * (y - 3))
    Tp = sp.expand(y**2 - 6 * v * Up)
    Sp = sp.expand(y**3 - 9 * v * Up * y)
    Lp = sp.expand(v**2 * (8 * v * Up - y**2))
    packet = [Lp, Tp, Up, Sp]
    minors = [
        sp.expand(jac(packet[i], packet[j], v, y))
        for i in range(4)
        for j in range(i + 1, 4)
    ]
    require(sp.factor(Sp**2 - Tp**3 - 27 * Lp * Up**2) == 0,
            "explicit packet identity")
    minor_ideal = sp.groebner(minors, v, y, order="lex")
    minor_unit = minor_ideal.reduce(sp.Integer(1))[1] == 0
    require(minor_unit, "immersive packet unit minor ideal")

    weights = sp.symbols("weight0:6")
    linear_expr = sp.Poly(
        sum(weight * minor for weight, minor in zip(weights, minors)) - 1,
        v, y,
    )
    linear_solution = sp.linsolve(
        [coefficient for _, coefficient in linear_expr.terms()], weights
    )
    require(linear_solution is sp.EmptySet or linear_solution == sp.EmptySet,
            "constant-linear projection no-go")

    print("CUSP-SQUARE PACKET -- MARKED-ROOT/KUMMER OWNER AUDIT")
    print("L=v^2*(8*v*U-y^2), T=y^2-6*v*U, S=y^3-9*v*U*y")
    print(f"S^2-T^3-27*L*U^2 = {identity}")
    print(f"E(X)=L*X^3+T*X+2*U factors as (v*X+1)*({sp.factor(quadratic)})")
    print(f"quadratic discriminant = {quadratic_disc}")
    print(f"linear/quadratic resultant = {resultant}")
    print(f"cubic discriminant = {cubic_disc}")
    print(f"Jac(T,S) = {jac_TS}")
    print("Jac(L,T), Jac(L,U), Jac(L,S) are divisible by v")
    print(f"Jac(T,U) = {jac_TU}")
    print(f"Jac(U,S)|_(v=0) = {jac_US_axis}")
    print(f"explicit packet minor Groebner basis contains 1: {minor_unit}")
    print(f"constant linear minor-combination solution set = {linear_solution}")
    print("VERDICT: the packet is one marked escaping root plus a quadratic Kummer pair.")
    print("SCOPE: natural/constant-linear projections fail; nonlinear projections remain open.")


if __name__ == "__main__":
    main()
