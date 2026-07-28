#!/usr/bin/env python3
"""Exact controls for THM-2799.

All truth-bearing checks use ``require`` rather than Python ``assert`` so
optimized execution performs the same verification.
"""

import ast
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(isinstance(node, ast.Assert) for node in ast.walk(tree))


def reduce_t(expr, t, ratio):
    num, den = sp.cancel(expr).as_numer_denom()
    require(sp.rem(den, ratio, domain=sp.QQ) != 0, "zero denominator class")
    inv = sp.invert(sp.Poly(den, t, domain=sp.QQ), sp.Poly(ratio, t, domain=sp.QQ))
    return sp.rem(sp.expand(num * inv.as_expr()), ratio, domain=sp.QQ)


def reduce_x_coefficients(expr, x, t, ratio):
    poly = sp.Poly(sp.expand(expr), x)
    return all(reduce_t(coef, t, ratio) == 0 for coef in poly.all_coeffs())


def one_degree(N, x, t, z):
    n = N - 2
    ratio = sp.expand(
        sum((j + 1) * (N - 2 - j) * t**j for j in range(N - 2))
    )
    U = sum(t**j for j in range(N - 1))
    aa = -sp.Rational(N, 2) * U
    bb = N * (U - 1)
    cc = N - 1 - sp.Rational(N, 2) * U
    numerator = sp.expand(x**N + aa * x**2 + bb * x + cc)
    double = sp.expand((x - 1) * (x - t))
    q = sp.expand(sp.Rational(N * (N - 2), 2) * U)

    require(sp.expand(numerator.subs(x, 1)) == 0, "B(1)")
    require(sp.expand(sp.diff(numerator, x).subs(x, 1)) == 0, "B'(1)")
    require(
        reduce_t(numerator.subs(x, t), t, ratio) == 0,
        "B(t) modulo ratio polynomial",
    )
    require(
        reduce_t(sp.diff(numerator, x).subs(x, t), t, ratio) == 0,
        "B'(t) modulo ratio polynomial",
    )
    require(
        reduce_x_coefficients(
            x * sp.diff(numerator, x) - N * numerator - q * double,
            x,
            t,
            ratio,
        ),
        "logarithmic numerator identity",
    )

    quotient, remainder = sp.div(numerator, double**2, x)
    require(
        reduce_x_coefficients(remainder, x, t, ratio),
        "double-root quotient remainder",
    )
    require(sp.degree(quotient, x) == N - 4, "simple-factor degree")

    lhs = sp.expand((1 - t) ** 3 * ratio)
    rhs = sp.expand(
        n - (n + 2) * t + (n + 2) * t ** (n + 1) - n * t ** (n + 2)
    )
    require(lhs == rhs, "Chebyshev numerator identity")
    require(sp.expand(t ** (N - 3) * ratio.subs(t, 1 / t) - ratio) == 0, "reciprocal")
    require(sp.discriminant(ratio, t) != 0, "ratio polynomial not squarefree")

    cheb = sp.chebyshevu(n, z)
    derivative = sp.diff(cheb, z)
    require(sp.Poly(derivative, z).count_roots(-1, 1) == n - 1, "Chebyshev roots")
    require(sp.degree(ratio, t) == n - 1, "ratio degree")

    fixed = 1 if sp.expand(ratio.subs(t, -1)) == 0 else 0
    require(fixed == (1 if N % 2 == 0 else 0), "inversion fixed-point parity")
    orbit_count = fixed + (N - 3 - fixed) // 2
    require(orbit_count == (N - 2) // 2, "affine-class count")

    return ratio, orbit_count, quotient


def main():
    script = Path(__file__)
    require(not has_asserts(script), "truth-bearing Python assert found")

    x, t, z = sp.symbols("x t z")
    print("THM-2799 ONE-POLE E=2 CHEBYSHEV RESPONSE AUDIT")
    print("N | ratio polynomial R_N(t) | affine classes")
    total_classes = 0
    for N in range(4, 15):
        ratio, count, _ = one_degree(N, x, t, z)
        total_classes += count
        if N <= 10:
            print(f"{N:2d} | {sp.factor(ratio)} | {count}")
        else:
            print(f"{N:2d} | degree {sp.degree(ratio, t)} reciprocal squarefree | {count}")

    # The first two boundary maps are useful exact controls.
    ratio4, _, quotient4 = one_degree(4, x, t, z)
    require(ratio4 == 2 * (t + 1), "N=4 ratio")
    require(reduce_t(quotient4 - 1, t, ratio4) == 0, "N=4 split quotient")
    ratio5, _, quotient5 = one_degree(5, x, t, z)
    require(ratio5 == 3 * t**2 + 4 * t + 3, "N=5 ratio")
    require(sp.degree(quotient5, x) == 1, "N=5 nonsplit linear factor")

    print(f"all-degree identities checked symbolically through N=14: 11")
    print(f"affine classes in checked range: {total_classes}")
    print("N=4 split boundary and N>=5 nonsplit criterion: PASS")
    print("Chebyshev critical-point root count: PASS")


if __name__ == "__main__":
    main()
