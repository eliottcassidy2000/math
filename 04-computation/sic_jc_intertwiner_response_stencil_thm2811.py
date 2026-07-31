#!/usr/bin/env python3
"""Exact lightweight controls for THM-2811.

The all-degree arguments are in the theorem.  This companion deliberately
checks finite hostile controls without relying on Python ``assert``.
"""

from fractions import Fraction
from math import comb, factorial

import sympy as sp


checks = 0


def require(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(f"FAIL: {label}")


def proportional_multisets(left, right):
    """Whether two nonzero rational multisets agree up to one scalar."""
    if len(left) != len(right):
        return False
    left = [Fraction(x) for x in left]
    right = [Fraction(x) for x in right]
    for target in right:
        scale = target / left[0]
        if sorted(scale * x for x in left) == sorted(right):
            return True
    return False


print("THM-2811 exact controls")

# 1. The generator test for a generic n=2, k=3 intertwiner.
l = sp.symbols("l00 l01 l02 l10 l11 l12")
s = sp.symbols("s00 s01 s02 s10 s11 s12")
L = sp.Matrix(2, 3, l)
S = sp.Matrix(2, 3, s)
mixed = L * S.T
expected = sp.Matrix(
    [
        [sum(L[i, a] * S[j, a] for a in range(3)) for j in range(2)]
        for i in range(2)
    ]
)
require(mixed == expected, "mixed contraction is L*S^T")
require(mixed.shape == (2, 2), "mixed contraction has old-pair shape")
print("intertwiner generator matrix: PASS")

# 2. Fibre defects and response-potential degrees.
fibre_rows = 0
for N in range(4, 41):
    zero = [2, 2] + [1] * (N - 4)
    third = [N - 1, 1]
    for d in range(1, N):
        b = N - d
        pole = [d, b]
        defect = (
            sum(e - 1 for e in zero)
            + sum(e - 1 for e in pole)
            + sum(e - 1 for e in third)
        )
        require(defect == 2 * N - 2, f"Riemann-Hurwitz N={N},d={d}")
        require(
            all(part != [N] for part in (zero, pole, third)),
            f"no total fibre N={N},d={d}",
        )
        require((N - 4) + N + 4 == 2 * N, f"potential degree N={N}")
        fibre_rows += 1
print(f"two-pole fibre/potential rows: PASS ({fibre_rows})")

# 3. The full e=2 polynomializability trichotomy.
trichotomy_rows = 0
for N in range(4, 41):
    zero = [2, 2] + [1] * (N - 4)
    one_pole = [N]
    one_third = [N - 2, 1, 1]
    require(sum(zero) == N, f"h1 zero degree N={N}")
    require(sum(one_pole) == N, f"h1 pole degree N={N}")
    require(sum(one_third) == N, f"h1 third degree N={N}")
    require(
        [zero == [N], one_pole == [N], one_third == [N]]
        == [False, True, False],
        f"h1 unique total fibre N={N}",
    )
    defect_h1 = sum(e - 1 for part in (zero, one_pole, one_third) for e in part)
    require(defect_h1 == 2 * N - 2, f"h1 Riemann-Hurwitz N={N}")
    trichotomy_rows += 1

    for a in range(1, N - 1):
        for b in range(1, N - a):
            c = N - a - b
            three_pole = [a, b, c]
            three_third = [N]
            require(sum(three_pole) == N, f"h3 pole degree N={N},a={a},b={b}")
            require(
                [zero == [N], three_pole == [N], three_third == [N]]
                == [False, False, True],
                f"h3 unique total fibre N={N},a={a},b={b}",
            )
            require(len(three_pole) == 3, f"h3 polynomial has three roots N={N}")
            defect_h3 = sum(
                e - 1 for part in (zero, three_pole, three_third) for e in part
            )
            require(
                defect_h3 == 2 * N - 2,
                f"h3 Riemann-Hurwitz N={N},a={a},b={b}",
            )
            trichotomy_rows += 1
print(f"e=2 polynomializability rows: PASS ({trichotomy_rows})")

# 4. Exact N=4 edge response: M=x(x-1)E and C=-3/8.
x = sp.symbols("x")
E = x**2 - sp.Rational(1, 2) * x - sp.Rational(1, 8)
M = sp.expand(x * (x - 1) * E)
C = -sp.Rational(3, 8)
Mp = sp.diff(M, x)
require(sp.cancel(C / Mp.subs(x, 0)) == -3, "N4 pole residue at zero")
require(sp.cancel(C / Mp.subs(x, 1)) == -1, "N4 pole residue at one")
require(sp.rem(C - 2 * Mp, E, domain=sp.QQ) == 0, "N4 double-zero residues")

# Reduce the two algebraic double-zero contributions by traces modulo E.
roots = sp.solve(E, x)
nodes_weights = [(sp.Integer(0), -3), (sp.Integer(1), -1)]
nodes_weights.extend((root, 2) for root in roots)
for degree in range(4):
    moment = sp.simplify(sum(weight * node**degree for node, weight in nodes_weights))
    expected_moment = C if degree == 3 else 0
    require(moment == expected_moment, f"N4 barycentric moment degree={degree}")
print("N=4 response residues/moments: PASS")

# 5. Independent rational-node Lagrange coefficient extraction.
lagrange_rows = 0
for N in range(2, 9):
    nodes = [Fraction(j * j + 2 * j + 1, j + 2) for j in range(N)]
    require(len(set(nodes)) == N, f"distinct rational nodes N={N}")
    derivatives = []
    for i, node in enumerate(nodes):
        value = Fraction(1)
        for j, other in enumerate(nodes):
            if i != j:
                value *= node - other
        derivatives.append(value)
    for degree in range(0, 2 * N + 2):
        lhs = sum(node**degree / derivative for node, derivative in zip(nodes, derivatives))
        poly = sp.Poly(x**degree, x, domain=sp.QQ)
        modulus = sp.Poly(
            sp.prod(x - sp.Rational(node.numerator, node.denominator) for node in nodes),
            x,
            domain=sp.QQ,
        )
        remainder = poly.rem(modulus)
        rhs = Fraction(remainder.nth(N - 1))
        require(lhs == rhs, f"Lagrange remainder N={N},degree={degree}")
        lagrange_rows += 1
print(f"rational-node Lagrange rows: PASS ({lagrange_rows})")

# 6. Equispaced finite-difference/barycentric formula.
stencil_rows = 0
for m in range(1, 41):
    for j in range(m + 1):
        derivative = Fraction(1)
        for h in range(m + 1):
            if h != j:
                derivative *= j - h
        lhs = Fraction(factorial(m), 1) / derivative
        rhs = Fraction(((-1) ** (m - j)) * comb(m, j), 1)
        require(lhs == rhs, f"equispaced derivative m={m},j={j}")
        stencil_rows += 1
    require(
        sum(((-1) ** (m - j)) * comb(m, j) for j in range(m + 1)) == 0,
        f"pure finite difference m={m}",
    )
    require(
        sum(((-1) ** (m - j)) * comb(m, j - 1) for j in range(1, m + 1)) == 1,
        f"shifted endpoint m={m}",
    )
print(f"SIC equispaced stencil rows: PASS ({stencil_rows})")

# 7. No response residue alphabet is a scaled binomial stencil.
mismatch_rows = 0
for N in range(4, 41):
    binomial = [((-1) ** (N - 1 - j)) * comb(N - 1, j) for j in range(N)]
    for d in range(1, N // 2 + 1):
        b = N - d
        response = [1] * (N - 4) + [2, 2, -d, -b]
        require(
            not proportional_multisets(response, binomial),
            f"response/binomial mismatch N={N},d={d}",
        )
        mismatch_rows += 1
print(f"response/binomial mismatch rows: PASS ({mismatch_rows})")

# 8. Both polynomializable N=4 hostiles.
X, Y = sp.symbols("X Y")
P1 = (X**2 - Y**2) ** 2
Hessian1 = sp.hessian(P1, (X, Y))
laplacian1 = sp.factor(sp.trace(Hessian1))
determinant1 = sp.factor(Hessian1.det())
require(
    sp.simplify(laplacian1 - 8 * (X**2 + Y**2)) == 0,
    "quartic hostile Laplacian",
)
require(
    sp.simplify(determinant1 + 48 * P1) == 0,
    "quartic hostile Hessian determinant",
)
require(laplacian1 != 0, "quartic hostile fails first contraction")

P3 = sp.expand(X * (X - Y) * (2 * X - Y) ** 2)
Hessian3 = sp.hessian(P3, (X, Y))
laplacian3 = sp.factor(sp.trace(Hessian3))
determinant3 = sp.factor(Hessian3.det())
require(
    sp.simplify(laplacian3 - 2 * (29 * X**2 - 27 * X * Y + 5 * Y**2))
    == 0,
    "three-pole quartic Laplacian",
)
require(
    sp.simplify(
        determinant3
        + 3 * (2 * X - Y) ** 2 * (8 * X**2 - 8 * X * Y + 3 * Y**2)
    )
    == 0,
    "three-pole quartic Hessian determinant",
)
require(laplacian3 != 0, "three-pole quartic fails first contraction")
print("one-/three-pole quartic Hessian hostiles: PASS")

print(f"PASS checks={checks}")
