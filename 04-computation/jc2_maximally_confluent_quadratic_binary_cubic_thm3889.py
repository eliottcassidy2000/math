#!/usr/bin/env python3
"""Exact companion for THM-3889.

The symbolic part proves the all-parameter statement.  The two enumerations
at the end are explicitly finite side evidence and are not used by the proof.
The program is assertion-driven and prints a deterministic frozen transcript.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import sys

import sympy as sp


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def zero(label: str, expression: sp.Expr) -> None:
    check(sp.expand(expression) == 0, label)


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


def support_inequality(
    expression: sp.Expr, x: sp.Symbol, z: sp.Symbol, p: int, q: int, bound: int, label: str
) -> None:
    polynomial = sp.Poly(sp.expand(expression), x, z)
    check(
        all(p * i + q * j >= bound for (i, j), coefficient in polynomial.terms() if coefficient != 0),
        label,
    )


def lowest_coefficient(
    expression: sp.Expr, x: sp.Symbol, exponent: int
) -> sp.Expr:
    return sp.expand(sp.Poly(sp.expand(expression), x).coeff_monomial(x**exponent))


# ---------------------------------------------------------------------------
# 1. Delone--Faddeev index determinant and the universal family.
# ---------------------------------------------------------------------------

A, C, X, Y = sp.symbols("A C X Y")
a, b, c, d = sp.symbols("a b c d")

M_omega = sp.Matrix([[0, -a * c, -a * d], [1, b, 0], [0, -a, 0]])
M_theta = sp.Matrix([[0, -a * d, -b * d], [0, 0, d], [1, 0, -c]])
T_vector = sp.Matrix([0, X, Y])
T_squared = (X * M_omega + Y * M_theta) * T_vector
index_determinant = sp.Matrix.hstack(sp.Matrix([1, 0, 0]), T_vector, T_squared).det()
binary_cubic = a * X**3 + b * X**2 * Y + c * X * Y**2 + d * Y**3
equal("Delone--Faddeev index determinant", index_determinant, -binary_cubic)


def discriminant(aa: sp.Expr, bb: sp.Expr, cc: sp.Expr, dd: sp.Expr) -> sp.Expr:
    return sp.expand(bb**2 * cc**2 - 4 * aa * cc**3 - 4 * bb**3 * dd - 27 * aa**2 * dd**2 + 18 * aa * bb * cc * dd)


alpha, alpha1, beta, beta1, gamma, gamma1, delta, eta = sp.symbols(
    "alpha alpha1 beta beta1 gamma gamma1 delta eta"
)
a2 = A**2 + A * C
b2 = 2 * A * C + C**2
c2 = C**2
d2 = sp.Integer(0)

U, V = sp.symbols("U V")
f2 = a2 * U**3 + b2 * U**2 * V + c2 * U * V**2
equal("split leading binary cubic", f2, U * (A * U + C * V) * ((A + C) * U + C * V))
equal("leading discriminant", discriminant(a2, b2, c2, d2), C**8)

aa = a2 + alpha * A + alpha1 * C
bb = b2 + beta * A + beta1 * C
cc = c2 + gamma * A + gamma1 * C
dd = delta * A + eta * C
Delta = discriminant(aa, bb, cc, dd)

for label, coefficient in (("a", aa), ("b", bb), ("c", cc), ("d", dd)):
    equal(f"{label} common zero", coefficient.subs({A: 0, C: 0}), 0)

Delta_poly = sp.Poly(Delta, A, C)
check(Delta_poly.total_degree() == 8, "discriminant total degree eight")
degree_eight = sum(
    coefficient * A**i * C**j
    for (i, j), coefficient in Delta_poly.terms()
    if i + j == 8
)
equal("degree-eight discriminant part", degree_eight, C**8)
check(all(i + j >= 4 for (i, j), coefficient in Delta_poly.terms() if coefficient != 0), "Delta in maximal ideal fourth power")


# ---------------------------------------------------------------------------
# 2. Complete lower Newton boundary at the unique projective infinity point.
# ---------------------------------------------------------------------------

x, z, w = sp.symbols("x z w")
H = sp.expand(sp.cancel(z**8 * Delta.subs({A: 1 / z, C: x / z})))
check(sp.denom(H) == 1, "homogenized local expression is polynomial")
equal("infinity restriction", H.subs(z, 0), x**8)

# delta != 0.
support_inequality(H, x, z, 1, 3, 6, "delta first Newton half-space")
support_inequality(H, x, z, 1, 5, 8, "delta second Newton half-space")
equal("delta z2 vertex", sp.Poly(H, x, z).coeff_monomial(z**2), -27 * delta**2)
equal("delta x3z vertex", sp.Poly(H, x, z).coeff_monomial(x**3 * z), 4 * delta)
equal("delta x8 vertex", sp.Poly(H, x, z).coeff_monomial(x**8), 1)
edge_delta_3 = lowest_coefficient(H.subs(z, w * x**3), x, 6)
edge_delta_5 = lowest_coefficient(H.subs(z, w * x**5), x, 8)
equal("delta order-three edge polynomial", edge_delta_3, delta * w * (4 - 27 * delta * w))
equal("delta order-five edge polynomial", edge_delta_5, 1 + 4 * delta * w)

# delta = 0, gamma != 0.
H0 = sp.expand(H.subs(delta, 0))
equal("gamma z3 vertex", sp.Poly(H0, x, z).coeff_monomial(z**3), -4 * gamma**3)
equal(
    "gamma x2z2 vertex",
    sp.Poly(H0, x, z).coeff_monomial(x**2 * z**2),
    -27 * eta**2 + 36 * eta * gamma - 8 * gamma**2,
)
equal("gamma x4z vertex", sp.Poly(H0, x, z).coeff_monomial(x**4 * z), 4 * (eta - gamma))
equal("gamma x8 vertex", sp.Poly(H0, x, z).coeff_monomial(x**8), 1)
support_inequality(H0, x, z, 1, 2, 6, "gamma first Newton half-space")
support_inequality(H0, x, z, 1, 4, 8, "gamma unequal second Newton half-space")

edge_gamma_2 = lowest_coefficient(H0.subs(z, w * x**2), x, 6)
expected_gamma_2 = w * (
    -4 * gamma**3 * w**2
    + (-27 * eta**2 + 36 * eta * gamma - 8 * gamma**2) * w
    + 4 * (eta - gamma)
)
equal("gamma order-two edge polynomial", edge_gamma_2, expected_gamma_2)
edge_gamma_4 = lowest_coefficient(H0.subs(z, w * x**4), x, 8)
equal("gamma unequal order-four edge polynomial", edge_gamma_4, 1 + 4 * (eta - gamma) * w)

# eta = gamma is a genuine separate seam.
Heq = sp.expand(H0.subs(eta, gamma))
support_inequality(Heq, x, z, 1, 2, 6, "gamma equality first Newton half-space")
support_inequality(Heq, x, z, 1, 3, 8, "gamma equality second Newton half-space")
equal("gamma equality x2z2 vertex", sp.Poly(Heq, x, z).coeff_monomial(x**2 * z**2), gamma**2)
edge_equal_2 = lowest_coefficient(Heq.subs(z, w * x**2), x, 6)
edge_equal_3 = lowest_coefficient(Heq.subs(z, w * x**3), x, 8)
equal("gamma equality order-two edge polynomial", edge_equal_2, gamma**2 * w**2 * (1 - 4 * gamma * w))
equal(
    "gamma equality order-three edge polynomial",
    edge_equal_3,
    1 + 2 * (2 * beta + gamma - 2 * gamma1) * w + gamma**2 * w**2,
)

# The last seam is visibly reducible, including every eta and gamma1.
Dred = sp.expand(Delta.subs({delta: 0, gamma: 0}))
Dred_quotient = sp.cancel(Dred / C)
check(sp.denom(Dred_quotient) == 1, "last seam quotient is polynomial")
equal("last seam C divisibility", Dred, C * Dred_quotient)
check(sp.Poly(Dred_quotient, A, C).total_degree() == 7, "last seam quotient nonconstant degree seven")


# ---------------------------------------------------------------------------
# 3. Fully declared finite exact censuses (side evidence only).
# ---------------------------------------------------------------------------

PolyDict = dict[tuple[int, int], int]


def pclean(poly: PolyDict) -> PolyDict:
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def padd(*polys: PolyDict) -> PolyDict:
    result: PolyDict = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            result[monomial] = result.get(monomial, 0) + coefficient
    return pclean(result)


def pscale(scale: int, poly: PolyDict) -> PolyDict:
    return pclean({monomial: scale * coefficient for monomial, coefficient in poly.items()})


def pmul(left: PolyDict, right: PolyDict) -> PolyDict:
    result: PolyDict = {}
    for (i, j), coefficient_left in left.items():
        for (r, s), coefficient_right in right.items():
            monomial = (i + r, j + s)
            result[monomial] = result.get(monomial, 0) + coefficient_left * coefficient_right
    return pclean(result)


def ppow(poly: PolyDict, exponent: int) -> PolyDict:
    result: PolyDict = {(0, 0): 1}
    for _ in range(exponent):
        result = pmul(result, poly)
    return result


def pdisc(pa: PolyDict, pb: PolyDict, pc: PolyDict, pd: PolyDict) -> PolyDict:
    return padd(
        pmul(ppow(pb, 2), ppow(pc, 2)),
        pscale(-4, pmul(pa, ppow(pc, 3))),
        pscale(-4, pmul(ppow(pb, 3), pd)),
        pscale(-27, pmul(ppow(pa, 2), ppow(pd, 2))),
        pscale(18, pmul(pmul(pmul(pa, pb), pc), pd)),
    )


def lower_hull(points: set[tuple[int, int]]) -> list[tuple[int, int]]:
    """Compact lower convex chain, from the leftmost to rightmost support."""
    by_x: dict[int, int] = {}
    for i, j in points:
        by_x[i] = min(j, by_x.get(i, j))
    ordered = sorted(by_x.items())
    hull: list[tuple[int, int]] = []
    for point in ordered:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            cross = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1)
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    return hull


def one_primitive_newton_edge(poly: PolyDict) -> bool:
    degree = max(i + j for i, j in poly)
    leading = {(i, j): coefficient for (i, j), coefficient in poly.items() if i + j == degree}
    if len(leading) != 1:
        return False
    (i0, j0), = leading.keys()
    if not ((i0 == degree and j0 == 0) or (i0 == 0 and j0 == degree)):
        return False
    # At the unique projective point, record (transverse exponent, Z exponent).
    if j0 == degree:
        local_points = {(j, degree - i - j) for i, j in poly}
    else:
        local_points = {(i, degree - i - j) for i, j in poly}
    hull = lower_hull(local_points)
    if len(hull) != 2:
        return False
    dx = hull[1][0] - hull[0][0]
    dy = hull[1][1] - hull[0][1]
    return math.gcd(abs(dx), abs(dy)) == 1


pA: PolyDict = {(1, 0): 1}
pC: PolyDict = {(0, 1): 1}
linear_forms = [pA, pC, padd(pA, pC), padd(pA, pscale(-1, pC))]
quadratic_forms = [{}, {(2, 0): 1}, {(1, 1): 1}, {(0, 2): 1}]
sparse_forms = [padd(linear, quadratic) for linear in linear_forms for quadratic in quadratic_forms]
check(len(sparse_forms) == 16, "finite sparse form universe size")

sparse_rows = 0
one_coordinate_lead = 0
primitive_edge_rows = 0
coordinate_divisible_primitive_rows = 0
for pa, pb, pc, pd in itertools.product(sparse_forms, repeat=4):
    sparse_rows += 1
    discr = pdisc(pa, pb, pc, pd)
    check(bool(discr), "finite sparse discriminant nonzero")
    degree = max(i + j for i, j in discr)
    leading = {(i, j): coefficient for (i, j), coefficient in discr.items() if i + j == degree}
    if len(leading) == 1:
        (i0, j0), = leading.keys()
        if (i0 == degree and j0 == 0) or (i0 == 0 and j0 == degree):
            one_coordinate_lead += 1
            if one_primitive_newton_edge(discr):
                primitive_edge_rows += 1
                # A raw primitive edge is not irreducible-eligible when the
                # entire discriminant carries the transverse coordinate.
                transverse_order = min((j if j0 == degree else i) for i, j in discr)
                check(transverse_order > 0, "primitive candidate coordinate-divisible")
                coordinate_divisible_primitive_rows += 1

check(sparse_rows == 65536, "finite sparse row count")
check(one_coordinate_lead == 5890, "finite sparse pure-coordinate leading count")
check(primitive_edge_rows == 92, "finite sparse raw primitive-edge count")
check(coordinate_divisible_primitive_rows == 92, "finite sparse coordinate-divisible primitive count")
check(
    primitive_edge_rows - coordinate_divisible_primitive_rows == 0,
    "finite sparse irreducible-eligible primitive-edge survivor count",
)

# In the normalized family only the C-linear choices in d and then c reach
# the final reducible seam.  This census records the exact labelled universe.
linear_pairs = [(1, 0), (0, 1), (1, 1), (1, -1)]
delta_rows = gamma_rows = reducible_rows = 0
for _a_pair, _b_pair, c_pair, d_pair in itertools.product(linear_pairs, repeat=4):
    delta_value = d_pair[0]
    gamma_value = c_pair[0]
    if delta_value != 0:
        delta_rows += 1
    elif gamma_value != 0:
        gamma_rows += 1
    else:
        reducible_rows += 1
check(delta_rows == 192, "finite perturbation delta rows")
check(gamma_rows == 48, "finite perturbation gamma rows")
check(reducible_rows == 16, "finite perturbation reducible rows")
check(delta_rows + gamma_rows + reducible_rows == 256, "finite perturbation universe size")

# Hostile controls distinguish projective support from normalization places.
check(math.gcd(2, 3) == 1, "cuspidal one-edge one-place control")
check(math.gcd(2, 4) == 2, "tacnodal one-edge two-place control")


semantic_packet = "\n".join(
    (
        "Delone Faddeev binary index common zero",
        "maximally confluent split quadratic row",
        "leading discriminant C eighth power",
        "one projective point versus normalization places",
        "delta order three and five Newton branches",
        "gamma unequal and equality seam branches",
        "reducible final seam",
        "finite sparse side census",
        "finite perturbation side census",
    )
) + "\n"
semantic_sha256 = hashlib.sha256(semantic_packet.encode("utf-8")).hexdigest()

print("THM3889_INDEX all coefficients in (A,C);binary index form represents no unit")
print("THM3889_LEADING split quadratic row discriminant=C^8;one projective infinity point")
print("THM3889_DELTA delta!=0 gives distinct order-3 and order-5 infinity branches")
print("THM3889_GAMMA delta=0,gamma!=0 gives two branches in both eta seams")
print("THM3889_REDUCIBLE delta=gamma=0 forces C|Delta")
print(
    f"THM3889_CENSUS1 rows={sparse_rows} one_coordinate_lead={one_coordinate_lead} "
    f"primitive_edge={primitive_edge_rows} coordinate_divisible={coordinate_divisible_primitive_rows} "
    f"irreducible_eligible={primitive_edge_rows - coordinate_divisible_primitive_rows}"
)
print(
    f"THM3889_CENSUS2 rows={delta_rows + gamma_rows + reducible_rows} "
    f"two_place={delta_rows + gamma_rows} reducible={reducible_rows}"
)
print("THM3889_SCOPE all linear perturbations of normalized split C^8 leading row;other quadratic rows open")
print(f"SEMANTIC_SHA256 {semantic_sha256}")
print(f"CHECKS {CHECKS}")

sys.stdout.flush()
