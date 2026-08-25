#!/usr/bin/env python3
"""Exact controls for THM-4122 and THM-4124.

The proofs of those theorems use cited results of Nguyen, Jelonek--Lason,
and Lang.  This companion checks the finite arithmetic, target-shear Newton
polygons, coefficient sidecar, normalization-inflation bookkeeping, and the
non-Keller 2:3 resonant-observable firewall.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def cross(origin, left, right):
    return ((left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0]))


def convex_hull(points):
    points = sorted(set(points))
    if len(points) <= 1:
        return points
    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def exponent_coefficients(poly, x, y):
    return {monomial: coefficient
            for monomial, coefficient in sp.Poly(sp.expand(poly), x, y).terms()
            if coefficient}


def newton_vertices(poly, x, y):
    support = set(exponent_coefficients(poly, x, y)) | {(0, 0)}
    return convex_hull(support)


def scaled_vertices(vertices, factor):
    return sorted((factor * x, factor * y) for x, y in vertices)


def jacobian(left, right, x, y):
    return sp.expand(sp.diff(left, x) * sp.diff(right, y)
                     - sp.diff(left, y) * sp.diff(right, x))


def width_bound(m, n, width):
    common = gcd(m, n)
    d, e = m // common, n // common
    cap = min(max(m, n) - 1, width)
    return d, e, tuple(
        rho for rho in range(1, common + 1)
        if rho * max(d, e) <= cap
    )


def main() -> None:
    x, y = sp.symbols("x y")

    # THM-4122's exact (72,108) conditional width arithmetic.
    d, e, rho_width6 = width_bound(72, 108, 6)
    _, _, rho_width5 = width_bound(72, 108, 5)
    _, _, rho_width30 = width_bound(72, 108, 30)
    require((d, e) == (2, 3), "reduced degree pair changed")
    require(rho_width6 == (1, 2), "exact width-(4,6) atlas changed")
    require(rho_width5 == (1,), "width-five sharpening changed")
    require(rho_width30 == tuple(range(1, 11)),
            "unconditional-width hostile changed")

    # Parametrization inflation: intrinsic (rho*d,rho*e), followed by a
    # polynomial cover of degree h, yields Nguyen-style (M*d,M*e) with
    # M=h*rho.  Raw M is not intrinsic.
    inflation_rows = []
    for rho in (1, 2, 5):
        for h_degree in (1, 3, 4):
            intrinsic = (rho * d, rho * e)
            inflated = (h_degree * intrinsic[0], h_degree * intrinsic[1])
            M = h_degree * rho
            require(inflated == (M * d, M * e), "inflation factor changed")
            inflation_rows.append((rho, h_degree, intrinsic, inflated))

    # Triangular-automorphism boundary controls for integral ratios. The
    # ratio-one row is a boundary case, not a nonautomorphic model.
    positive_rows = []
    P = x + y**2
    vertices_P = newton_vertices(P, x, y)
    coefficients_P = exponent_coefficients(P, x, y)
    for ratio in (1, 2, 3, 4):
        Q = y + P**ratio
        vertices_Q = newton_vertices(Q, x, y)
        require(jacobian(P, Q, x, y) == 1, "positive Keller control changed")
        require(sorted(vertices_Q) == scaled_vertices(vertices_P, ratio),
                f"Newton scaling changed at ratio {ratio}")
        coefficients_Q = exponent_coefficients(Q, x, y)
        ratios = []
        for vertex in vertices_P:
            if vertex == (0, 0):
                continue
            scaled = (ratio * vertex[0], ratio * vertex[1])
            ratios.append(sp.simplify(
                coefficients_Q[scaled] / coefficients_P[vertex] ** ratio
            ))
        require(set(ratios) == {1}, "all-vertex ratio lost synchronization")
        contracted = sp.expand(Q - P**ratio)
        require(contracted == y, "target shear did not contract to y")
        positive_rows.append((ratio, vertices_Q, tuple(ratios), sp.total_degree(contracted)))

    # Same scaled polygon but unequal vertex coefficients: support alone is
    # insufficient, and the map is deliberately non-Keller.
    Q_hostile = y + x**2 + 2 * y**4
    hostile_vertices = newton_vertices(Q_hostile, x, y)
    require(sorted(hostile_vertices) == scaled_vertices(vertices_P, 2),
            "scaled-polygon hostile disappeared")
    hostile_coefficients = exponent_coefficients(Q_hostile, x, y)
    hostile_ratios = (
        hostile_coefficients[(2, 0)] / coefficients_P[(1, 0)]**2,
        hostile_coefficients[(0, 4)] / coefficients_P[(0, 2)]**2,
    )
    hostile_jacobian = jacobian(P, Q_hostile, x, y)
    require(hostile_ratios == (1, 2), "coefficient hostile changed")
    require(not hostile_jacobian.is_constant(), "hostile accidentally became Keller")

    # The tempting 2:3 cancellation Q^2-cP^3 is an observable, not a target
    # automorphism.  The chain-rule identity is checked on a generic-enough
    # exact polynomial pair and then compared literally.
    c = sp.symbols("c")
    P_probe = x + 2 * y**2 + x * y
    Q_probe = y + 3 * x**2 + x * y**2
    resonant = sp.expand(Q_probe**2 - c * P_probe**3)
    resonant_jacobian = jacobian(P_probe, resonant, x, y)
    chain_rule = sp.expand(2 * Q_probe * jacobian(P_probe, Q_probe, x, y))
    require(sp.expand(resonant_jacobian - chain_rule) == 0,
            "2:3 resonant-observable firewall changed")
    require(not resonant_jacobian.is_constant(),
            "2:3 observable falsely preserved the Keller condition")

    # Radial-vertex lemma on exact rational polygons: rv is outside lambda P
    # for lambda<r.  Barycentric factor lambda/r lies strictly between 0,1.
    radial_rows = []
    for ratio in (2, 3, 5):
        for numerator in range(1, ratio):
            lam = Fraction(numerator, 1)
            barycentric_factor = lam / ratio
            require(0 < barycentric_factor < 1, "radial convexity changed")
            radial_rows.append((ratio, lam, barycentric_factor))

    print("THM4122/4124 ASYMPTOTIC WIDTH AND SHEAR CONTROL CERTIFICATE")
    print("scope=planar_nonproper_Keller_necessary_conditions;JC2=OPEN")
    print(f"degree_pair_72_108_reduced={(d,e)}")
    print(f"rho_exact_width6={rho_width6};rho_width5={rho_width5}")
    print(f"unconditional_width30_hostile={rho_width30}")
    print(f"inflation_rows={inflation_rows}")
    print(f"integral_shear_positive_rows={positive_rows}")
    print(f"hostile_scaled_vertices={hostile_vertices}")
    print(f"hostile_vertex_ratios={hostile_ratios};hostile_J={hostile_jacobian}")
    print("resonant_2_3_identity=J(P,Q^2-cP^3)=2Q*J(P,Q)")
    print(f"radial_convexity_rows={radial_rows}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
