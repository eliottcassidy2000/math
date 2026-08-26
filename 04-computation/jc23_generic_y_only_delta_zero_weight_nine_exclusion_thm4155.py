#!/usr/bin/env python3
"""Exact primary certificate for THM-4155.

This script works on the Y-only, Delta=0 exact-weight-nine wall of the
normalized live (2,3) planar-Jacobian seam.  It computes the general symbolic
critical resultants in both source coordinates (s,p) and normalized
coordinates (X,T), reconstructs the valued Newton polygon and all labelled
boundary places, and checks the finite/full monodromy ledgers.

No ``assert`` is used, so every check remains active under ``python -O``.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction as F
from itertools import combinations
from math import gcd

import sympy as sp


CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(poly: sp.Expr, variable: sp.Symbol) -> int:
    terms = sp.Poly(poly, variable).terms()
    require(bool(terms), "zero polynomial has no valuation")
    return min(monomial[0] for monomial, coefficient in terms if coefficient != 0)


def remainder_at(
    poly: sp.Expr,
    t_value: sp.Rational,
    modulus: sp.Expr,
    X: sp.Symbol,
    T: sp.Symbol,
) -> sp.Expr:
    specialized = sp.cancel(poly.subs(T, t_value))
    require(sp.denom(specialized) == 1, "unexpected specialized denominator")
    return sp.factor(
        sp.rem(sp.Poly(specialized, X), sp.Poly(modulus, X)).as_expr()
    )


def convex_hull(points: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    points = sorted(set(points))

    def cross(
        origin: tuple[int, int],
        left: tuple[int, int],
        right: tuple[int, int],
    ) -> int:
        return (
            (left[0] - origin[0]) * (right[1] - origin[1])
            - (right[0] - origin[0]) * (left[1] - origin[1])
        )

    lower: list[tuple[int, int]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(
    vertices: tuple[tuple[int, int], ...],
) -> tuple[int, int, int]:
    area_twice = abs(
        sum(
            vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
            - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
            for index in range(len(vertices))
        )
    )
    boundary = sum(
        gcd(
            abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]),
        )
        for index in range(len(vertices))
    )
    require((area_twice - boundary + 2) % 2 == 0, "Pick parity failed")
    return area_twice, boundary, (area_twice - boundary + 2) // 2


def edge_packet(
    vertices: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, ...], tuple[tuple[object, ...], ...]]:
    packet: list[int] = []
    rows: list[tuple[object, ...]] = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        level = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - level
        if start[0] == end[0] == 0:
            require((length, index) == (3, 1), "affine vertical edge changed")
            continue
        packet.extend([index] * length)
        rows.append((start, end, length, inward, level, index))
    return tuple(sorted(packet, reverse=True)), tuple(rows)


def valued_support(
    coefficients: dict[tuple[int, int], sp.Expr],
) -> dict[tuple[int, int, int], sp.Expr]:
    raw: list[tuple[int, int, int, sp.Expr]] = [
        (2, 0, 0, sp.Integer(1)),
        (0, 1, 0, sp.Integer(-1)),
    ]
    for (i, j), coefficient in coefficients.items():
        raw.append((j + 2, i + j, 1, -coefficient))
        raw.append((j, i + j + 1, 1, coefficient))
    collapsed: dict[tuple[int, int, int], sp.Expr] = {}
    for x_coordinate, y_coordinate, height, coefficient in raw:
        key = (x_coordinate, y_coordinate, height)
        collapsed[key] = sp.expand(collapsed.get(key, sp.Integer(0)) + coefficient)
    return {key: value for key, value in collapsed.items() if value != 0}


def lower_planes(
    support: dict[tuple[int, int, int], sp.Expr],
) -> tuple[tuple[F, F, F], ...]:
    points = tuple(sorted(support))
    planes: set[tuple[F, F, F]] = set()
    for first, second, third in combinations(points, 3):
        x0, y0, z0 = first
        x1, y1, z1 = second
        x2, y2, z2 = third
        determinant = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        if determinant == 0:
            continue
        slope_x = F(
            (z1 - z0) * (y2 - y0) - (z2 - z0) * (y1 - y0),
            determinant,
        )
        slope_y = F(
            (x1 - x0) * (z2 - z0) - (x2 - x0) * (z1 - z0),
            determinant,
        )
        constant = F(z0) - slope_x * x0 - slope_y * y0
        gaps = [
            F(z) - slope_x * x - slope_y * y - constant
            for x, y, z in points
        ]
        if min(gaps) >= 0:
            planes.add((slope_x, slope_y, constant))
    return tuple(sorted(planes))


def main() -> None:
    s, p = sp.symbols("s p")
    X, T = sp.symbols("X T")
    Theta, Phi, Zeta = sp.symbols("Theta Phi Zeta")
    w, q = sp.symbols("w q")

    K0 = sp.Rational(2848, 45)
    t = p - s**2
    H = sp.expand(
        -3 * p
        + sp.Rational(8, 3) * p**2
        - sp.Rational(1376, 135) * p**3
        + K0 * s**2 * p**2
        + Phi * s * p**3
        + Theta * s**2 * p**3
        + Zeta * s**3 * p**3
    )
    G_source = -s**2 / (2 * t) + H

    P = T + X**2 * T**2
    Y = X * T * P
    G = sp.expand(
        -X**2 * T / 2
        - 3 * P
        + sp.Rational(8, 3) * P**2
        - sp.Rational(1376, 135) * P**3
        + K0 * Y**2
        + Phi * P**2 * Y
        + Theta * P * Y**2
        + Zeta * Y**3
    )
    require(
        sp.factor(
            sp.cancel(G_source.subs({s: X * T, p: P}) - G)
        ) == 0,
        "source-to-normalized reconstruction changed",
    )

    # General symbolic source critical projection.
    A = sp.factor(sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p))
    C_critical = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    B = sp.factor(sp.cancel((C_critical + s * A) / t**2))
    require(sp.denom(A).has(s, p) is False, "A acquired a source denominator")
    require(sp.denom(B).has(s, p) is False, "B acquired a source denominator")
    require(
        sp.factor(sp.cancel(t**2 * sp.diff(G_source, s) - p * A)) == 0,
        "first source critical identity changed",
    )
    require(
        sp.factor(
            sp.cancel(
                2 * t**2 * sp.diff(G_source, p) - (t**2 * B - s * A)
            )
        ) == 0,
        "second source critical identity changed",
    )
    require((sp.degree(A, s), sp.degree(B, s)) == (6, 3),
            "source critical degrees changed")
    require(sp.factor(sp.Poly(A, s).LC() - 3 * Zeta * p**2) == 0,
            "A leading row changed")
    require(sp.factor(sp.Poly(B, s).LC() - 9 * Zeta * p**2) == 0,
            "B leading row changed")

    top_face = Zeta * w**3 + Theta * w**2 + Phi * w - sp.Rational(1376, 135)
    top_discriminant = sp.factor(sp.discriminant(top_face, w))
    inner_wall = sp.factor(4 * Theta * K0**2 - 27 * Zeta**2)

    source_resultant = sp.factor(sp.resultant(A, B, s))
    require(valuation(source_resultant, p) == 6, "source p-artifact changed")
    source_residual_expr = sp.cancel(source_resultant / p**6)
    require(sp.denom(source_residual_expr) == 1,
            "source critical quotient is not polynomial")
    source_residual = sp.Poly(source_residual_expr, p)
    require(source_residual.degree() == 18, "source residual is not R18")
    require(
        sp.factor(source_residual.TC() + 46656 * Zeta * inner_wall) == 0,
        "source residual constant endpoint changed",
    )
    require(
        sp.factor(
            source_residual.LC()
            + 236196 * Zeta**5 * top_discriminant
        ) == 0,
        "source residual leading endpoint changed",
    )

    # General symbolic normalized critical projection.  This is intentionally
    # recomputed rather than inferred from the source projection.
    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    require(sp.denom(f) == 1, "G_X/T is not polynomial")
    require((sp.degree(f, X), sp.degree(h, X)) == (8, 9),
            "normalized critical X-degrees changed")
    require(sp.factor(sp.Poly(f, X).LC() - 9 * Zeta * T**8) == 0,
            "f common leading row changed")
    require(sp.factor(sp.Poly(h, X).LC() - 9 * Zeta * T**8) == 0,
            "h common leading row changed")

    normalized_resultant = sp.resultant(f, h, X)
    require(valuation(normalized_resultant, T) == 56,
            "normalized T-artifact changed")
    normalized_residual_expr = sp.cancel(
        normalized_resultant / (T**56 * (6 * T + 1)**2)
    )
    require(sp.denom(normalized_residual_expr) == 1,
            "normalized critical quotient is not polynomial")
    normalized_residual = sp.Poly(normalized_residual_expr, T)
    require(normalized_residual.degree() == 18,
            "normalized residual is not Q18")
    require(
        sp.factor(
            normalized_residual.TC()
            + sp.Rational(3**15, 2**7) * Zeta**7
        ) == 0,
        "normalized residual constant endpoint changed",
    )
    require(
        sp.factor(
            normalized_residual.LC()
            + sp.Rational(6561, 4)
            * Zeta**3
            * inner_wall**2
            * top_discriminant
        ) == 0,
        "normalized residual leading endpoint changed",
    )

    # The two universal pairs restore four points omitted by the residual
    # projection.  Their Hessian signs keep them visible as hostile controls.
    hessian = sp.factor(sp.det(sp.hessian(G, (X, T))))
    require(sp.factor(f.subs(T, 0)) == -X, "wrong f at T=0")
    require(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0,
            "wrong h at T=0")
    require(remainder_at(G, sp.Rational(0), X**2 + 6, X, T) == 0,
            "wrong zero critical value")
    require(remainder_at(hessian, sp.Rational(0), X**2 + 6, X, T) == 6,
            "zero critical pair is not Morse")
    require(remainder_at(f, -sp.Rational(1, 6), X**2 - 6, X, T) == 0,
            "wrong f at T=-1/6")
    require(remainder_at(h, -sp.Rational(1, 6), X**2 - 6, X, T) == 0,
            "wrong h at T=-1/6")
    require(
        remainder_at(
            G - sp.Rational(1, 2),
            -sp.Rational(1, 6),
            X**2 - 6,
            X,
            T,
        ) == 0,
        "wrong half critical value",
    )
    require(
        remainder_at(hessian, -sp.Rational(1, 6), X**2 - 6, X, T) == -6,
        "half critical pair is not Morse",
    )
    critical_length = normalized_residual.degree() + 4
    require(critical_length == 22, "affine critical length changed")

    # Nonempty exact control.  It also checks that no residual-discriminant
    # gate is being smuggled into the theorem.
    control = {Theta: sp.Integer(2), Phi: sp.Integer(3), Zeta: sp.Integer(5)}
    require(inner_wall.subs(control) == sp.Rational(63521957, 2025),
            "control inner wall changed")
    require(top_discriminant.subs(control) == -sp.Rational(10233928, 135),
            "control top discriminant changed")
    control_residual = sp.Poly(normalized_residual.as_expr().subs(control), T)
    require(sp.gcd(control_residual, control_residual.diff()).degree() == 0,
            "control residual is not squarefree")

    # Valued Newton polygon and labelled boundary packet.
    coefficients = {
        (1, 0): -sp.Integer(3),
        (2, 0): sp.Rational(8, 3),
        (3, 0): -sp.Rational(1376, 135),
        (0, 2): K0,
        (2, 1): Phi,
        (1, 2): Theta,
        (0, 3): Zeta,
    }
    support = valued_support(coefficients)
    polygon = convex_hull([(x_coordinate, y_coordinate)
                           for x_coordinate, y_coordinate, height in support])
    expected_polygon = ((0, 1), (2, 0), (5, 3), (3, 4), (0, 4))
    require(polygon == expected_polygon, "Delta-zero Newton polygon changed")
    require(polygon_ledger(polygon) == (27, 11, 9), "Pick ledger changed")
    require(
        lower_planes(support)
        == ((F(0), F(1, 3), F(-1, 3)), (F(1, 9), F(2, 9), F(-2, 9))),
        "lower planes changed",
    )
    packet, edge_rows = edge_packet(polygon)
    require(packet == (8, 3, 3, 3, 2, 2, 2, 1),
            "labelled packet changed")
    require(sum(index - 1 for index in packet) == 16,
            "packet defect changed")
    genus = 9
    require(sum(index - 1 for index in packet) == 2 * genus - 2,
            "Riemann--Hurwitz genus ledger changed")

    # The top edge gives three rational index-three places exactly off its
    # discriminant.  The index-two edge remains one prime cubic carrier.
    require(sp.degree(top_face, w) == 3, "top face is not cubic")
    carrier = Zeta * w**3 + K0 * w**2 - (q - sp.Rational(1, 2))
    carrier_discriminant = sp.factor(sp.discriminant(carrier, w))
    expected_carrier_discriminant = sp.factor(
        (q - sp.Rational(1, 2))
        * (4 * K0**3 - 27 * Zeta**2 * (q - sp.Rational(1, 2)))
    )
    require(carrier_discriminant == expected_carrier_discriminant,
            "cubic carrier discriminant changed")
    require(edge_rows[1][-1] == 2 and edge_rows[1][2] == 3,
            "cubic carrier edge changed")
    require(edge_rows[3][-1] == 3 and edge_rows[3][2] == 3,
            "rational top edge changed")

    # Complete response ledgers and strict monodromy contradictions.
    full_degree = sum(packet)
    full_required_support = full_degree - packet.count(1)
    full_overlap = full_degree - critical_length
    full_commutator_bound = 3 * full_overlap
    require((full_degree, full_overlap, full_commutator_bound,
             full_required_support) == (24, 2, 6, 23),
            "full response ledger changed")
    require(full_commutator_bound < full_required_support,
            "full response is no longer excluded")

    finite_origin_packet = (8, 3, 3, 3, 1)
    finite_degree = sum(finite_origin_packet)
    carrier_index = 3
    both_nonidentity = 2 * finite_degree - critical_length - 2 + carrier_index
    one_identity = 2 * finite_degree - critical_length - 1 + carrier_index
    both_identity = carrier_index
    finite_required = finite_degree - 1
    require(
        (finite_degree, carrier_index, both_nonidentity, one_identity,
         both_identity, finite_required) == (18, 3, 15, 16, 3, 17),
        "finite response ledger changed",
    )
    require(max(both_nonidentity, one_identity, both_identity) < finite_required,
            "finite response is no longer excluded")

    semantic_rows = (
        "wall=eta=Delta=0;K0=2848/45;zeta!=0",
        "gate=zeta*(4ThetaK0^2-27zeta^2)*Disc(C)!=0",
        "C=zeta*w^3+Theta*w^2+Phi*w-1376/135",
        "source=p^6*R18;TC=-46656*zeta*inner;LC=-236196*zeta^5*Disc(C)",
        "normalized=T^56*(6T+1)^2*Q18",
        "Q18_TC=-(3^15/2^7)*zeta^7",
        "Q18_LC=-(6561/4)*zeta^3*inner^2*Disc(C)",
        "critical_length=22",
        "polygon=(0,1),(2,0),(5,3),(3,4),(0,4);A2=27;B=11;g=9",
        "packet=8,3,3,3,2,2,2,1;defect=16",
        "labels=rational(8,3,3,3,1)+cubic(2,2,2)",
        "responses=full24;finite18;beta3",
        "monodromy=full6<23;finite15,16,3<17",
        "remaining=Disc(C)=0|inner=0|zeta=0",
    )
    semantic_digest = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("THM-4155 GENERIC Y-ONLY DELTA-ZERO EXACT CERTIFICATE")
    print(f"checks={CHECKS}")
    print("wall=eta=0;Delta=0;K0=2848/45")
    print("top_face=C(w)=zeta*w^3+Theta*w^2+Phi*w-1376/135")
    print("gate=zeta*(4*Theta*K0^2-27*zeta^2)*Disc(C)!=0")
    print("source_resultant=p^6*R18")
    print("source_R18_constant=-46656*zeta*(4*Theta*K0^2-27*zeta^2)")
    print("source_R18_leading=-236196*zeta^5*Disc(C)")
    print("normalized_resultant=T^56*(6*T+1)^2*Q18")
    print("Q18_constant=-(3^15/2^7)*zeta^7")
    print("Q18_leading=-(6561/4)*zeta^3*(4*Theta*K0^2-27*zeta^2)^2*Disc(C)")
    print("critical_length=22")
    print("control=Theta2;Phi3;zeta5;inner=63521957/2025;Disc(C)=-10233928/135;squarefree=YES")
    print("polygon=" + repr(polygon) + ";Pick=(27,11,9)")
    print("lower_planes=" + repr(lower_planes(support)))
    print("packet=8,3,3,3,2,2,2,1;defect=16;genus=9")
    print("labels=rational(8,3,3,3,1)+prime_cubic_carrier(2,2,2)")
    print("carrier=q-1/2=K0*w^2+zeta*w^3;degree=3;separable=YES")
    print("full_response=n24;overlap2;commutator_support<=6<23")
    print("finite_response=n18;beta3;capacities=15,16,3<17")
    print("remaining_walls=zeta=0|4*Theta*K0^2-27*zeta^2=0|Disc(C)=0")
    print("semantic_sha256=" + semantic_digest)
    print("verdict=PASS")


if __name__ == "__main__":
    main()
