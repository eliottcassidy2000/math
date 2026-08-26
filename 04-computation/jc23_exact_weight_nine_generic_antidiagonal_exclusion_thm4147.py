#!/usr/bin/env python3
"""Primary exact audit extending THM-4147 across its anti-diagonal wall.

The script is standalone.  It reconstructs the complete normalized M=9
source, checks THM-4147's three generic critical projections at an exact
positive control, and treats the formerly excluded ``zeta=-eta`` wall
symbolically over QQ[Delta,Theta,Phi,eta].  The latter calculation uses a
source (s,p) critical projection independent of the control calculation in
(X,T).  It also normalizes the repeated boundary point, freezes the labelled
packet, and checks the finite/full monodromy budgets.

No ``assert`` is used, so every gate remains live under ``python -O``.
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
    return min(monomial[0] for monomial, _ in terms)


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
            - (left[1] - origin[1]) * (right[0] - origin[0])
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


def valued_support(
    coefficients: dict[tuple[int, int], sp.Rational],
) -> dict[tuple[int, int, int], sp.Rational]:
    raw = [
        (2, 0, 0, sp.Rational(1)),
        (0, 1, 0, sp.Rational(-1)),
    ]
    for (i, j), coefficient in coefficients.items():
        raw.append((j + 2, i + j, 1, -coefficient))
        raw.append((j, i + j + 1, 1, coefficient))
    collapsed: dict[tuple[int, int, int], sp.Rational] = {}
    for x_coordinate, y_coordinate, height, coefficient in raw:
        key = (x_coordinate, y_coordinate, height)
        collapsed[key] = collapsed.get(key, sp.Rational(0)) + coefficient
    return {key: value for key, value in collapsed.items() if value != 0}


def lower_planes(
    points: dict[tuple[int, int, int], sp.Rational],
) -> tuple[tuple[F, F, F], ...]:
    items = tuple(sorted(points))
    planes: set[tuple[F, F, F]] = set()
    for first, second, third in combinations(items, 3):
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
            for x, y, z in items
        ]
        if min(gaps) >= 0:
            planes.add((slope_x, slope_y, constant))
    return tuple(sorted(planes))


def raw_edge_packet(
    vertices: tuple[tuple[int, int], ...],
) -> tuple[int, ...]:
    packet: list[int] = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward_x, inward_y = -dy // length, dx // length
        constant = inward_x * start[0] + inward_y * start[1]
        index = inward_x + inward_y - constant
        # The length-four s=0 edge consists of affine points, not punctures.
        if start[0] == end[0] == 0:
            require((length, index) == (4, 1), "vertical affine edge changed")
            continue
        packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def make_normalized_G(
    X: sp.Symbol,
    T: sp.Symbol,
    Delta: sp.Expr,
    Theta: sp.Expr,
    Phi: sp.Expr,
    eta: sp.Expr,
    zeta: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    P = T + X**2 * T**2
    Y = X * T * P
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    G = sp.expand(
        -X**2 * T / 2
        - 3 * P
        + sp.Rational(8, 3) * P**2
        - sp.Rational(1376, 135) * P**3
        + K * Y**2
        + Phi * P**2 * Y
        + Delta * P**4
        + Theta * P * Y**2
        + eta * P**3 * Y
        + zeta * Y**3
    )
    return G, P, Y, K


def critical_xt(
    G: sp.Expr,
    X: sp.Symbol,
    T: sp.Symbol,
    artifact_power: int,
) -> tuple[sp.Expr, sp.Expr, sp.Poly]:
    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    require(sp.denom(f) == 1, "G_X/T is not polynomial")
    resultant = sp.factor(sp.resultant(f, h, X))
    require(valuation(resultant, T) == artifact_power, "T artifact changed")
    quotient = sp.cancel(resultant / (T**artifact_power * (6 * T + 1) ** 2))
    require(sp.denom(quotient) == 1, "critical quotient is not polynomial")
    return f, h, sp.Poly(quotient, T)


def q_degree_on_wall(
    X: sp.Symbol,
    T: sp.Symbol,
    Delta: sp.Rational,
    Theta: sp.Rational,
    Phi: sp.Rational,
    eta: sp.Rational,
    zeta: sp.Rational,
) -> tuple[int, int, int]:
    G, _, _, _ = make_normalized_G(X, T, Delta, Theta, Phi, eta, zeta)
    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    resultant = sp.factor(sp.resultant(f, h, X))
    t_order = valuation(resultant, T)
    universal_order = 0
    remainder = resultant
    while sp.rem(sp.Poly(remainder, T), sp.Poly(6 * T + 1, T)) == 0:
        remainder = sp.cancel(remainder / (6 * T + 1))
        universal_order += 1
    quotient = sp.Poly(sp.cancel(remainder / T**t_order), T)
    return t_order, universal_order, quotient.degree()


def main() -> None:
    X, T = sp.symbols("X T")
    Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta

    # Complete weight-nine support and its exact anti-diagonal contraction.
    weighted = {
        (i, j)
        for i in range(5)
        for j in range(4)
        if 0 < 2 * i + 3 * j <= 9
    } - {(0, 1), (1, 1)}
    require(
        weighted
        == {
            (1, 0), (2, 0), (3, 0), (0, 2), (2, 1),
            (4, 0), (1, 2), (3, 1), (0, 3),
        },
        "complete M=9 support changed",
    )
    G_symbolic, P, Y, _ = make_normalized_G(
        X, T, Delta, Theta, Phi, eta, -eta
    )
    require(
        sp.factor(P**3 * Y - Y**3 - X * T**2 * P**3) == 0,
        "anti-diagonal top contraction changed",
    )

    # Universal critical pairs, checked with every anti-diagonal coefficient
    # still symbolic.  The target Hessian signs make both pairs Morse.
    f_symbolic = sp.cancel(sp.diff(G_symbolic, X) / T)
    h_symbolic = sp.diff(G_symbolic, T)
    hessian = sp.factor(sp.det(sp.hessian(G_symbolic, (X, T))))
    require(sp.factor(f_symbolic.subs(T, 0)) == -X, "wrong f at T=0")
    require(
        sp.factor(h_symbolic.subs(T, 0) + (X**2 + 6) / 2) == 0,
        "wrong h at T=0",
    )
    require(remainder_at(G_symbolic, sp.Rational(0), X**2 + 6, X, T) == 0,
            "wrong T=0 critical value")
    require(remainder_at(hessian, sp.Rational(0), X**2 + 6, X, T) == 6,
            "T=0 pair is not Morse")
    require(remainder_at(f_symbolic, -sp.Rational(1, 6), X**2 - 6, X, T) == 0,
            "universal f pair changed")
    require(remainder_at(h_symbolic, -sp.Rational(1, 6), X**2 - 6, X, T) == 0,
            "universal h pair changed")
    require(
        remainder_at(G_symbolic - sp.Rational(1, 2), -sp.Rational(1, 6), X**2 - 6, X, T) == 0,
        "wrong universal value",
    )
    require(remainder_at(hessian, -sp.Rational(1, 6), X**2 - 6, X, T) == -6,
            "universal pair is not Morse")

    # Symbolic source-coordinate critical projection on zeta=-eta.
    s_source, p_source = sp.symbols("s_source p_source")
    t_source = p_source - s_source**2
    H_source = sp.expand(
        -3 * p_source
        + sp.Rational(8, 3) * p_source**2
        - sp.Rational(1376, 135) * p_source**3
        + K * s_source**2 * p_source**2
        + Phi * s_source * p_source**3
        + Delta * p_source**4
        + Theta * s_source**2 * p_source**3
        + eta * s_source * p_source**4
        - eta * s_source**3 * p_source**3
    )
    G_source = -s_source**2 / (2 * t_source) + H_source
    require(
        sp.factor(
            sp.cancel(
                G_source.subs({s_source: X * T, p_source: P}) - G_symbolic
            )
        ) == 0,
        "source reconstruction changed",
    )
    A = sp.factor(
        sp.cancel(
            (-s_source * p_source + t_source**2 * sp.diff(H_source, s_source))
            / p_source
        )
    )
    C = sp.expand(s_source**2 + 2 * t_source**2 * sp.diff(H_source, p_source))
    B = sp.factor(sp.cancel((C + s_source * A) / t_source**2))
    require(sp.denom(A).has(s_source, p_source) is False,
            "A acquired a source denominator")
    require(sp.denom(B).has(s_source, p_source) is False,
            "B acquired a source denominator")
    require(
        sp.factor(sp.cancel(t_source**2 * sp.diff(G_source, s_source) - p_source * A)) == 0,
        "first source critical identity changed",
    )
    require(
        sp.factor(
            sp.cancel(2 * t_source**2 * sp.diff(G_source, p_source)
                      - (t_source**2 * B - s_source * A))
        ) == 0,
        "second source critical identity changed",
    )
    source_resultant = sp.factor(sp.resultant(A, B, s_source))
    require(valuation(source_resultant, p_source) == 6,
            "anti source p-artifact changed")
    source_residual = sp.Poly(sp.cancel(source_resultant / p_source**6), p_source)
    source_wall = 4 * Theta * K**2 - 27 * eta**2
    require(source_residual.degree() == 19, "anti source residual is not Q19")
    require(
        sp.factor(source_residual.TC() - 46656 * eta * source_wall) == 0,
        "anti source residual constant changed",
    )
    require(
        sp.factor(source_residual.LC() - 1327104 * eta**5 * (Delta + Theta)**4) == 0,
        "anti source residual leading row changed",
    )

    # Direct (X,T) resultants at the upstream rational control.  The first
    # three rows cross-check THM-4147; the fourth is the new wall result.
    delta0 = sp.Rational(1)
    theta0 = sp.Rational(19, 11)
    phi0 = sp.Rational(11, 7)
    eta0 = sp.Rational(23, 13)
    zeta0 = sp.Rational(29, 17)
    k0 = sp.Rational(5591, 90)
    require(k0 == sp.Rational(2848, 45) - sp.Rational(7, 6) * delta0,
            "control K changed")
    controls = (
        ("P", eta0, sp.Rational(0), 20, 56,
         -sp.Rational(3**15, 2**7) * eta0**7,
         72900 * eta0**5 * k0**4 * theta0**4),
        ("Y", sp.Rational(0), zeta0, 20, 56,
         -sp.Rational(3**15, 2**7) * zeta0**7,
         78732 * delta0**2 * zeta0**5
         * (4 * theta0 * k0**2 - 27 * zeta0**2) ** 2),
        ("B", eta0, zeta0, 21, 56,
         -sp.Rational(3**15, 2**7) * (eta0 + zeta0) ** 7,
         sp.Rational(59049, 4) * eta0**3 * zeta0**2
         * (eta0 + zeta0) ** 2
         * (4 * theta0 * k0**2 - 27 * zeta0**2) ** 2),
        ("A", eta0, -eta0, 19, 42,
         -12288 * (delta0 + theta0) ** 6,
         -1458 * (delta0 + theta0) * eta0**4
         * (4 * theta0 * k0**2 - 27 * eta0**2) ** 2),
    )
    control_rows: list[tuple[str, int, int, str]] = []
    for name, eta_case, zeta_case, degree, artifact, constant, leading in controls:
        G_case, _, _, _ = make_normalized_G(
            X, T, delta0, theta0, phi0, eta_case, zeta_case
        )
        f_case, h_case, residual = critical_xt(G_case, X, T, artifact)
        require(residual.degree() == degree, f"{name}: residual degree changed")
        require(sp.factor(residual.TC() - constant) == 0,
                f"{name}: residual constant changed")
        require(sp.factor(residual.LC() - leading) == 0,
                f"{name}: residual leading row changed")
        require(sp.gcd(residual, residual.diff()).degree() == 0,
                f"{name}: positive control residual is not squarefree")
        require(residual.eval(-sp.Rational(1, 6)) != 0,
                f"{name}: residual meets the universal T=-1/6 fibre")
        if name == "A":
            require((sp.degree(f_case, X), sp.degree(h_case, X)) == (7, 8),
                    "anti critical X-degrees changed")
            require(
                sp.factor(sp.Poly(f_case, X).LC() - 8 * (delta0 + theta0) * T**7) == 0,
                "anti f leading row changed",
            )
            require(
                sp.factor(sp.Poly(h_case, X).LC() - 8 * (delta0 + theta0) * T**7) == 0,
                "anti h leading row changed",
            )
        digest_payload = ",".join(str(coefficient) for coefficient in residual.monic().all_coeffs())
        control_rows.append(
            (name, degree, degree + 4, hashlib.sha256(digest_payload.encode()).hexdigest())
        )

    # Raw B-polygon, repeated top-edge node, and normalized anti packet.
    coefficients = {
        (1, 0): -sp.Rational(3),
        (2, 0): sp.Rational(8, 3),
        (3, 0): -sp.Rational(1376, 135),
        (0, 2): k0,
        (2, 1): phi0,
        (4, 0): delta0,
        (1, 2): theta0,
        (3, 1): eta0,
        (0, 3): -eta0,
    }
    support = valued_support(coefficients)
    polygon = convex_hull([(x, y) for x, y, _ in support])
    require(polygon == ((0, 1), (2, 0), (5, 3), (1, 5), (0, 5)),
            "anti raw polygon changed")
    require(polygon_ledger(polygon) == (31, 11, 11),
            "anti raw Pick ledger changed")
    require(
        lower_planes(support)
        == ((F(0), F(1, 4), F(-1, 4)), (F(1, 9), F(2, 9), F(-2, 9))),
        "anti lower planes changed",
    )
    require(raw_edge_packet(polygon) == (8, 8, 4, 2, 2, 2, 1),
            "anti raw packet changed")

    Q_parameter, z_boundary, a_boundary = sp.symbols(
        "Q_parameter z_boundary a_boundary"
    )
    source_fibre = sp.expand(
        (s_source**2 - p_source) * (1 - Q_parameter * H_source)
        - Q_parameter * s_source**2 / 2
    )
    r_boundary = 1 - a_boundary
    local_equation = sp.factor(
        sp.cancel(
            z_boundary**11
            * source_fibre.subs(
                {s_source: 1 / z_boundary, p_source: r_boundary / z_boundary**2}
            )
        )
    )
    local_poly = sp.Poly(sp.expand(local_equation), a_boundary, z_boundary)
    tangent = sp.factor(
        sum(
            coefficient * a_boundary**monomial[0] * z_boundary**monomial[1]
            for monomial, coefficient in local_poly.terms()
            if sum(monomial) == 2
        )
    )
    require(
        sp.factor(tangent - Q_parameter * a_boundary
                  * (eta * a_boundary - (Delta + Theta) * z_boundary)) == 0,
        "anti boundary tangent cone changed",
    )
    require(
        all(sum(monomial) >= 2 for monomial, _ in local_poly.terms()),
        "anti boundary acquired a lower-order term",
    )
    local_a_derivative = sp.diff(local_equation, a_boundary)
    branch_slopes = (sp.Rational(0), (Delta + Theta) / eta)
    derivative_leads = tuple(
        sp.factor(
            sp.limit(
                local_a_derivative.subs(a_boundary, slope * z_boundary)
                / z_boundary,
                z_boundary,
                0,
            )
        )
        for slope in branch_slopes
    )
    require(
        derivative_leads
        == (-Q_parameter * (Delta + Theta), Q_parameter * (Delta + Theta)),
        "anti node branch derivatives changed",
    )
    # omega=Q*z^7 dz/L_a has order six on each branch, hence index seven.
    normalized_packet = (7, 7, 4, 2, 2, 2, 1)
    normalized_genus = 10
    require(sum(index - 1 for index in normalized_packet) == 2 * normalized_genus - 2,
            "anti normalized RH ledger changed")

    # The nonrational edge is one separable cubic closed point.
    w, q = sp.symbols("w q")
    carrier_polynomial = -eta * w**3 + K * w**2 - (q - sp.Rational(1, 2))
    carrier_discriminant = sp.factor(sp.discriminant(carrier_polynomial, w))
    expected_discriminant = sp.factor(
        (q - sp.Rational(1, 2))
        * (4 * K**3 - 27 * eta**2 * (q - sp.Rational(1, 2)))
    )
    require(sp.factor(carrier_discriminant - expected_discriminant) == 0,
            "anti cubic carrier discriminant changed")

    # Exact monodromy budgets.  Full: n=25,L=23 and packet defect 18.
    # Finite cubic carrier: n=19,L=23,beta=3.
    critical_length = 23
    full_degree = sum(normalized_packet)
    packet_defect = sum(index - 1 for index in normalized_packet)
    full_overlap = full_degree - critical_length
    require((full_degree, packet_defect, full_overlap) == (25, 18, 2),
            "anti full ledger changed")
    require(2 * full_overlap < packet_defect,
            "anti full commutator bound no longer contradicts packet")
    finite_degree = full_degree - 6
    carrier_index = 3
    finite_capacity = 2 * finite_degree - critical_length - 1 + carrier_index
    require((finite_degree, carrier_index, finite_capacity) == (19, 3, 17),
            "anti finite ledger changed")
    require(finite_capacity < finite_degree - 1,
            "anti finite carrier can now merge all sheets")

    # Non-theorem telemetry: first exact critical-resultant drop strata.
    delta_drop = sp.Rational(2)
    k_drop = sp.Rational(2743, 45)
    eta_drop = sp.Rational(7)
    phi_drop = sp.Rational(5)
    theta_critical = sp.factor(27 * eta_drop**2 / (4 * k_drop**2))
    drop_rows = (
        ("anti_sum_zero",)
        + q_degree_on_wall(X, T, delta_drop, -delta_drop, phi_drop, eta_drop, -eta_drop),
        ("anti_critical_infinity",)
        + q_degree_on_wall(X, T, delta_drop, theta_critical, phi_drop, eta_drop, -eta_drop),
    )
    require(drop_rows[0][1:] == (30, 2, 17), "anti sum-zero drop changed")
    require(drop_rows[1][1:] == (42, 2, 18), "anti endpoint drop changed")

    semantic_rows = (
        "normalized=G_base+eta*P^3Y+zeta*Y^3",
        "anti=zeta=-eta;top=eta*X*T^2*P^3",
        "source_resultant=p^6*R19",
        "source_endpoints=46656*eta*(4ThetaK^2-27eta^2);1327104*eta^5*(Delta+Theta)^4",
        "xt_resultant=T^42*(6T+1)^2*Q19",
        "xt_endpoints=-12288*(Delta+Theta)^6;-1458*(Delta+Theta)*eta^4*(4ThetaK^2-27eta^2)^2",
        "critical_length=23",
        "open=eta*Delta*(Delta+Theta)*(4ThetaK^2-27eta^2)*Disc(Q19)*Q19(-1/6)!=0",
        "raw=genus11;packet=8,8,4,2,2,2,1",
        "normalized=genus10;packet=7,7,4,2,2,2,1",
        "responses=full25;finite19;carrier_index3",
        "monodromy=full:4<18;finite:17<18",
    )
    semantic_digest = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("THM-4147 GENERIC CROSS-CHECK + ANTI-DIAGONAL EXTENSION")
    print(f"checks={CHECKS}")
    print("complete_weight9_support=" + repr(tuple(sorted(weighted))))
    print("forced_K=2848/45-(7/6)Delta")
    for name, degree, length, digest in control_rows:
        print(f"control={name};residual_degree={degree};critical_length={length};sha256={digest}")
    print("anti_source_resultant=p^6*R19")
    print("anti_source_R19_constant=46656*eta*(4*Theta*K^2-27*eta^2)")
    print("anti_source_R19_leading=1327104*eta^5*(Delta+Theta)^4")
    print("anti_xt_resultant=T^42*(6T+1)^2*Q19")
    print("anti_xt_Q19_constant=-12288*(Delta+Theta)^6")
    print("anti_xt_Q19_leading=-1458*(Delta+Theta)*eta^4*(4*Theta*K^2-27*eta^2)^2")
    print("anti_open_gates=eta*Delta*(Delta+Theta)*(4*Theta*K^2-27*eta^2)!=0")
    print("anti_critical_open=Disc_T(Q19)*Q19(-1/6)!=0")
    print("universal_pairs=T0:X^2=-6,value0,hess6;T-1/6:X^2=6,value1/2,hess-6")
    print("anti_raw_polygon=" + repr(polygon))
    print("anti_lower_planes=" + repr(lower_planes(support)))
    print("anti_tangent=Q*a*(eta*a-(Delta+Theta)*z)")
    print("anti_node_branches=a~0*z,a~((Delta+Theta)/eta)*z;indices=7,7")
    print("anti_normalized=genus10;packet=7,7,4,2,2,2,1;defect=18")
    print("anti_responses=full:n25,finite:n19,beta3")
    print("anti_full=commutator_index<=4<18")
    print("anti_finite=merger_capacity17<n-1=18")
    print("drop_telemetry=" + repr(drop_rows))
    print("generic_crosscheck=THM4147_PYB_endpoints_and_lengths:PASS")
    print("semantic_sha256=" + semantic_digest)
    print("verdict=PASS")


if __name__ == "__main__":
    main()
