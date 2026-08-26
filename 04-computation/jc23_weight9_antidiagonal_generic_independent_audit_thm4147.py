#!/usr/bin/env python3
"""Clean-room exact audit of THM-4147's generic anti-diagonal chamber.

This script imports no THM-4147 implementation.  It starts with the rational
source in independent ``(s,p)`` coordinates, uses the disjoint exact control

    Delta=-64/105, K=64, Phi=1, Theta=2, eta=3, zeta=-3,

and recomputes both critical eliminants.  It then constructs the boundary
chart directly from the cleared generic fibre, verifies the ordinary node,
and obtains the two normalized e=7 branches from the order of the residue
differential.  Finally it reconstructs the Newton/Pick ledger, the cubic
carrier, and the full/finite permutation budgets.

No ``assert`` occurs, so all gates remain active under ``python -O``.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd

import sympy as sy


CHECKS = 0


def require(condition: bool, label: str) -> None:
    """Count and enforce an audit gate even in optimized Python."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def order_at_zero(expression: sy.Expr, variable: sy.Symbol) -> int:
    """Exact order of a nonzero polynomial at variable=0."""
    polynomial = sy.Poly(sy.expand(expression), variable)
    terms = [(monomial[0], coefficient)
             for monomial, coefficient in polynomial.terms()
             if coefficient != 0]
    require(bool(terms), f"zero polynomial has no {variable}-order")
    return min(exponent for exponent, _ in terms)


def factor_multiplicity(
    expression: sy.Expr,
    factor: sy.Expr,
    variable: sy.Symbol,
) -> tuple[int, sy.Expr]:
    """Return the exact multiplicity of factor and the remaining quotient."""
    multiplicity = 0
    remainder = sy.Poly(sy.expand(expression), variable)
    divisor = sy.Poly(factor, variable)
    while True:
        quotient, residue = sy.div(remainder, divisor)
        if not residue.is_zero:
            return multiplicity, sy.factor(remainder.as_expr())
        multiplicity += 1
        remainder = quotient


def monic_digest(polynomial: sy.Poly) -> str:
    """Stable digest of a univariate rational polynomial."""
    payload = ",".join(str(coefficient)
                       for coefficient in polynomial.monic().all_coeffs())
    return sha256(payload.encode("utf-8")).hexdigest()


def convex_hull(points: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    """Counterclockwise strict convex hull of lattice points."""
    ordered = sorted(set(points))

    def cross(
        origin: tuple[int, int],
        left: tuple[int, int],
        right: tuple[int, int],
    ) -> int:
        return ((left[0] - origin[0]) * (right[1] - origin[1])
                - (left[1] - origin[1]) * (right[0] - origin[0]))

    lower: list[tuple[int, int]] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(
    vertices: tuple[tuple[int, int], ...],
) -> tuple[int, int, int, tuple[int, ...]]:
    """Return twice-area, boundary count, Pick genus, and raw packet."""
    signed_area_twice = sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    )
    require(signed_area_twice > 0, "Newton polygon orientation changed")
    boundary = 0
    packet: list[int] = []
    for start, stop in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = stop[0] - start[0], stop[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        require(length > 0, "zero Newton edge")
        boundary += length
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        ramification_index = inward[0] + inward[1] - constant
        if start[0] == stop[0] == 0:
            require((length, ramification_index) == (4, 1),
                    "affine vertical edge changed")
            continue
        packet.extend([ramification_index] * length)
    require((signed_area_twice - boundary + 2) % 2 == 0,
            "Pick parity changed")
    genus = (signed_area_twice - boundary + 2) // 2
    return signed_area_twice, boundary, genus, tuple(sorted(packet, reverse=True))


def valued_support(
    fibre: sy.Expr,
    s: sy.Symbol,
    p: sy.Symbol,
    q_inverse: sy.Symbol,
) -> tuple[tuple[int, int, int, sy.Expr], ...]:
    """Combine coefficients before taking their Q-adic order."""
    rows: list[tuple[int, int, int, sy.Expr]] = []
    for (s_power, p_power), coefficient in sy.Poly(
        sy.expand(fibre), s, p
    ).terms():
        coefficient_poly = sy.Poly(coefficient, q_inverse)
        nonzero = [(power[0], value)
                   for power, value in coefficient_poly.terms()
                   if value != 0]
        require(bool(nonzero), "zero projected coefficient survived")
        height = min(power for power, _ in nonzero)
        rows.append((s_power, p_power, height,
                     sy.factor(coefficient_poly.nth(height))))
    return tuple(sorted(rows, key=str))


def make_source(
    s: sy.Symbol,
    p: sy.Symbol,
    delta: sy.Expr,
    theta: sy.Expr,
    phi: sy.Expr,
    eta: sy.Expr,
) -> tuple[sy.Expr, sy.Expr, sy.Expr, sy.Expr]:
    """Construct H and G on zeta=-eta without using the primary audit."""
    t = p - s**2
    kappa = sy.Rational(2848, 45) - sy.Rational(7, 6) * delta
    h = sy.expand(
        -3 * p
        + sy.Rational(8, 3) * p**2
        - sy.Rational(1376, 135) * p**3
        + kappa * s**2 * p**2
        + phi * s * p**3
        + delta * p**4
        + theta * s**2 * p**3
        + eta * s * p**4
        - eta * s**3 * p**3
    )
    g = -s**2 / (2 * t) + h
    return g, h, t, kappa


def main() -> None:
    # A rational point disjoint from the primary anti-diagonal control.
    delta = -sy.Rational(64, 105)
    theta = sy.Rational(2)
    phi = sy.Rational(1)
    eta = sy.Rational(3)

    s, p, q_inverse = sy.symbols("s p Q")
    source_g, source_h, t, kappa = make_source(
        s, p, delta, theta, phi, eta
    )
    require(kappa == 64, "forced K/Delta relation changed")
    carrier_gate = sy.factor(4 * theta * kappa**2 - 27 * eta**2)
    require(delta != 0, "Delta gate failed")
    require(eta != 0, "eta gate failed")
    require(delta + theta != 0, "ordinary-node gate failed")
    require(carrier_gate != 0, "cubic endpoint gate failed")

    # Source-coordinate critical equations are derived from the rational
    # function, then saturated only by the displayed coordinate factors.
    first_numerator = -s * p + t**2 * sy.diff(source_h, s)
    first_sat = sy.cancel(first_numerator / p)
    second_numerator = s**2 + 2 * t**2 * sy.diff(source_h, p)
    second_sat = sy.cancel((second_numerator + s * first_sat) / t**2)
    require(sy.denom(first_sat) == 1,
            "source first critical equation is not polynomial")
    require(not sy.denom(second_sat).has(s, p),
            "source second critical equation has a variable denominator")
    require(
        sy.factor(sy.cancel(t**2 * sy.diff(source_g, s)
                            - p * first_sat)) == 0,
        "source first derivative identity failed",
    )
    require(
        sy.factor(sy.cancel(2 * t**2 * sy.diff(source_g, p)
                            - (t**2 * second_sat - s * first_sat))) == 0,
        "source second derivative identity failed",
    )
    source_resultant = sy.factor(sy.resultant(first_sat, second_sat, s))
    source_artifact = order_at_zero(source_resultant, p)
    require(source_artifact == 6, "source p-artifact is not six")
    source_residual_expr = sy.cancel(source_resultant / p**source_artifact)
    require(sy.denom(source_residual_expr) == 1,
            "source residual has a variable denominator")
    source_residual = sy.Poly(source_residual_expr, p)
    require(source_residual.degree() == 19, "source residual is not R19")
    require(sy.gcd(source_residual, source_residual.diff()).degree() == 0,
            "source R19 is not squarefree")
    require(source_residual.TC() != 0 and source_residual.LC() != 0,
            "source R19 endpoint failed")
    require(
        sy.factor(source_residual.TC()
                  - 46656 * eta * carrier_gate) == 0,
        "source R19 constant changed",
    )
    require(
        sy.factor(source_residual.LC()
                  - 1327104 * eta**5 * (delta + theta)**4) == 0,
        "source R19 leading row changed",
    )

    # The collapsed p=0 stratum restores the universal pair at T=-1/6.
    require(sy.factor(sy.diff(source_g, s).subs(p, 0)) == 0,
            "collapsed p=0 first derivative changed")
    require(
        sy.factor(sy.diff(source_g, p).subs(p, 0)
                  - (1 - 6 * s**2) / (2 * s**2)) == 0,
        "collapsed p=0 second derivative changed",
    )
    require(
        sy.factor(source_g.subs({p: 0, s**2: sy.Rational(1, 6)})
                  - sy.Rational(1, 2)) == 0,
        "collapsed p=0 critical value changed",
    )

    # Independently rebuild the polynomial in (X,T) by substituting into
    # the rational source.  This does not reuse an expanded primary formula.
    X, T = sy.symbols("X T")
    P = T + X**2 * T**2
    xt_g = sy.cancel(source_g.subs({s: X * T, p: P}))
    require(sy.denom(xt_g) == 1, "source substitution did not polynomialize")
    xt_g = sy.expand(xt_g)
    gx_over_t = sy.cancel(sy.diff(xt_g, X) / T)
    gt = sy.diff(xt_g, T)
    require(sy.denom(gx_over_t) == 1, "G_X/T is not polynomial")
    require((sy.degree(gx_over_t, X), sy.degree(gt, X)) == (7, 8),
            "anti-diagonal critical X degrees changed")
    expected_x_lead = 8 * (delta + theta) * T**7
    require(sy.factor(sy.Poly(gx_over_t, X).LC() - expected_x_lead) == 0,
            "G_X/T leading X row changed")
    require(sy.factor(sy.Poly(gt, X).LC() - expected_x_lead) == 0,
            "G_T leading X row changed")

    xt_resultant = sy.factor(sy.resultant(gx_over_t, gt, X))
    xt_artifact = order_at_zero(xt_resultant, T)
    require(xt_artifact == 42, "XT T-artifact is not forty-two")
    after_t = sy.cancel(xt_resultant / T**xt_artifact)
    require(sy.denom(after_t) == 1, "XT T-saturation failed")
    universal_multiplicity, xt_residual_expr = factor_multiplicity(
        after_t, 6 * T + 1, T
    )
    require(universal_multiplicity == 2,
            "universal (6T+1) multiplicity is not two")
    xt_residual = sy.Poly(xt_residual_expr, T)
    require(xt_residual.degree() == 19, "XT residual is not Q19")
    require(sy.gcd(xt_residual, xt_residual.diff()).degree() == 0,
            "XT Q19 is not squarefree")
    require(xt_residual.eval(-sy.Rational(1, 6)) != 0,
            "XT Q19 meets T=-1/6")
    require(xt_residual.TC() != 0 and xt_residual.LC() != 0,
            "XT Q19 endpoint failed")
    require(
        sy.factor(xt_residual.TC()
                  + 12288 * (delta + theta)**6) == 0,
        "XT Q19 constant changed",
    )
    require(
        sy.factor(xt_residual.LC()
                  + 1458 * (delta + theta) * eta**4 * carrier_gate**2) == 0,
        "XT Q19 leading row changed",
    )

    # Restore both pairs discarded by the eliminant saturations and check
    # they are Morse critical points with the advertised values.
    hessian = sy.factor(sy.det(sy.hessian(xt_g, (X, T))))
    require(sy.factor(gx_over_t.subs(T, 0) + X) == 0,
            "T=0 first critical row changed")
    require(sy.factor(gt.subs(T, 0) + (X**2 + 6) / 2) == 0,
            "T=0 second critical row changed")
    require(sy.rem(sy.Poly(xt_g.subs(T, 0), X),
                   sy.Poly(X**2 + 6, X)).is_zero,
            "T=0 pair has wrong value")
    require(
        sy.rem(sy.Poly(hessian.subs(T, 0) - 6, X),
               sy.Poly(X**2 + 6, X)).is_zero,
        "T=0 pair is not Morse",
    )
    special_first = sy.Poly(gx_over_t.subs(T, -sy.Rational(1, 6)), X)
    special_second = sy.Poly(gt.subs(T, -sy.Rational(1, 6)), X)
    special_gcd = sy.gcd(special_first, special_second).monic()
    require(special_gcd.as_expr() == X**2 - 6,
            "T=-1/6 universal pair changed")
    require(
        sy.rem(sy.Poly(xt_g.subs(T, -sy.Rational(1, 6))
                       - sy.Rational(1, 2), X), special_gcd).is_zero,
        "T=-1/6 pair has wrong value",
    )
    require(
        sy.rem(sy.Poly(hessian.subs(T, -sy.Rational(1, 6)) + 6, X),
               special_gcd).is_zero,
        "T=-1/6 pair is not Morse",
    )
    critical_length_source = source_residual.degree() + 2 + 2
    critical_length_xt = xt_residual.degree() + 2 + 2
    require((critical_length_source, critical_length_xt) == (23, 23),
            "source/XT critical lengths disagree")

    # Build the generic cleared fibre and its valued Newton support directly.
    cleared_fibre = sy.expand(
        (s**2 - p) * (1 - q_inverse * source_h)
        - q_inverse * s**2 / 2
    )
    support = valued_support(cleared_fibre, s, p, q_inverse)
    polygon = convex_hull(tuple((row[0], row[1]) for row in support))
    expected_polygon = ((0, 1), (2, 0), (5, 3), (1, 5), (0, 5))
    require(polygon == expected_polygon, "raw anti Newton polygon changed")
    area_twice, boundary, raw_genus, raw_packet = polygon_ledger(polygon)
    require((area_twice, boundary, raw_genus) == (31, 11, 11),
            "raw anti Pick ledger changed")
    require(raw_packet == (8, 8, 4, 2, 2, 2, 1),
            "raw anti boundary packet changed")

    # The top length-two edge has a square initial form; this is the collision
    # whose local normalization the next chart resolves.
    cleared_poly = sy.Poly(cleared_fibre, s, p)
    top_coefficients = tuple(
        sy.factor(cleared_poly.coeff_monomial(s**(5 - 2 * index)
                                              * p**(3 + index)))
        for index in range(3)
    )
    require(top_coefficients == (
        q_inverse * eta, -2 * q_inverse * eta, q_inverse * eta
    ), "top edge coefficients are not a square")
    edge_coordinate = sy.symbols("edge_coordinate")
    top_initial = sy.factor(sum(
        top_coefficients[index] * edge_coordinate**index
        for index in range(3)
    ))
    require(top_initial == eta * q_inverse * (edge_coordinate - 1)**2,
            "top edge repeated point changed")

    # Reconstruct the boundary chart with symbolic live coefficients.  The
    # tangent calculation is generic, while the disjoint control witnesses
    # the nonempty open set.
    delta_s, theta_s, phi_s, eta_s = sy.symbols(
        "Delta_s Theta_s Phi_s eta_s"
    )
    ss, pp, z, a, Q = sy.symbols("ss pp z a Q")
    _, h_symbolic, _, _ = make_source(
        ss, pp, delta_s, theta_s, phi_s, eta_s
    )
    symbolic_fibre = sy.expand(
        (ss**2 - pp) * (1 - Q * h_symbolic) - Q * ss**2 / 2
    )
    local_equation = sy.cancel(
        z**11 * symbolic_fibre.subs({ss: 1 / z, pp: (1 - a) / z**2})
    )
    require(sy.denom(local_equation) == 1,
            "boundary chart retained a denominator")
    local_equation = sy.expand(local_equation)
    local_terms = sy.Poly(local_equation, a, z).terms()
    require(all(sum(monomial) >= 2 for monomial, _ in local_terms),
            "boundary point has order below two")
    tangent = sy.factor(sum(
        coefficient * a**monomial[0] * z**monomial[1]
        for monomial, coefficient in local_terms
        if sum(monomial) == 2
    ))
    expected_tangent = Q * a * (
        eta_s * a - (delta_s + theta_s) * z
    )
    require(sy.factor(tangent - expected_tangent) == 0,
            "ordinary-node tangent cone changed")

    # Derive the residue differential rather than reading e from the raw
    # polygon.  Since L=z^11 F and p=(1-a)/z^2,
    #       L_a=-z^9 F_p,  ds/F_p=z^7 dz/L_a.
    local_a_derivative = sy.diff(local_equation, a)
    substituted_fibre_p = sy.diff(symbolic_fibre, pp).subs(
        {ss: 1 / z, pp: (1 - a) / z**2}
    )
    require(
        sy.factor(sy.cancel(local_a_derivative
                            + z**9 * substituted_fibre_p)) == 0,
        "local residue-differential identity changed",
    )
    branch_slopes = (sy.Rational(0),
                     (delta_s + theta_s) / eta_s)
    branch_derivative_leads = tuple(
        sy.factor(sy.limit(
            local_a_derivative.subs(a, slope * z) / z,
            z,
            0,
        ))
        for slope in branch_slopes
    )
    require(branch_derivative_leads == (
        -Q * (delta_s + theta_s), Q * (delta_s + theta_s)
    ), "the two node branches have wrong differential leading terms")
    differential_orders = tuple(7 - 1 for _ in branch_slopes)
    branch_indices = tuple(order + 1 for order in differential_orders)
    require((differential_orders, branch_indices) == ((6, 6), (7, 7)),
            "normalized node branch indices changed")

    node_delta = 1
    normalized_genus = raw_genus - node_delta
    normalized_packet = tuple(sorted(
        branch_indices + tuple(index for index in raw_packet if index != 8),
        reverse=True,
    ))
    require(normalized_genus == 10, "normalized genus is not ten")
    require(normalized_packet == (7, 7, 4, 2, 2, 2, 1),
            "normalized packet changed")
    packet_defect = sum(index - 1 for index in normalized_packet)
    require(packet_defect == 2 * normalized_genus - 2 == 18,
            "normalized Riemann-Hurwitz ledger changed")

    # The only nonrational boundary point is a degree-three carrier.
    w, q = sy.symbols("w q")
    carrier = -eta * w**3 + kappa * w**2 - (q - sy.Rational(1, 2))
    carrier_discriminant = sy.factor(sy.discriminant(carrier, w))
    expected_carrier_discriminant = sy.factor(
        (q - sy.Rational(1, 2))
        * (4 * kappa**3 - 27 * eta**2 * (q - sy.Rational(1, 2)))
    )
    require(carrier_discriminant == expected_carrier_discriminant,
            "cubic carrier discriminant changed")
    carrier_over_q = sy.Poly(carrier, w, domain=sy.QQ.frac_field(q))
    require(carrier_over_q.degree() == 3 and carrier_over_q.is_irreducible,
            "carrier is not one cubic closed point over QQ(q)")

    # Permutation arithmetic relative to the inherited fixed-sheet,
    # transitivity, support-index, and commutator-overlap inequalities.
    critical_length = critical_length_xt
    full_degree = sum(normalized_packet)
    full_overlap_bound = full_degree - critical_length
    full_commutator_bound = 2 * full_overlap_bound
    require((full_degree, critical_length, full_overlap_bound,
             full_commutator_bound, packet_defect) == (25, 23, 2, 4, 18),
            "full n=25 permutation ledger changed")
    require(full_commutator_bound < packet_defect,
            "full commutator budget no longer contradicts the packet")

    carrier_meridians = 3
    carrier_removed_degree = 2 * carrier_meridians
    finite_degree = full_degree - carrier_removed_degree
    finite_capacity = (
        2 * finite_degree - critical_length - 1 + carrier_meridians
    )
    require((finite_degree, carrier_meridians, finite_capacity,
             finite_degree - 1) == (19, 3, 17, 18),
            "finite n=19 permutation ledger changed")
    require(finite_capacity < finite_degree - 1,
            "finite carrier budget can now be transitive")

    semantic_rows = (
        "control=Delta:-64/105,K:64,Phi:1,Theta:2,eta:3,zeta:-3",
        "open=Delta*eta*(Delta+Theta)*(4ThetaK^2-27eta^2)!=0",
        "source=p^6*R19;squarefree;length23",
        "xt=T^42*(6T+1)^2*Q19;squarefree;Q19(-1/6)!=0;length23",
        "raw=polygon31,11,11;packet8,8,4,2,2,2,1",
        "node=tangent:Q*a*(eta*a-(Delta+Theta)z);delta1",
        "differential=z^7*dz/L_a;orders6,6;indices7,7",
        "normalized=genus10;packet7,7,4,2,2,2,1;defect18",
        "carrier=irreducible cubic;beta3",
        "full=n25,L23,commutator4<18",
        "finite=n19,L23,beta3,capacity17<18",
    )
    semantic_digest = sha256("\n".join(semantic_rows).encode("utf-8")).hexdigest()

    print("THM-4147 ANTI-DIAGONAL GENERIC CLEAN-ROOM AUDIT")
    print(f"checks={CHECKS}")
    print("control=Delta:-64/105,K:64,Phi:1,Theta:2,eta:3,zeta:-3")
    print("open_gates=Delta:-64/105;eta:3;Delta+Theta:146/105;"
          f"4ThetaK^2-27eta^2:{carrier_gate};PASS")
    print("source_critical=p^6*R19;degree19;squarefree;length23;PASS")
    print(f"source_R19_sha256={monic_digest(source_residual)}")
    print("xt_critical=T^42*(6T+1)^2*Q19;degree19;"
          "squarefree;Q19(-1/6)!=0;length23;PASS")
    print(f"xt_Q19_sha256={monic_digest(xt_residual)}")
    print("universal_pairs=T0:X^2=-6,value0,hessian6;"
          "T-1/6:X^2=6,value1/2,hessian-6;PASS")
    print(f"raw_polygon={polygon};Pick={(area_twice, boundary, raw_genus)};"
          f"packet={raw_packet}")
    print("top_edge_initial=3*Q*(edge_coordinate-1)^2")
    print("ordinary_node_tangent=Q*a*(eta*a-(Delta+Theta)*z);delta=1;PASS")
    print("branch_differential=Q*z^7*dz/L_a;orders=(6,6);indices=(7,7);PASS")
    print(f"normalized=genus:{normalized_genus};packet:{normalized_packet};"
          f"defect:{packet_defect}")
    print("cubic_carrier=-3*w^3+64*w^2-(q-1/2);"
          "irreducible_over_QQ(q);separable;beta=3;PASS")
    print("full_budget=n25,L23,overlap<=2,commutator_index<=4<defect18;PASS")
    print("finite_budget=n19,L23,beta3,merger_capacity17<n-1=18;PASS")
    print(f"semantic_sha256={semantic_digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
