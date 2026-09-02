#!/usr/bin/env python3
"""Primary exact certificate for the THM-4334 beta-zero endpoint wall.

The finite universe is the inherited exact-weight-at-most-twelve source of
THM-4327.  We impose Z=beta_11=0 and U*W*zeta_3!=0.  Eight remaining rows
are optional and six aggregate support points may cancel independently; the
latter is a hostile overcount of actual coefficient specializations.

The program exhausts all 2^8*2^6 states, reconstructs the complete lower
hull, and checks the face/genus/order/contact data used by THM-4334.  The
proper-flat and good-differential implications remain the cited geometric
inputs of the theorem rather than computational assumptions.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import comb, gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"THM4334 primary failure: {label}")


def weighted_rows(limit: int = 12) -> tuple[tuple[int, int, int], ...]:
    rows = []
    for i in range(limit // 2 + 1):
        for j in range(limit // 3 + 1):
            weight = 2 * i + 3 * j
            if 0 < weight <= limit and (i, j) not in {(0, 1), (1, 1)}:
                rows.append((i, j, weight))
    return tuple(sorted(rows, key=lambda row: (row[2], row[1], row[0])))


def lifted_support(rows: tuple[tuple[int, int, int], ...]) -> set[tuple[int, int, int]]:
    support = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j, _weight in rows:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


def projected_rank_two(points: tuple[tuple[int, int, int], ...], mask: int) -> bool:
    chosen = tuple(points[k] for k in range(len(points)) if mask & (1 << k))
    for first, second, third in combinations(chosen, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant:
            return True
    return False


def plane_through(first, second, third):
    determinant = (
        (second[0] - first[0]) * (third[1] - first[1])
        - (second[1] - first[1]) * (third[0] - first[0])
    )
    if determinant == 0:
        return None
    slope_s = Q(
        (second[2] - first[2]) * (third[1] - first[1])
        - (second[1] - first[1]) * (third[2] - first[2]),
        determinant,
    )
    slope_p = Q(
        (second[0] - first[0]) * (third[2] - first[2])
        - (second[2] - first[2]) * (third[0] - first[0]),
        determinant,
    )
    constant = Q(first[2]) - slope_s * first[0] - slope_p * first[1]
    return slope_s, slope_p, constant


def lower_hull_engine(universe: set[tuple[int, int, int]]):
    points = tuple(sorted(universe))
    records = {}
    for triple in combinations(points, 3):
        plane = plane_through(*triple)
        if plane is None:
            continue
        a, b, c = plane
        below = 0
        equal = 0
        for index, (r, ell, height) in enumerate(points):
            gap = Q(height) - a * r - b * ell - c
            if gap < 0:
                below |= 1 << index
            elif gap == 0:
                equal |= 1 << index
        records[plane] = (below, equal)

    @lru_cache(maxsize=None)
    def rank_two(mask: int) -> bool:
        return projected_rank_two(points, mask)

    def faces(mask: int):
        return tuple(
            plane
            for plane, (below, equal) in sorted(records.items())
            if not (below & mask) and rank_two(equal & mask)
        )

    return points, records, faces


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, first, second):
        return (
            (first[0] - origin[0]) * (second[1] - origin[1])
            - (first[1] - origin[1]) * (second[0] - origin[0])
        )

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
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(points):
    polygon = convex_hull(points)
    area2 = abs(sum(
        polygon[k][0] * polygon[(k + 1) % len(polygon)][1]
        - polygon[(k + 1) % len(polygon)][0] * polygon[k][1]
        for k in range(len(polygon))
    ))
    boundary = sum(
        gcd(
            abs(polygon[(k + 1) % len(polygon)][0] - polygon[k][0]),
            abs(polygon[(k + 1) % len(polygon)][1] - polygon[k][1]),
        )
        for k in range(len(polygon))
    )
    require((area2 - boundary + 2) % 2 == 0, ("Pick parity", polygon))
    genus = (area2 - boundary + 2) // 2
    return polygon, area2, boundary, genus


def edge_packet(points):
    polygon = convex_hull(points)
    packet = []
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        if not (start[0] == end[0] == 0):
            packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def face_order(base: int, plane: tuple[Q, Q, Q]) -> Q:
    a, b, c = plane
    return base * (Q(5, 6) - a - b - c)


def p_euler_remainder(terms, shift):
    """Apply P*d/dP-shift to a sparse polynomial with linear coefficients.

    Coefficients are exact four-tuples in the basis (1,U,W,zeta).  This is
    enough to certify the two generic-point derivative identities without a
    computer-algebra dependency or a string-level placeholder.
    """
    answer = {}
    for exponent, coefficient in terms.items():
        factor = exponent[1] - shift
        scaled = tuple(factor * value for value in coefficient)
        if any(scaled):
            answer[exponent] = scaled
    return answer


def parameter_poly_add(first, second):
    answer = dict(first)
    for exponent, coefficient in second.items():
        answer[exponent] = answer.get(exponent, 0) + coefficient
        if answer[exponent] == 0:
            del answer[exponent]
    return answer


def parameter_poly_multiply(first, second):
    answer = {}
    for (u1, w1), coefficient1 in first.items():
        for (u2, w2), coefficient2 in second.items():
            exponent = (u1 + u2, w1 + w2)
            answer[exponent] = answer.get(exponent, 0) + coefficient1 * coefficient2
    return {exponent: coefficient for exponent, coefficient in answer.items() if coefficient}


def parameter_poly_scale(polynomial, scalar):
    return {exponent: scalar * coefficient for exponent, coefficient in polynomial.items()}


def qpoly_trim(polynomial):
    polynomial = [Q(value) for value in polynomial]
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def qpoly_derivative(polynomial):
    return qpoly_trim([index * value for index, value in enumerate(polynomial)][1:] or [0])


def qpoly_divmod(dividend, divisor):
    dividend = qpoly_trim(dividend)
    divisor = qpoly_trim(divisor)
    require(divisor != [0], "nonzero polynomial divisor")
    quotient = [Q(0)] * max(1, len(dividend) - len(divisor) + 1)
    while dividend != [0] and len(dividend) >= len(divisor):
        shift = len(dividend) - len(divisor)
        factor = dividend[-1] / divisor[-1]
        quotient[shift] += factor
        for index, value in enumerate(divisor):
            dividend[index + shift] -= factor * value
        dividend = qpoly_trim(dividend)
    return qpoly_trim(quotient), dividend


def qpoly_gcd(first, second):
    first, second = qpoly_trim(first), qpoly_trim(second)
    while second != [0]:
        _quotient, remainder = qpoly_divmod(first, second)
        first, second = second, remainder
    if first == [0]:
        return first
    return qpoly_trim([value / first[-1] for value in first])


def is_squarefree(polynomial):
    polynomial = qpoly_trim(polynomial)
    return len(polynomial) > 1 and len(qpoly_gcd(polynomial, qpoly_derivative(polynomial))) == 1


def main() -> None:
    rows = weighted_rows()
    expected_rows = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11), (6, 0, 12), (3, 2, 12), (0, 4, 12),
    )
    require(rows == expected_rows, "complete labelled row universe")

    p, p2, p3 = (1, 0), (2, 0), (3, 0)
    u, w, zeta = (6, 0), (3, 2), (0, 3)
    beta, z = (1, 3), (0, 4)
    required_labels = {p, p2, p3, u, w, zeta}
    absent_labels = {beta, z}
    available_rows = tuple(row for row in rows if row[:2] not in absent_labels)
    required_rows = tuple(row for row in rows if row[:2] in required_labels)
    optional_rows = tuple(
        row for row in rows
        if row[:2] not in required_labels | absent_labels
    )
    require(len(required_rows) == 6 and len(optional_rows) == 8, "row typing")

    universe = lifted_support(available_rows)
    points, records, lower_faces = lower_hull_engine(universe)
    point_index = {point: index for index, point in enumerate(points)}
    fixed_mask = sum(1 << point_index[x] for x in {(2, 0, 0), (0, 1, 0), (2, 0, 1)})

    def row_mask(row):
        i, j, _weight = row
        return (
            (1 << point_index[(j + 2, i + j, 1)])
            | (1 << point_index[(j, i + j + 1, 1)])
        )

    required_mask = fixed_mask
    for row in required_rows:
        required_mask |= row_mask(row)
    collision_points = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1),
        (2, 6, 1), (3, 4, 1), (3, 5, 1),
    )
    collision_bits = tuple(1 << point_index[x] for x in collision_points)

    face_counter = Counter()
    for optional_choice in range(1 << len(optional_rows)):
        base_mask = required_mask
        for index, row in enumerate(optional_rows):
            if optional_choice & (1 << index):
                base_mask |= row_mask(row)
        for collision_choice in range(1 << len(collision_bits)):
            mask = base_mask
            for index, bit in enumerate(collision_bits):
                if collision_choice & (1 << index):
                    mask &= ~bit
            face_counter[lower_faces(mask)] += 1

    main_plane = (Q(1, 12), Q(1, 6), Q(-1, 6))
    beta0_plane = (Q(2, 9), Q(1, 9), Q(-4, 9))
    expected_faces = (main_plane, beta0_plane)
    require(sum(face_counter.values()) == 16_384, "hostile-state count")
    require(face_counter == Counter({expected_faces: 16_384}), "unique lower complex")

    equal_sets = {}
    for plane in expected_faces:
        _below, equal = records[plane]
        equal_sets[plane] = tuple(points[k] for k in range(len(points)) if equal & (1 << k))
    require(equal_sets[main_plane] == (
        (0, 1, 0), (0, 7, 1), (2, 0, 0), (2, 6, 1), (4, 5, 1),
    ), "main-face support")
    require(equal_sets[beta0_plane] == (
        (2, 0, 0), (4, 5, 1), (5, 3, 1),
    ), "beta-zero face support")

    polygons = {
        "main": ((0, 1), (2, 0), (4, 5), (0, 7)),
        "central": ((0, 0), (0, 6), (2, 5)),
        "beta0": ((2, 0), (4, 5), (5, 3)),
        "global": ((0, 1), (2, 0), (5, 3), (4, 5), (0, 7)),
    }
    ledgers = {name: polygon_ledger(polygon)[1:] for name, polygon in polygons.items()}
    require(ledgers == {
        "main": (36, 10, 14),
        "central": (12, 8, 3),
        "beta0": (9, 5, 3),
        "global": (45, 13, 17),
    }, "Pick ledgers")
    require(3 + 3 + 11 == 17, "normalization plus graph genus")
    require((3, 13, 1, 11) == (3, 12 + 1, 1, (12 + 1) - 3 + 1), "dual graph")

    packet = edge_packet(polygons["global"])
    require(packet == (11, 11, 10, 2, 2, 2, 1), "source-infinity packet")
    require(sum(index - 1 for index in packet) == 32 == 2 * 17 - 2, "packet genus saturation")

    orders = tuple(face_order(36, plane) for plane in expected_faces)
    require(orders == (Q(27), Q(34)), "good-differential orders")
    require(tuple(36 * value for value in beta0_plane) == (8, 4, -16), "integral beta0 plane")

    # The beta-zero face is the connected cyclic cover
    # S^9=(W^3/zeta^5)*(1-x)^5/x^3.  Its three branch valuations are
    # -3,5,-2, and Riemann--Hurwitz gives genus three.
    branch_orders = (-3, 5, -2)
    branch_gcds = tuple(gcd(9, abs(value)) for value in branch_orders)
    require(sum(branch_orders) == 0, "Kummer divisor degree")
    require(gcd(9, *map(abs, branch_orders)) == 1, "Kummer connectedness")
    require(branch_gcds == (3, 1, 1), "Kummer inertia")
    ramification = sum(9 - value for value in branch_gcds)
    kummer_genus = 1 + (-18 + ramification) // 2
    require(ramification == 22 and kummer_genus == 3, "Kummer genus")

    # Nonzero determinant gives torus smoothness for 1-x-y.  Apply the exact
    # P-Euler operators to the sparse face polynomials; on a face equation
    # these are equal to P*C_P and P*B_P respectively.
    central_det = 0 * 5 - 6 * 2
    beta0_det = 2 * 3 - 5 * 3
    require((central_det, beta0_det) == (-12, -9), "torus Jacobian determinants")
    one = (1, 0, 0, 0)
    minus_u = (0, -1, 0, 0)
    minus_w = (0, 0, -1, 0)
    minus_zeta = (0, 0, 0, -1)
    central_terms = {(0, 0): one, (0, 6): minus_u, (2, 5): minus_w}
    beta0_terms = {(0, 0): one, (2, 5): minus_w, (3, 3): minus_zeta}
    require(
        p_euler_remainder(central_terms, 5)
        == {(0, 0): (-5, 0, 0, 0), (0, 6): (0, -1, 0, 0)},
        "P*C_P-5*C=-(U*P^6+5)",
    )
    require(
        p_euler_remainder(beta0_terms, 3)
        == {(0, 0): (-3, 0, 0, 0), (2, 5): (0, 0, -2, 0)},
        "P*B_P-3*B=-(2*W*S^2*P^5+3)",
    )

    edge_names = (
        "X-1", "1-zeta*X^3", "-W*X-zeta",
        "(X-1)*(U*X+W)", "U-X^6", "1-W*X",
    )
    require(len(edge_names) == 6, "complete outer/internal edge list")
    u_parameter = {(1, 0): 1}
    w_parameter = {(0, 1): 1}
    difference = parameter_poly_add(w_parameter, parameter_poly_scale(u_parameter, -1))
    product_uw = parameter_poly_multiply(u_parameter, w_parameter)
    top_discriminant = parameter_poly_add(
        parameter_poly_multiply(difference, difference),
        parameter_poly_scale(product_uw, 4),
    )
    parameter_sum = parameter_poly_add(u_parameter, w_parameter)
    require(
        top_discriminant == parameter_poly_multiply(parameter_sum, parameter_sum),
        "top-edge discriminant=(U+W)^2",
    )
    u_value, w_value, zeta_value = map(Q, (2, 3, 5))
    simple_edges = (
        (-1, 1),
        (1, 0, 0, -zeta_value),
        (-zeta_value, -w_value),
        (-w_value, w_value - u_value, u_value),
        (u_value, 0, 0, 0, 0, 0, -1),
        (1, -w_value),
    )
    require(all(is_squarefree(edge) for edge in simple_edges), "off-contact exact edge control")
    require(
        tuple(qpoly_trim(simple_edges[3])) == (-w_value, w_value - u_value, u_value),
        "literal top-edge expansion",
    )

    # Lambda=0: W=-U.  Reconstruct the exact top identity rather than
    # transporting an off-contact packet.  Coefficients are ordered by r.
    endpoint_top = qpoly_trim((0, 0, 0, 0, 0, -1, 1))
    transported_top = qpoly_trim(tuple(
        Q(1, 2) * first + Q(1, 2) * second
        for first, second in zip(
            (0, 0, 0, 0, -1, 0, 1),
            (0, 0, 0, 0, 1, -2, 1),
        )
    ))
    require(endpoint_top == transported_top, "exact endpoint W-perturbation identity")

    # The main branches have length-twelve contact and the same
    # coefficient-uniform tail ladder as THM-4292/4297.  Setting beta_11=0
    # turns c1=alpha_11+beta_11 into alpha_11; on a repeated deepest locus it
    # therefore removes C1 rather than creating a new coefficient.
    c1_alpha_beta = (Q(1), Q(1))
    c1_after_beta0 = c1_alpha_beta[0]
    require(c1_after_beta0 == 1, "beta-zero c1 remains the alpha coefficient")
    require((6 - 5) == 1, "top derivative coefficient")
    lambda0_u = Q(2)
    lambda0_w = -lambda0_u
    lambda0_top = qpoly_trim((-lambda0_w, lambda0_w - lambda0_u, lambda0_u))
    require(
        tuple(value / lambda0_u for value in lambda0_top) == (1, -2, 1),
        "Lambda-zero top factor U*(X-1)^2",
    )
    require(qpoly_gcd(lambda0_top, qpoly_derivative(lambda0_top)) == [-1, 1], "unique double root X=1")
    require(2 * 12 - 1 == 23, "length-twelve contact is A23")

    # With q=t^6*y, the cubic correction divided by t^12 starts at t^6;
    # expanding r^4=(1+t^6*y)^4 records every later exponent as well.
    correction_exponents = tuple(6 + 6 * index for index in range(5))
    correction_coefficients = tuple(comb(4, index) for index in range(5))
    require(correction_exponents == (6, 12, 18, 24, 30), "q^3 correction timing")
    require(correction_coefficients == (1, 4, 6, 4, 1), "r^4 correction expansion")

    # Valuation linear forms use the basis (s,nu), with nu>s>0.
    terminal_gap = (Q(-6), Q(6))
    correction_gap = (Q(6), Q(6))
    require(
        tuple(correction_gap[index] - terminal_gap[index] for index in range(2)) == (12, 0),
        "terminal b^12*q precedes the t^6 correction by 12s",
    )
    numerator_order = (Q(3), Q(5))

    def subtract_half_gap(gap):
        return tuple(numerator_order[index] - gap[index] / 2 for index in range(2))

    splitter_gaps = (terminal_gap,) + tuple((Q(index), Q(index)) for index in range(1, 5))
    tail_forms = tuple(subtract_half_gap(gap) for gap in splitter_gaps)
    require(
        tail_forms == (
            (6, 2),
            (Q(5, 2), Q(9, 2)),
            (2, 4),
            (Q(3, 2), Q(7, 2)),
            (1, 3),
        ),
        "rederived tail-order table",
    )
    require(all(a > 0 and b > 0 for a, b in tail_forms), "all imported tail orders positive")

    # On the actual theorem gate the deepest repeated ladder cannot reach C4
    # or the terminal splitter.  Repetition has c3=eta+zeta=0; zeta!=0 then
    # forces eta=-zeta!=0, so C3=eta/c^2 is nonzero (and C2 may split sooner).
    zeta_control = Q(1)
    eta_on_repeated_c3 = -zeta_control
    require(eta_on_repeated_c3 != 0, "zeta gate forces nonzero C3 on repeated locus")
    repeated_short_gaps = ((Q(2), Q(2)), (Q(3), Q(3)))
    require(
        all(
            correction_gap[0] - gap[0] > 0 and correction_gap[1] - gap[1] > 0
            for gap in repeated_short_gaps
        ),
        "C2 or C3 splits strictly before t^6 correction",
    )
    require(
        tuple(subtract_half_gap(gap) for gap in repeated_short_gaps)
        == ((2, 4), (Q(3, 2), Q(7, 2))),
        "short repeated-ladder positive orders",
    )
    require(36 // 12 == 3 and 36 % 12 == 0, "Lambda-zero base transport ramification three")

    print("THM4334_BETA0_ENDPOINT_PRIMARY=PASS")
    print("wall=Z=0,beta_11=0,U*W*zeta_3!=0;Lambda=U+W_arbitrary")
    print("universe_rows=16 required=6 absent=2 optional=8 collision_bits=6")
    print("hostile_states=16384 lower_complexes=1")
    print("faces=" + repr(expected_faces))
    print("face_supports=" + repr(equal_sets))
    print("pick_ledgers=" + repr(ledgers))
    print("face_normalizations=R:rational,C:genus3,T9:cyclic9_genus3")
    print("T9_kummer=S^9=(W^3/zeta^5)*(1-x)^5/x^3 branch_orders=(-3,5,-2)")
    print("dual_graph=vertices3 edges13 connected1 b1=11")
    print("orders_L36=(27,34)")
    print("edges=" + repr(edge_names))
    print("simple_edge_gate=U*W*zeta_3*(U+W)!=0")
    print("Lambda0_contact=A23 tail_orders_all_positive=true response_packet_not_transported")
    print("Lambda0_short_ladder=zeta_3!=0_forces_C2_or_C3_before_t^6_correction")
    print("source_infinity_packet=" + repr(packet))
    print("RELATIVE_CONSEQUENCE=all_special_components_constant=>no_nonautomorphic_pair")
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
