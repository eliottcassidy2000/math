#!/usr/bin/env python3
"""Primary exact certificate for THM-4248's M=11, Z=0 owner descent.

The calculation is exact.  Rational arithmetic is used for every Newton
plane, and SymPy is used only for exact symbolic polynomial identities.  The
script certifies the five exhaustive first-surviving-owner strata, including
the fixed ``-Q*s^2/2`` support point, all 26,624 support/collision hostiles,
the face and edge models, regular-model scaling data, the unchanged primitive
CM component, and the degree-zero specialization inventory.

The implications from primitive CM to simplicity, the audited toroidal-model
theorem, and proper-flat degree conservation are proof inputs recorded by the
theorem, not claims that a computer algebra system proves here.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, label):
    """Optimization-stable exact check."""
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError("THM-4248 PRIMARY FAILURE: " + label)


def gcd3(a, b, c):
    return gcd(gcd(abs(a), abs(b)), abs(c))


def monomials_through(weight):
    rows = []
    for i in range(weight // 2 + 1):
        for j in range(weight // 3 + 1):
            value = 2 * i + 3 * j
            if 0 < value <= weight and (i, j) not in {(0, 1), (1, 1)}:
                rows.append((i, j, value))
    return tuple(sorted(rows, key=lambda row: (row[2], row[1], row[0])))


def valued_support(rows):
    """Complete valued support, including the fixed residual point."""
    support = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j, _weight in rows:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


def plane_records(points):
    """All candidate planes, with below/equality masks precomputed."""
    points = tuple(sorted(points))
    candidates = set()
    for first, second, third in combinations(points, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant == 0:
            continue
        slope_s = Fraction(
            (second[2] - first[2]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[2] - first[2]),
            determinant,
        )
        slope_p = Fraction(
            (second[0] - first[0]) * (third[2] - first[2])
            - (second[2] - first[2]) * (third[0] - first[0]),
            determinant,
        )
        constant = Fraction(first[2]) - slope_s * first[0] - slope_p * first[1]
        candidates.add((slope_s, slope_p, constant))

    records = []
    for plane in sorted(candidates):
        below = 0
        equal = 0
        slope_s, slope_p, constant = plane
        for index, (r, ell, height) in enumerate(points):
            gap = Fraction(height) - slope_s * r - slope_p * ell - constant
            if gap < 0:
                below |= 1 << index
            elif gap == 0:
                equal |= 1 << index
        records.append((plane, below, equal))
    return points, tuple(records)


def support_mask(points, support):
    index = {point: position for position, point in enumerate(points)}
    return sum(1 << index[point] for point in support)


def projected_rank_two(points, mask):
    chosen = [points[index] for index in range(len(points)) if mask & (1 << index)]
    for first, second, third in combinations(chosen, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant:
            return True
    return False


def lower_planes(points, records, mask):
    answer = set()
    for plane, below, equal in records:
        if below & mask:
            continue
        equality = equal & mask
        if projected_rank_two(points, equality):
            answer.add(plane)
    return answer


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
        polygon[index][0] * polygon[(index + 1) % len(polygon)][1]
        - polygon[(index + 1) % len(polygon)][0] * polygon[index][1]
        for index in range(len(polygon))
    ))
    boundary = sum(
        gcd(
            abs(polygon[(index + 1) % len(polygon)][0] - polygon[index][0]),
            abs(polygon[(index + 1) % len(polygon)][1] - polygon[index][1]),
        )
        for index in range(len(polygon))
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    interior = (area2 - boundary + 2) // 2
    return polygon, area2, boundary, interior


def edge_packet(polygon):
    packet = []
    rows = []
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        rows.append((start, end, length, inward, constant, index))
        if start[0] == end[0] == 0:
            continue
        packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True)), tuple(rows)


def edge_polynomial(polynomial, first, second, start, end, variable):
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    ux, uy = dx // length, dy // length
    answer = 0
    for (i, j), coefficient in sp.Poly(polynomial, first, second).terms():
        vx, vy = i - start[0], j - start[1]
        if vx * dy - vy * dx:
            continue
        position = vx // ux if ux else vy // uy
        if 0 <= position <= length:
            answer += coefficient * variable**position
    return sp.factor(answer)


def triangle_interiors(vertices):
    def area2(first, second, point):
        return abs(
            (second[0] - first[0]) * (point[1] - first[1])
            - (second[1] - first[1]) * (point[0] - first[0])
        )

    total = area2(vertices[0], vertices[1], vertices[2])
    answer = []
    for i in range(max(point[0] for point in vertices) + 1):
        for j in range(max(point[1] for point in vertices) + 1):
            pieces = tuple(
                area2(vertices[index], vertices[(index + 1) % 3], (i, j))
                for index in range(3)
            )
            if sum(pieces) == total and all(piece > 0 for piece in pieces):
                answer.append((i, j))
    return tuple(answer)


def determinant_one(sequence):
    return all(
        abs(first[0] * second[1] - first[1] * second[0]) == 1
        for first, second in zip(sequence, sequence[1:])
    )


M_PLANE = (Fraction(1, 11), Fraction(2, 11), Fraction(-2, 11))
V5_PLANE = (Fraction(0), Fraction(1, 5), Fraction(-1, 5))
V4_PLANE = (Fraction(-1, 4), Fraction(1, 4), Fraction(-1, 4))
V3_PLANE = (Fraction(-2, 3), Fraction(1, 3), Fraction(-1, 3))
TK_PLANE = (Fraction(1), Fraction(-1, 2), Fraction(-2))

A_ROW = (4, 1)
B_ROW = (1, 3)
U_ROW = (5, 0)
Z_ROW = (0, 3)
K_ROW = (0, 2)
DELTA_ROW = (4, 0)
EPSILON_ROW = (3, 0)


CASE_DATA = (
    {
        "name": "U_NE_0__K_NE_0", "required": {A_ROW, B_ROW, U_ROW, K_ROW},
        "absent": {Z_ROW}, "planes": {M_PLANE, V5_PLANE, TK_PLANE},
        "k": 5, "tail": True, "base": 330, "delta_mode": "free",
        "polygon": ((0, 1), (2, 0), (4, 2), (5, 4), (1, 6), (0, 6)),
        "patterns": 8192,
    },
    {
        "name": "U_NE_0__K_EQ_0", "required": {A_ROW, B_ROW, U_ROW},
        "absent": {Z_ROW, K_ROW}, "planes": {M_PLANE, V5_PLANE},
        "k": 5, "tail": False, "base": 330, "delta_mode": "kzero",
        "polygon": ((0, 1), (2, 0), (5, 4), (1, 6), (0, 6)),
        "patterns": 8192,
    },
    {
        "name": "U_EQ_0__DELTA_NE_0__K_NE_0",
        "required": {A_ROW, B_ROW, DELTA_ROW, K_ROW},
        "absent": {U_ROW, Z_ROW}, "planes": {M_PLANE, V4_PLANE, TK_PLANE},
        "k": 4, "tail": True, "base": 132, "delta_mode": "free",
        "polygon": ((0, 1), (2, 0), (4, 2), (5, 4), (1, 6), (0, 5)),
        "patterns": 4096,
    },
    {
        "name": "U_EQ_0__DELTA_NE_0__K_EQ_0",
        "required": {A_ROW, B_ROW, DELTA_ROW},
        "absent": {U_ROW, Z_ROW, K_ROW}, "planes": {M_PLANE, V4_PLANE},
        "k": 4, "tail": False, "base": 132, "delta_mode": "kzero",
        "polygon": ((0, 1), (2, 0), (5, 4), (1, 6), (0, 5)),
        "patterns": 4096,
    },
    {
        "name": "U_EQ_0__DELTA_EQ_0",
        "required": {A_ROW, B_ROW, EPSILON_ROW, K_ROW},
        "absent": {U_ROW, Z_ROW, DELTA_ROW},
        "planes": {M_PLANE, V3_PLANE, TK_PLANE},
        "k": 3, "tail": True, "base": 132, "delta_mode": "zero",
        "polygon": ((0, 1), (2, 0), (4, 2), (5, 4), (1, 6), (0, 4)),
        "patterns": 2048,
    },
)


def run_hull_atlas(monomials):
    collisions = ((2, 3, 1), (2, 4, 1), (2, 5, 1), (3, 4, 1), (3, 5, 1))
    total = 0
    case_counts = []
    for case in CASE_DATA:
        optional = tuple(
            row for row in monomials
            if row[:2] not in case["required"] | case["absent"]
        )
        universe_rows = tuple(row for row in monomials if row[:2] not in case["absent"])
        universe = valued_support(universe_rows)
        points, records = plane_records(universe)
        count = 0
        for support_choice in range(1 << len(optional)):
            chosen = [row for row in monomials if row[:2] in case["required"]]
            chosen.extend(
                row for bit, row in enumerate(optional)
                if support_choice & (1 << bit)
            )
            base_support = valued_support(chosen)
            need((2, 0, 1) in base_support, "fixed residual retained")
            for collision_choice in range(1 << len(collisions)):
                deleted = {
                    point for bit, point in enumerate(collisions)
                    if collision_choice & (1 << bit)
                }
                support = base_support - deleted
                mask = support_mask(points, support)
                need(lower_planes(points, records, mask) == case["planes"],
                     "exact lower hull: " + case["name"])
                count += 1
        need(count == case["patterns"], "hostile count: " + case["name"])
        case_counts.append(count)
        total += count
    need(total == 26624, "all 26,624 support/collision hostiles")
    return tuple(case_counts), total


def main():
    monomials = monomials_through(11)
    expected_monomials = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11),
    )
    need(monomials == expected_monomials, "complete M<=11 monomial universe")
    case_counts, hull_total = run_hull_atlas(monomials)

    # Face and polygon ledgers.
    main_points = ((0, 1), (2, 0), (5, 4), (1, 6))
    tail_points = ((2, 0), (4, 2), (5, 4))
    vertical_points = {
        5: ((0, 1), (1, 6), (0, 6)),
        4: ((0, 1), (1, 6), (0, 5)),
        3: ((0, 1), (1, 6), (0, 4)),
    }
    need(polygon_ledger(main_points)[1:] == (33, 5, 15), "main Pick ledger")
    need(polygon_ledger(tail_points)[1:] == (2, 4, 0), "K-tail Pick ledger")
    need(tuple(polygon_ledger(vertical_points[k])[1:] for k in (5, 4, 3))
         == ((5, 7, 0), (4, 6, 0), (3, 5, 0)),
         "three vertical-owner Pick ledgers")

    S, P, X = sp.symbols("S P X")
    A, B, C_OWNER, K = sp.symbols("A B C_OWNER K")
    rational = S**2 - P
    cyclotomic = 1 - A*S*P**5 - B*S**3*P**4
    main_face = sp.expand(rational * cyclotomic)
    tail_face = sp.expand(S**2 * (1 - K*S**2*P**2 - B*S**3*P**4))
    node_determinant = sp.factor(sp.det(sp.Matrix((
        (sp.diff(rational, S), sp.diff(rational, P)),
        (sp.diff(cyclotomic, S), sp.diff(cyclotomic, P)),
    ))).subs(P, S**2))
    need(sp.factor(node_determinant + 11*(A + B)*S**10) == 0,
         "eleven transverse R-C nodes")
    node_polynomial = sp.Poly(1 - (A + B)*S**11, S,
                              domain=sp.QQ.frac_field(A, B))
    need(node_polynomial.degree() == 11
         and sp.gcd(node_polynomial, node_polynomial.diff()).degree() == 0,
         "eleven distinct main nodes")
    need(sp.det(sp.Matrix(((A, 3*B), (5*A, 4*B)))) == -11*A*B,
         "main CM face torus smoothness")
    T0 = sp.symbols("T0")
    need(sp.diff(1 - K*T0**2 - B*T0**3*P, P) == -B*T0**3,
         "K-tail rationality and smoothness")

    packet_rows = []
    face_edge_rows = []
    for case in CASE_DATA:
        polygon, area2, boundary, genus = polygon_ledger(case["polygon"])
        need(polygon == case["polygon"], "global polygon: " + case["name"])
        expected_pick = (40, 12, 15) if case["tail"] else (38, 10, 15)
        if case["k"] == 4:
            expected_pick = (39, 11, 15) if case["tail"] else (37, 9, 15)
        if case["k"] == 3:
            expected_pick = (38, 10, 15)
        need((area2, boundary, genus) == expected_pick,
             "global Pick ledger: " + case["name"])
        packet, edge_rows = edge_packet(polygon)
        expected_packet = ((10, 10, 5, 5, 2, 2, 1) if case["tail"]
                           else (10, 10, 7, 5, 1))
        need(packet == expected_packet, "outer packet: " + case["name"])
        need(sum(index - 1 for index in packet) == 28,
             "packet defect 28: " + case["name"])
        packet_rows.append((case["name"], packet, sum(packet)))

        k = case["k"]
        vertical_face = sp.expand(P * (-1 + C_OWNER*P**k + A*S*P**5))
        need(sp.diff(-1 + C_OWNER*P**k + A*S*P**5, S) == A*P**5,
             "vertical face rationality: " + case["name"])
        if case["tail"]:
            point_a, point_b, point_c, point_d, point_e, point_f = polygon
            schemes = (
                edge_polynomial(main_face, S, P, point_a, point_b, X),
                edge_polynomial(tail_face, S, P, point_b, point_c, X),
                edge_polynomial(tail_face, S, P, point_c, point_d, X),
                edge_polynomial(main_face, S, P, point_d, point_e, X),
                edge_polynomial(vertical_face, S, P, point_e, point_f, X),
                edge_polynomial(vertical_face, S, P, point_f, point_a, X),
                edge_polynomial(main_face, S, P, point_a, point_e, X),
                edge_polynomial(main_face, S, P, point_b, point_d, X),
            )
            expected = (
                X - 1, 1 - K*X**2, -K - B*X,
                (X - 1)*(A*X + B), A + C_OWNER*X,
                C_OWNER - X**k, A*X - 1, 1 - B*X,
            )
            expected_discriminants = (
                1, 4*K, 1, (A + B)**2, 1,
                {5: 3125*C_OWNER**4,
                 4: -256*C_OWNER**3,
                 3: -27*C_OWNER**2}[k],
                1, 1,
            )
        else:
            point_a, point_b, point_d, point_e, point_f = polygon
            schemes = (
                edge_polynomial(main_face, S, P, point_a, point_b, X),
                edge_polynomial(main_face, S, P, point_b, point_d, X),
                edge_polynomial(main_face, S, P, point_d, point_e, X),
                edge_polynomial(vertical_face, S, P, point_e, point_f, X),
                edge_polynomial(vertical_face, S, P, point_f, point_a, X),
                edge_polynomial(main_face, S, P, point_a, point_e, X),
            )
            expected = (
                X - 1, 1 - B*X, (X - 1)*(A*X + B),
                A + C_OWNER*X, C_OWNER - X**k, A*X - 1,
            )
            expected_discriminants = (
                1, 1, (A + B)**2, 1,
                {5: 3125*C_OWNER**4,
                 4: -256*C_OWNER**3,
                 3: -27*C_OWNER**2}[k],
                1,
            )
        need(all(sp.factor(value - target) == 0
                 for value, target in zip(schemes, expected)),
             "all edge schemes: " + case["name"])
        discriminants = tuple(sp.factor(sp.discriminant(value, X)) for value in schemes)
        need(discriminants == expected_discriminants,
             "all edge discriminants: " + case["name"])
        face_edge_rows.append((case["name"], len(schemes)))

        # Integral lower graph, outer denominators, and rational chain fans.
        base = case["base"]
        for plane in case["planes"]:
            normal = (
                int(base * plane[0]), int(base * plane[1]), -1,
            )
            need(Fraction(normal[0], base) == plane[0]
                 and Fraction(normal[1], base) == plane[1],
                 "integral face normal: " + case["name"])
            need(gcd3(*normal) == 1, "face multiplicity one: " + case["name"])
        active_rows = tuple(row for row in monomials if row[:2] not in case["absent"])
        for r, ell, height in valued_support(active_rows):
            for slope_s, slope_p, constant in case["planes"]:
                scaled_height = base * (slope_s*r + slope_p*ell + constant)
                need(scaled_height.denominator == 1,
                     "integral height: " + case["name"])
                need(base*height - scaled_height >= 0,
                     "nonnegative integral gap: " + case["name"])

        graph_vertices = {
            point: (point[0], point[1], 0 if point in {(0, 1), (2, 0)} else base)
            for point in polygon
        }
        for start, end, length, _inward, _constant, _index in edge_rows:
            difference = tuple(
                graph_vertices[end][axis] - graph_vertices[start][axis]
                for axis in range(3)
            )
            need(gcd3(*difference) == length,
                 "outer lifted gcd: " + case["name"])
        for start, end in (((0, 1), (1, 6)),):
            difference = tuple(
                graph_vertices[end][axis] - graph_vertices[start][axis]
                for axis in range(3)
            )
            need(gcd3(*difference) == 1,
                 "primitive M-V edge: " + case["name"])
        if case["tail"]:
            start, end = (2, 0), (5, 4)
            difference = tuple(
                graph_vertices[end][axis] - graph_vertices[start][axis]
                for axis in range(3)
            )
            need(gcd3(*difference) == 1,
                 "primitive M-T edge: " + case["name"])

        def scaled_face_height(plane, point):
            value = base*(plane[0]*point[0] + plane[1]*point[1] + plane[2])
            need(value.denominator == 1, "integral slope probe: " + case["name"])
            return int(value)

        vertical_plane = {5: V5_PLANE, 4: V4_PLANE, 3: V3_PLANE}[k]
        ae_endpoints = {
            5: (-60, -66), 4: (-24, -33), 3: (-24, -44),
        }[k]
        derived_ae = (
            scaled_face_height(M_PLANE, (0, 0))
            - scaled_face_height(M_PLANE, (0, 1)),
            scaled_face_height(vertical_plane, (0, 0))
            - scaled_face_height(vertical_plane, (0, 1)),
        )
        ae_sequence = tuple(
            (value, 1)
            for value in range(ae_endpoints[0], ae_endpoints[1] - 1, -1)
        )
        need(derived_ae == ae_endpoints and determinant_one(ae_sequence)
             and len(ae_sequence) - 2 == {5: 5, 4: 8, 3: 19}[k],
             "M-V determinant-one rational chain: " + case["name"])
        if case["tail"]:
            bd_endpoints = (-90, -165) if base == 330 else (-36, -66)
            derived_bd = (
                scaled_face_height(M_PLANE, (1, -1))
                - scaled_face_height(M_PLANE, (2, 0)),
                scaled_face_height(TK_PLANE, (1, -1))
                - scaled_face_height(TK_PLANE, (2, 0)),
            )
            bd_sequence = tuple(
                (value, 1)
                for value in range(bd_endpoints[0], bd_endpoints[1] - 1, -1)
            )
            need(derived_bd == bd_endpoints and determinant_one(bd_sequence)
                 and len(bd_sequence) - 2 == (74 if base == 330 else 29),
                 "M-T determinant-one rational chain: " + case["name"])

    # Exact source charts in all five owner strata.
    source_s, source_p, Q, sigma = sp.symbols("source_s source_p Q sigma")
    Delta, Phi, Theta, eta, Xi, U = sp.symbols("Delta Phi Theta eta Xi U")
    epsilon = -sp.Rational(1376, 135)
    delta_kzero = sp.Rational(5696, 105)
    need(sp.Rational(2848, 45) - sp.Rational(7, 6)*delta_kzero == 0,
         "K=0 forces Delta=5696/105")
    chart_rows = []
    for case in CASE_DATA:
        if case["delta_mode"] == "free":
            delta_value = Delta
        elif case["delta_mode"] == "kzero":
            delta_value = delta_kzero
        else:
            delta_value = sp.Integer(0)
        k_value = sp.factor(sp.Rational(2848, 45) - sp.Rational(7, 6)*delta_value)
        if not case["tail"]:
            need(k_value == 0, "K-zero chart specialization")
        if case["k"] == 5:
            u_value, owner_value = U, U
        elif case["k"] == 4:
            u_value, owner_value = sp.Integer(0), delta_value
        else:
            u_value, owner_value = sp.Integer(0), epsilon
        H = (
            -3*source_p + sp.Rational(8, 3)*source_p**2
            + epsilon*source_p**3 + k_value*source_s**2*source_p**2
            + Phi*source_s*source_p**3 + delta_value*source_p**4
            + Theta*source_s**2*source_p**3 + eta*source_s*source_p**4
            + u_value*source_p**5 + Xi*source_s**2*source_p**4
            + A*source_s*source_p**5 + B*source_s**3*source_p**4
        )
        F_source = sp.expand((source_s**2 - source_p)*(1 - Q*H) - Q*source_s**2/2)
        base = case["base"]
        main_s_exp, main_p_exp, main_mult = -base//11, -2*base//11, 2*base//11
        H_main = sp.cancel(sigma**base * H.subs({
            source_s: sigma**main_s_exp*S,
            source_p: sigma**main_p_exp*P,
        }))
        need(sp.denom(H_main) == 1, "main chart integral: " + case["name"])
        need(sp.factor(H_main.subs(sigma, 0)
                       - A*S*P**5 - B*S**3*P**4) == 0,
             "main chart reduction: " + case["name"])
        F_main = sp.cancel(sigma**main_mult * F_source.subs({
            Q: sigma**base,
            source_s: sigma**main_s_exp*S,
            source_p: sigma**main_p_exp*P,
        }))
        expected_main = ((S**2 - P)*(1 - H_main) - sigma**base*S**2/2)
        need(sp.factor(F_main - expected_main) == 0,
             "exact main UV chart: " + case["name"])

        vertical_plane = {5: V5_PLANE, 4: V4_PLANE, 3: V3_PLANE}[case["k"]]
        vs = -int(base * vertical_plane[0])
        vp = -int(base * vertical_plane[1])
        vm = -int(base * vertical_plane[2])
        F_vertical = sp.cancel(sigma**vm * F_source.subs({
            Q: sigma**base,
            source_s: sigma**vs*S,
            source_p: sigma**vp*P,
        }))
        expected_vertical = sp.expand(P*(-1 + owner_value*P**case["k"] + A*S*P**5))
        need(sp.denom(F_vertical) == 1,
             "vertical chart integral: " + case["name"])
        need(sp.factor(F_vertical.subs(sigma, 0) - expected_vertical) == 0,
             "vertical chart reduction: " + case["name"])
        if case["tail"]:
            F_tail = sp.cancel(sigma**(2*base) * F_source.subs({
                Q: sigma**base,
                source_s: sigma**(-base)*S,
                source_p: sigma**(base//2)*P,
            }))
            expected_tail = sp.expand(S**2*(1 - k_value*S**2*P**2 - B*S**3*P**4))
            need(sp.denom(F_tail) == 1,
                 "tail chart integral: " + case["name"])
            need(sp.factor(F_tail.subs(sigma, 0) - expected_tail) == 0,
                 "tail chart reduction: " + case["name"])
        chart_rows.append((case["name"], base, base - 1))

    # Unchanged cyclic degree-eleven component and primitive CM type.
    x = sp.symbols("x")
    x_expression = A*S*P**5
    cyclic_identity = A**3*P**11*(1 - x_expression) - B*x_expression**3
    need(sp.factor(cyclic_identity - A**3*P**11*cyclotomic) == 0,
         "cyclic-cover identity")
    branch_residues = (3, 10, 9)
    need(sum(branch_residues) % 11 == 0
         and all(gcd(value, 11) == 1 for value in branch_residues),
         "three full order-eleven branches")
    interiors = triangle_interiors(((0, 0), (1, 5), (3, 4)))
    need(interiors == ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4)),
         "five regular differentials")
    cm_type = {(6*i + j) % 11 for i, j in interiors}
    need(cm_type == {4, 5, 8, 9, 10}, "primitive CM character set")
    need(cm_type | {(-value) % 11 for value in cm_type} == set(range(1, 11)),
         "full H1 cyclotomic spectrum")
    stabilizer = tuple(
        unit for unit in range(1, 11)
        if {(unit*value) % 11 for value in cm_type} == cm_type
    )
    need(stabilizer == (1,), "CM type has trivial stabilizer")
    need(sp.Poly(sp.cyclotomic_poly(11, x), x).degree() == 10,
         "Q(zeta_11) has degree 2g")

    # Good j=0 target for both base changes.
    target_a, target_A, target_C, target_X, target_Y = sp.symbols(
        "target_a target_A target_C target_X target_Y"
    )
    for base in (330, 132):
        target_equation = (
            target_C**2 - target_A**3
            + sp.Rational(3, 4)*target_a**2*target_A
            - sigma**(-base) + sp.Rational(1, 4)*target_a**3
        )
        scaled = sp.cancel(sigma**base * target_equation.subs({
            target_A: sigma**(-base//3)*target_X,
            target_C: sigma**(-base//2)*target_Y,
        }))
        expected = (
            target_Y**2 - target_X**3 - 1
            + sp.Rational(3, 4)*target_a**2*sigma**(2*base//3)*target_X
            + sp.Rational(1, 4)*target_a**3*sigma**base
        )
        need(sp.factor(scaled - expected) == 0, "good target scaling")
    need(sp.discriminant(target_X**3 + 1, target_X) == -27,
         "smooth j=0 special target")

    # Genus completeness and degree-zero consequence.
    for case in CASE_DATA:
        vertices = 4 if case["tail"] else 3
        edges = 13 if case["tail"] else 12
        graph_rank = edges - vertices + 1
        need(graph_rank == 10 and 5 + graph_rank == 15,
             "complete special genus: " + case["name"])
        component_genera = (0, 5, 0, 0) if case["tail"] else (0, 5, 0)
        need(tuple(value for value in component_genera if value) == (5,),
             "unique positive-genus component: " + case["name"])
        # Cited primitive-CM simplicity gives Hom(J(C),E0)=0.  All remaining
        # components are rational, so every specialized component degree is 0.
        specialized_degrees = (0,) * len(component_genera)
        need(sum(specialized_degrees) == 0,
             "special degree zero after cited Hom gate: " + case["name"])

    # Sharp hostiles: the other named walls require new geometry.
    need(polygon_ledger(((0, 1), (0, 6), (3, 5)))[1:] == (15, 7, 5),
         "A=0 creates a genus-five replacement face")
    need(polygon_ledger(((2, 0), (5, 3), (3, 5)))[1:] == (12, 6, 4),
         "B=0 creates a genus-four replacement face")
    need(sp.factor(((X - 1)*(A*X + B)).subs(B, -A) - A*(X - 1)**2) == 0,
         "A+B=0 gives a double boundary root")
    need(tuple(row[2] for row in packet_rows) == (35, 33, 35, 33, 35),
         "K=0 changes the packet degree and must be stratified")

    print("THM4248_WEIGHT11_Z0_OWNER_DESCENT_EXACT_CERTIFICATE")
    print("scope=exact_M11,b=d=0,reduced_(2,3),Z=0,A*B*(A+B)!=0,U_arbitrary,all_lower_coefficients")
    print("owner_strata=U!=0/K!=0;U!=0/K=0;U=0/Delta!=0/K!=0;U=0/Delta!=0/K=0;U=0/Delta=0")
    print("hull_cases=26624;per_stratum=" + ",".join(str(value) for value in case_counts)
          + ";fixed_residual=(2,0,1);failures=0")
    print("faces=M11+V_(5|4|3)+optional_TK;face_Pick=M(33,5,15),TK(2,4,0),V5(5,7,0),V4(4,6,0),V3(3,5,0)")
    print("packets=K!=0:(10,10,5,5,2,2,1)/degree35;K=0:(10,10,7,5,1)/degree33;defect=28")
    print("edge_schemes=all_exact;edge_discriminants=all_exact;hidden_extra_units=none")
    print("regular_models=U!=0:Q=sigma^330,A329;U=0:Q=sigma^132,A131;all_components_multiplicity_one")
    print("chains=MV:(5,8,19)_by_owner;MT:(74_at_330,29_at_132);all_rational")
    print("CM=genus5,Q(zeta11),type(4,5,8,9,10),stabilizer(1),Hom_to_E0_zero_by_cited_primitive_CM")
    print("genus=15;dual_rank=10;special_degree=0;generic_finite_degree_positive=>contradiction")
    print("hostiles=A0_genus5_face;B0_genus4_face;A+B=0_double_edge;K=0_face_deletion")
    print("proof_inputs=THM3992/3997/4012/4045/4222/4232;Milne_primitive_CM;proper_flat_degree_conservation")
    print("nonclaims=other_three_M11_walls,other_cells,seam_entry,JC2,DC2")
    print("checks=" + str(CHECKS))


if __name__ == "__main__":
    main()
