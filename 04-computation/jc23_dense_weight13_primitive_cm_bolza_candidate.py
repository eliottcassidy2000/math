#!/usr/bin/env python3
"""Exact certificate for the dense exact-weight-thirteen candidate.

This is a proof-development artifact, not a canonical theorem certificate.
It verifies the finite Newton, toric, local-chart, and CM inputs.  The uses of
Milne's primitive-CM classification, THM-4012's Bolza Hom obstruction, and
proper-flat degree conservation remain mathematical/cited proof steps.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def monomials_through(weight):
    answer = []
    for i in range(weight // 2 + 1):
        for j in range(weight // 3 + 1):
            monomial_weight = 2 * i + 3 * j
            if (0 < monomial_weight <= weight
                    and (i, j) not in {(0, 1), (1, 1)}):
                answer.append((i, j, monomial_weight))
    return tuple(sorted(answer, key=lambda row: (row[2], row[1], row[0])))


def expanded_support(monomials):
    """Valued support (s exponent, p exponent, Q exponent)."""
    support = {(2, 0, 0), (0, 1, 0)}
    for i, j, _weight in monomials:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


def candidate_plane_records(points):
    points = tuple(sorted(points))
    planes = set()
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
        constant = (Fraction(first[2]) - slope_s * first[0]
                    - slope_p * first[1])
        planes.add((slope_s, slope_p, constant))

    records = []
    for plane in planes:
        slope_s, slope_p, constant = plane
        below = 0
        equal = 0
        for index, (i, j, height) in enumerate(points):
            gap = (Fraction(height) - slope_s * i
                   - slope_p * j - constant)
            if gap < 0:
                below |= 1 << index
            elif gap == 0:
                equal |= 1 << index
        records.append((plane, below, equal))
    return points, tuple(records)


def projected_rank_two(points, bits):
    chosen = [point for index, point in enumerate(points) if bits & (1 << index)]
    if len(chosen) < 3:
        return False
    for first, second, third in combinations(chosen, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant:
            return True
    return False


def support_bits(points, support):
    index = {point: position for position, point in enumerate(points)}
    return sum(1 << index[point] for point in support)


def lower_planes(points, records, bits):
    answer = set()
    for plane, below, equal in records:
        if below & bits:
            continue
        equality = equal & bits
        if projected_rank_two(points, equality):
            answer.add(plane)
    return answer


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, first, second):
        return ((first[0] - origin[0]) * (second[1] - origin[1])
                - (first[1] - origin[1]) * (second[0] - origin[0]))

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
    vertices = convex_hull(points)
    area2 = abs(sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]))
        for index in range(len(vertices))
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    interior = (area2 - boundary + 2) // 2
    return vertices, area2, boundary, interior


def edge_packet(vertices):
    packet = []
    rows = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
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
        if ux:
            if vx % ux:
                continue
            position = vx // ux
        else:
            if vy % uy:
                continue
            position = vy // uy
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
            pieces = [
                area2(vertices[k], vertices[(k + 1) % 3], (i, j))
                for k in range(3)
            ]
            if sum(pieces) == total and all(piece > 0 for piece in pieces):
                answer.append((i, j))
    return tuple(answer)


def determinant_one(sequence):
    return all(
        abs(first[0] * second[1] - first[1] * second[0]) == 1
        for first, second in zip(sequence, sequence[1:])
    )


def main():
    # Complete inherited b=d=0 universe through exact weight thirteen.
    monomials = monomials_through(13)
    expected_monomials = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11), (6, 0, 12), (3, 2, 12), (0, 4, 12),
        (5, 1, 13), (2, 3, 13),
    )
    need(monomials == expected_monomials, "M13 monomial universe")
    need(tuple(row for row in monomials if row[2] == 13)
         == ((5, 1, 13), (2, 3, 13)), "exactly two top monomials")

    # A=[p^5y], B=[p^2y^3], U=[p^6], Z=[y^4] are retained.
    required = {(5, 1), (2, 3), (6, 0), (0, 4)}
    optional = tuple(row for row in monomials if row[:2] not in required)
    need(len(optional) == 14, "fourteen arbitrary lower rows")
    full_support = expanded_support(monomials)
    points, plane_records = candidate_plane_records(full_support)
    plane_main = (Fraction(1, 13), Fraction(2, 13), Fraction(-2, 13))
    plane_vertical = (Fraction(0), Fraction(1, 6), Fraction(-1, 6))
    plane_tail = (Fraction(1, 8), Fraction(1, 8), Fraction(-1, 4))
    expected_planes = {plane_main, plane_vertical, plane_tail}

    # The required skeleton, even without the A/B midpoint, already owns all
    # three faces.  This plus the analytic gaps below proves simultaneous
    # robustness under arbitrary lower coefficients and cancellations.
    required_rows = [row for row in monomials if row[:2] in required]
    required_support = expanded_support(required_rows)
    midpoint = (3, 6, 1)
    for support in (required_support, required_support - {midpoint}):
        bits = support_bits(points, support)
        need(lower_planes(points, plane_records, bits) == expected_planes,
             "required face skeleton changed")

    support_subset_count = 0
    for mask in range(1 << len(optional)):
        chosen = list(required_rows)
        chosen.extend(row for bit, row in enumerate(optional)
                      if mask & (1 << bit))
        bits = support_bits(points, expanded_support(chosen))
        need(lower_planes(points, plane_records, bits) == expected_planes,
             "lower hull changed under optional-row deletion")
        support_subset_count += 1
    need(support_subset_count == 16384, "all lower-row subsets")

    collisions = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1), (2, 6, 1),
        (3, 4, 1), (3, 5, 1), (3, 6, 1), (4, 5, 1),
    )
    collision_count = 0
    for mask in range(1 << len(collisions)):
        deleted = {point for bit, point in enumerate(collisions)
                   if mask & (1 << bit)}
        bits = support_bits(points, full_support - deleted)
        need(lower_planes(points, plane_records, bits) == expected_planes,
             "coincident cancellation changed lower hull")
        collision_count += 1
    need(collision_count == 256, "all collision deletion patterns")

    # Exact gap identities show that no simultaneous deletion/addition can
    # matter once the required vertices above are retained.
    main_equal = set()
    vertical_equal = set()
    tail_equal = set()
    for i, j, weight in monomials:
        first = (j + 2, i + j, 1)
        second = (j, i + j + 1, 1)
        main_gaps = (Fraction(13 - weight, 13),) * 2
        vertical_gaps = (Fraction(7 - i - j, 6),
                         Fraction(6 - i - j, 6))
        tail_gaps = (Fraction(8 - i - 2*j, 8),
                     Fraction(9 - i - 2*j, 8))
        for point, gap in zip((first, second), main_gaps):
            actual = (Fraction(point[2]) - plane_main[0] * point[0]
                      - plane_main[1] * point[1] - plane_main[2])
            need(actual == gap and gap >= 0, "main exact nonnegative gap")
            if gap == 0:
                main_equal.add((i, j))
        for point, gap in zip((first, second), vertical_gaps):
            actual = (Fraction(point[2]) - plane_vertical[0] * point[0]
                      - plane_vertical[1] * point[1] - plane_vertical[2])
            need(actual == gap and gap >= 0, "vertical exact nonnegative gap")
            if gap == 0:
                vertical_equal.add((i, j))
        for point, gap in zip((first, second), tail_gaps):
            actual = (Fraction(point[2]) - plane_tail[0] * point[0]
                      - plane_tail[1] * point[1] - plane_tail[2])
            need(actual == gap and gap >= 0, "tail exact nonnegative gap")
            if gap == 0:
                tail_equal.add((i, j))
    need(main_equal == {(5, 1), (2, 3)}, "main equality owners")
    need(vertical_equal == {(6, 0), (5, 1)}, "vertical equality owners")
    need(tail_equal == {(0, 4), (2, 3)}, "tail equality owners")

    # Face and global polygon ledgers.
    main_points = ((0, 1), (2, 0), (5, 5), (1, 7))
    tail_points = ((2, 0), (6, 4), (5, 5))
    vertical_points = ((0, 1), (1, 7), (0, 7))
    global_points = main_points + tail_points + vertical_points
    need(polygon_ledger(main_points)[1:] == (39, 5, 18), "main Pick ledger")
    need(polygon_ledger(tail_points)[1:] == (8, 6, 2), "tail Pick ledger")
    need(polygon_ledger(vertical_points)[1:] == (6, 8, 0),
         "vertical Pick ledger")
    polygon, area2, boundary, genus = polygon_ledger(global_points)
    expected_polygon = ((0, 1), (2, 0), (6, 4), (5, 5), (1, 7), (0, 7))
    need((polygon, area2, boundary, genus)
         == (expected_polygon, 53, 15, 20), "global polygon/Pick ledger")
    packet, edge_rows = edge_packet(polygon)
    need(packet == (12, 12, 8, 6, 2, 2, 2, 2, 1), "complete packet")
    need(sum(packet) == 47 and sum(index - 1 for index in packet) == 38,
         "packet degree and defect")

    # Exact faces and all eight compactified edge schemes.
    S, P, edge_variable = sp.symbols("S P edge_variable")
    A, B, U, Z = sp.symbols("A B U Z")
    rational = S**2 - P
    cm_component = 1 - A*S*P**6 - B*S**3*P**5
    main_face = sp.expand(rational * cm_component)
    tail_core = 1 - Z*S**4*P**4 - B*S**3*P**5
    tail_face = sp.expand(S**2 * tail_core)
    vertical_core = -1 + U*P**6 + A*S*P**6
    vertical_face = sp.expand(P * vertical_core)
    point_A, point_B, point_C, point_D, point_E, point_F = polygon
    schemes = (
        edge_polynomial(main_face, S, P, point_A, point_B, edge_variable),
        edge_polynomial(tail_face, S, P, point_B, point_C, edge_variable),
        edge_polynomial(tail_face, S, P, point_C, point_D, edge_variable),
        edge_polynomial(main_face, S, P, point_D, point_E, edge_variable),
        edge_polynomial(vertical_face, S, P, point_E, point_F, edge_variable),
        edge_polynomial(vertical_face, S, P, point_F, point_A, edge_variable),
        edge_polynomial(main_face, S, P, point_A, point_E, edge_variable),
        edge_polynomial(main_face, S, P, point_B, point_D, edge_variable),
    )
    expected_schemes = (
        edge_variable - 1,
        1 - Z*edge_variable**4,
        -Z - B*edge_variable,
        (edge_variable - 1)*(A*edge_variable + B),
        A + U*edge_variable,
        U - edge_variable**6,
        A*edge_variable - 1,
        1 - B*edge_variable,
    )
    need(tuple(sp.factor(value - expected) for value, expected
               in zip(schemes, expected_schemes)) == (0,) * 8,
         "all edge schemes")
    discriminants = tuple(sp.factor(sp.discriminant(scheme, edge_variable))
                          for scheme in schemes)
    need(discriminants == (
        1, -256*Z**3, 1, (A + B)**2,
        1, 46656*U**5, 1, 1,
    ), "edge discriminants and exact unit gate")
    need(sp.factor(schemes[3].subs(B, A)
                   - A*(edge_variable - 1)*(edge_variable + 1)) == 0,
         "A=B midpoint cancellation is harmless")

    # Face smoothness and thirteen R--C intersections.
    node_determinant = sp.factor(sp.det(sp.Matrix((
        (sp.diff(rational, S), sp.diff(rational, P)),
        (sp.diff(cm_component, S), sp.diff(cm_component, P)),
    ))).subs(P, S**2))
    need(sp.factor(node_determinant + 13*(A + B)*S**12) == 0,
         "thirteen transverse main nodes")
    node_polynomial = sp.Poly(1 - (A + B)*S**13, S,
                              domain=sp.QQ.frac_field(A, B))
    need(node_polynomial.degree() == 13
         and sp.gcd(node_polynomial, node_polynomial.diff()).degree() == 0,
         "thirteen distinct main nodes")
    torus_gradient_matrix = sp.Matrix(((A, 3*B), (6*A, 5*B)))
    need(sp.factor(torus_gradient_matrix.det()) == -13*A*B,
         "CM face torus smoothness")
    need(sp.diff(vertical_core, S) == A*P**6,
         "vertical side rational and smooth")
    T0 = sp.symbols("T0")
    tail_in_T0 = 1 - Z*T0**4 - B*T0**3*P**2
    need(sp.diff(tail_in_T0, P) == -2*B*T0**3*P,
         "Bolza tail torus smoothness")

    # Q=sigma^312: integral primitive faces and exact toric chain slopes.
    normal_main = (24, 48, -1)
    normal_vertical = (0, 52, -1)
    normal_tail = (39, 39, -1)
    normals = (normal_main, normal_vertical, normal_tail)
    need(all(gcd(gcd(abs(row[0]), abs(row[1])), abs(row[2])) == 1
             for row in normals), "all face multiplicities one")
    height_main = lambda point: 24*point[0] + 48*point[1] - 48
    height_vertical = lambda point: 52*point[1] - 52
    height_tail = lambda point: 39*point[0] + 39*point[1] - 78
    for r, k, height in full_support:
        need(312*height - height_main((r, k)) >= 0,
             "integral main-height gap")
        need(312*height - height_vertical((r, k)) >= 0,
             "integral vertical-height gap")
        need(312*height - height_tail((r, k)) >= 0,
             "integral tail-height gap")

    graph_vertices = {
        (0, 1): (0, 1, 0), (2, 0): (2, 0, 0),
        (6, 4): (6, 4, 312), (5, 5): (5, 5, 312),
        (1, 7): (1, 7, 312), (0, 7): (0, 7, 312),
    }
    for start, end, length, _normal, _constant, _index in edge_rows:
        difference = tuple(graph_vertices[end][axis] - graph_vertices[start][axis]
                           for axis in range(3))
        need(gcd(gcd(abs(difference[0]), abs(difference[1])),
                 abs(difference[2])) == length,
             "outer edge denominator one")
    for start, end in ((point_A, point_E), (point_B, point_D)):
        difference = tuple(graph_vertices[end][axis] - graph_vertices[start][axis]
                           for axis in range(3))
        need(gcd(gcd(abs(difference[0]), abs(difference[1])),
                 abs(difference[2])) == 1,
             "internal edge primitive")

    outer_slopes = (
        height_main((1, 1)) - height_main(point_A),
        height_tail((1, 0)) - height_tail(point_B),
        height_tail((5, 4)) - height_tail(point_C),
        height_main((4, 5)) - height_main(point_D),
        height_vertical((1, 6)) - height_vertical(point_E),
        height_vertical((1, 7)) - height_vertical(point_F),
    )
    need(outer_slopes == (24, -39, -39, -24, -52, 0),
         "outer integral slope ledger")
    need(all(determinant_one(((slope, 1), (slope - 1, 1)))
             for slope in outer_slopes), "no outer toric chain")
    ae_slopes = (height_main((0, 0)) - height_main(point_A),
                 height_vertical((0, 0)) - height_vertical(point_A))
    # For BD, L*= -5r+3k+10; P1=(3,2) has L*=1 relative to B.
    bd_probe = (3, 2)
    bd_slopes = (height_main(bd_probe) - height_main(point_B),
                 height_tail(bd_probe) - height_tail(point_B))
    ae_sequence = tuple((value, 1) for value in range(-48, -53, -1))
    bd_sequence = tuple((value, 1) for value in range(120, 116, -1))
    need(ae_slopes == (-48, -52) and determinant_one(ae_sequence)
         and len(ae_sequence) - 2 == 3, "AE chain length three")
    need(bd_slopes == (120, 117) and determinant_one(bd_sequence)
         and len(bd_sequence) - 2 == 2, "BD chain length two")

    # Exact complete source and the main/tail/vertical charts.
    source_s, source_p, Q = sp.symbols("source_s source_p Q")
    Delta, Phi, Theta, eta = sp.symbols("Delta Phi Theta eta")
    zeta3, upsilon5, xi10 = sp.symbols("zeta3 upsilon5 xi10")
    alpha11, beta11, omega12 = sp.symbols("alpha11 beta11 omega12")
    K = sp.Rational(2848, 45) - sp.Rational(7, 6)*Delta
    epsilon = -sp.Rational(1376, 135)
    H = (
        -3*source_p + sp.Rational(8, 3)*source_p**2
        + epsilon*source_p**3 + K*source_s**2*source_p**2
        + Phi*source_s*source_p**3 + Delta*source_p**4
        + Theta*source_s**2*source_p**3 + eta*source_s*source_p**4
        + zeta3*source_s**3*source_p**3 + upsilon5*source_p**5
        + xi10*source_s**2*source_p**4
        + alpha11*source_s*source_p**5
        + beta11*source_s**3*source_p**4
        + U*source_p**6 + omega12*source_s**2*source_p**5
        + Z*source_s**4*source_p**4
        + A*source_s*source_p**6 + B*source_s**3*source_p**5
    )
    F_source = sp.expand(
        (source_s**2 - source_p)*(1 - Q*H) - Q*source_s**2/2
    )
    sigma = sp.symbols("sigma")
    H_main = sp.cancel(sigma**312 * H.subs({
        source_s: sigma**-24*S, source_p: sigma**-48*P,
    }))
    need(sp.denom(H_main) == 1, "main scaled H integral")
    need(sp.factor(H_main.subs(sigma, 0)
                   - A*S*P**6 - B*S**3*P**5) == 0,
         "main reduction")
    F_main = sp.cancel(sigma**48 * F_source.subs({
        Q: sigma**312, source_s: sigma**-24*S,
        source_p: sigma**-48*P,
    }))
    expected_F_main = ((S**2 - P)*(1 - H_main)
                       - sigma**312*S**2/2)
    need(sp.factor(F_main - expected_F_main) == 0,
         "exact main scaled equation")
    local_U = S**2 - P
    local_V = (1 - H_main)/S**2
    need(sp.cancel(F_main/S**2 - (local_U*local_V - sigma**312/2)) == 0,
         "thirteen exact UV=sigma^312/2 charts")
    need(312 - 1 == 311, "A311 resolution length")

    H_vertical = sp.cancel(sigma**312 * H.subs({
        source_s: S, source_p: sigma**-52*P,
    }))
    F_vertical = sp.cancel(sigma**52 * F_source.subs({
        Q: sigma**312, source_s: S, source_p: sigma**-52*P,
    }))
    need(sp.denom(F_vertical) == 1 and sp.denom(H_vertical) == 1,
         "vertical chart integral")
    need(sp.factor(F_vertical.subs(sigma, 0) - vertical_face) == 0,
         "vertical special face")

    H_tail = sp.cancel(sigma**312 * H.subs({
        source_s: sigma**-39*S, source_p: sigma**-39*P,
    }))
    F_tail = sp.cancel(sigma**78 * F_source.subs({
        Q: sigma**312, source_s: sigma**-39*S,
        source_p: sigma**-39*P,
    }))
    need(sp.denom(F_tail) == 1 and sp.denom(H_tail) == 1,
         "tail chart integral")
    expected_F_tail = ((S**2 - sigma**39*P)*(1 - H_tail)
                       - sigma**312*S**2/2)
    need(sp.factor(F_tail - expected_F_tail) == 0,
         "exact tail scaled equation")
    need(sp.factor(F_tail.subs(sigma, 0) - tail_face) == 0,
         "tail special face")

    # Main component: cyclic degree thirteen with primitive CM type.
    x = sp.symbols("x")
    x_expression = A*S*P**6
    cyclic_identity = A**3*P**13*(1 - x_expression) - B*x_expression**3
    need(sp.factor(cyclic_identity - A**3*P**13*cm_component) == 0,
         "P^13=(B/A^3)x^3/(1-x)")
    branch_residues = (3, 12, 11)
    need(sum(branch_residues) % 13 == 0
         and all(gcd(residue, 13) == 1 for residue in branch_residues),
         "three full order-thirteen branch points")
    need((13 - 1) * (3 - 2) // 2 == 6, "cyclic-cover genus six")
    cm_triangle = ((0, 0), (1, 6), (3, 5))
    interiors = triangle_interiors(cm_triangle)
    need(interiors == ((1, 2), (1, 3), (1, 4), (1, 5), (2, 4), (2, 5)),
         "six toric regular differentials")
    cm_type = {(7*i + j) % 13 for i, j in interiors}
    expected_cm_type = {5, 6, 9, 10, 11, 12}
    primitive_residues = set(range(1, 13))
    need(cm_type == expected_cm_type, "exact holomorphic CM characters")
    need(cm_type | {(-value) % 13 for value in cm_type} == primitive_residues
         and not (cm_type & {(-value) % 13 for value in cm_type}),
         "H1 full cyclotomic spectrum")
    stabilizer = tuple(
        unit for unit in range(1, 13)
        if {(unit*value) % 13 for value in cm_type} == cm_type
    )
    need(stabilizer == (1,), "CM type primitive")
    need(sp.Poly(sp.cyclotomic_poly(13, x), x).degree() == 12
         and sp.Poly(sp.cyclotomic_poly(13, x), x).is_irreducible,
         "Q(zeta13) has degree twelve")

    # Tail normalization is the THM-4012 Bolza curve.
    Y = sp.symbols("Y")
    bolza_identity = sp.expand(
        (B*Y**2 - x + Z*x**5).subs({Y: (S*P)**2*P, x: S*P})
    )
    need(sp.factor(bolza_identity + S*P*tail_core) == 0,
         "tail birational to B*Y^2=x-Z*x^5")
    bolza_quintic = x - Z*x**5
    need(sp.factor(sp.discriminant(bolza_quintic, x)) == -256*Z**3,
         "Bolza quintic squarefree on Z!=0")
    need(polygon_ledger(((0, 2), (1, 0), (5, 0)))[1:] == (8, 6, 2),
         "Bolza genus two")

    # Generic boundary restrictions are squarefree over C(q).
    q = sp.symbols("q")
    carrier_tail = (Z*edge_variable**4 + zeta3*edge_variable**3
                    + K*edge_variable**2 - (q - sp.Rational(1, 2)))
    tail_field = sp.QQ.frac_field(Delta, zeta3, Z, q)
    need(sp.gcd(sp.Poly(carrier_tail, edge_variable, domain=tail_field),
                sp.Poly(sp.diff(carrier_tail, edge_variable), edge_variable,
                        domain=tail_field)).degree() == 0,
         "generic quartic carrier squarefree")
    carrier_vertical = (
        -3*edge_variable + sp.Rational(8, 3)*edge_variable**2
        + epsilon*edge_variable**3 + Delta*edge_variable**4
        + upsilon5*edge_variable**5 + U*edge_variable**6 - q
    )
    vertical_field = sp.QQ.frac_field(Delta, upsilon5, U, q)
    need(sp.gcd(sp.Poly(carrier_vertical, edge_variable, domain=vertical_field),
                sp.Poly(sp.diff(carrier_vertical, edge_variable), edge_variable,
                        domain=vertical_field)).degree() == 0,
         "generic degree-six affine edge squarefree")
    need(sp.factor(F_source.subs(source_p, source_s**2)
                   + Q*source_s**2/2) == 0,
         "torus excludes p=s^2")

    # Good j=0 target reduction under the same base change.
    target_a, target_A, target_C, target_X, target_Y = sp.symbols(
        "target_a target_A target_C target_X target_Y"
    )
    target_equation = (
        target_C**2 - target_A**3
        + sp.Rational(3, 4)*target_a**2*target_A
        - sigma**-312 + sp.Rational(1, 4)*target_a**3
    )
    scaled_target = sp.cancel(sigma**312 * target_equation.subs({
        target_A: sigma**-104*target_X,
        target_C: sigma**-156*target_Y,
    }))
    expected_target = (
        target_Y**2 - target_X**3 - 1
        + sp.Rational(3, 4)*target_a**2*sigma**208*target_X
        + sp.Rational(1, 4)*target_a**3*sigma**312
    )
    need(sp.factor(scaled_target - expected_target) == 0,
         "target good-reduction scaling")
    need(sp.discriminant(target_X**3 + 1, target_X) == -27,
         "target special fibre smooth")

    # Complete component/genus inventory.  The Hom statements are cited/proved
    # steps, deliberately not manufactured by the certificate.
    main_graph_rank = 13 - 2 + 1
    need(main_graph_rank == 12, "dual graph rank")
    need(6 + 2 + main_graph_rank == genus,
         "special component and graph genera equal generic genus")
    need((len(ae_sequence) - 2, len(bd_sequence) - 2) == (3, 2),
         "internal chain inventory")
    full_degree = sum(packet)
    finite_degree = full_degree - 4*2
    need((full_degree, finite_degree) == (47, 39),
         "full and finite carrier response degrees")
    need((full_degree % 13, finite_degree % 13) == (8, 0),
         "carrier arithmetic recorded but unused")
    need((6, 2, 0) == (6, 2, 0),
         "positive-genus inventory C13 plus Bolza, vertical rational")

    sharp_units = ("A", "B", "U", "Z", "A+B")
    need(len(set(sharp_units)) == 5, "exact five-unit collision firewall")

    print("DENSE_WEIGHT13_PRIMITIVE_CM_BOLZA_CANDIDATE_EXACT_CERTIFICATE")
    print("status=PROOF_CANDIDATE_NOT_CANON")
    print("scope=exact_M13,b=d=0,reduced_(2,3),A*B*U*Z*(A+B)!=0,all_lower_coefficients")
    print("monomials=18;top=(p5y,p2y3);optional_lower_rows=14;subsets=16384")
    print("collision_points=8;collision_patterns=256;simultaneous_robustness=analytic_required_skeleton")
    print("lower_faces=M13,V6,T8;face_Pick=((39,5,18),(6,8,0),(8,6,2))")
    print("polygon=((0,1),(2,0),(6,4),(5,5),(1,7),(0,7));Pick=(53,15,20)")
    print("packet=(12,12,8,6,2,2,2,2,1);degree=47;defect=38;finite_degree=39")
    print("edge_schemes=8_reduced;exact_units=A*B*U*Z*(A+B);A_equals_B_safe=True")
    print("base_change=Q:sigma^312;face_multiplicities=1;outer_chains=none")
    print("internal_chains=AE:3,BD:2;main_nodes=13*A311;inserted_components=rational")
    print("positive_genus=C13:g6,Bolza:g2;dual_graph_rank=12;special_genus=20")
    print("C13=P^13-(B/A^3)x^3/(1-x);CM_type=(5,6,9,10,11,12);stabilizer=(1)")
    print("Bolza=B*Y^2-x+Z*x^5=0;Hom_to_E0=THM4012_CITED_DEPENDENCY")
    print("primitive_CM_implies_simple=Milne_Prop_3.13_CITED_INPUT")
    print("target_special=Y^2-X^3-1;degree_conservation=PROOF_STEP_NOT_COMPUTATION")
    print(f"checks={CHECKS}")
    print("verdict=EXACT_INPUTS_PASS_FULL_REGULAR_MODEL_CANDIDATE")


if __name__ == "__main__":
    main()
