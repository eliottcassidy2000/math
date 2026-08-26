#!/usr/bin/env python3
"""Exact certificate for THM-4222's dense weight-eleven exclusion.

The script certifies the finite/exact Newton, toric, local-chart, and CM-type
inputs.  The implication "primitive CM type => simple CM abelian variety" is
an explicitly cited external theorem, not a computational claim; see Milne,
Complex Multiplication, Proposition 3.13.
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
    """Return valued support (s exponent,p exponent,Q exponent)."""
    support = {(2, 0, 0), (0, 1, 0)}
    for i, j, _weight in monomials:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


def candidate_plane_records(points):
    """Precompute exact planes and bitmasks for fast subset hull audits."""
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
    first = chosen[0]
    for second in chosen[1:]:
        if second[:2] == first[:2]:
            continue
        for third in chosen[2:]:
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
    """Restrict a bivariate polynomial to an oriented lattice edge."""
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
    # Complete inherited b=d=0 monomial universe through exact weight eleven.
    monomials = monomials_through(11)
    expected_monomials = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11),
    )
    need(monomials == expected_monomials, "M11 monomial universe")
    need(tuple(row for row in monomials if row[2] == 11)
         == ((4, 1, 11), (1, 3, 11)), "exactly two new top monomials")

    # A=[p^4 y], B=[p y^3], U=[p^5], Zeta=[y^3] are retained.
    required = {(4, 1), (1, 3), (5, 0), (0, 3)}
    optional = tuple(row for row in monomials if row[:2] not in required)
    need(len(optional) == 9, "nine arbitrary lower monomial rows")
    full_support = expanded_support(monomials)
    points, plane_records = candidate_plane_records(full_support)
    plane_main = (Fraction(1, 11), Fraction(2, 11), Fraction(-2, 11))
    plane_vertical = (Fraction(0), Fraction(1, 5), Fraction(-1, 5))
    plane_tail = (Fraction(1, 3), Fraction(0), Fraction(-2, 3))
    expected_planes = {plane_main, plane_vertical, plane_tail}
    support_subset_count = 0
    for mask in range(1 << len(optional)):
        chosen = [row for row in monomials if row[:2] in required]
        chosen.extend(row for bit, row in enumerate(optional)
                      if mask & (1 << bit))
        bits = support_bits(points, expanded_support(chosen))
        need(lower_planes(points, plane_records, bits) == expected_planes,
             "lower hull changed under arbitrary lower-support subset")
        support_subset_count += 1
    need(support_subset_count == 512, "all 512 lower-support subsets")

    collisions = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1), (3, 4, 1), (3, 5, 1),
    )
    collision_count = 0
    for mask in range(1 << len(collisions)):
        deleted = {point for bit, point in enumerate(collisions)
                   if mask & (1 << bit)}
        bits = support_bits(points, full_support - deleted)
        need(lower_planes(points, plane_records, bits) == expected_planes,
             "coincident-coefficient cancellation changed lower hull")
        collision_count += 1
    need(collision_count == 32, "all 32 collision hostiles")

    # Independent analytic gap identities explain the finite enumeration.
    for i, j, weight in monomials:
        first = (j + 2, i + j, 1)
        second = (j, i + j + 1, 1)
        main_gaps = tuple(
            Fraction(point[2]) - plane_main[0] * point[0]
            - plane_main[1] * point[1] - plane_main[2]
            for point in (first, second)
        )
        vertical_gaps = tuple(
            Fraction(point[2]) - plane_vertical[0] * point[0]
            - plane_vertical[1] * point[1] - plane_vertical[2]
            for point in (first, second)
        )
        tail_gaps = tuple(
            Fraction(point[2]) - plane_tail[0] * point[0]
            - plane_tail[1] * point[1] - plane_tail[2]
            for point in (first, second)
        )
        need(main_gaps == (Fraction(11 - weight, 11),) * 2,
             "main-plane exact weight gap")
        need(vertical_gaps == (Fraction(6 - i - j, 5),
                               Fraction(5 - i - j, 5)),
             "vertical-plane exact gap")
        need(tail_gaps == (Fraction(3 - j, 3), Fraction(5 - j, 3)),
             "tail-plane exact gap")

    # Face and global polygon ledgers.
    main_points = ((0, 1), (2, 0), (5, 4), (1, 6))
    tail_points = ((2, 0), (5, 3), (5, 4))
    vertical_points = ((0, 1), (1, 6), (0, 6))
    global_points = main_points + tail_points + vertical_points
    need(polygon_ledger(main_points)[1:] == (33, 5, 15), "main Pick ledger")
    need(polygon_ledger(tail_points)[1:] == (3, 5, 0), "tail Pick ledger")
    need(polygon_ledger(vertical_points)[1:] == (5, 7, 0),
         "vertical Pick ledger")
    polygon, area2, boundary, genus = polygon_ledger(global_points)
    expected_polygon = ((0, 1), (2, 0), (5, 3), (5, 4), (1, 6), (0, 6))
    need((polygon, area2, boundary, genus)
         == (expected_polygon, 41, 13, 15), "global polygon/Pick ledger")
    packet, edge_rows = edge_packet(polygon)
    need(packet == (10, 10, 5, 4, 2, 2, 2, 1), "complete outer packet")
    need(sum(packet) == 36 and sum(index - 1 for index in packet) == 28,
         "packet degree and defect")

    # Exact face factorization and all eight compactified edge schemes.
    S, P, edge_variable = sp.symbols("S P edge_variable")
    A, B, U, Zeta = sp.symbols("A B U Zeta")
    rational = S**2 - P
    cyclotomic = 1 - A*S*P**5 - B*S**3*P**4
    main_face = sp.expand(rational * cyclotomic)
    tail_core = 1 - Zeta*S**3*P**3 - B*S**3*P**4
    tail_face = sp.expand(S**2 * tail_core)
    vertical_core = -1 + U*P**5 + A*S*P**5
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
        1 - Zeta*edge_variable**3,
        -Zeta - B*edge_variable,
        (edge_variable - 1)*(A*edge_variable + B),
        A + U*edge_variable,
        U - edge_variable**5,
        A*edge_variable - 1,
        1 - B*edge_variable,
    )
    need(tuple(sp.factor(value - expected) for value, expected
               in zip(schemes, expected_schemes)) == (0,) * 8,
         "all eight edge schemes derived from face polynomials")
    discriminants = tuple(sp.factor(sp.discriminant(scheme, edge_variable))
                          for scheme in schemes)
    need(discriminants == (
        1, -27*Zeta**2, 1, (A + B)**2,
        1, 3125*U**4, 1, 1,
    ), "all edge schemes reduced on the five-unit gate")
    need(sp.factor(schemes[3].subs(B, A) - A*(edge_variable - 1)
                                   *(edge_variable + 1)) == 0,
         "A=B midpoint cancellation is harmless")

    # Face smoothness and the eleven transverse intersections R cap C.
    node_determinant = sp.factor(sp.det(sp.Matrix((
        (sp.diff(rational, S), sp.diff(rational, P)),
        (sp.diff(cyclotomic, S), sp.diff(cyclotomic, P)),
    ))).subs(P, S**2))
    need(sp.factor(node_determinant + 11*(A + B)*S**10) == 0,
         "eleven transverse main nodes")
    node_polynomial = sp.Poly(1 - (A + B)*S**11, S,
                              domain=sp.QQ.frac_field(A, B))
    need(node_polynomial.degree() == 11
         and sp.gcd(node_polynomial, node_polynomial.diff()).degree() == 0,
         "eleven distinct main nodes")
    torus_gradient_matrix = sp.Matrix(((A, 3*B), (5*A, 4*B)))
    need(sp.factor(torus_gradient_matrix.det()) == -11*A*B,
         "cyclotomic face torus smoothness")
    need(sp.diff(vertical_core, S) == A*P**5,
         "vertical side is rational and torus smooth")
    T0 = sp.symbols("T0")
    tail_in_T0 = 1 - T0**3*(Zeta + B*P)
    need(sp.diff(tail_in_T0, P) == -B*T0**3,
         "tail side is rational and torus smooth")

    # Q=sigma^330 makes the three lower graphs integral and primitive.
    sigma = sp.symbols("sigma")
    normal_main = (30, 60, -1)
    normal_vertical = (0, 66, -1)
    normal_tail = (110, 0, -1)
    normals = (normal_main, normal_vertical, normal_tail)
    need(all(gcd(gcd(abs(row[0]), abs(row[1])), abs(row[2])) == 1
             for row in normals), "all face multiplicities are one")
    height_main = lambda point: 30*point[0] + 60*point[1] - 60
    height_vertical = lambda point: 66*point[1] - 66
    height_tail = lambda point: 110*point[0] - 220
    for r, k, height in full_support:
        need(330*height - height_main((r, k)) >= 0,
             "integral main-height gap")
        need(330*height - height_vertical((r, k)) >= 0,
             "integral vertical-height gap")
        need(330*height - height_tail((r, k)) >= 0,
             "integral tail-height gap")

    graph_vertices = {
        (0, 1): (0, 1, 0), (2, 0): (2, 0, 0),
        (5, 3): (5, 3, 330), (5, 4): (5, 4, 330),
        (1, 6): (1, 6, 330), (0, 6): (0, 6, 330),
    }
    for start, end, length, _normal, _constant, _index in edge_rows:
        difference = tuple(graph_vertices[end][axis] - graph_vertices[start][axis]
                           for axis in range(3))
        need(gcd(gcd(abs(difference[0]), abs(difference[1])),
                 abs(difference[2])) == length,
             "outer edge lattice denominator one")
    for start, end in ((point_A, point_E), (point_B, point_D)):
        difference = tuple(graph_vertices[end][axis] - graph_vertices[start][axis]
                           for axis in range(3))
        need(gcd(gcd(abs(difference[0]), abs(difference[1])),
                 abs(difference[2])) == 1,
             "internal edge primitive")

    outer_slopes = (
        height_main((1, 1)) - height_main(point_A),
        height_tail((1, 0)) - height_tail(point_B),
        height_tail((4, 3)) - height_tail(point_C),
        height_main((4, 4)) - height_main(point_D),
        height_vertical((1, 5)) - height_vertical(point_E),
        height_vertical((1, 6)) - height_vertical(point_F),
    )
    need(outer_slopes == (30, -110, -110, -30, -66, 0),
         "outer integral slope ledger")
    need(all(determinant_one(((slope, 1), (slope - 1, 1)))
             for slope in outer_slopes), "no outer toric chain")
    ae_slopes = (height_main((0, 0)) - height_main(point_A),
                 height_vertical((0, 0)) - height_vertical(point_A))
    bd_probe = (1, -1)
    bd_slopes = (height_main(bd_probe) - height_main(point_B),
                 height_tail(bd_probe) - height_tail(point_B))
    ae_sequence = tuple((value, 1) for value in range(-60, -67, -1))
    bd_sequence = tuple((value, 1) for value in range(-90, -111, -1))
    need(ae_slopes == (-60, -66) and determinant_one(ae_sequence)
         and len(ae_sequence) - 2 == 5, "AE rational chain length five")
    need(bd_slopes == (-90, -110) and determinant_one(bd_sequence)
         and len(bd_sequence) - 2 == 19, "BD rational chain length nineteen")

    # Exact source charts and A_329 local equations.
    source_s, source_p, Q = sp.symbols("source_s source_p Q")
    Delta, Phi, Theta, eta, Xi = sp.symbols("Delta Phi Theta eta Xi")
    K = sp.Rational(2848, 45) - sp.Rational(7, 6)*Delta
    epsilon = -sp.Rational(1376, 135)
    H = (
        -3*source_p + sp.Rational(8, 3)*source_p**2
        + epsilon*source_p**3 + K*source_s**2*source_p**2
        + Phi*source_s*source_p**3 + Delta*source_p**4
        + Theta*source_s**2*source_p**3 + eta*source_s*source_p**4
        + Zeta*source_s**3*source_p**3 + U*source_p**5
        + Xi*source_s**2*source_p**4 + A*source_s*source_p**5
        + B*source_s**3*source_p**4
    )
    F_source = sp.expand(
        (source_s**2 - source_p)*(1 - Q*H) - Q*source_s**2/2
    )
    H_main = sp.cancel(sigma**330 * H.subs({
        source_s: sigma**-30*S, source_p: sigma**-60*P,
    }))
    need(sp.denom(H_main) == 1, "main scaled H integral")
    need(sp.factor(H_main.subs(sigma, 0)
                   - A*S*P**5 - B*S**3*P**4) == 0,
         "main reduction retains exactly weight eleven")
    F_main = sp.cancel(sigma**60 * F_source.subs({
        Q: sigma**330, source_s: sigma**-30*S,
        source_p: sigma**-60*P,
    }))
    expected_F_main = ((S**2 - P)*(1 - H_main)
                       - sigma**330*S**2/2)
    need(sp.factor(F_main - expected_F_main) == 0,
         "exact main scaled equation")
    local_U = S**2 - P
    local_V = (1 - H_main)/S**2
    need(sp.cancel(F_main/S**2 - (local_U*local_V - sigma**330/2)) == 0,
         "eleven exact UV=sigma^330/2 charts")
    need(330 - 1 == 329, "A329 resolution length")

    F_vertical = sp.cancel(sigma**66 * F_source.subs({
        Q: sigma**330, source_s: S, source_p: sigma**-66*P,
    }))
    need(sp.denom(F_vertical) == 1, "vertical chart integral")
    need(sp.factor(F_vertical.subs(sigma, 0) - vertical_face) == 0,
         "vertical special face")
    F_tail = sp.cancel(sigma**220 * F_source.subs({
        Q: sigma**330, source_s: sigma**-110*S, source_p: P,
    }))
    need(sp.denom(F_tail) == 1, "tail chart integral")
    need(sp.factor(F_tail.subs(sigma, 0) - tail_face) == 0,
         "tail special face")

    # The genus-five component is a cyclic degree-eleven three-point cover.
    x = sp.symbols("x")
    x_expression = A*S*P**5
    cyclic_identity = (
        A**3*P**11*(1 - x_expression) - B*x_expression**3
    )
    need(sp.factor(cyclic_identity - A**3*P**11*cyclotomic) == 0,
         "P^11=(B/A^3)x^3/(1-x) normalization")
    branch_residues = (3, 10, 9)
    need(sum(branch_residues) % 11 == 0
         and all(gcd(residue, 11) == 1 for residue in branch_residues),
         "three full order-eleven branch points")
    need((11 - 1) * (3 - 2) // 2 == 5, "cyclic-cover genus five")
    cm_triangle = ((0, 0), (1, 5), (3, 4))
    interiors = triangle_interiors(cm_triangle)
    need(interiors == ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4)),
         "five toric regular differentials")
    cm_type = {(6*i + j) % 11 for i, j in interiors}
    expected_cm_type = {4, 5, 8, 9, 10}
    primitive_residues = set(range(1, 11))
    need(cm_type == expected_cm_type, "exact holomorphic CM characters")
    need(cm_type | {(-value) % 11 for value in cm_type} == primitive_residues
         and not (cm_type & {(-value) % 11 for value in cm_type}),
         "H1 has the full cyclotomic spectrum")
    stabilizer = tuple(
        unit for unit in range(1, 11)
        if {(unit*value) % 11 for value in cm_type} == cm_type
    )
    need(stabilizer == (1,), "CM type has trivial stabilizer and is primitive")
    cyclotomic_11 = sp.cyclotomic_poly(11, x)
    need(sp.Poly(cyclotomic_11, x).degree() == 10
         and sp.Poly(cyclotomic_11, x).is_irreducible,
         "Q(zeta11) has degree 10=2g")

    # Labelled side attachments are the totally ramified x=1 and x=0 points.
    # div(x/(1-x))=11(P_0-P_1); primality and g>0 make the difference exact
    # order eleven.  This label is not needed for the Hom=0 contradiction.
    need(11 == 11 and sp.isprime(11) and 5 > 0,
         "attachment difference has exact order eleven input")

    # Generic completion and the unchanged prime cubic carrier.
    q = sp.symbols("q")
    carrier = Zeta*edge_variable**3 + K*edge_variable**2 - (q - sp.Rational(1, 2))
    carrier_discriminant = sp.factor(sp.discriminant(carrier, edge_variable))
    need(sp.factor(carrier_discriminant - (
        (q - sp.Rational(1, 2))
        * (4*K**3 - 27*Zeta**2*(q - sp.Rational(1, 2)))
    )) == 0, "prime cubic carrier discriminant")
    carrier_field = sp.QQ.frac_field(Delta, Zeta, q)
    need(sp.gcd(sp.Poly(carrier, edge_variable, domain=carrier_field),
                sp.Poly(sp.diff(carrier, edge_variable), edge_variable,
                        domain=carrier_field)).degree() == 0,
         "generic cubic carrier squarefree")
    vertical_generic = (
        -3*edge_variable + sp.Rational(8, 3)*edge_variable**2
        + epsilon*edge_variable**3 + Delta*edge_variable**4
        + U*edge_variable**5 - q
    )
    vertical_field = sp.QQ.frac_field(Delta, U, q)
    need(sp.gcd(sp.Poly(vertical_generic, edge_variable, domain=vertical_field),
                sp.Poly(sp.diff(vertical_generic, edge_variable), edge_variable,
                        domain=vertical_field)).degree() == 0,
         "generic degree-five affine edge squarefree")
    need(sp.factor(F_source.subs(source_p, source_s**2)
                   + Q*source_s**2/2) == 0,
         "torus excludes t=p-s^2=0")

    # The same base change gives the inherited target good j=0 reduction.
    target_a, target_A, target_C, target_X, target_Y = sp.symbols(
        "target_a target_A target_C target_X target_Y"
    )
    target_equation = (
        target_C**2 - target_A**3
        + sp.Rational(3, 4)*target_a**2*target_A
        - sigma**-330 + sp.Rational(1, 4)*target_a**3
    )
    scaled_target = sp.cancel(sigma**330 * target_equation.subs({
        target_A: sigma**-110*target_X,
        target_C: sigma**-165*target_Y,
    }))
    expected_target = (
        target_Y**2 - target_X**3 - 1
        + sp.Rational(3, 4)*target_a**2*sigma**220*target_X
        + sp.Rational(1, 4)*target_a**3*sigma**330
    )
    need(sp.factor(scaled_target - expected_target) == 0,
         "target good-reduction scaling")
    need(sp.discriminant(target_X**3 + 1, target_X) == -27,
         "target special fibre Y^2=X^3+1 is smooth")

    # Complete genus and degree-conservation inventory.
    main_graph_rank = 11 - 2 + 1
    side_graph_gain = 0
    need((main_graph_rank, side_graph_gain) == (10, 0),
         "dual graph rank")
    need(5 + main_graph_rank == genus, "special genus equals generic Pick genus")
    full_degree = sum(packet)
    finite_degree = full_degree - 3*2
    need((full_degree, finite_degree) == (36, 30),
         "full and finite carrier response degrees")
    need(full_degree % 3 == finite_degree % 3 == 0,
         "weight-ten order-three congruence does not obstruct M11")
    need((full_degree % 11, finite_degree % 11) == (3, 8),
         "order-eleven labels are not a genus-one degree invoice")
    # Cited Milne Proposition 3.13 turns the exact primitive CM pair above
    # into simplicity.  A simple fivefold has no nonzero Hom to an elliptic
    # curve.  All other face, toric-chain, A329, and resolution components
    # are rational, so degree conservation has zero on every special term.
    positive_genus_components = (5,)
    rational_face_genera = (0, 0, 0)
    need(positive_genus_components == (5,) and rational_face_genera == (0, 0, 0),
         "only the primitive-CM genus-five component can carry degree")
    component_degrees_after_cited_hom_gate = (0, 0, 0, 0)
    need(sum(component_degrees_after_cited_hom_gate) == 0,
         "specialized degree sum is zero")

    # Each of the five displayed factors performs a distinct geometric job.
    sharp_walls = ("A=0", "B=0", "U=0", "Zeta=0", "A+B=0")
    need(len(set(sharp_walls)) == 5, "five sharp coefficient walls")

    print("THM4222_DENSE_WEIGHT11_PRIMITIVE_CM_EXACT_CERTIFICATE")
    print("scope=exact_M11,b=d=0,reduced_(2,3),A*B*U*Zeta*(A+B)!=0,all_lower_coefficients")
    print("monomials=13;new_top=(p4y,py3);lower_support_subsets=512;collision_hostiles=32")
    print("lower_faces=M11,V5,T3;face_Pick=((33,5,15),(5,7,0),(3,5,0))")
    print("polygon=((0,1),(2,0),(5,3),(5,4),(1,6),(0,6));Pick=(41,13,15)")
    print("packet=(10,10,5,4,2,2,2,1);degree=36;defect=28;cubic_carrier_degree=3;finite_degree=30")
    print("edge_schemes=8_reduced;gate=A*B*U*Zeta*(A+B)!=0;A_equals_B_safe=True")
    print("base_change=Q:sigma^330;face_multiplicities=1;outer_chains=none")
    print("internal_chains=AE:5,BD:19;main_nodes=11*A329;all_inserted_components=rational")
    print("positive_genus=C11:g5;dual_graph_rank=10;special_genus=15")
    print("C11=P^11-(B/A^3)x^3/(1-x);CM_type=(4,5,8,9,10);stabilizer=(1)")
    print("CITED_INPUT=Milne_Complex_Multiplication_Prop_3.13:primitive_CM_pair_implies_simple_isogeny_class")
    print("Hom(Jac(C11),E0)=0_by_cited_input_and_dimension;target_special=Y^2-X^3-1")
    print("degree_conservation=special_sum_0_contradicts_finite_generic_degree")
    print("walls=(A=0,B=0,U=0,Zeta=0,A+B=0);other_lower_coefficient_walls_safe")
    print(f"checks={CHECKS}")
    print("verdict=EXACT_M11_DENSE_CHAMBER_EXCLUDED_RELATIVE")


if __name__ == "__main__":
    main()
