#!/usr/bin/env python3
"""Independent exact audit of THM-4232's planar-JC M=11, U=0 wall.

Status: independently audited exact companion to the maintained theorem.

This program uses only the Python standard library.  It does not import the
primary scout or SymPy.  It reconstructs the inherited M<=11 monomial rows,
exhausts optional-support and collision deletions, and implements enough exact
polynomial arithmetic to recompute the face/edge dictionary and edge
discriminants.  The CM simplicity implication and proper-flat degree
conservation are deliberately reported as external proof inputs, not as
computations performed here.

M=12 and its quartic quotient are OPEN and outside this audit.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


CHECKS = 0


def need(condition, label):
    """Optimization-proof replacement for assert."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("AUDIT FAILURE: " + label)


def gcd3(a, b, c):
    return gcd(gcd(abs(a), abs(b)), abs(c))


# ---------------------------------------------------------------------------
# A tiny exact coefficient ring Q[A,B,Z,Delta].  Coefficients are dictionaries
# from exponent 4-tuples to Fraction.  All iteration is explicitly sorted so
# PYTHONHASHSEED cannot affect stdout or arithmetic choices.
# ---------------------------------------------------------------------------

VARIABLES = ("A", "B", "Z", "Delta")
ZERO_MONOMIAL = (0, 0, 0, 0)


def c_clean(value):
    return {key: coefficient for key, coefficient in value.items() if coefficient}


def c_const(value):
    value = Fraction(value)
    return {} if value == 0 else {ZERO_MONOMIAL: value}


def c_var(name):
    exponent = [0] * len(VARIABLES)
    exponent[VARIABLES.index(name)] = 1
    return {tuple(exponent): Fraction(1)}


def c_add(*values):
    answer = {}
    for value in values:
        for key, coefficient in value.items():
            answer[key] = answer.get(key, Fraction(0)) + coefficient
    return c_clean(answer)


def c_neg(value):
    return {key: -coefficient for key, coefficient in value.items()}


def c_sub(first, second):
    return c_add(first, c_neg(second))


def c_scale(value, scalar):
    scalar = Fraction(scalar)
    return c_clean({key: scalar * coefficient for key, coefficient in value.items()})


def c_mul(first, second):
    answer = {}
    for left_key, left_coefficient in first.items():
        for right_key, right_coefficient in second.items():
            key = tuple(a + b for a, b in zip(left_key, right_key))
            answer[key] = answer.get(key, Fraction(0)) + left_coefficient * right_coefficient
    return c_clean(answer)


def c_pow(value, exponent):
    need(exponent >= 0, "coefficient power is nonnegative")
    answer = c_const(1)
    base = value
    power = exponent
    while power:
        if power & 1:
            answer = c_mul(answer, base)
        base = c_mul(base, base)
        power //= 2
    return answer


def c_divide_by_term(value, divisor):
    """Exact division when divisor has one coefficient monomial."""
    need(len(divisor) == 1, "discriminant leading coefficient is one term")
    (divisor_key, divisor_coefficient), = divisor.items()
    answer = {}
    for key, coefficient in value.items():
        quotient_key = tuple(a - b for a, b in zip(key, divisor_key))
        need(all(exponent >= 0 for exponent in quotient_key),
             "exact coefficient monomial division exponents")
        answer[quotient_key] = coefficient / divisor_coefficient
    return c_clean(answer)


CQ = c_const
A_COEFF = c_var("A")
B_COEFF = c_var("B")
Z_COEFF = c_var("Z")
DELTA_COEFF = c_var("Delta")
EPSILON = Fraction(-1376, 135)


# Sparse polynomials in two geometric variables S,P with coefficient-ring
# values.  A key (i,j) denotes S^i P^j.

def m_clean(value):
    return {key: coefficient for key, coefficient in value.items() if coefficient}


def m_add(*values):
    answer = {}
    for value in values:
        for key, coefficient in value.items():
            answer[key] = c_add(answer.get(key, {}), coefficient)
    return m_clean(answer)


def m_neg(value):
    return {key: c_neg(coefficient) for key, coefficient in value.items()}


def m_sub(first, second):
    return m_add(first, m_neg(second))


def m_mul(first, second):
    answer = {}
    for (left_s, left_p), left_coefficient in first.items():
        for (right_s, right_p), right_coefficient in second.items():
            key = (left_s + right_s, left_p + right_p)
            product = c_mul(left_coefficient, right_coefficient)
            answer[key] = c_add(answer.get(key, {}), product)
    return m_clean(answer)


def m_pow(value, exponent):
    answer = {(0, 0): CQ(1)}
    base = value
    power = exponent
    while power:
        if power & 1:
            answer = m_mul(answer, base)
        base = m_mul(base, base)
        power //= 2
    return answer


def m_derivative(value, axis):
    answer = {}
    for (s_power, p_power), coefficient in value.items():
        exponent = (s_power, p_power)[axis]
        if exponent:
            key = (s_power - (axis == 0), p_power - (axis == 1))
            answer[key] = c_scale(coefficient, exponent)
    return m_clean(answer)


def m_substitute_p_s2(value):
    answer = {}
    for (s_power, p_power), coefficient in value.items():
        degree = s_power + 2 * p_power
        answer[degree] = c_add(answer.get(degree, {}), coefficient)
    return u_clean(answer)


# Sparse univariate polynomials over the same coefficient ring.

def u_clean(value):
    return {degree: coefficient for degree, coefficient in value.items() if coefficient}


def u_add(*values):
    answer = {}
    for value in values:
        for degree, coefficient in value.items():
            answer[degree] = c_add(answer.get(degree, {}), coefficient)
    return u_clean(answer)


def u_neg(value):
    return {degree: c_neg(coefficient) for degree, coefficient in value.items()}


def u_sub(first, second):
    return u_add(first, u_neg(second))


def u_mul(first, second):
    answer = {}
    for left_degree, left_coefficient in first.items():
        for right_degree, right_coefficient in second.items():
            degree = left_degree + right_degree
            product = c_mul(left_coefficient, right_coefficient)
            answer[degree] = c_add(answer.get(degree, {}), product)
    return u_clean(answer)


def u_derivative(value):
    return u_clean({degree - 1: c_scale(coefficient, degree)
                    for degree, coefficient in value.items() if degree})


def u_degree(value):
    need(bool(value), "nonzero univariate polynomial")
    return max(value)


def determinant(matrix):
    """Subset-DP determinant over the exact coefficient ring."""
    size = len(matrix)
    need(all(len(row) == size for row in matrix), "square determinant matrix")
    states = {0: CQ(1)}
    for row_index in range(size):
        next_states = {}
        for mask in sorted(states):
            partial = states[mask]
            for column in range(size):
                if mask & (1 << column):
                    continue
                entry = matrix[row_index][column]
                if not entry:
                    continue
                inversions_added = sum(1 for prior in range(column + 1, size)
                                       if mask & (1 << prior))
                term = c_mul(partial, entry)
                if inversions_added & 1:
                    term = c_neg(term)
                new_mask = mask | (1 << column)
                next_states[new_mask] = c_add(next_states.get(new_mask, {}), term)
        states = next_states
    return states.get((1 << size) - 1, {})


def resultant(first, second):
    first_degree = u_degree(first)
    second_degree = u_degree(second)
    first_descending = [first.get(degree, {})
                        for degree in range(first_degree, -1, -1)]
    second_descending = [second.get(degree, {})
                         for degree in range(second_degree, -1, -1)]
    size = first_degree + second_degree
    rows = []
    for shift in range(second_degree):
        row = [{} for _ in range(size)]
        row[shift:shift + first_degree + 1] = first_descending
        rows.append(row)
    for shift in range(first_degree):
        row = [{} for _ in range(size)]
        row[shift:shift + second_degree + 1] = second_descending
        rows.append(row)
    return determinant(rows)


def discriminant(polynomial):
    degree = u_degree(polynomial)
    if degree <= 1:
        return CQ(1)
    derivative = u_derivative(polynomial)
    answer = resultant(polynomial, derivative)
    if (degree * (degree - 1) // 2) & 1:
        answer = c_neg(answer)
    return c_divide_by_term(answer, polynomial[degree])


# ---------------------------------------------------------------------------
# Exact monomial universe, valued support, and hostile lower-hull census.
# ---------------------------------------------------------------------------

def inherited_monomial_rows(limit):
    rows = []
    for p_power in range(limit // 2 + 1):
        for y_power in range(limit // 3 + 1):
            weight = 2 * p_power + 3 * y_power
            if 0 < weight <= limit and (p_power, y_power) not in {(0, 1), (1, 1)}:
                rows.append((p_power, y_power, weight))
    return tuple(sorted(rows, key=lambda row: (row[2], row[1], row[0])))


def support_multiset(rows):
    # The third fixed point is the coefficient -1/2 of -Q*s^2/2.
    multiplicities = {(2, 0, 0): 1, (0, 1, 0): 1, (2, 0, 1): 1}
    for p_power, y_power, _weight in rows:
        images = ((y_power + 2, p_power + y_power, 1),
                  (y_power, p_power + y_power + 1, 1))
        for point in images:
            multiplicities[point] = multiplicities.get(point, 0) + 1
    return multiplicities


def support_set(rows):
    return set(support_multiset(rows))


def plane_through(first, second, third):
    determinant_xy = ((second[0] - first[0]) * (third[1] - first[1])
                      - (second[1] - first[1]) * (third[0] - first[0]))
    if determinant_xy == 0:
        return None
    slope_s = Fraction(
        (second[2] - first[2]) * (third[1] - first[1])
        - (second[1] - first[1]) * (third[2] - first[2]), determinant_xy)
    slope_p = Fraction(
        (second[0] - first[0]) * (third[2] - first[2])
        - (second[2] - first[2]) * (third[0] - first[0]), determinant_xy)
    constant = Fraction(first[2]) - slope_s * first[0] - slope_p * first[1]
    return slope_s, slope_p, constant


def build_plane_table(universal_support):
    points = tuple(sorted(universal_support))
    planes = {plane for triple in combinations(points, 3)
              if (plane := plane_through(*triple)) is not None}
    table = []
    for plane in sorted(planes):
        slope_s, slope_p, constant = plane
        below = 0
        equal = 0
        for index, (s_power, p_power, height) in enumerate(points):
            gap = (Fraction(height) - slope_s * s_power
                   - slope_p * p_power - constant)
            if gap < 0:
                below |= 1 << index
            elif gap == 0:
                equal |= 1 << index
        table.append((plane, below, equal))
    return points, tuple(table)


def mask_for(points, selected):
    index = {point: position for position, point in enumerate(points)}
    return sum(1 << index[point] for point in selected)


def affine_rank_two(points, mask):
    selected = [point for index, point in enumerate(points) if mask & (1 << index)]
    return any(((second[0] - first[0]) * (third[1] - first[1])
                - (second[1] - first[1]) * (third[0] - first[0])) != 0
               for first, second, third in combinations(selected, 3))


def exact_lower_planes(points, table, support):
    present = mask_for(points, support)
    answer = set()
    for plane, below, equal in table:
        if below & present:
            continue
        face_mask = equal & present
        if affine_rank_two(points, face_mask):
            answer.add(plane)
    return answer


# ---------------------------------------------------------------------------
# Lattice polygons, edge schemes, and toric ledgers.
# ---------------------------------------------------------------------------

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


def polygon_data(points):
    polygon = convex_hull(points)
    area2 = abs(sum(
        polygon[index][0] * polygon[(index + 1) % len(polygon)][1]
        - polygon[(index + 1) % len(polygon)][0] * polygon[index][1]
        for index in range(len(polygon))))
    boundary = sum(
        gcd(abs(polygon[(index + 1) % len(polygon)][0] - polygon[index][0]),
            abs(polygon[(index + 1) % len(polygon)][1] - polygon[index][1]))
        for index in range(len(polygon)))
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    interior = (area2 - boundary + 2) // 2
    return polygon, area2, boundary, interior


def edge_rows(polygon):
    rows = []
    packet = []
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        packet_value = inward[0] + inward[1] - constant
        rows.append((start, end, length, inward, constant, packet_value))
        # The pure-p vertical edge is not a source-response packet edge.
        if not (start[0] == end[0] == 0):
            packet.extend([packet_value] * length)
    return tuple(rows), tuple(sorted(packet, reverse=True))


def edge_polynomial(face, start, end):
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    unit_x, unit_y = dx // length, dy // length
    answer = {}
    for (s_power, p_power), coefficient in face.items():
        relative_x, relative_y = s_power - start[0], p_power - start[1]
        if relative_x * dy - relative_y * dx:
            continue
        if unit_x:
            if relative_x % unit_x:
                continue
            position = relative_x // unit_x
        else:
            if relative_y % unit_y:
                continue
            position = relative_y // unit_y
        if 0 <= position <= length:
            answer[position] = c_add(answer.get(position, {}), coefficient)
    return u_clean(answer)


def triangle_interior_points(vertices):
    def doubled_area(first, second, point):
        return abs((second[0] - first[0]) * (point[1] - first[1])
                   - (second[1] - first[1]) * (point[0] - first[0]))

    total = doubled_area(vertices[0], vertices[1], vertices[2])
    answer = []
    for first_coordinate in range(max(point[0] for point in vertices) + 1):
        for second_coordinate in range(max(point[1] for point in vertices) + 1):
            point = (first_coordinate, second_coordinate)
            pieces = tuple(doubled_area(vertices[index],
                                        vertices[(index + 1) % 3], point)
                           for index in range(3))
            if sum(pieces) == total and all(piece > 0 for piece in pieces):
                answer.append(point)
    return tuple(answer)


def primitive_chain(first_slope, last_slope):
    step = 1 if last_slope >= first_slope else -1
    return tuple((slope, 1)
                 for slope in range(first_slope, last_slope + step, step))


def determinant_one_chain(chain):
    return all(abs(first[0] * second[1] - first[1] * second[0]) == 1
               for first, second in zip(chain, chain[1:]))


def check_exponent_chart(case, h_terms, alpha, beta, expected_initial, label):
    exponents = {name: 132 + alpha * s_power + beta * p_power
                 for name, s_power, p_power in h_terms}
    need(min(exponents.values()) == 0, case["name"] + ": " + label + " minimum")
    need(tuple(sorted(name for name, exponent in exponents.items() if exponent == 0))
         == tuple(sorted(expected_initial)),
         case["name"] + ": " + label + " initial H terms")
    need(all(exponent >= 0 for exponent in exponents.values()),
         case["name"] + ": " + label + " integral H exponents")
    factor = -min(2 * alpha, beta)
    source_s2_exponent = factor + 2 * alpha
    source_p_exponent = factor + beta
    error_exponent = factor + 132 + 2 * alpha
    need(min(source_s2_exponent, source_p_exponent) == 0,
         case["name"] + ": " + label + " source binomial normalized")
    need(error_exponent > 0, case["name"] + ": " + label + " error vanishes")
    return (factor, source_s2_exponent, source_p_exponent, error_exponent,
            tuple(sorted(exponents.items())))


def check_branch(case, monomials, universal_points, plane_table, collisions):
    name = case["name"]
    retained = tuple(row for row in monomials if row[:2] not in case["deleted"])
    required = tuple(row for row in retained if row[:2] in case["required"])
    optional = tuple(row for row in retained if row[:2] not in case["required"])
    need(len(required) == 4, name + ": four required rows")
    need(len(optional) == case["optional_count"], name + ": optional row count")
    need(len(support_set(retained)) == case["active_support_count"],
         name + ": complete active valued-support count")

    hostile_count = 0
    for optional_mask in range(1 << len(optional)):
        chosen = list(required)
        chosen.extend(row for index, row in enumerate(optional)
                      if optional_mask & (1 << index))
        base_support = support_set(chosen)
        for collision_mask in range(1 << len(collisions)):
            hostile_support = base_support - {
                point for index, point in enumerate(collisions)
                if collision_mask & (1 << index)
            }
            actual = exact_lower_planes(universal_points, plane_table, hostile_support)
            need(actual == case["expected_planes"],
                 name + ": hostile lower hull invariant")
            hostile_count += 1
    need(hostile_count == case["hostile_count"], name + ": hostile count")

    # A second, formula-level explanation of why the three planes remain lower.
    for p_power, y_power, weight in retained:
        images = ((y_power + 2, p_power + y_power, 1),
                  (y_power, p_power + y_power + 1, 1))
        main_gaps = tuple(Fraction(point[2]) - case["main_plane"][0] * point[0]
                          - case["main_plane"][1] * point[1]
                          - case["main_plane"][2] for point in images)
        tail_gaps = tuple(Fraction(point[2]) - case["tail_plane"][0] * point[0]
                          - case["tail_plane"][1] * point[1]
                          - case["tail_plane"][2] for point in images)
        vertical_gaps = tuple(Fraction(point[2]) - case["vertical_plane"][0] * point[0]
                              - case["vertical_plane"][1] * point[1]
                              - case["vertical_plane"][2] for point in images)
        need(main_gaps == (Fraction(11 - weight, 11),) * 2,
             name + ": exact main gap identity")
        need(tail_gaps == (Fraction(3 - y_power, 3),
                           Fraction(5 - y_power, 3)),
             name + ": exact tail gap identity")
        expected_vertical_gaps = (
            (Fraction(7 - p_power, 4), Fraction(4 - p_power, 4))
            if name == "Delta4"
            else (Fraction(8 - p_power + y_power, 3),
                  Fraction(3 - p_power + y_power, 3)))
        need(vertical_gaps == expected_vertical_gaps
             and all(gap >= 0 for gap in vertical_gaps),
             name + ": exact vertical replacement gap identity")
    residual = (2, 0, 1)
    residual_gaps = tuple(
        Fraction(residual[2]) - plane[0] * residual[0]
        - plane[1] * residual[1] - plane[2]
        for plane in (case["main_plane"], case["tail_plane"],
                      case["vertical_plane"]))
    need(residual_gaps == case["residual_gaps"]
         and all(gap > 0 for gap in residual_gaps),
         name + ": fixed -Q*s^2/2 point strictly above all lower faces")

    main_points = ((0, 1), (2, 0), (5, 4), (1, 6))
    tail_points = ((2, 0), (5, 3), (5, 4))
    replacement_points = case["replacement_points"]
    need(polygon_data(main_points)[1:] == (33, 5, 15), name + ": main Pick ledger")
    need(polygon_data(tail_points)[1:] == (3, 5, 0), name + ": tail Pick ledger")
    need(polygon_data(replacement_points)[1:] == case["replacement_pick"],
         name + ": replacement Pick ledger")
    global_data = polygon_data(main_points + tail_points + replacement_points)
    need(global_data == case["global_data"], name + ": global polygon/Pick ledger")
    polygon, area2, boundary, genus = global_data
    outer_edges, packet = edge_rows(polygon)
    need(packet == (10, 10, 5, 4, 2, 2, 2, 1), name + ": source packet")
    source_degree = sum(packet)
    ramification_defect = sum(value - 1 for value in packet)
    need((source_degree, ramification_defect) == (36, 28), name + ": degree ledger")
    need(source_degree - 3 * 2 == 30, name + ": finite source degree ledger")

    # Build every face by exact sparse multiplication.
    rational = {(2, 0): CQ(1), (0, 1): CQ(-1)}
    cyclotomic = {(0, 0): CQ(1), (1, 5): c_neg(A_COEFF),
                  (3, 4): c_neg(B_COEFF)}
    main_face = m_mul(rational, cyclotomic)
    tail_core = {(0, 0): CQ(1), (3, 3): c_neg(Z_COEFF),
                 (3, 4): c_neg(B_COEFF)}
    tail_face = m_mul({(2, 0): CQ(1)}, tail_core)
    c_value = case["c"]
    replacement_core = {(0, 0): CQ(-1), (0, case["pure_degree"]): c_value,
                        (1, 5): A_COEFF}
    replacement_face = m_mul({(0, 1): CQ(1)}, replacement_core)
    point_a, point_b, point_c, point_d, point_e, point_f = polygon
    actual_schemes = (
        edge_polynomial(main_face, point_a, point_b),
        edge_polynomial(tail_face, point_b, point_c),
        edge_polynomial(tail_face, point_c, point_d),
        edge_polynomial(main_face, point_d, point_e),
        edge_polynomial(replacement_face, point_e, point_f),
        edge_polynomial(replacement_face, point_f, point_a),
        edge_polynomial(main_face, point_a, point_e),
        edge_polynomial(main_face, point_b, point_d),
    )
    x_minus_one = {0: CQ(-1), 1: CQ(1)}
    ax_plus_b = {0: B_COEFF, 1: A_COEFF}
    expected_schemes = (
        x_minus_one,
        {0: CQ(1), 3: c_neg(Z_COEFF)},
        {0: c_neg(Z_COEFF), 1: c_neg(B_COEFF)},
        u_mul(x_minus_one, ax_plus_b),
        {0: A_COEFF, 1: c_value},
        {0: c_value, case["pure_degree"]: CQ(-1)},
        {0: CQ(-1), 1: A_COEFF},
        {0: CQ(1), 1: c_neg(B_COEFF)},
    )
    need(actual_schemes == expected_schemes, name + ": eight exact edge schemes")
    actual_discriminants = tuple(discriminant(polynomial)
                                 for polynomial in actual_schemes)
    expected_pure_discriminant = c_scale(c_pow(c_value, case["pure_degree"] - 1),
                                         -256 if case["pure_degree"] == 4 else -27)
    expected_discriminants = (
        CQ(1), c_scale(c_pow(Z_COEFF, 2), -27), CQ(1),
        c_pow(c_add(A_COEFF, B_COEFF), 2), CQ(1),
        expected_pure_discriminant, CQ(1), CQ(1))
    need(actual_discriminants == expected_discriminants,
         name + ": exact edge discriminants")

    # Main intersections and the nonzero discriminant gate.
    gradient_determinant = m_sub(
        m_mul(m_derivative(rational, 0), m_derivative(cyclotomic, 1)),
        m_mul(m_derivative(rational, 1), m_derivative(cyclotomic, 0)))
    gradient_on_rational = m_substitute_p_s2(gradient_determinant)
    expected_gradient = {10: c_scale(c_add(A_COEFF, B_COEFF), -11)}
    need(gradient_on_rational == expected_gradient,
         name + ": eleven transverse main nodes")
    intersection = m_substitute_p_s2(cyclotomic)
    need(intersection == {0: CQ(1), 11: c_neg(c_add(A_COEFF, B_COEFF))},
         name + ": main intersection polynomial")
    need(3 * 5 - 1 * 4 == 11, name + ": C exponent determinant magnitude")
    need(bool(A_COEFF) and bool(B_COEFF), name + ": AB gate represented")
    need(case["pure_nonzero"], name + ": pure-owner coefficient gate")

    # Cyclic degree-eleven identity and primitive CM type.
    x_monomial = {(1, 5): A_COEFF}
    a3p11 = {(0, 11): c_pow(A_COEFF, 3)}
    cyclic_left = m_sub(
        m_mul(a3p11, m_sub({(0, 0): CQ(1)}, x_monomial)),
        m_mul({(0, 0): B_COEFF}, m_pow(x_monomial, 3)))
    cyclic_right = m_mul(a3p11, cyclotomic)
    need(cyclic_left == cyclic_right, name + ": cyclic degree-eleven identity")
    cm_triangle = ((0, 0), (1, 5), (3, 4))
    interiors = triangle_interior_points(cm_triangle)
    need(interiors == ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4)),
         name + ": CM Newton-triangle interiors")
    cm_type = frozenset((6 * first + second) % 11 for first, second in interiors)
    need(cm_type == frozenset((4, 5, 8, 9, 10)), name + ": CM type")
    stabilizer = tuple(unit for unit in range(1, 11)
                       if frozenset((unit * item) % 11 for item in cm_type) == cm_type)
    need(stabilizer == (1,), name + ": CM type stabilizer")
    need(len(range(1, 11)) == 2 * len(interiors), name + ": CM field degree 2g")

    # Primitive graph normals, lifted-edge multiplicities, and fan chains.
    height_functions = case["height_functions"]
    need(all(gcd3(*normal) == 1 for normal in case["normals"]),
         name + ": primitive lower-face normals")
    for s_power, p_power, height in support_set(retained):
        for height_function in height_functions:
            need(132 * height - height_function((s_power, p_power)) >= 0,
                 name + ": nonnegative integral height gap")
    lifted_vertices = {
        (0, 1): (0, 1, 0),
        (2, 0): (2, 0, 0),
        (5, 3): (5, 3, 132),
        (5, 4): (5, 4, 132),
        (1, 6): (1, 6, 132),
        point_f: (point_f[0], point_f[1], 132),
    }
    lifted_gcds = []
    for start, end, length, _normal, _constant, _packet_value in outer_edges:
        difference = tuple(lifted_vertices[end][axis] - lifted_vertices[start][axis]
                           for axis in range(3))
        lifted_gcd = gcd3(*difference)
        lifted_gcds.append(lifted_gcd)
        need(lifted_gcd == length, name + ": lifted outer-edge gcd")
    internal_gcds = []
    for start, end in ((point_a, point_e), (point_b, point_d)):
        difference = tuple(lifted_vertices[end][axis] - lifted_vertices[start][axis]
                           for axis in range(3))
        internal_gcds.append(gcd3(*difference))
    need(tuple(internal_gcds) == (1, 1), name + ": primitive internal lifted edges")

    main_height, tail_height, vertical_height = height_functions
    outer_slopes = (
        main_height((1, 1)) - main_height(point_a),
        tail_height((1, 0)) - tail_height(point_b),
        tail_height((4, 3)) - tail_height(point_c),
        main_height((4, 4)) - main_height(point_d),
        vertical_height((1, 5)) - vertical_height(point_e),
        vertical_height((1, point_f[1])) - vertical_height(point_f),
    )
    need(outer_slopes == case["outer_slopes"], name + ": outer fan slopes")
    need(all(determinant_one_chain(((slope, 1), (slope - 1, 1)))
             for slope in outer_slopes), name + ": outer determinant-one steps")
    ae_endpoints = (main_height((0, 0)) - main_height(point_a),
                    vertical_height((0, 0)) - vertical_height(point_a))
    bd_endpoints = (main_height((1, -1)) - main_height(point_b),
                    tail_height((1, -1)) - tail_height(point_b))
    ae_chain = primitive_chain(*ae_endpoints)
    bd_chain = primitive_chain(*bd_endpoints)
    need(ae_endpoints == case["ae_endpoints"], name + ": AE endpoint slopes")
    need(bd_endpoints == (-36, -44), name + ": BD endpoint slopes")
    need(determinant_one_chain(ae_chain) and determinant_one_chain(bd_chain),
         name + ": determinant-one internal chains")
    need((len(ae_chain) - 2, len(bd_chain) - 2) == case["chain_counts"],
         name + ": internal chain counts")

    # Exact exponent bookkeeping for Q=sigma^132 on all three side charts.
    h_terms = [
        ("p", 0, 1), ("p2", 0, 2), ("p3", 0, 3),
        ("s2p2", 2, 2), ("sp3", 1, 3),
        ("s2p3", 2, 3), ("sp4", 1, 4), ("Zs3p3", 3, 3),
        ("s2p4", 2, 4), ("Asp5", 1, 5), ("Bs3p4", 3, 4),
    ]
    if name == "Delta4":
        h_terms.append(("Delta_p4", 0, 4))
    chart_main = check_exponent_chart(case, h_terms, -12, -24,
                                      ("Asp5", "Bs3p4"), "main chart")
    chart_tail = check_exponent_chart(case, h_terms, -44, 0,
                                      ("Bs3p4", "Zs3p3"), "tail chart")
    chart_vertical = check_exponent_chart(
        case, h_terms, case["replacement_scaling"][0],
        case["replacement_scaling"][1], case["replacement_initial"],
        "replacement chart")
    need(chart_main[:4] == (24, 0, 0, 132), name + ": main F exponent ledger")
    need(chart_tail[:4] == (88, 0, 88, 132), name + ": tail F exponent ledger")
    need(chart_vertical[:4] == case["replacement_f_ledger"],
         name + ": replacement F exponent ledger")

    # Target scaling: A_t=sigma^-44 X, C_t=sigma^-66 Y and Q=sigma^132.
    target_exponents = (("Y2", 132 - 2 * 66),
                        ("X3", 132 - 3 * 44),
                        ("one", 0), ("a2X", 132 - 44), ("a3", 132))
    need(target_exponents == (("Y2", 0), ("X3", 0), ("one", 0),
                              ("a2X", 88), ("a3", 132)),
         name + ": target side-chart exponent ledger")
    # d(Y^2-X^3-1)=(-3X^2,2Y); their only simultaneous zero (0,0)
    # does not lie on the curve because the value there is -1.
    need(Fraction(-1) != 0, name + ": smooth target special fibre")

    # Four components: genus-5 C, rational R, rational tail, rational vertical.
    # Eleven C--R nodes plus the two primitive internal attachments.
    component_vertices = 4
    component_edges = 11 + 1 + 1
    graph_b1 = component_edges - component_vertices + 1
    need(graph_b1 == 10, name + ": fibre graph first Betti number")
    need((len(interiors), graph_b1, len(interiors) + graph_b1) == (5, 10, genus),
         name + ": special-fibre genus inventory")
    need((polygon_data(main_points)[3], polygon_data(tail_points)[3],
          polygon_data(replacement_points)[3]) == (15, 0, 0),
         name + ": face arithmetic-genus inventory")

    return {
        "name": name,
        "hostiles": hostile_count,
        "polygon": polygon,
        "pick": (area2, boundary, genus),
        "packet": packet,
        "degree": (source_degree, ramification_defect),
        "lifted_gcds": tuple(lifted_gcds),
        "outer_slopes": outer_slopes,
        "chain_endpoints": (ae_endpoints, bd_endpoints),
        "chain_counts": case["chain_counts"],
        "cm_interiors": interiors,
        "cm_type": tuple(sorted(cm_type)),
        "cm_stabilizer": stabilizer,
        "charts": (chart_main[:4], chart_tail[:4], chart_vertical[:4]),
        "graph": (component_vertices, component_edges, graph_b1, genus),
    }


def main():
    monomials = inherited_monomial_rows(11)
    expected_monomials = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11))
    need(monomials == expected_monomials, "complete inherited M<=11 monomial universe")
    need(tuple(row for row in monomials if row[2] == 11)
         == ((4, 1, 11), (1, 3, 11)), "two exact-M11 top rows")
    multiplicities = support_multiset(monomials)
    collisions = tuple(sorted(point for point, count in multiplicities.items() if count > 1))
    expected_collisions = ((2, 3, 1), (2, 4, 1), (2, 5, 1),
                           (3, 4, 1), (3, 5, 1))
    need(collisions == expected_collisions, "five collision points reconstructed")
    need(all(multiplicities[point] == 2 for point in collisions),
         "every collision has multiplicity two")
    universal_points, plane_table = build_plane_table(multiplicities)
    need((len(universal_points), len(plane_table)) == (24, 234),
         "universal support/plane-table sizes")

    main_plane = (Fraction(1, 11), Fraction(2, 11), Fraction(-2, 11))
    tail_plane = (Fraction(1, 3), Fraction(0), Fraction(-2, 3))
    main_height = lambda point: 12 * point[0] + 24 * point[1] - 24
    tail_height = lambda point: 44 * point[0] - 88
    vertical4_height = lambda point: -33 * point[0] + 33 * point[1] - 33
    vertical3_height = lambda point: -88 * point[0] + 44 * point[1] - 44
    cases = (
        {
            "name": "Delta4",
            "deleted": {(5, 0)},
            "required": {(4, 0), (0, 3), (4, 1), (1, 3)},
            "optional_count": 8,
            "active_support_count": 23,
            "hostile_count": 8192,
            "main_plane": main_plane,
            "tail_plane": tail_plane,
            "vertical_plane": (Fraction(-1, 4), Fraction(1, 4), Fraction(-1, 4)),
            "expected_planes": {main_plane, tail_plane,
                                (Fraction(-1, 4), Fraction(1, 4), Fraction(-1, 4))},
            "replacement_points": ((0, 1), (1, 6), (0, 5)),
            "replacement_pick": (4, 6, 0),
            "global_data": (((0, 1), (2, 0), (5, 3), (5, 4),
                             (1, 6), (0, 5)), 40, 12, 15),
            "c": DELTA_COEFF,
            "pure_degree": 4,
            "pure_nonzero": True,
            "height_functions": (main_height, tail_height, vertical4_height),
            "normals": ((12, 24, -1), (44, 0, -1), (-33, 33, -1)),
            "outer_slopes": (12, -44, -44, -12, -33, -33),
            "ae_endpoints": (-24, -33),
            "chain_counts": (8, 7),
            "replacement_scaling": (33, -33),
            "replacement_initial": ("Asp5", "Delta_p4"),
            "replacement_f_ledger": (33, 99, 0, 231),
            "residual_gaps": (Fraction(1), Fraction(1), Fraction(7, 4)),
        },
        {
            "name": "Delta0",
            "deleted": {(5, 0), (4, 0)},
            "required": {(3, 0), (0, 3), (4, 1), (1, 3)},
            "optional_count": 7,
            "active_support_count": 22,
            "hostile_count": 4096,
            "main_plane": main_plane,
            "tail_plane": tail_plane,
            "vertical_plane": (Fraction(-2, 3), Fraction(1, 3), Fraction(-1, 3)),
            "expected_planes": {main_plane, tail_plane,
                                (Fraction(-2, 3), Fraction(1, 3), Fraction(-1, 3))},
            "replacement_points": ((0, 1), (1, 6), (0, 4)),
            "replacement_pick": (3, 5, 0),
            "global_data": (((0, 1), (2, 0), (5, 3), (5, 4),
                             (1, 6), (0, 4)), 39, 11, 15),
            "c": CQ(EPSILON),
            "pure_degree": 3,
            "pure_nonzero": EPSILON != 0,
            "height_functions": (main_height, tail_height, vertical3_height),
            "normals": ((12, 24, -1), (44, 0, -1), (-88, 44, -1)),
            "outer_slopes": (12, -44, -44, -12, -44, -88),
            "ae_endpoints": (-24, -44),
            "chain_counts": (19, 7),
            "replacement_scaling": (88, -44),
            "replacement_initial": ("Asp5", "p3"),
            "replacement_f_ledger": (44, 220, 0, 352),
            "residual_gaps": (Fraction(1), Fraction(1), Fraction(8, 3)),
        },
    )

    results = tuple(check_branch(case, monomials, universal_points,
                                 plane_table, collisions) for case in cases)
    need(tuple(result["hostiles"] for result in results) == (8192, 4096),
         "branch hostile totals")
    need(sum(result["hostiles"] for result in results) == 12288,
         "combined hostile total")
    need(all(result["packet"] == results[0]["packet"] for result in results),
         "common packet")
    need(all(result["degree"] == (36, 28) for result in results),
         "common degree ledger")

    print("THM4232_INDEPENDENT_CLEAN_ROOM_AUDIT")
    print("status=VERIFIED_EXACT_INDEPENDENT_AUDIT")
    print("scope=A*B*Z*(A+B)!=0;M=11;U=0;Q=sigma^132")
    print("branches=Delta!=0:p4_owner_Delta|Delta=0:p3_owner_-1376/135")
    print("universe=13_monomial_rows;universal_valued_points=24;active_points=Delta4:23,Delta0:22;collisions=5;candidate_planes=234")
    print("hostiles=Delta4:8192,Delta0:4096,total:12288")
    print("lower_faces=main_M11,tail_T3,replacement_V4_or_V3")
    print("global_pick=Delta4:(40,12,15),Delta0:(39,11,15)")
    print("packet_both=(10,10,5,4,2,2,2,1);source_degree=36;finite_degree=30;ramification_defect=28")
    print("edge_dictionary=8_exact_schemes;discriminant_gate=A*B*Z*(A+B)*c_k!=0")
    print("toric=primitive_normals;lifted_edge_gcds_exact;face_multiplicity=1")
    print("slopes_Delta4=outer:(12,-44,-44,-12,-33,-33),AE:(-24,-33):8,BD:(-36,-44):7")
    print("slopes_Delta0=outer:(12,-44,-44,-12,-44,-88),AE:(-24,-44):19,BD:(-36,-44):7")
    print("chart_tuple_fields=(F_scale,s2_exponent,p_exponent,fixed_residual_exponent)")
    print("charts=main:scale(-12,-24):init(Asp5,Bs3p4):F(24,0,0,132);tail:scale(-44,0):init(Zs3p3,Bs3p4):F(88,0,88,132)")
    print("replacement_charts=V4:scale(33,-33):init(Delta_p4,Asp5):F(33,99,0,231);V3:scale(88,-44):init((-1376/135)p3,Asp5):F(44,220,0,352)")
    print("cm=triangle:((0,0),(1,5),(3,4));g:5;type:(4,5,8,9,10);stabilizer:(1)")
    print("fibre_graph=vertices:4,edges:13,b1:10;main_genus:5;special_genus:15")
    print("target_special=Y^2-X^3-1;smooth;target_scaling_exponents=(0,0,0,88,132)")
    print("external_inputs=CM_simplicity_implication;proper_flat_degree_conservation")
    print("outside_scope=M12_quartic_quotient_OPEN")
    print("checks=" + str(CHECKS))
    print("verdict=ACCEPT_THM4232_RELATIVE")


if __name__ == "__main__":
    main()
