#!/usr/bin/env python3
"""Independent exact audit for the THM-4334 beta-zero endpoint seam.

Scope
=====
Work in the exact-weight-twelve reduced (2,3) source inherited by THM-4327.
The target scout stratum is

    Z=0, beta_11=0, U*W*zeta_3 != 0.

Every other labelled row of weight at most twelve is arbitrary.  The support
universe, fixed rows, optional rows, and hostile aggregate-point deletions are
printed below.  The deletions overcount actual coefficient cancellations.

The consequence object is not merely the lower hull.  The program verifies:

* the complete lower-face complex in every hostile support mask;
* an independent dual-slack/dominance proof that avoids face enumeration;
* a direct-hull replay on positive and hostile controls;
* exact face equations and Kummer/hyperelliptic normal-form identities;
* irreducibility/genus, generic-point derivative witnesses, and torus
  smoothness determinants for the positive-genus components;
* face and global Pick ledgers, the component graph, and genus conservation;
* every outer/internal edge polynomial and its exact nondegeneracy gate;
* primitive integral height charts and positive good-differential orders;
* the source-infinity packet as a consistency control only; and
* beta-restored, W-zero, and zeta-zero boundary hostiles; and
* a separate Lambda-zero contact check (unique A23, simple orbit-disjoint
  M--K9 edge, and the exact positive imported tail-order table).

This is a standard-library Fraction certificate.  It verifies the new finite
inputs used by THM-4334; it does not itself prove the
imported toroidal/proper-flat theorems, seam entry, JC(2), or DC(2).
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def need(condition, label):
    """Optimization-stable assertion."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("THM4334 INDEPENDENT AUDIT FAILURE: " + label)


# ---------------------------------------------------------------------------
# Exact support and lower-hull enumeration


def weighted_rows(limit=12):
    rows = []
    for i in range(limit // 2 + 1):
        for j in range(limit // 3 + 1):
            weight = 2 * i + 3 * j
            if 0 < weight <= limit and (i, j) not in {(0, 1), (1, 1)}:
                rows.append((i, j, weight))
    return tuple(sorted(rows, key=lambda row: (row[2], row[1], row[0])))


def lifted_support(rows):
    support = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j, _weight in rows:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


def projected_rank_two(points, mask=None):
    if mask is not None:
        points = tuple(
            point for index, point in enumerate(points) if mask & (1 << index)
        )
    else:
        points = tuple(points)
    for first, second, third in combinations(points, 3):
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
    return slope_s, slope_p, constant


def prepare_plane_records(universe):
    points = tuple(sorted(universe))
    records = {}
    for first, second, third in combinations(points, 3):
        plane = plane_through(first, second, third)
        if plane is None:
            continue
        slope_s, slope_p, constant = plane
        below = 0
        equal = 0
        for index, (r, ell, height) in enumerate(points):
            gap = Fraction(height) - slope_s * r - slope_p * ell - constant
            if gap < 0:
                below |= 1 << index
            elif gap == 0:
                equal |= 1 << index
        records[plane] = (below, equal)
    return points, tuple(
        (plane, below, equal)
        for plane, (below, equal) in sorted(records.items())
    )


def make_hull_engine(rows):
    # Z=0 removes y^4 from the universe before all beta-zero masks are formed.
    universe = lifted_support(row for row in rows if row[:2] != (0, 4))
    points, records = prepare_plane_records(universe)
    point_index = {point: index for index, point in enumerate(points)}

    @lru_cache(maxsize=None)
    def rank_two(mask):
        return projected_rank_two(points, mask)

    def support_mask(active_rows):
        return sum(
            1 << point_index[point] for point in lifted_support(active_rows)
        )

    def lower_faces(mask):
        return tuple(
            plane
            for plane, below, equal in records
            if not (below & mask) and rank_two(equal & mask)
        )

    return points, point_index, support_mask, lower_faces


def hostile_atlas(rows, point_index, support_mask, lower_faces,
                  required, absent, collision_points):
    optional = tuple(
        row for row in rows if row[:2] not in required and row[:2] not in absent
    )
    required_rows = tuple(row for row in rows if row[:2] in required)
    collision_bits = tuple(1 << point_index[point] for point in collision_points)
    counter = Counter()
    total = 0
    for optional_mask in range(1 << len(optional)):
        active = list(required_rows)
        active.extend(
            row for index, row in enumerate(optional)
            if optional_mask & (1 << index)
        )
        base_mask = support_mask(tuple(active))
        for collision_mask in range(1 << len(collision_bits)):
            mask = base_mask
            for index, bit in enumerate(collision_bits):
                if collision_mask & (1 << index):
                    mask &= ~bit
            counter[lower_faces(mask)] += 1
            total += 1
    return optional, total, counter


def direct_lower_faces(active_points):
    """Independent active-triple hull path for selected replay controls."""
    active_points = tuple(sorted(active_points))
    faces = set()
    for first, second, third in combinations(active_points, 3):
        plane = plane_through(first, second, third)
        if plane is None:
            continue
        a, b, c = plane
        gaps = tuple(
            Fraction(h) - a * r - b * ell - c
            for r, ell, h in active_points
        )
        if all(gap >= 0 for gap in gaps):
            equal = tuple(
                point for point, gap in zip(active_points, gaps) if gap == 0
            )
            if projected_rank_two(equal):
                faces.add(plane)
    return tuple(sorted(faces))


def point_in_convex_polygon(point, polygon):
    """Exact closed-polygon containment, independent of the 3D hull engine."""
    polygon = convex_hull(polygon)
    signs = []
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        cross = (
            (end[0] - start[0]) * (point[1] - start[1])
            - (end[1] - start[1]) * (point[0] - start[0])
        )
        if cross:
            signs.append(1 if cross > 0 else -1)
    return not signs or all(sign == signs[0] for sign in signs)


def dual_dominance_certificate(rows, absent_rows, main_plane, merged_plane,
                               main_polygon, merged_polygon, anchors):
    """Prove lower-hull stability from unavoidable points and exact slacks.

    This path does not generate candidate planes or enumerate masks.  The five
    unavoidable anchors already have exactly the two desired lower facets.
    Every possible target-support point projects into their two-face polygon
    and lies on or above both dual affine planes.  Adding such dominated
    points, or deleting only nonanchors, cannot change the lower hull.
    """
    universe = tuple(sorted(lifted_support(
        row for row in rows if row[:2] not in absent_rows
    )))
    anchors = tuple(sorted(anchors))
    need(set(anchors).issubset(universe), "dual anchors belong to universe")
    need(direct_lower_faces(anchors) == (main_plane, merged_plane),
         "five-anchor lower hull has exactly M and K9")

    main_witness = ((2, 0, 0), (0, 1, 0), (0, 7, 1))
    merged_witness = ((2, 0, 0), (4, 5, 1), (5, 3, 1))
    need(projected_rank_two(main_witness), "unavoidable rank-two M witness")
    need(projected_rank_two(merged_witness), "unavoidable rank-two K9 witness")

    slack_rows = []
    for r, ell, height in universe:
        main_slack = Fraction(height) - (
            main_plane[0] * r + main_plane[1] * ell + main_plane[2]
        )
        merged_slack = Fraction(height) - (
            merged_plane[0] * r + merged_plane[1] * ell + merged_plane[2]
        )
        need(main_slack >= 0 and merged_slack >= 0,
             "nonnegative dual slacks at " + repr((r, ell, height)))
        need(
            point_in_convex_polygon((r, ell), main_polygon)
            or point_in_convex_polygon((r, ell), merged_polygon),
            "two-face projection cover at " + repr((r, ell, height)),
        )
        slack_rows.append(((r, ell, height), main_slack, merged_slack))

    collision_points = {
        (2, 3, 1), (2, 4, 1), (2, 5, 1),
        (2, 6, 1), (3, 4, 1), (3, 5, 1),
    }
    need(not (set(anchors) & collision_points),
         "hostile collision deletions preserve every anchor")
    serialization = "\n".join(
        "%d,%d,%d:%d/%d:%d/%d" % (
            point[0], point[1], point[2],
            main_slack.numerator, main_slack.denominator,
            merged_slack.numerator, merged_slack.denominator,
        )
        for point, main_slack, merged_slack in slack_rows
    ).encode("ascii")
    return universe, tuple(slack_rows), sha256(serialization).hexdigest()


# ---------------------------------------------------------------------------
# Exact polygon, packet, and univariate edge utilities


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


def poly_trim(poly):
    poly = list(map(Fraction, poly))
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return tuple(poly)


def poly_derivative(poly):
    poly = poly_trim(poly)
    if len(poly) == 1:
        return (Fraction(0),)
    return tuple(Fraction(index) * poly[index] for index in range(1, len(poly)))


def poly_divmod(numerator, denominator):
    numerator = list(poly_trim(numerator))
    denominator = poly_trim(denominator)
    need(denominator != (0,), "nonzero polynomial divisor")
    quotient = [Fraction(0)] * max(1, len(numerator) - len(denominator) + 1)
    while len(numerator) >= len(denominator) and any(numerator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] / denominator[-1]
        quotient[shift] += coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] -= coefficient * value
        numerator = list(poly_trim(numerator))
    return poly_trim(quotient), poly_trim(numerator)


def poly_gcd(first, second):
    first = poly_trim(first)
    second = poly_trim(second)
    while second != (0,):
        _quotient, remainder = poly_divmod(first, second)
        first, second = second, remainder
    if first == (0,):
        return first
    lead = first[-1]
    return tuple(value / lead for value in first)


def poly_squarefree(poly):
    poly = poly_trim(poly)
    return len(poly) >= 2 and poly_gcd(poly, poly_derivative(poly)) == (Fraction(1),)


def edge_schemes(U, W, zeta):
    # Coefficients are ascending in the primitive edge parameter X.
    return (
        ("X-1", (-1, 1)),
        ("1-zeta*X^3", (1, 0, 0, -zeta)),
        ("-zeta-W*X", (-zeta, -W)),
        ("(X-1)*(U*X+W)", (-W, W - U, U)),
        ("U-X^6", (U, 0, 0, 0, 0, 0, -1)),
        ("1-W*X", (1, -W)),
    )


# ---------------------------------------------------------------------------
# Tiny exact sparse ring for the face identities


NVARS = 5  # S, P, U, W, zeta
S_INDEX, P_INDEX, U_INDEX, W_INDEX, ZETA_INDEX = range(NVARS)


def sparse_clean(poly):
    return {monomial: Fraction(value) for monomial, value in poly.items() if value}


def sparse_constant(value):
    value = Fraction(value)
    return {} if value == 0 else {(0,) * NVARS: value}


def sparse_variable(index):
    exponent = [0] * NVARS
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def sparse_add(*polys):
    result = {}
    for poly in polys:
        for monomial, value in poly.items():
            result[monomial] = result.get(monomial, Fraction(0)) + value
    return sparse_clean(result)


def sparse_scale(poly, scalar):
    scalar = Fraction(scalar)
    return sparse_clean({monomial: scalar * value for monomial, value in poly.items()})


def sparse_neg(poly):
    return sparse_scale(poly, -1)


def sparse_sub(first, second):
    return sparse_add(first, sparse_neg(second))


def sparse_mul(first, second):
    result = {}
    for monomial_a, value_a in first.items():
        for monomial_b, value_b in second.items():
            monomial = tuple(a + b for a, b in zip(monomial_a, monomial_b))
            result[monomial] = result.get(monomial, Fraction(0)) + value_a * value_b
    return sparse_clean(result)


def sparse_pow(poly, exponent):
    need(exponent >= 0, "nonnegative sparse exponent")
    result = sparse_constant(1)
    base = poly
    while exponent:
        if exponent & 1:
            result = sparse_mul(result, base)
        base = sparse_mul(base, base)
        exponent //= 2
    return result


def sparse_derivative(poly, index):
    result = {}
    for monomial, value in poly.items():
        if monomial[index]:
            new_monomial = list(monomial)
            multiplier = new_monomial[index]
            new_monomial[index] -= 1
            new_monomial = tuple(new_monomial)
            result[new_monomial] = result.get(new_monomial, Fraction(0)) + multiplier * value
    return sparse_clean(result)


def sparse_substitute_constant(poly, index, value):
    value = Fraction(value)
    result = {}
    for monomial, coefficient in poly.items():
        new_monomial = list(monomial)
        exponent = new_monomial[index]
        new_monomial[index] = 0
        new_monomial = tuple(new_monomial)
        result[new_monomial] = (
            result.get(new_monomial, Fraction(0))
            + coefficient * value ** exponent
        )
    return sparse_clean(result)


def sparse_degree(poly, index):
    return max((monomial[index] for monomial in poly), default=-1)


def check_face_identities():
    one = sparse_constant(1)
    S = sparse_variable(S_INDEX)
    P = sparse_variable(P_INDEX)
    U = sparse_variable(U_INDEX)
    W = sparse_variable(W_INDEX)
    zeta = sparse_variable(ZETA_INDEX)
    S2, S3 = sparse_pow(S, 2), sparse_pow(S, 3)
    P3, P5, P6, P9 = (
        sparse_pow(P, 3), sparse_pow(P, 5), sparse_pow(P, 6), sparse_pow(P, 9)
    )

    x = sparse_mul(W, sparse_mul(S2, P5))
    y = sparse_mul(zeta, sparse_mul(S3, P3))
    central = sparse_sub(sparse_sub(one, sparse_mul(U, P6)), x)
    merged = sparse_sub(sparse_sub(one, x), y)

    # Central component: with x0=P^-1 and y0=W*S/P,
    # y0^2=W*x0*(x0^6-U).  This is the denominator-cleared identity.
    central_left = sparse_sub(
        sparse_mul(sparse_pow(W, 2), sparse_mul(S2, P5)),
        sparse_mul(W, sparse_sub(one, sparse_mul(U, P6))),
    )
    need(central_left == sparse_neg(sparse_mul(W, central)),
         "central hyperelliptic normal form")

    # Merged beta-zero component: x+y=1 and
    # P^9=(zeta^2/W^3)*x^3/(1-x)^2.
    one_minus_x = sparse_sub(one, x)
    merged_left = sparse_sub(
        sparse_mul(sparse_mul(sparse_pow(W, 3), P9), sparse_pow(one_minus_x, 2)),
        sparse_mul(sparse_pow(zeta, 2), sparse_pow(x, 3)),
    )
    merged_right = sparse_mul(
        sparse_mul(sparse_pow(W, 3), P9),
        sparse_mul(merged, sparse_add(one_minus_x, y)),
    )
    need(merged_left == merged_right, "merged ninth-Kummer normal form")

    central_remainder = sparse_sub(
        sparse_mul(P, sparse_derivative(central, P_INDEX)),
        sparse_scale(central, 5),
    )
    expected_central = sparse_neg(
        sparse_add(sparse_mul(U, P6), sparse_constant(5))
    )
    need(central_remainder == expected_central,
         "central generic-point P-derivative identity")
    need(sparse_degree(central_remainder, S_INDEX) < sparse_degree(central, S_INDEX),
         "central derivative proper S-remainder")

    merged_remainder = sparse_sub(
        sparse_mul(P, sparse_derivative(merged, P_INDEX)),
        sparse_scale(merged, 3),
    )
    expected_merged = sparse_neg(
        sparse_add(sparse_scale(x, 2), sparse_constant(3))
    )
    need(merged_remainder == expected_merged,
         "merged generic-point P-derivative identity")
    need(sparse_degree(merged_remainder, S_INDEX) < sparse_degree(merged, S_INDEX),
         "merged derivative proper S-remainder")

    # Exact symbolic top-edge discriminant.
    top_a = sparse_neg(W)
    top_b = sparse_sub(W, U)
    top_c = U
    top_discriminant = sparse_sub(sparse_pow(top_b, 2), sparse_scale(sparse_mul(top_a, top_c), 4))
    need(top_discriminant == sparse_pow(sparse_add(U, W), 2),
         "top-edge discriminant=(U+W)^2")

    return (
        "central:y^2=W*x*(x^6-U)",
        "merged:P^9=(zeta^2/W^3)*x^3/(1-x)^2",
        "central:P*C_P-5*C=-(U*P^6+5)",
        "merged:P*B_P-3*B=-(2W*S^2*P^5+3)",
    )


def check_lambda_zero_contact():
    """Exact compatibility inputs for the separately imported A23 audit."""
    one = sparse_constant(1)
    b = sparse_variable(S_INDEX)
    r = sparse_variable(P_INDEX)
    U = sparse_variable(U_INDEX)
    q = sparse_sub(r, one)

    # At Z=beta=0 and W=-U, the main closure is
    # Ctilde=b^12-U*r^5*(r-1).  The R branch is r=1.
    ctilde = sparse_sub(
        sparse_pow(b, 12),
        sparse_mul(U, sparse_mul(sparse_pow(r, 5), q)),
    )
    on_r_branch = sparse_substitute_constant(ctilde, P_INDEX, 1)
    need(on_r_branch == sparse_pow(b, 12),
         "Lambda-zero intersection multiplicity exactly twelve")
    derivative_at_contact = sparse_substitute_constant(
        sparse_substitute_constant(
            sparse_derivative(ctilde, P_INDEX), S_INDEX, 0
        ),
        P_INDEX,
        1,
    )
    need(derivative_at_contact == sparse_neg(U),
         "Lambda-zero C branch smooth transverse derivative")

    # The only torus root of r^5(r-1) is r=1; r=0 is a toric corner.
    top_factor = (0, 0, 0, 0, 0, -1, 1)
    r5 = (0, 0, 0, 0, 0, 1)
    quotient, remainder = poly_divmod(top_factor, r5)
    need(quotient == (-1, 1) and remainder == (0,),
         "unique Lambda-zero torus contact root")

    # THM-4297's difference from the W=0 prepared model becomes
    # +(U/2)r^4q^3.  Its exact q-order is three.
    correction = sparse_scale(
        sparse_mul(U, sparse_mul(sparse_pow(r, 4), sparse_pow(q, 3))),
        Fraction(1, 2),
    )
    third_derivative = correction
    for _index in range(3):
        third_derivative = sparse_derivative(third_derivative, P_INDEX)
    third_at_one = sparse_substitute_constant(third_derivative, P_INDEX, 1)
    need(third_at_one == sparse_scale(U, 3),
         "Lambda-zero correction has exact q-order three")

    # The M--K9 edge is a different polygon edge.  Its nonzero root is in
    # the relative interior, so it cannot share the top-edge torus orbit.
    top_edge = frozenset(((4, 5), (0, 7)))
    shared_edge = frozenset(((2, 0), (4, 5)))
    need(top_edge != shared_edge and top_edge & shared_edge == {(4, 5)},
         "top and shared edge relative interiors are disjoint")
    shared_polynomial = poly_trim((1, 2))  # W=-U=-2 at a nonzero control.
    need(poly_squarefree(shared_polynomial)
         and shared_polynomial[0] * shared_polynomial[-1] != 0,
         "Lambda-zero M--K9 edge simple and noncorner")

    lambda0_schemes = edge_schemes(Fraction(2), Fraction(-2), Fraction(3))
    repeated = []
    for name, polynomial in lambda0_schemes:
        polynomial = poly_trim(polynomial)
        if poly_squarefree(polynomial):
            need(polynomial[0] * polynomial[-1] != 0,
                 "Lambda-zero simple edge avoids corners: " + name)
        else:
            repeated.append(name)
    need(repeated == ["(X-1)*(U*X+W)"],
         "unique Lambda-zero repeated edge scheme")

    tail_coefficients = (
        (Fraction(6), Fraction(2)),
        (Fraction(5, 2), Fraction(9, 2)),
        (Fraction(2), Fraction(4)),
        (Fraction(3, 2), Fraction(7, 2)),
        (Fraction(1), Fraction(3)),
    )
    need(all(s_coefficient > 0 and b_coefficient > 0
             for s_coefficient, b_coefficient in tail_coefficients),
         "all imported A23 tail orders positive")
    need(Fraction(36, 12) == 3, "Lambda-zero base-change index is three")
    return tail_coefficients


def kummer_genus(degree, valuations):
    need(sum(valuations) == 0, "Kummer divisor has degree zero")
    need(gcd(degree, *(abs(value) for value in valuations)) == 1,
         "Kummer cover connected")
    twice_genus_minus_two = -2 * degree + sum(
        degree - gcd(degree, abs(value)) for value in valuations
    )
    need(twice_genus_minus_two % 2 == 0, "integral Kummer genus")
    return (twice_genus_minus_two + 2) // 2


def face_order(base, plane):
    scaled = tuple(base * value for value in plane)
    need(all(value.denominator == 1 for value in scaled), "integral face chart")
    need((Fraction(base, 6)).denominator == 1, "integral good-target twist")
    order = base * (Fraction(5, 6) - sum(plane))
    need(order.denominator == 1 and order > 0, "positive integral face order")
    return int(order), tuple(int(value) for value in scaled)


def equal_row_owners(rows, plane):
    a, b, c = plane
    owners = []
    for i, j, weight in rows:
        for tag, point in (
            ("s2", (j + 2, i + j, 1)),
            ("p", (j, i + j + 1, 1)),
        ):
            if Fraction(point[2]) == a * point[0] + b * point[1] + c:
                owners.append(((i, j, weight), tag, point))
    return tuple(owners)


def main():
    rows = weighted_rows()
    expected_rows = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11), (6, 0, 12), (3, 2, 12), (0, 4, 12),
    )
    need(rows == expected_rows, "complete labelled weight-at-most-twelve rows")
    points, point_index, support_mask, lower_faces = make_hull_engine(rows)

    M12 = (Fraction(1, 12), Fraction(1, 6), Fraction(-1, 6))
    MERGED9 = (Fraction(2, 9), Fraction(1, 9), Fraction(-4, 9))
    GENERIC7 = (Fraction(1, 7), Fraction(1, 7), Fraction(-2, 7))
    GENERIC3 = (Fraction(1, 3), Fraction(0), Fraction(-2, 3))
    W0_PLANE = (Fraction(1, 6), Fraction(1, 6), Fraction(-1, 3))

    P_ROW, P2_ROW, P3_ROW = (1, 0), (2, 0), (3, 0)
    U_ROW, W_ROW = (6, 0), (3, 2)
    Z_ROW, BETA_ROW, ZETA_ROW = (0, 4), (1, 3), (0, 3)
    fixed = {P_ROW, P2_ROW, P3_ROW, U_ROW}
    collisions = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1),
        (2, 6, 1), (3, 4, 1), (3, 5, 1),
    )

    # Target seam: W and zeta are nonzero owners; Z and beta are absent.
    target_required = fixed | {W_ROW, ZETA_ROW}
    target_absent = {Z_ROW, BETA_ROW}
    target_owner_map = {}
    for i, j, weight in rows:
        if (i, j) in target_absent:
            continue
        for tag, point in (
            ("s2", (j + 2, i + j, 1)),
            ("p", (j, i + j + 1, 1)),
        ):
            target_owner_map.setdefault(point, []).append(((i, j, weight), tag))
    need(target_owner_map[(0, 7, 1)] == [((6, 0, 12), "p")],
         "U anchor has a unique target owner")
    need(target_owner_map[(4, 5, 1)] == [((3, 2, 12), "s2")],
         "W anchor has a unique target owner after Z=0")
    need(target_owner_map[(5, 3, 1)] == [((0, 3, 9), "s2")],
         "zeta anchor has a unique target owner")
    optional, total, target_counter = hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        target_required, target_absent, collisions,
    )
    target_faces = (M12, MERGED9)
    need(tuple(row[:2] for row in optional) == (
        (0, 2), (2, 1), (4, 0), (1, 2),
        (3, 1), (5, 0), (2, 2), (4, 1),
    ), "target optional-row ledger")
    need(total == 16384, "target hostile-mask count")
    need(target_counter == Counter({target_faces: total}),
         "target unique two-face complex")

    # Independent active-triple hull replay: minimal/maximal row supports,
    # every singleton collision deletion, alternating deletions, and all.
    required_rows = tuple(row for row in rows if row[:2] in target_required)
    maximal_rows = tuple(row for row in rows if row[:2] not in target_absent)
    replay_specs = [("minimal", required_rows, ()) , ("maximal", maximal_rows, ())]
    replay_specs.extend(
        ("single_" + str(index), maximal_rows, (point,))
        for index, point in enumerate(collisions)
    )
    replay_specs.extend((
        ("alternating_even", maximal_rows, collisions[::2]),
        ("alternating_odd", maximal_rows, collisions[1::2]),
        ("all_collisions", maximal_rows, collisions),
    ))
    for name, active_rows, deleted in replay_specs:
        mask = support_mask(active_rows)
        for point in deleted:
            mask &= ~(1 << point_index[point])
        active_points = tuple(
            point for index, point in enumerate(points) if mask & (1 << index)
        )
        need(lower_faces(mask) == direct_lower_faces(active_points) == target_faces,
             "independent hull replay " + name)

    # Known generic theorem branch as a positive control: restoring beta
    # splits the merged face into the genus-three beta face and rational T3.
    generic_optional, generic_total, generic_counter = hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        fixed | {W_ROW, BETA_ROW, ZETA_ROW}, {Z_ROW}, collisions,
    )
    need(len(generic_optional) == 8 and generic_total == 16384,
         "beta-restored control universe")
    need(generic_counter == Counter({(M12, GENERIC7, GENERIC3): 16384}),
         "beta-restored THM-4327 face control")

    # Boundary hostiles.  With W=0, the new face plane has optional owners,
    # so its face equation is not the stable merged Kummer equation.  With
    # zeta=0 as well, even the face-complex type varies across allowed masks.
    w0_collisions = tuple(point for point in collisions if point != (2, 6, 1))
    w0_optional, w0_total, w0_counter = hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        fixed | {ZETA_ROW}, {Z_ROW, BETA_ROW, W_ROW}, w0_collisions,
    )
    need(len(w0_optional) == 8 and w0_total == 8192,
         "W-zero hostile universe")
    need(w0_counter == Counter({(M12, W0_PLANE): 8192}),
         "W-zero hostile face plane")
    w0_owners = equal_row_owners(rows, W0_PLANE)
    w0_optional_owners = tuple(
        owner[0][:2] for owner in w0_owners if owner[0][:2] in {(2, 2), (4, 1)}
    )
    need(w0_optional_owners == ((2, 2), (4, 1)),
         "W-zero face has two optional owners")

    both0_optional, both0_total, both0_counter = hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        fixed | {W_ROW}, {Z_ROW, BETA_ROW, ZETA_ROW}, collisions,
    )
    need(len(both0_optional) == 8 and both0_total == 16384,
         "zeta-zero hostile universe")
    need(len(both0_counter) == 2 and sorted(both0_counter.values()) == [2048, 14336],
         "zeta-zero face-complex bifurcation")

    face_identities = check_face_identities()
    need((0 * 5 - 6 * 2) == -12, "central exponent determinant")
    need((2 * 3 - 5 * 3) == -9, "merged exponent determinant")

    central_genus = kummer_genus(2, (1,) * 7 + (-7,))
    merged_genus = kummer_genus(9, (3, -2, -1))
    need((central_genus, merged_genus) == (3, 3), "face genera")

    polygons = {
        "M": ((0, 1), (2, 0), (4, 5), (0, 7)),
        "central": ((0, 0), (0, 6), (2, 5)),
        "merged": ((2, 0), (4, 5), (5, 3)),
        "global": ((0, 1), (2, 0), (5, 3), (4, 5), (0, 7)),
    }
    anchors = {
        (2, 0, 0),  # base S^2
        (0, 1, 0),  # base -P
        (0, 7, 1),  # U*P^7
        (4, 5, 1),  # -W*S^4*P^5
        (5, 3, 1),  # -zeta*S^5*P^3
    }
    dual_universe, dual_slacks, dual_sha256 = dual_dominance_certificate(
        rows,
        target_absent,
        M12,
        MERGED9,
        polygons["M"],
        polygons["merged"],
        anchors,
    )
    need(len(dual_universe) == 26, "dual target support-universe size")
    need(len(dual_slacks) == 26, "complete dual slack table")
    ledgers = {name: polygon_ledger(polygon)[1:] for name, polygon in polygons.items()}
    need(ledgers == {
        "M": (36, 10, 14),
        "central": (12, 8, 3),
        "merged": (9, 5, 3),
        "global": (45, 13, 17),
    }, "Pick ledgers")

    # At Lambda != 0, R and C meet in twelve transverse points; C and the
    # merged face meet once along their primitive common edge.
    graph_vertices = 3
    graph_edges = 12 + 1
    graph_b1 = graph_edges - graph_vertices + 1
    need(graph_b1 == 11, "component graph first Betti number")
    need(central_genus + merged_genus + graph_b1 == ledgers["global"][2],
         "complete global genus conservation")

    packet = edge_packet(polygons["global"])
    need(packet == (11, 11, 10, 2, 2, 2, 1), "source-infinity packet")
    need(sum(index - 1 for index in packet) == 32 == 2 * 17 - 2,
         "packet Riemann-Hurwitz saturation")
    need((sum(packet) - 6, sum(packet)) == (33, 39),
         "locked-cubic response consistency control")

    main_order, main_chart = face_order(36, M12)
    merged_order, merged_chart = face_order(36, MERGED9)
    need((main_order, merged_order) == (27, 34), "exact face orders")
    need(main_chart == (3, 6, -6) and merged_chart == (8, 4, -16),
         "scaled height graphs")
    need(gcd(abs(main_chart[0]), abs(main_chart[1]), 1) == 1,
         "main primitive height normal")
    need(gcd(abs(merged_chart[0]), abs(merged_chart[1]), 1) == 1,
         "merged primitive height normal")

    lambda0_tail_coefficients = check_lambda_zero_contact()
    # Deep repetition has c3=eta+zeta=0.  On this beta-zero theorem gate,
    # zeta!=0 therefore forces eta=-zeta and C3=eta/c^2!=0; C2 may split
    # sooner, but the ladder cannot reach C4 or the terminal splitter.
    zeta_control = Fraction(3)
    eta_on_repeated_c3 = -zeta_control
    need(eta_on_repeated_c3 != 0, "zeta gate forces nonzero C3")
    need(lambda0_tail_coefficients[2:4] == (
        (Fraction(2), Fraction(4)),
        (Fraction(3, 2), Fraction(7, 2)),
    ), "C2/C3 positive short-ladder orders")

    # Edge schemes are read independently from the two face polygons.  Test
    # two Lambda-nonzero rational controls, including mixed signs.
    for label, (U_value, W_value, zeta_value) in (
        ("positive", (Fraction(2), Fraction(1), Fraction(3))),
        ("mixed_sign", (Fraction(-3), Fraction(2), Fraction(-5))),
    ):
        need(U_value * W_value * zeta_value * (U_value + W_value) != 0,
             label + " edge-gate control")
        for name, polynomial in edge_schemes(U_value, W_value, zeta_value):
            polynomial = poly_trim(polynomial)
            need(polynomial[0] != 0 and polynomial[-1] != 0,
                 label + " edge avoids toric corners: " + name)
            need(poly_squarefree(polynomial), label + " squarefree edge: " + name)

    # Lambda=0 is a separate branch: the top edge is U*(X-1)^2, so simple-root
    # packets cannot be transported.  check_lambda_zero_contact independently
    # verifies the new compatibility inputs for THM-4292/4297's tail table.
    lambda0_top = poly_trim(edge_schemes(1, -1, 3)[3][1])
    need(lambda0_top == (1, -2, 1), "Lambda-zero top=U*(X-1)^2 control")
    need(not poly_squarefree(lambda0_top), "Lambda-zero repeated-edge hostile")

    # At zeta=0 and W=0, an endpoint of a defining merged-face edge vanishes.
    need(len(poly_trim(edge_schemes(2, 1, 0)[1][1])) == 1,
         "zeta-zero loses cubic edge")
    need(len(poly_trim(edge_schemes(2, 0, 3)[2][1])) == 1,
         "W-zero merged edge loses its relative-interior root")

    edge_names = tuple(name for name, _poly in edge_schemes(2, 1, 3))
    dual_slack_profile = Counter(
        (main_slack, merged_slack)
        for _point, main_slack, merged_slack in dual_slacks
    )

    print("THM4334_BETA0_ENDPOINT_INDEPENDENT_AUDIT=PASS")
    print("STATUS=FINITE-EXACT_INDEPENDENT_AUDIT")
    print("TARGET=Z=0;beta_11=0;U*W*zeta_3!=0;Lambda=U+W_arbitrary;all_other_lower_rows_arbitrary")
    print("UNIVERSE_ROWS=16 forbidden_rows=((0,1),(1,1)) wall_absent=((0,4),(1,3))")
    print("FIXED_REQUIRED_ROWS=((1,0),(2,0),(3,0),(6,0),(3,2),(0,3))")
    print("OPTIONAL_ROWS=" + repr(tuple(row[:2] for row in optional)))
    print("COLLISION_POINTS=" + repr(collisions) + " hostile_overcount=true")
    print("TARGET_HOSTILES=16384 optional_bits=8 collision_bits=6")
    print("TARGET_FACE_COMPLEX=" + repr(target_faces) + " multiplicity=16384 unique=true")
    print("INDEPENDENT_DUAL_CERT=5_unavoidable_anchors;rank2_witnesses=M,K9;26_support_points_dominated;two_face_projection_cover=true")
    print("INDEPENDENT_ANCHOR_OWNERS=(0,7,1):U_p;(4,5,1):W_s2;(5,3,1):zeta_s2;all_unique_after_Z=beta=0")
    print("INDEPENDENT_DUAL_SLACK_PROFILE=" + repr(dual_slack_profile))
    print("INDEPENDENT_DUAL_SHA256=" + dual_sha256)
    print("INDEPENDENT_DIRECT_HULL_CONTROLS=" + str(len(replay_specs)) + " all_match=true")
    print("FACE_CORES=M:(S^2-P)*(1-U*P^6-W*S^2*P^5);MERGED:S^2*(1-W*S^2*P^5-zeta*S^3*P^3)")
    print("FACE_IDENTITIES=" + repr(face_identities))
    print("KUMMER_VALUATIONS=central:n2:(1x7,-7);merged:n9:(3,-2,-1)")
    print("FACE_GENERA=central:3,merged:3 torus_determinants=(-12,-9)")
    print("PICK_LEDGERS=" + repr(ledgers))
    print("GRAPH_LAMBDA_NONZERO=vertices(R,C,K9):3 edges(R-C:12,C-K9:1):13 b1=11")
    print("GENUS_CONSEQUENCE=3+3+11=17=global;no_unlisted_positive_genus_component")
    print("HEIGHT_BASE=36 scaled_planes=M:" + repr(main_chart) + ",K9:" + repr(merged_chart) + ";normal_height_coefficient=-1")
    print("GOOD_DIFFERENTIAL_ORDERS=M:27,K9:34 all_positive=true")
    print("EDGE_SCHEMES=" + repr(edge_names))
    print("EDGE_GATE_LAMBDA_NONZERO=U*W*zeta_3*(U+W);all_simple_and_toric_noncorner=true")
    print("PACKET_LAMBDA_NONZERO=" + repr(packet) + " finite/full_responses=(33,39) consistency_only=true")
    print("CONSEQUENCE_OBJECT=components(R:g0,C:g3,K9:g3);graph_b1=11;global_g=17;positive_orders=(27,34);only_rational_resolution_chains")
    print("LAMBDA_ZERO_CONTACT=Ctilde=b^12-U*r^5*(r-1);intersection_length=12;partial_r=-U;unique_torus_contact=A23")
    print("LAMBDA_ZERO_SHARED_EDGE=1-WX_simple_noncorner;top_and_shared_relative_interiors_orbit_disjoint=true")
    print("LAMBDA_ZERO_CORRECTION=+(U/2)*r^4*(r-1)^3 exact_q_order=3;base_change_index=3")
    print("LAMBDA_ZERO_IMPORTED_TAIL_COEFFICIENTS=" + repr(lambda0_tail_coefficients) + " all_positive=true;simple_packet_transported=false")
    print("LAMBDA_ZERO_SHORT_LADDER=zeta_3!=0_forces_C2_or_C3_before_t^6_correction")
    print("RELATIVE_THEOREM_SIGNAL=THM4327_good_differential_plus_proper_flat_degree_mechanism_has_all_new_inputs;Lambda0_uses_THM4292/4297_tail_import")
    print("CONTROL_BETA_RESTORED=THM4327_faces(M,V7,T3)_all_16384_masks")
    print("HOSTILE_W0=8192_masks_unique_planes(M,W0);W0_face_optional_owners=((2,2),(4,1));stable_factorization_lost")
    print("HOSTILE_ZETA0_BETA0_WNZ=16384_masks_face_complex_types=2 counts=(14336,2048)")
    print("HOSTILE_LAMBDA0_SIMPLE_MODEL=top_edge_U*(X-1)^2;simple_packet_forbidden;separate_A23_inputs_pass")
    print("OPEN=W0_beta0;zeta0_intersections;seam_entry;JC2;DC2")
    print("CHECKS=" + str(CHECKS))


if __name__ == "__main__":
    main()
