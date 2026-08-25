#!/usr/bin/env python3
"""Independent exact audit for THM-4103's theta-boundary response sieve.

This audit deliberately does not use SymPy or the primary certificate's hull
implementation.  It reconstructs the theta-only expansion term by term,
certifies the proposed polygons by their primitive supporting inequalities,
uses shoelace plus boundary gcds for Pick's theorem, and enumerates all 2^7
labelled response subsets literally.  It checks both live optional atoms
phi and, as robustness hostiles, the counterfactual lambda=0 or alpha=0 patterns;
none changes the Newton polygon or any boundary datum.

No ``assert`` statements are used, so every gate survives ``python -O``.
The mathematical input imported from THM-4053 is that the six polygon vertex
coefficients are nonzero, the outer edge polynomials are squarefree on the
nonresonant theta-only stratum, the generic source has the toric smooth model,
and its degree is an Eisenstein norm.
"""

from hashlib import sha256
from itertools import product
from json import dumps
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add_coefficient(table, exponent, atom, scalar):
    coefficient = table.setdefault(exponent, {})
    coefficient[atom] = coefficient.get(atom, 0) + scalar
    if coefficient[atom] == 0:
        del coefficient[atom]


def theta_support(optional_nonzero):
    """Expand (s^2-p)(1-QH)+gamma*Q*s^2 over a formal coefficient basis."""
    table = {}
    add_coefficient(table, (2, 0), "one", 1)
    add_coefficient(table, (0, 1), "one", -1)
    add_coefficient(table, (2, 0), "gammaQ", 1)

    # H monomials after y=sp.  epsilon, kappa, and theta are inherited
    # nonzero.  lambda and alpha are allowed to vanish here as an extra
    # robustness check; phi genuinely may vanish under Delta_V != 0.
    h_terms = (
        ((0, 1), "lambdaQ", optional_nonzero[0]),
        ((0, 2), "alphaQ", optional_nonzero[1]),
        ((0, 3), "epsilonQ", True),
        ((2, 2), "kappaQ", True),
        ((1, 3), "phiQ", optional_nonzero[2]),
        ((2, 3), "thetaQ", True),
    )
    for (s_power, p_power), atom, active in h_terms:
        if not active:
            continue
        # -Q*s^2*H + Q*p*H
        add_coefficient(table, (s_power + 2, p_power), atom, -1)
        add_coefficient(table, (s_power, p_power + 1), atom, 1)
    return {exponent: coefficient for exponent, coefficient in table.items()
            if coefficient}


def polygon_certificate(support, vertices, expected_edges):
    """Certify a convex polygon directly from primitive supporting lines."""
    require(all(vertex in support for vertex in vertices),
            "a claimed polygon vertex is absent")
    rows = []
    area_twice = 0
    boundary = 0
    for index, (left, right) in enumerate(zip(vertices, vertices[1:] + vertices[:1])):
        dx = right[0] - left[0]
        dy = right[1] - left[1]
        length = gcd(abs(dx), abs(dy))
        require(length > 0, "repeated polygon vertex")
        normal = (-dy // length, dx // length)
        height = normal[0] * left[0] + normal[1] * left[1]
        require((length, normal, height) == expected_edges[index],
                "edge primitive data changed")
        require(normal[0] * right[0] + normal[1] * right[1] == height,
                "right endpoint left its supporting line")
        require(all(normal[0] * point[0] + normal[1] * point[1] >= height
                    for point in support),
                "support point violates a polygon half-space")
        equality_points = tuple(sorted(
            point for point in support
            if normal[0] * point[0] + normal[1] * point[1] == height))
        require(left in equality_points and right in equality_points,
                "supporting edge lost an endpoint")
        rows.append((left, right, length, normal, height, equality_points))
        area_twice += left[0] * right[1] - right[0] * left[1]
        boundary += length
    require(area_twice > 0, "polygon is not counterclockwise")
    require((area_twice - boundary) % 2 == 0, "Pick numerator is odd")
    genus = (area_twice - boundary + 2) // 2
    return rows, area_twice, boundary, genus


def boundary_response(rows, residue_point, affine_normal):
    infinity = []
    affine = []
    edge_data = []
    for left, right, length, normal, height, equality_points in rows:
        distance = (normal[0] * residue_point[0]
                    + normal[1] * residue_point[1] - height)
        require(distance >= 1, "residue point is not strictly interior")
        order = distance - 1
        packet = [distance] * length
        if normal == affine_normal:
            affine.extend(packet)
        else:
            infinity.extend(packet)
        edge_data.append((left, right, length, normal, height, distance,
                          order, length * order, equality_points))
    return tuple(edge_data), tuple(infinity), tuple(affine)


def prime_factorization(value):
    factors = []
    remaining = value
    divisor = 2
    while divisor * divisor <= remaining:
        valuation = 0
        while remaining % divisor == 0:
            remaining //= divisor
            valuation += 1
        if valuation:
            factors.append((divisor, valuation))
        divisor += 1
    if remaining > 1:
        factors.append((remaining, 1))
    return tuple(factors)


def local_eisenstein_norm(value):
    return all(prime % 3 != 2 or valuation % 2 == 0
               for prime, valuation in prime_factorization(value))


def multiply_generating_factor(coefficients, weight):
    result = list(coefficients) + [0] * weight
    for degree, coefficient in enumerate(coefficients):
        result[degree + weight] += coefficient
    return result


def main():
    vertices = ((0, 1), (2, 0), (4, 2), (4, 3), (2, 4), (0, 4))
    expected_edges = (
        (1, (1, 2), 2),
        (2, (-1, 1), -2),
        (1, (-1, 0), -4),
        (1, (-1, -2), -10),
        (2, (0, -1), -4),
        (3, (1, 0), 0),
    )

    branch_ledgers = []
    reference_rows = None
    # Audit all eight presence patterns.  Only the last bit (phi) varies on
    # the live seam; lambda/alpha deletions are counterfactual hostiles.
    for optional_nonzero in product((False, True), repeat=3):
        coefficient_table = theta_support(optional_nonzero)
        support = frozenset(coefficient_table)
        rows, area2, boundary, genus = polygon_certificate(
            support, vertices, expected_edges)
        require((area2, boundary, genus) == (24, 10, 8),
                "theta polygon Pick ledger changed")
        # The shared (2,3) coefficient is kappa-epsilon and is inherited
        # nonzero.  Its formal two-atom signature must survive reconstruction.
        require(coefficient_table[(2, 3)]
                == {"epsilonQ": -1, "kappaQ": 1},
                "shared kappa-epsilon coefficient changed")
        if optional_nonzero[2]:
            require((3, 3) in support and (1, 4) in support,
                    "phi support points missing")
        else:
            require((3, 3) not in support and (1, 4) not in support,
                    "phi=0 support points survived")
        row_core = tuple((row[0], row[1], row[2], row[3], row[4])
                         for row in rows)
        if reference_rows is None:
            reference_rows = row_core
        require(row_core == reference_rows,
                "optional support changed polygon edge data")
        branch_ledgers.append((optional_nonzero, tuple(sorted(support))))

    # Use the maximal support only to print edge equality-point ledgers; the
    # preceding loop proves that optional deletions do not alter the polygon.
    maximal_support = frozenset(theta_support((True, True, True)))
    rows, area2, boundary, genus = polygon_certificate(
        maximal_support, vertices, expected_edges)
    edge_data, infinity_packet, affine_packet = boundary_response(
        rows, (1, 1), (1, 0))
    require(tuple(row[5] for row in edge_data) == (1, 2, 3, 7, 3, 1),
            "toric residue distances changed")
    require(infinity_packet == (1, 2, 2, 3, 7, 3, 3),
            "ordered source-infinity packet changed")
    require(tuple(sorted(infinity_packet)) == (1, 2, 2, 3, 3, 3, 7),
            "source-infinity multiset changed")
    require(affine_packet == (1, 1, 1),
            "three affine s=0 branches changed")
    ramification = sum(index - 1 for index in infinity_packet)
    pole_weight = sum(infinity_packet)
    require((ramification, pole_weight) == (14, 21),
            "boundary ramification/pole-weight ledger changed")
    require(ramification == 2 * genus - 2,
            "boundary ledger no longer saturates Riemann-Hurwitz")

    # Independent arithmetic routes to the Eisenstein degree intersection.
    represented = set()
    witnesses = {}
    for a in range(-24, 25):
        for b in range(-24, 25):
            norm = a * a - a * b + b * b
            if 7 <= norm <= pole_weight:
                represented.add(norm)
                witnesses.setdefault(norm, (a, b))
    local = {value for value in range(7, pole_weight + 1)
             if local_eisenstein_norm(value)}
    require(represented == local,
            "representation and local Eisenstein criteria disagree")
    allowed_degrees = tuple(sorted(represented))
    require(allowed_degrees == (7, 9, 12, 13, 16, 19, 21),
            "Eisenstein degree intersection changed")

    # Give every puncture an edge label.  This makes the 23 edge-refined
    # profiles independent of the coarser grouping by equal index.
    punctures = (
        ("AB", 0, 1),
        ("BC", 0, 2), ("BC", 1, 2),
        ("CD", 0, 3),
        ("DE", 0, 7),
        ("EF", 0, 3), ("EF", 1, 3),
    )
    require(tuple(item[2] for item in punctures) == infinity_packet,
            "labelled punctures do not match the boundary packet")
    labelled_by_degree = {degree: 0 for degree in allowed_degrees}
    coarse_profiles = set()
    edge_profiles = set()
    for mask in range(1 << len(punctures)):
        finite = tuple(punctures[index] for index in range(len(punctures))
                       if mask & (1 << index))
        debt = sum(item[2] for item in finite)
        degree = pole_weight - debt
        if degree not in labelled_by_degree:
            continue
        labelled_by_degree[degree] += 1
        counts_by_index = tuple(sum(item[2] == index for item in finite)
                                for index in (1, 2, 3, 7))
        coarse_profiles.add((degree,) + counts_by_index)
        counts_by_edge = tuple(sum(item[0] == edge for item in finite)
                               for edge in ("AB", "BC", "CD", "EF", "DE"))
        edge_profiles.add((degree,) + counts_by_edge)
    expected_labelled = {7: 7, 9: 9, 12: 9, 13: 10,
                         16: 7, 19: 2, 21: 1}
    require(labelled_by_degree == expected_labelled,
            "labelled response census changed")
    require((len(coarse_profiles), len(edge_profiles),
             sum(labelled_by_degree.values())) == (16, 23, 45),
            "16/23/45 response census changed")

    # AB is forced into the O-fibre, so its finite-image bit is zero.  The two
    # BC punctures form one quadratic closed point over k(q), so their bits
    # agree.  Re-enumerate masks literally instead of filtering grouped rows.
    sharp_labelled = {degree: 0 for degree in allowed_degrees}
    sharp_coarse = set()
    sharp_edges = set()
    for mask in range(1 << len(punctures)):
        if mask & 1:
            continue
        if bool(mask & (1 << 1)) != bool(mask & (1 << 2)):
            continue
        finite = tuple(punctures[index] for index in range(len(punctures))
                       if mask & (1 << index))
        degree = pole_weight - sum(item[2] for item in finite)
        if degree not in sharp_labelled:
            continue
        sharp_labelled[degree] += 1
        by_index = tuple(sum(item[2] == index for item in finite)
                         for index in (1, 2, 3, 7))
        sharp_coarse.add((degree,) + by_index)
        by_edge = tuple(sum(item[0] == edge for item in finite)
                        for edge in ("AB", "BC", "CD", "EF", "DE"))
        sharp_edges.add((degree,) + by_edge)
    sharp_labelled = {degree: count for degree, count in sharp_labelled.items()
                      if count}
    require(sharp_labelled == {7: 3, 12: 1, 21: 1},
            f"sharpened labelled responses changed: {sharp_labelled}")
    require(sharp_coarse == {
        (7, 0, 2, 1, 1),
        (12, 0, 0, 3, 0),
        (21, 0, 0, 0, 0),
    }, f"sharpened coarse profiles changed: {sharp_coarse}")
    require((len(sharp_coarse), len(sharp_edges), sum(sharp_labelled.values()))
            == (3, 4, 5), "sharpened 3/4/5 census changed")

    # DE support floors from the Laurent caps.  Apart from x^2 in A and x^3
    # in C, a t^n coefficient has x-degree at most n+1 or n+2.  Audit a
    # deliberately redundant box and the two sharp boundary pairs.
    a_competitors = []
    c_competitors = []
    for n in range(1, 120):
        for m in range(n + 2):
            if 7 * m - 6 * n >= 14:
                a_competitors.append((n, m, n + m))
        for m in range(n + 3):
            if 7 * m - 6 * n >= 21:
                c_competitors.append((n, m, n + m))
    require(min(n for n, _, _ in a_competitors) == 7
            and min(total for _, _, total in a_competitors) == 15
            and (7, 8, 15) in a_competitors,
            "independent DE A support floor changed")
    require(min(n for n, _, _ in c_competitors) == 7
            and min(total for _, _, total in c_competitors) == 16
            and (7, 9, 16) in c_competitors,
            "independent DE C support floor changed")

    # Independent generating-function checksum for the literal 2^7 census.
    generating = [1]
    for index in infinity_packet:
        generating = multiply_generating_factor(generating, index)
    require(sum(generating) == 1 << len(infinity_packet),
            "response generating polynomial does not count 2^7 subsets")
    require({degree: generating[pole_weight - degree]
             for degree in allowed_degrees} == expected_labelled,
            "generating-function and literal subset counts disagree")

    # Counterfactual observer hostile: delete the inherited kappa vertex.
    # This is not another live stratum; it tests the loss from retaining only
    # a highest face.  Direct supporting lines replace the old hull routine.
    hostile_support = maximal_support - {(4, 2)}
    hostile_vertices = ((0, 1), (2, 0), (4, 3), (2, 4), (0, 4))
    hostile_edges = (
        (1, (1, 2), 2),
        (1, (-3, 2), -6),
        (1, (-1, -2), -10),
        (2, (0, -1), -4),
        (3, (1, 0), 0),
    )
    hostile_rows, hostile_area2, hostile_boundary, hostile_genus = (
        polygon_certificate(hostile_support, hostile_vertices, hostile_edges))
    hostile_edge_data, hostile_infinity, hostile_affine = boundary_response(
        hostile_rows, (1, 1), (1, 0))
    require((hostile_area2, hostile_boundary, hostile_genus)
            == (22, 8, 8), "support-deletion Pick ledger changed")
    require(hostile_infinity == (1, 5, 7, 3, 3),
            "support-deletion infinity packet changed")
    require(hostile_affine == (1, 1, 1),
            "support-deletion affine packet changed")
    require(sum(index - 1 for index in hostile_infinity) == 14,
            "support-deletion RH checksum changed")
    require(sum(hostile_infinity) == 19,
            "support-deletion pole weight changed")

    semantic_object = {
        "vertices": vertices,
        "branch_ledgers": branch_ledgers,
        "edge_data": edge_data,
        "pick": (area2, boundary, genus),
        "infinity_packet": infinity_packet,
        "affine_packet": affine_packet,
        "ramification": ramification,
        "pole_weight": pole_weight,
        "allowed_degrees": allowed_degrees,
        "norm_witnesses": tuple((degree, witnesses[degree])
                                for degree in allowed_degrees),
        "labelled_by_degree": labelled_by_degree,
        "coarse_profiles": tuple(sorted(coarse_profiles)),
        "edge_profiles": tuple(sorted(edge_profiles)),
        "sharp_labelled": sharp_labelled,
        "sharp_coarse": tuple(sorted(sharp_coarse)),
        "sharp_edges": tuple(sorted(sharp_edges)),
        "de_support_floors": ((7, 15), (7, 16)),
        "generating": tuple(generating),
        "hostile": (hostile_vertices, hostile_edge_data,
                    hostile_infinity, hostile_affine,
                    hostile_area2, hostile_boundary, hostile_genus),
    }
    semantic_digest = sha256(
        dumps(semantic_object, sort_keys=True, separators=(",", ":"),
              default=list).encode("utf-8")).hexdigest()

    print("THM-4103 INDEPENDENT THETA-BOUNDARY RESPONSE AUDIT")
    print("status=PASS;FINITE-EXACT;JC2=OPEN;no_sympy=1")
    print("support_pattern_count=8;live_patterns=2;live_optional=phi;"
          "lambda_or_alpha_zero_hostiles=6")
    print("phi_zero_gate=PASS;optional_phi_points=(3,3),(1,4);polygon_unchanged=1")
    print(f"polygon={vertices}")
    print(f"pick=area2:{area2},boundary:{boundary},genus:{genus}")
    for row in edge_data:
        print("edge=" + repr(row))
    print(f"source_infinity_packet={infinity_packet};pole_weight={pole_weight}")
    print(f"affine_s0_packet={affine_packet}")
    print(f"riemann_hurwitz={ramification};expected={2 * genus - 2}")
    print(f"allowed_eisenstein_degrees={allowed_degrees}")
    print("norm_witnesses=" + repr(tuple((degree, witnesses[degree])
                                         for degree in allowed_degrees)))
    print("labelled_by_degree=" + repr(tuple(sorted(labelled_by_degree.items()))))
    print("response_counts=coarse:16,edge_refined:23,labelled:45")
    print("AB_target_infinity=1;BC_quadratic_orbit=joint_response")
    print("sharpened_labelled=" + repr(tuple(sorted(sharp_labelled.items()))))
    print("sharpened_counts=degrees:3,coarse:3,edge_refined:4,labelled:5")
    print("DE_support_floors=A:(7,15),C:(7,16)")
    print("generating_coefficients=" + repr(tuple(generating)))
    print(f"hostile_polygon={hostile_vertices}")
    print(f"hostile_pick=area2:{hostile_area2},boundary:{hostile_boundary},genus:{hostile_genus}")
    print(f"hostile_infinity_packet={hostile_infinity};pole_weight={sum(hostile_infinity)}")
    print(f"hostile_affine_packet={hostile_affine};riemann_hurwitz={sum(index - 1 for index in hostile_infinity)}")
    print(f"semantic_sha256={semantic_digest}")


if __name__ == "__main__":
    main()
