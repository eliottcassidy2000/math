#!/usr/bin/env python3
"""Exact certificate for THM-4103's theta-only boundary response.

Universe: the complete generic Newton polygon of

    F_Q=(s^2-p)(1-Q H(p,sp))+gamma*Q*s^2

on THM-4053's live ``delta=0, theta!=0, Delta_V!=0`` survivor.  The six
polygon vertices and the lambda/alpha pairs are mandatory; only phi is
optional.  All eight lambda/alpha/phi presence patterns are checked, so the
six patterns deleting lambda or alpha serve as robustness hostiles.

This verifies the lattice, ramification, degree, response-profile, and
support-deletion ledgers used in the proof.  No Python assertions are used,
so ``python`` and ``python -O`` execute the same gates.
"""

from hashlib import sha256
from itertools import product
from math import comb, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cross(origin, left, right):
    return ((left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0]))


def convex_hull(points):
    points = sorted(set(points))
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
    return lower[:-1] + upper[:-1]


def edge_ledger(polygon, residue_point):
    rows = []
    area_twice = 0
    boundary_length = 0
    canonical_degree = 0
    for left, right in zip(polygon, polygon[1:] + polygon[:1]):
        dx = right[0] - left[0]
        dy = right[1] - left[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * left[0] + inward[1] * left[1]
        distance = (inward[0] * residue_point[0]
                    + inward[1] * residue_point[1] - constant)
        order = distance - 1
        rows.append((left, right, length, inward, distance,
                     order, length * order))
        area_twice += left[0] * right[1] - right[0] * left[1]
        boundary_length += length
        canonical_degree += length * order
    genus = (area_twice - boundary_length + 2) // 2
    return rows, area_twice, boundary_length, genus, canonical_degree


def eisenstein_local_criterion(value):
    remaining = value
    prime = 2
    while prime * prime <= remaining:
        valuation = 0
        while remaining % prime == 0:
            remaining //= prime
            valuation += 1
        if prime % 3 == 2 and valuation % 2:
            return False
        prime += 1
    return not (remaining > 1 and remaining % 3 == 2)


def main():
    mandatory = {
        (0, 1), (2, 0), (4, 2), (4, 3), (2, 4), (0, 4),
        # Forced epsilon/kappa support which is not a polygon vertex.
        (2, 3),
    }
    optional_pairs = (
        {(2, 1), (0, 2)},       # lambda (forced live; deletion is hostile)
        {(2, 2), (0, 3)},       # alpha (forced live; deletion is hostile)
        {(3, 3), (1, 4)},       # phi; phi=0 is permitted
    )
    expected_polygon = [(0, 1), (2, 0), (4, 2),
                        (4, 3), (2, 4), (0, 4)]
    support_branches = []
    for flags in product((0, 1), repeat=3):
        support = set(mandatory)
        for flag, pair in zip(flags, optional_pairs):
            if flag:
                support.update(pair)
        hull = convex_hull(support)
        require(hull == expected_polygon,
                f"Newton polygon changed for optional flags {flags}: {hull}")
        support_branches.append((flags, tuple(sorted(support))))

    residue_point = (1, 1)
    rows, area2, boundary, genus, ramification_total = edge_ledger(
        expected_polygon, residue_point)
    require((area2, boundary, genus) == (24, 10, 8), "Pick ledger changed")
    require([row[3] for row in rows]
            == [(1, 2), (-1, 1), (-1, 0),
                (-1, -2), (0, -1), (1, 0)],
            "primitive inward normals changed")
    require([row[4] for row in rows] == [1, 2, 3, 7, 3, 1],
            "residue distances changed")
    for row in rows:
        normal = row[3]
        constant = normal[0] * row[0][0] + normal[1] * row[0][1]
        require(all(normal[0] * point[0] + normal[1] * point[1] >= constant
                    for point in mandatory),
                f"mandatory support escaped half-space {row[0:2]}")
    require(ramification_total == 14 == 2 * genus - 2,
            "canonical/Riemann-Hurwitz checksum changed")

    # The last edge has s=0 and three nonzero p roots, hence (x,t)=(0,p)
    # in the affine source.  Every point on the first five edges is at source
    # infinity.  A residue order e-1 is a ramification index e.
    affine_indices = [rows[-1][4]] * rows[-1][2]
    infinity_indices = []
    for row in rows[:-1]:
        infinity_indices.extend([row[4]] * row[2])
    require(affine_indices == [1, 1, 1], "affine edge packet changed")
    require(sorted(infinity_indices) == [1, 2, 2, 3, 3, 3, 7],
            "source-infinity packet changed")
    require(sum(index - 1 for index in infinity_indices) == 14,
            "source-infinity ramification total changed")
    require(sum(infinity_indices) == 21, "source-infinity pole weight changed")

    norms = {a * a - a * b + b * b
             for a in range(-30, 31) for b in range(-30, 31)}
    allowed_degrees = sorted(value for value in norms if 7 <= value <= 21)
    require(allowed_degrees == [7, 9, 12, 13, 16, 19, 21],
            f"Eisenstein degree list changed: {allowed_degrees}")
    require(all((value in norms) == eisenstein_local_criterion(value)
                for value in range(1, 22)),
            "representation and local norm tests disagree")

    # A profile (a,b,c,d) counts finite-image punctures of index 1,2,3,7.
    profiles = {}
    for degree in allowed_degrees:
        debt = 21 - degree
        degree_profiles = []
        labelled_count = 0
        for count1 in range(2):
            for count2 in range(3):
                for count3 in range(4):
                    for count7 in range(2):
                        if count1 + 2 * count2 + 3 * count3 + 7 * count7 != debt:
                            continue
                        multiplicity = comb(2, count2) * comb(3, count3)
                        degree_profiles.append(
                            (count1, count2, count3, count7, multiplicity))
                        labelled_count += multiplicity
        profiles[degree] = (debt, tuple(degree_profiles), labelled_count)
    require(sum(len(row[1]) for row in profiles.values()) == 16,
            "coarse profile count changed")
    require(sum(row[2] for row in profiles.values()) == 45,
            "labelled profile count changed")

    labelled_by_degree = {degree: 0 for degree in allowed_degrees}
    for mask in range(1 << len(infinity_indices)):
        debt = sum(index for bit, index in enumerate(infinity_indices)
                   if (mask >> bit) & 1)
        degree = 21 - debt
        if degree in labelled_by_degree:
            labelled_by_degree[degree] += 1
    require([labelled_by_degree[degree] for degree in allowed_degrees]
            == [7, 9, 9, 10, 7, 2, 1],
            f"literal labelled census changed: {labelled_by_degree}")
    require(sum(labelled_by_degree.values()) == 45,
            "literal labelled total changed")

    edge_profile_count = 0
    for degree in allowed_degrees:
        debt = 21 - degree
        for count1 in range(2):
            for count2 in range(3):
                for count3_vertical in range(2):
                    for count3_horizontal in range(3):
                        for count7 in range(2):
                            weight = (count1 + 2 * count2
                                      + 3 * count3_vertical
                                      + 3 * count3_horizontal + 7 * count7)
                            edge_profile_count += (weight == debt)
    require(edge_profile_count == 23, "edge-refined profile count changed")

    # Algebraic sidecars proved in THM-4103 sharpen this combinatorial census:
    # the AB index-one puncture maps to target infinity, while the two BC
    # index-two punctures are one quadratic closed point over k(q) and respond
    # together.  Thus count1=0 and count2 is 0 or 2.
    sharpened_profiles = {}
    sharpened_edge_count = 0
    sharpened_labelled_count = 0
    for degree in allowed_degrees:
        kept = tuple(profile for profile in profiles[degree][1]
                     if profile[0] == 0 and profile[1] in (0, 2))
        if not kept:
            continue
        labelled = sum(profile[4] for profile in kept)
        edge_count = 0
        for profile in kept:
            count3 = profile[2]
            edge_count += sum(vertical + horizontal == count3
                              for vertical in range(2)
                              for horizontal in range(3))
        sharpened_profiles[degree] = (kept, labelled, edge_count)
        sharpened_labelled_count += labelled
        sharpened_edge_count += edge_count
    require(sharpened_profiles == {
        7: (((0, 2, 1, 1, 3),), 3, 2),
        12: (((0, 0, 3, 0, 1),), 1, 1),
        21: (((0, 0, 0, 0, 1),), 1, 1),
    }, f"sharpened profile census changed: {sharpened_profiles}")
    require((len(sharpened_profiles), sharpened_edge_count,
             sharpened_labelled_count) == (3, 4, 5),
            "sharpened 3/4/5 profile totals changed")

    # DE hostile dichotomy.  Apart from the unique deleted-line monomials
    # x^2 and x^3, the Laurent depth caps give m<=n+1 in A and m<=n+2 in C.
    # Any term able to compete with their DE orders -14 and -21 therefore
    # has the exact sharp support floors below.  If no such term exists, the
    # unique deleted-line leading coefficient is incompatible with either a
    # finite image or the index-seven target pole coefficient.
    a_competitors = [(n, m, m + n, -7 * m + 6 * n)
                     for n in range(1, 80) for m in range(n + 2)
                     if -7 * m + 6 * n <= -14]
    c_competitors = [(n, m, m + n, -7 * m + 6 * n)
                     for n in range(1, 80) for m in range(n + 3)
                     if -7 * m + 6 * n <= -21]
    a_floor = (min(row[0] for row in a_competitors),
               min(row[2] for row in a_competitors))
    c_floor = (min(row[0] for row in c_competitors),
               min(row[2] for row in c_competitors))
    require(a_floor == (7, 15) and (7, 8, 15, -14) in a_competitors,
            f"DE A support floor changed: {a_floor}")
    require(c_floor == (7, 16) and (7, 9, 16, -21) in c_competitors,
            f"DE C support floor changed: {c_floor}")

    # Hostile: deleting the forced kappa vertex preserves genus and total
    # canonical degree, but changes the boundary packet.
    hostile_support = mandatory - {(4, 2)}
    hostile_polygon = convex_hull(hostile_support)
    hostile_rows, hostile_area2, hostile_boundary, hostile_genus, hostile_ram = (
        edge_ledger(hostile_polygon, residue_point))
    hostile_indices = []
    for row in hostile_rows:
        hostile_indices.extend([row[4]] * row[2])
    require((hostile_area2, hostile_boundary, hostile_genus, hostile_ram)
            == (22, 8, 8, 14), "support-deletion hostile checksum changed")
    require(sorted(hostile_indices) == [1, 1, 1, 1, 3, 3, 5, 7],
            "support-deletion hostile packet changed")

    semantics = repr((support_branches, expected_polygon, rows,
                      affine_indices, infinity_indices, allowed_degrees,
                      profiles, labelled_by_degree, edge_profile_count,
                      sharpened_profiles, sharpened_edge_count,
                      sharpened_labelled_count, a_floor, c_floor,
                      hostile_polygon, hostile_indices))
    digest = sha256(semantics.encode("utf-8")).hexdigest()

    print("JC23 THETA-ONLY BOUNDARY RESPONSE CERTIFICATE")
    print("status=FINITE-EXACT;THM4103=PROVED_RELATIVE_TO_THM3992_4053;JC2=OPEN")
    print(f"support_patterns={len(support_branches)};live_patterns=2;"
          "live_optional=phi;lambda_or_alpha_zero_hostiles=6")
    print(f"polygon={expected_polygon}")
    print(f"pick=area2:{area2},boundary:{boundary},genus:{genus}")
    for row in rows:
        print("edge=" + repr(row))
    print(f"ramification_total={ramification_total}")
    print(f"affine_x0_indices={affine_indices}")
    print(f"source_infinity_indices={infinity_indices};weight={sum(infinity_indices)}")
    print(f"allowed_eisenstein_degrees={allowed_degrees}")
    for degree in allowed_degrees:
        debt, degree_profiles, labelled_count = profiles[degree]
        print(f"degree={degree};finite_debt={debt};profiles={degree_profiles};"
              f"labelled={labelled_count}")
    print("profile_totals=coarse:16,edge_refined:23,labelled:45")
    print("AB_target_infinity=1;BC_quadratic_orbit=joint_response")
    print(f"sharpened_profiles={sharpened_profiles}")
    print("sharpened_totals=degrees:3,coarse:3,edge_refined:4,labelled:5")
    print(f"DE_support_floors=A:{a_floor},C:{c_floor}")
    print(f"support_deletion_hostile_polygon={hostile_polygon}")
    print(f"support_deletion_hostile_indices={hostile_indices}")
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()
