#!/usr/bin/env python3
"""Exact certificate for the generic exact-M=9 JC(2,3) gate in THM-4147.

It checks three nonempty coefficient controls, the complete weight-nine
monomial universe, valued lower supports and Newton polygons, boundary
packets, two independent critical projections, and the finite/full
monodromy budgets.  The canonical theorem supplies the geometric transport
and permutation proofs.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd

import sympy as sp


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valued_lower_support(coefficients):
    """Collapsed Q-adic lower support of the source-fibre equation.

    The dominated Q*S^2 contribution has the same (S,P) exponent as the
    valuation-zero S^2 term, so its Q-adic valuation is zero and it does not
    appear as a separate height-one point.
    """
    raw = [(2, 0, 0, F(1)), (0, 1, 0, F(-1))]
    for (i, j), coefficient in coefficients.items():
        raw.append((j + 2, i + j, 1, -coefficient))
        raw.append((j, i + j + 1, 1, coefficient))
    answer = {}
    for x, y, z, coefficient in raw:
        key = (x, y, z)
        answer[key] = answer.get(key, F(0)) + coefficient
    return {key: value for key, value in answer.items() if value}


def lower_faces(points):
    items = tuple(sorted(points))
    answer = {}
    for triple in combinations(range(len(items)), 3):
        q0, q1, q2 = (items[index] for index in triple)
        x0, y0, z0 = q0
        x1, y1, z1 = q1
        x2, y2, z2 = q2
        determinant = ((x1 - x0) * (y2 - y0)
                       - (x2 - x0) * (y1 - y0))
        if determinant == 0:
            continue
        ax = F((z1 - z0) * (y2 - y0)
               - (z2 - z0) * (y1 - y0), determinant)
        ay = F((x1 - x0) * (z2 - z0)
               - (x2 - x0) * (z1 - z0), determinant)
        constant = F(z0) - ax * x0 - ay * y0
        gaps = {
            point: F(point[2]) - ax * point[0] - ay * point[1] - constant
            for point in items
        }
        if min(gaps.values()) >= 0:
            face = frozenset(point for point, gap in gaps.items() if gap == 0)
            answer[face] = (ax, ay, constant)
    return answer


def convex_hull(points):
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

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
    vertices = convex_hull((x, y) for x, y, _ in points)
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]))
        for i in range(len(vertices))
    )
    require((area2 - boundary + 2) % 2 == 0, "Pick parity failed")
    genus = (area2 - boundary + 2) // 2
    return vertices, area2, boundary, genus


def edge_ledger(vertices):
    answer = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        answer.append((start, end, inward, length, constant, index))
    return tuple(answer)


def valuation(poly, variable):
    dictionary = sp.Poly(poly, variable).as_dict()
    require(bool(dictionary), "zero polynomial has no valuation")
    return min(key[0] for key in dictionary)


def polynomial_digest(poly):
    monic = poly.monic()
    coefficients = []
    for coefficient in monic.all_coeffs():
        coefficient = sp.Rational(coefficient)
        coefficients.append(f"{coefficient.p}/{coefficient.q}")
    payload = f"degree={monic.degree()};" + ",".join(coefficients)
    return sha256(payload.encode()).hexdigest()


def permutation_index(permutation):
    seen = set()
    cycles = 0
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycles += 1
        point = start
        while point not in seen:
            seen.add(point)
            point = permutation[point]
    return len(permutation) - cycles


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    answer = [0] * len(permutation)
    for source, target in enumerate(permutation):
        answer[target] = source
    return tuple(answer)


def make_coefficients(K, Phi, Delta, Theta, eta, zeta):
    answer = {
        (1, 0): F(-3),
        (2, 0): F(8, 3),
        (3, 0): F(-1376, 135),
        (0, 2): K,
        (2, 1): Phi,
        (4, 0): Delta,
        (1, 2): Theta,
    }
    if eta:
        answer[(3, 1)] = eta
    if zeta:
        answer[(0, 3)] = zeta
    return answer


def verify_source_monomial_universe():
    """Enumerate the complete b=d=0 source support through weight nine."""
    weighted = {
        (i, j)
        for i in range(5)
        for j in range(4)
        if 0 < 2 * i + 3 * j <= 9
    }
    # The live normal form removes [y]R and [p*y]R exactly.
    normalized = weighted - {(0, 1), (1, 1)}
    require(normalized == {
        (1, 0), (2, 0), (3, 0), (4, 0),
        (0, 2), (2, 1), (1, 2), (3, 1), (0, 3),
    }, "exact weight-nine source monomial universe changed")


def verify_cubic_carrier_discriminant():
    """Freeze every finite Hurwitz value of the cubic carrier."""
    w, q, kappa, cubic = sp.symbols("w q kappa cubic")
    polynomial = cubic * w ** 3 + kappa * w ** 2 - (q - sp.Rational(1, 2))
    expected = ((q - sp.Rational(1, 2))
                * (4 * kappa ** 3
                   - 27 * cubic ** 2 * (q - sp.Rational(1, 2))))
    require(sp.factor(sp.discriminant(polynomial, w) - expected) == 0,
            "cubic carrier Hurwitz discriminant changed")


def main():
    verify_source_monomial_universe()
    verify_cubic_carrier_discriminant()
    X, T, s, p = sp.symbols("X T s p")
    P = T + X ** 2 * T ** 2
    Y = X * T * P

    Delta = F(1)
    K = F(2848, 45) - F(7, 6) * Delta
    Phi = F(11, 7)
    Theta = F(19, 11)
    eta = F(23, 13)
    zeta = F(29, 17)
    require(K == F(5591, 90), "forced K row changed")

    cases = (
        {
            "name": "p3y_only",
            "eta": eta,
            "zeta": F(0),
            "polygon": ((0, 1), (2, 0), (4, 2), (4, 3),
                        (3, 4), (1, 5), (0, 5)),
            "planes": {
                (F(0), F(1, 4), F(-1, 4)),
                (F(1, 9), F(2, 9), F(-2, 9)),
                (F(1, 5), F(1, 5), F(-2, 5)),
                (F(1, 2), F(0), F(-1)),
            },
            "area2": 29,
            "boundary": 11,
            "genus": 10,
            "packet": (8, 5, 4, 3, 2, 2, 1),
            "residual_degree": 20,
            "independent_p_valuation": 4,
            "finite_degree": 21,
            "full_degree": 25,
            "carrier_index": 2,
            "q0": -F(3 ** 15, 2 ** 7) * eta ** 7,
            "lc": F(72900) * eta ** 5 * K ** 4 * Theta ** 4,
        },
        {
            "name": "y3_only",
            "eta": F(0),
            "zeta": zeta,
            "polygon": ((0, 1), (2, 0), (5, 3), (3, 4), (0, 5)),
            "planes": {
                (F(1, 12), F(1, 4), F(-1, 4)),
                (F(1, 9), F(2, 9), F(-2, 9)),
            },
            "area2": 30,
            "boundary": 10,
            "genus": 11,
            "packet": (11, 8, 2, 2, 2, 1),
            "residual_degree": 20,
            "independent_p_valuation": 8,
            "finite_degree": 20,
            "full_degree": 26,
            "carrier_index": 3,
            "q0": -F(3 ** 15, 2 ** 7) * zeta ** 7,
            "lc": None,
        },
        {
            "name": "both_top_terms",
            "eta": eta,
            "zeta": zeta,
            "polygon": ((0, 1), (2, 0), (5, 3), (1, 5), (0, 5)),
            "planes": {
                (F(0), F(1, 4), F(-1, 4)),
                (F(1, 9), F(2, 9), F(-2, 9)),
            },
            "area2": 31,
            "boundary": 11,
            "genus": 11,
            "packet": (8, 8, 4, 2, 2, 2, 1),
            "residual_degree": 21,
            "independent_p_valuation": 8,
            "finite_degree": 21,
            "full_degree": 27,
            "carrier_index": 3,
            "q0": -F(3 ** 15, 2 ** 7) * (eta + zeta) ** 7,
            "lc": None,
        },
    )

    semantic_rows = []
    printed_rows = []
    for case in cases:
        name = case["name"]
        eta_case = case["eta"]
        zeta_case = case["zeta"]
        coefficients = make_coefficients(
            K, Phi, Delta, Theta, eta_case, zeta_case
        )
        support = valued_lower_support(coefficients)
        faces = lower_faces(support)
        vertices, area2, boundary, genus = polygon_ledger(support)
        require(vertices == case["polygon"], f"{name}: polygon changed")
        require(set(faces.values()) == case["planes"],
                f"{name}: lower-face planes changed")
        require((area2, boundary, genus)
                == (case["area2"], case["boundary"], case["genus"]),
                f"{name}: Pick ledger changed")

        edges = edge_ledger(vertices)
        toric_indices = []
        for start, end, inward, length, constant, index in edges:
            # The x=0 vertical edge consists of four affine s=0 points.
            if start[0] == end[0] == 0:
                require(length == 4 and index == 1,
                        f"{name}: affine vertical edge changed")
                continue
            toric_indices.extend([index] * length)
        packet = tuple(sorted(toric_indices, reverse=True))
        require(packet == case["packet"], f"{name}: packet changed")
        defect = sum(index - 1 for index in packet)
        require(defect == 2 * genus - 2,
                f"{name}: boundary no longer saturates Riemann--Hurwitz")

        G = (
            -X ** 2 * T / 2
            - 3 * P
            + sp.Rational(8, 3) * P ** 2
            - sp.Rational(1376, 135) * P ** 3
            + sp.Rational(K.numerator, K.denominator) * Y ** 2
            + sp.Rational(Phi.numerator, Phi.denominator) * P ** 2 * Y
            + sp.Rational(Delta.numerator, Delta.denominator) * P ** 4
            + sp.Rational(Theta.numerator, Theta.denominator) * P * Y ** 2
            + sp.Rational(eta_case.numerator, eta_case.denominator) * P ** 3 * Y
            + sp.Rational(zeta_case.numerator, zeta_case.denominator) * Y ** 3
        )
        f = sp.cancel(sp.diff(G, X) / T)
        h = sp.diff(G, T)
        resultant = sp.resultant(f, h, X)
        require(valuation(resultant, T) == 56,
                f"{name}: primary T-valuation changed")
        quotient = sp.cancel(resultant / (T ** 56 * (6 * T + 1) ** 2))
        require(sp.denom(quotient) == 1, f"{name}: quotient not polynomial")
        qpoly = sp.Poly(quotient, T)
        require(qpoly.degree() == case["residual_degree"],
                f"{name}: primary residual degree changed")
        require(qpoly.nth(0) == sp.Rational(case["q0"].numerator,
                                            case["q0"].denominator),
                f"{name}: primary residual constant changed")
        if case["lc"] is not None:
            expected_lc = sp.Rational(case["lc"].numerator,
                                      case["lc"].denominator)
            require(qpoly.LC() == expected_lc,
                    f"{name}: primary residual leading row changed")
        require(qpoly.LC() != 0 and sp.gcd(qpoly, qpoly.diff()).degree() == 0,
                f"{name}: control is not resultant-generic")
        require(qpoly.eval(-sp.Rational(1, 6)) != 0,
                f"{name}: residual meets the universal T=-1/6 fibre")

        # Independent rational (s,p) projection.
        t = p - s ** 2
        y = s * p
        H = (
            -3 * p
            + sp.Rational(8, 3) * p ** 2
            - sp.Rational(1376, 135) * p ** 3
            + sp.Rational(K.numerator, K.denominator) * y ** 2
            + sp.Rational(Phi.numerator, Phi.denominator) * p ** 2 * y
            + sp.Rational(Delta.numerator, Delta.denominator) * p ** 4
            + sp.Rational(Theta.numerator, Theta.denominator) * p * y ** 2
            + sp.Rational(eta_case.numerator, eta_case.denominator) * p ** 3 * y
            + sp.Rational(zeta_case.numerator, zeta_case.denominator) * y ** 3
        )
        Gsp = -s ** 2 / (2 * t) + H
        A = sp.cancel(t ** 2 * sp.diff(Gsp, s) / p)
        C = sp.cancel(2 * t ** 2 * sp.diff(Gsp, p))
        require(sp.denom(A) == sp.denom(C) == 1,
                f"{name}: independent critical pair is not polynomial")
        independent_resultant = sp.resultant(A, C, s)
        p_valuation = valuation(independent_resultant, p)
        require(p_valuation == case["independent_p_valuation"],
                f"{name}: independent p-valuation changed")
        independent_residual = sp.cancel(independent_resultant / p ** p_valuation)
        require(sp.denom(independent_residual) == 1,
                f"{name}: independent residual not polynomial")
        independent_poly = sp.Poly(independent_residual, p)
        require(independent_poly.degree() == case["residual_degree"],
                f"{name}: independent residual degree changed")
        require(independent_poly.nth(0) != 0
                and independent_poly.LC() != 0
                and sp.gcd(independent_poly, independent_poly.diff()).degree() == 0,
                f"{name}: independent residual control degenerated")

        critical_length = case["residual_degree"] + 4
        finite_degree = case["finite_degree"]
        full_degree = case["full_degree"]
        beta = case["carrier_index"]

        # Finite response: if one handle generator is nonidentity, its pair
        # contributes at most 2n-L-1 merger units; if both are identities,
        # only beta remains.  Both alternatives miss transitivity.
        finite_support_sum = 2 * finite_degree - critical_length
        finite_one_handle_capacity = finite_support_sum - 1 + beta
        finite_two_handle_capacity = finite_support_sum - 2 + beta
        require(beta < finite_degree - 1,
                f"{name}: both-identity finite hostile not excluded")
        require(finite_one_handle_capacity < finite_degree - 1,
                f"{name}: one-handle finite budget no longer excludes")
        require(finite_two_handle_capacity < finite_degree - 1,
                f"{name}: two-handle finite budget no longer excludes")
        require(critical_length > finite_degree + beta,
                f"{name}: compact finite-response inequality changed")

        # Full response: X,Y alone are transitive, so their supports cover all
        # sheets and k=|supp X intersect supp Y| <= n-L.  The proved theorem
        # lemma ind([X,Y])<=2k then conflicts with the boundary defect.
        overlap_cap = full_degree - critical_length
        commutator_index_cap = 2 * overlap_cap
        require(defect > commutator_index_cap,
                f"{name}: full-response commutator gate no longer excludes")

        row = (
            f"{name}:polygon={vertices};planes={sorted(faces.values())};"
            f"Pick=({area2},{boundary},{genus});packet={packet};defect={defect};"
            f"primary=T^56*(6T+1)^2*Q{qpoly.degree()};"
            f"independent=p^{p_valuation}*R{independent_poly.degree()};"
            "Tminus1over6=universal-only;"
            f"L={critical_length};finite=(n={finite_degree},beta={beta},"
            f"cap1={finite_one_handle_capacity}<{finite_degree - 1});"
            f"full=(n={full_degree},k<={overlap_cap},"
            f"comm_ind<={commutator_index_cap}<{defect})"
        )
        semantic_rows.append(row)
        printed_rows.append((name, vertices, sorted(faces.values()), area2,
                             boundary, genus, packet, defect, qpoly.degree(),
                             p_valuation, independent_poly.degree(),
                             critical_length, finite_degree, beta,
                             finite_one_handle_capacity, full_degree,
                             overlap_cap, commutator_index_cap,
                             qpoly.nth(0), qpoly.LC(),
                             polynomial_digest(qpoly),
                             polynomial_digest(independent_poly)))

    # The genuinely repeated weight-nine top edge eta+zeta=0 is not part of
    # the theorem chamber.  It changes the critical strict transform, unlike
    # the harmless midpoint support cancellation eta=zeta.
    collision_zeta = -eta
    G_collision = (
        -X ** 2 * T / 2 - 3 * P + sp.Rational(8, 3) * P ** 2
        - sp.Rational(1376, 135) * P ** 3
        + sp.Rational(K.numerator, K.denominator) * Y ** 2
        + sp.Rational(Phi.numerator, Phi.denominator) * P ** 2 * Y
        + sp.Rational(Delta.numerator, Delta.denominator) * P ** 4
        + sp.Rational(Theta.numerator, Theta.denominator) * P * Y ** 2
        + sp.Rational(eta.numerator, eta.denominator) * P ** 3 * Y
        + sp.Rational(collision_zeta.numerator,
                      collision_zeta.denominator) * Y ** 3
    )
    collision_resultant = sp.resultant(
        sp.cancel(sp.diff(G_collision, X) / T),
        sp.diff(G_collision, T), X
    )
    require(valuation(collision_resultant, T) == 42,
            "top-collision strict-transform valuation changed")
    collision_quotient = sp.cancel(
        collision_resultant / (T ** 42 * (6 * T + 1) ** 2)
    )
    require(sp.denom(collision_quotient) == 1
            and sp.degree(collision_quotient, T) == 19,
            "top-collision residual degree changed")

    equal_coefficients = make_coefficients(
        K, Phi, Delta, Theta, eta, eta
    )
    equal_support = valued_lower_support(equal_coefficients)
    require((3, 4, 1) not in equal_support,
            "eta=zeta did not cancel the lifted midpoint")
    require(convex_hull((x, y) for x, y, _ in equal_support)
            == cases[2]["polygon"],
            "eta=zeta changed the generic top polygon")
    G_equal = (
        -X ** 2 * T / 2 - 3 * P + sp.Rational(8, 3) * P ** 2
        - sp.Rational(1376, 135) * P ** 3
        + sp.Rational(K.numerator, K.denominator) * Y ** 2
        + sp.Rational(Phi.numerator, Phi.denominator) * P ** 2 * Y
        + sp.Rational(Delta.numerator, Delta.denominator) * P ** 4
        + sp.Rational(Theta.numerator, Theta.denominator) * P * Y ** 2
        + sp.Rational(eta.numerator, eta.denominator) * P ** 3 * Y
        + sp.Rational(eta.numerator, eta.denominator) * Y ** 3
    )
    equal_resultant = sp.resultant(
        sp.cancel(sp.diff(G_equal, X) / T), sp.diff(G_equal, T), X
    )
    require(valuation(equal_resultant, T) == 56,
            "harmless midpoint cancellation changed T-valuation")
    equal_quotient = sp.cancel(
        equal_resultant / (T ** 56 * (6 * T + 1) ** 2)
    )
    require(sp.denom(equal_quotient) == 1
            and sp.degree(equal_quotient, T) == 21,
            "harmless midpoint cancellation changed residual degree")

    # Exhaustive hostile check of the commutator-overlap inequality on S_n,
    # n<=5.  The report contains the general partial-permutation proof.
    for degree in range(1, 6):
        from itertools import permutations
        all_permutations = tuple(permutations(range(degree)))
        for left in all_permutations:
            left_inverse = inverse(left)
            left_support = {i for i, value in enumerate(left) if i != value}
            for right in all_permutations:
                right_inverse = inverse(right)
                right_support = {i for i, value in enumerate(right) if i != value}
                commutator = compose(
                    compose(compose(left, right), left_inverse), right_inverse
                )
                overlap = len(left_support & right_support)
                require(permutation_index(commutator) <= 2 * overlap,
                        "commutator-overlap hostile failed")

    semantic = "\n".join(semantic_rows) + "\n"
    print("THM4147 GENERIC EXACT-M=9 CERTIFICATE")
    print(f"checks={CHECKS}")
    print("forced_row: K=2848/45-(7/6)Delta; control Delta=1 gives K=5591/90")
    print("top_row: eta*p^3*y+zeta*y^3, (eta,zeta)!=(0,0)")
    print("cubic_carrier_discriminant="
          "(q-1/2)*(4K^3-27zeta^2*(q-1/2)):PASS")
    for row in printed_rows:
        (name, vertices, planes, area2, boundary, genus, packet, defect,
         primary_degree, p_valuation, independent_degree, critical_length,
         finite_degree, beta, finite_cap, full_degree, overlap_cap,
         commutator_cap, q0, leading, primary_digest,
         independent_digest) = row
        print(f"case={name}")
        print(f"  polygon={vertices}")
        print(f"  lower_planes={planes}")
        print(f"  Pick=area2:{area2},boundary:{boundary},genus:{genus}")
        print(f"  packet={packet}, defect={defect}")
        print(f"  primary=T^56*(6T+1)^2*Q{primary_degree}; "
              f"independent=p^{p_valuation}*R{independent_degree}; "
              f"critical_length={critical_length}")
        print(f"  Q_constant={q0}; Q_leading={leading}")
        print(f"  primary_residual_sha256={primary_digest}")
        print(f"  independent_residual_sha256={independent_digest}")
        print(f"  finite_response: n={finite_degree}, beta={beta}, "
              f"one_handle_capacity={finite_cap} < n-1={finite_degree - 1}")
        print(f"  full_response: n={full_degree}, overlap<={overlap_cap}, "
              f"commutator_index<={commutator_cap} < defect={defect}")
    print("top_collision_eta_plus_zeta_zero="
          "T^42*(6T+1)^2*Q19:EXCLUDED_FROM_CHAMBER")
    print("midpoint_cancellation_eta_equals_zeta="
          "support_midpoint_deleted,T^56*(6T+1)^2*Q21:HARMLESS")
    print("commutator_overlap_hostile=S_n exhaustive through n=5:PASS")
    print("critical_open_extra_gate=Q(-1/6)!=0:PASS")
    print(f"semantic_sha256={sha256(semantic.encode()).hexdigest()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
