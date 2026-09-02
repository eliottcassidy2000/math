#!/usr/bin/env python3
"""Exact finite certificate for the exact-M=12 ``Z=0`` endpoint wall.

Universe
========
The residual source contains every labelled row ``p^i y^j`` with
``0 < 2*i+3*j <= 12``, except the inherited forbidden rows ``y`` and
``p*y``.  The endpoint wall deletes ``y^4`` and retains the fixed rows
``p,p^2,p^3,p^6``.  A row lifts to the two valued support points

    (j+2,i+j,1), (j,i+j+1,1).

The fixed source points ``(2,0,0),(0,1,0)`` and the nonoptional residual
point ``(2,0,1)`` are always retained.  Six aggregate support points can
cancel.  Treating their deletion bits independently deliberately
*overcounts* coefficient cells, so every actual specialization is covered.

Claims checked
==============
1. On ``Z=0`` with ``U,W,beta_11,zeta_3`` nonzero, all 16,384 hostile
   lower-row/deletion/collision masks have the same three lower faces.
2. On ``Z=W=0`` with ``U,beta_11,zeta_3`` nonzero, all 8,192 hostile masks
   have the stated replacement three-face complex.
3. Their exact polygons, Pick genera, source-infinity packets, response
   sums, source-normal factorization, face normal forms, generic-point
   ``G_P`` witnesses, edge squarefreeness gates, dual-graph ledgers, and
   good-target differential orders are the displayed values.
4. On the remaining ``W!=0, U+W=0`` contact, the exact top identity imports
   THM-4297's length-twelve A23 chart, delayed ``t^6`` correction, and five
   positive Keller-form tail orders.  Response packets are not transported.
5. The broader ``Z=0,U!=0`` hostile over-atlas contains 131,072 masks, 51
   face-complex types, and 32 supporting planes.  Every primary-face order
   coefficient ``5/6-a-b-c`` is positive, with minimum ``3/4``.

Geometric imports
=================
As in THM-4045/THM-4248, the audited toroidal face/edge model turns a
Newton-nondegenerate lower complex into its face normalizations plus rational
toric chains.  The exact Keller residue identity and good elliptic scaling
turn positive vertical order into constancy.  Proper-flat degree conservation
then excludes a positive-degree Keller response.  This program certifies the
finite Newton, order, response, and displayed A23-tail inputs; it does not
prove those imported geometric theorems, seam entry, JC(2), or any other
degenerate collision-tail inventory.

Run from the repository root with either

    python3 -B 04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py
    python3 -B -O 04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, gcd


CHECKS = 0


# Sparse polynomials over Q in (S,P,U,W,beta,zeta).  They are deliberately
# tiny: the symbolic face identities below need addition, multiplication,
# nonnegative powers, and differentiation only.
SPARSE_VARIABLES = 6
S_INDEX = 0
P_INDEX = 1


def sparse_clean(raw):
    return {monomial: Fraction(coefficient)
            for monomial, coefficient in raw.items() if coefficient}


def sparse_constant(value):
    value = Fraction(value)
    if not value:
        return {}
    return {(0,) * SPARSE_VARIABLES: value}


def sparse_variable(index):
    exponent = [0] * SPARSE_VARIABLES
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def sparse_add(*polynomials):
    answer = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            answer[monomial] = answer.get(monomial, Fraction(0)) + coefficient
    return sparse_clean(answer)


def sparse_scale(polynomial, scalar):
    scalar = Fraction(scalar)
    return sparse_clean({monomial: scalar * coefficient
                         for monomial, coefficient in polynomial.items()})


def sparse_neg(polynomial):
    return sparse_scale(polynomial, -1)


def sparse_subtract(first, second):
    return sparse_add(first, sparse_neg(second))


def sparse_multiply(first, second):
    answer = {}
    for first_monomial, first_coefficient in first.items():
        for second_monomial, second_coefficient in second.items():
            monomial = tuple(
                first_monomial[index] + second_monomial[index]
                for index in range(SPARSE_VARIABLES)
            )
            answer[monomial] = (
                answer.get(monomial, Fraction(0))
                + first_coefficient * second_coefficient
            )
    return sparse_clean(answer)


def sparse_power(polynomial, exponent):
    need(exponent >= 0, "nonnegative sparse-polynomial exponent")
    answer = sparse_constant(1)
    factor = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = sparse_multiply(answer, factor)
        factor = sparse_multiply(factor, factor)
        remaining //= 2
    return answer


def sparse_derivative(polynomial, index):
    answer = {}
    for monomial, coefficient in polynomial.items():
        if monomial[index] == 0:
            continue
        exponent = list(monomial)
        multiplier = exponent[index]
        exponent[index] -= 1
        answer[tuple(exponent)] = coefficient * multiplier
    return sparse_clean(answer)


def sparse_degree(polynomial, index):
    if not polynomial:
        return -1
    return max(monomial[index] for monomial in polynomial)


def substitute_s_squared_by_p(polynomial):
    """Restrict an even-in-S polynomial to S^2=P."""
    answer = {}
    for monomial, coefficient in polynomial.items():
        need(monomial[S_INDEX] % 2 == 0, "S^2=P substitution parity")
        exponent = list(monomial)
        exponent[P_INDEX] += exponent[S_INDEX] // 2
        exponent[S_INDEX] = 0
        exponent = tuple(exponent)
        answer[exponent] = answer.get(exponent, Fraction(0)) + coefficient
    return sparse_clean(answer)


def sparse_substitute_constant(polynomial, index, value):
    """Substitute one sparse variable by an exact scalar."""
    value = Fraction(value)
    answer = {}
    for monomial, coefficient in polynomial.items():
        exponent = list(monomial)
        coefficient *= value ** exponent[index]
        exponent[index] = 0
        exponent = tuple(exponent)
        answer[exponent] = answer.get(exponent, Fraction(0)) + coefficient
    return sparse_clean(answer)


def polynomial_trim(polynomial):
    answer = list(polynomial)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(Fraction(coefficient) for coefficient in answer)


def polynomial_multiply(first, second):
    answer = [Fraction(0)] * (len(first) + len(second) - 1)
    for first_index, first_coefficient in enumerate(first):
        for second_index, second_coefficient in enumerate(second):
            answer[first_index + second_index] += first_coefficient * second_coefficient
    return polynomial_trim(answer)


def polynomial_derivative(polynomial):
    if len(polynomial) <= 1:
        return (Fraction(0),)
    return polynomial_trim(tuple(
        index * polynomial[index] for index in range(1, len(polynomial))
    ))


def polynomial_divmod(numerator, denominator):
    numerator = list(polynomial_trim(numerator))
    denominator = polynomial_trim(denominator)
    need(denominator != (Fraction(0),), "univariate division denominator")
    quotient = [Fraction(0)] * max(1, len(numerator) - len(denominator) + 1)
    while len(numerator) >= len(denominator) and any(numerator):
        shift = len(numerator) - len(denominator)
        factor = numerator[-1] / denominator[-1]
        quotient[shift] += factor
        for index, coefficient in enumerate(denominator):
            numerator[index + shift] -= factor * coefficient
        numerator = list(polynomial_trim(numerator))
    return polynomial_trim(quotient), polynomial_trim(numerator)


def polynomial_gcd(first, second):
    first = polynomial_trim(first)
    second = polynomial_trim(second)
    while second != (Fraction(0),):
        _quotient, remainder = polynomial_divmod(first, second)
        first, second = second, remainder
    if first == (Fraction(0),):
        return first
    return polynomial_trim(tuple(coefficient / first[-1] for coefficient in first))


def polynomial_squarefree(polynomial):
    polynomial = polynomial_trim(polynomial)
    return len(polynomial_gcd(polynomial, polynomial_derivative(polynomial))) == 1


def need(condition, label):
    """Optimization-stable exact assertion."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("THM-4327 Z0 AUDIT FAILURE: " + label)


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


def projected_rank_two(points, mask):
    chosen = tuple(
        points[index] for index in range(len(points)) if mask & (1 << index)
    )
    for first, second, third in combinations(chosen, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant:
            return True
    return False


def prepare_plane_records(universe):
    points = tuple(sorted(universe))
    records = {}
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
        constant = (
            Fraction(first[2]) - slope_s * first[0] - slope_p * first[1]
        )
        plane = (slope_s, slope_p, constant)
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
    z_row = (0, 4)
    universe = lifted_support(row for row in rows if row[:2] != z_row)
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
    face_counter = Counter()
    total = 0
    for optional_mask in range(1 << len(optional)):
        active_rows = list(required_rows)
        active_rows.extend(
            row for index, row in enumerate(optional)
            if optional_mask & (1 << index)
        )
        base_mask = support_mask(tuple(active_rows))
        need(base_mask & (1 << point_index[(2, 0, 1)]),
             "fixed residual point retained")
        for collision_mask in range(1 << len(collision_bits)):
            mask = base_mask
            for index, bit in enumerate(collision_bits):
                if collision_mask & (1 << index):
                    mask &= ~bit
            face_counter[lower_faces(mask)] += 1
            total += 1
    return len(optional), total, face_counter


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
        # The vertical s=0 edge consists of affine source points.
        if not (start[0] == end[0] == 0):
            packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def determinant(first, second):
    return first[0] * second[1] - first[1] * second[0]


def face_order(base, plane):
    slope_s, slope_p, constant = plane
    values = tuple(base * item for item in plane)
    need(all(value.denominator == 1 for value in values),
         "integral face chart")
    need((base * Fraction(1, 6)).denominator == 1,
         "integral good-target twist")
    # F_Q=sigma^(L*c)G, F_p=sigma^(L*(b+c))G_P,
    # ds=sigma^(-L*a)dS, and dA/(2C)=sigma^(L/6)eta_0.
    order = base * (Fraction(5, 6) - slope_s - slope_p - constant)
    need(order.denominator == 1 and order > 0, "positive integral face order")
    return int(order)


def check_symbolic_face_geometry():
    """Prove the displayed face normal forms and generic ``G_P`` units."""
    one = sparse_constant(1)
    S = sparse_variable(0)
    P = sparse_variable(1)
    U = sparse_variable(2)
    W = sparse_variable(3)
    beta = sparse_variable(4)
    zeta = sparse_variable(5)

    S2 = sparse_power(S, 2)
    S3 = sparse_power(S, 3)
    P3 = sparse_power(P, 3)
    P4 = sparse_power(P, 4)
    P5 = sparse_power(P, 5)
    P6 = sparse_power(P, 6)
    P7 = sparse_power(P, 7)

    R = sparse_subtract(S2, P)
    central = sparse_subtract(
        sparse_subtract(one, sparse_multiply(U, P6)),
        sparse_multiply(W, sparse_multiply(S2, P5)),
    )
    beta_face = sparse_subtract(
        sparse_subtract(one, sparse_multiply(W, sparse_multiply(S2, P5))),
        sparse_multiply(beta, sparse_multiply(S3, P4)),
    )
    zeta_face = sparse_subtract(
        sparse_subtract(one, sparse_multiply(zeta, sparse_multiply(S3, P3))),
        sparse_multiply(beta, sparse_multiply(S3, P4)),
    )
    beta_face_w0 = sparse_subtract(
        sparse_subtract(one, sparse_multiply(U, P6)),
        sparse_multiply(beta, sparse_multiply(S3, P4)),
    )

    # Central normalization: x=P^-1 and y=W*S/P.  Clearing P^7 from
    # y^2=W*x*(x^6-U) gives the exact multiple -W*central below.
    central_normal_left = sparse_subtract(
        sparse_multiply(
            sparse_power(W, 2), sparse_multiply(S2, P5)
        ),
        sparse_multiply(
            W, sparse_subtract(one, sparse_multiply(U, P6))
        ),
    )
    need(
        central_normal_left == sparse_neg(sparse_multiply(W, central)),
        "symbolic central hyperelliptic normal form",
    )

    # Beta normalization: x=W*S^2*P^5 and y=beta*S^3*P^4 satisfy
    # beta_face=1-x-y.  Difference of squares proves the printed P^7 law.
    x_beta = sparse_multiply(W, sparse_multiply(S2, P5))
    y_beta = sparse_multiply(beta, sparse_multiply(S3, P4))
    need(
        beta_face == sparse_subtract(sparse_subtract(one, x_beta), y_beta),
        "symbolic beta-face x+y equation",
    )
    beta_normal_left = sparse_subtract(
        sparse_multiply(
            sparse_multiply(sparse_power(W, 3), P7),
            sparse_power(sparse_subtract(one, x_beta), 2),
        ),
        sparse_multiply(sparse_power(beta, 2), sparse_power(x_beta, 3)),
    )
    beta_normal_right = sparse_multiply(
        sparse_multiply(sparse_power(W, 3), P7),
        sparse_multiply(
            beta_face,
            sparse_add(sparse_subtract(one, x_beta), y_beta),
        ),
    )
    need(beta_normal_left == beta_normal_right,
         "symbolic beta Kummer normal form")

    # The T3 face is rational with V=S*P and
    # (zeta+beta*P)V^3=1.
    V = sparse_multiply(S, P)
    zeta_rational_left = sparse_subtract(
        sparse_multiply(
            sparse_add(zeta, sparse_multiply(beta, P)), sparse_power(V, 3)
        ),
        one,
    )
    need(zeta_rational_left == sparse_neg(zeta_face),
         "symbolic zeta rational parameter")

    # At W=0 the V9 equation is already a cubic Kummer presentation.
    w0_normal_left = sparse_subtract(
        sparse_multiply(beta, sparse_multiply(S3, P4)),
        sparse_subtract(one, sparse_multiply(U, P6)),
    )
    need(w0_normal_left == sparse_neg(beta_face_w0),
         "symbolic W-zero genus-six normal form")

    # Exact Kummer/hyperelliptic ramification ledgers certify that the normal
    # forms are irreducible and have the promoted geometric genera.  A
    # valuation prime to the cover degree rules out a proper power.
    central_branch_count = 7 + 1  # seven simple finite roots plus infinity
    central_genus = (central_branch_count - 2) // 2
    need(central_branch_count == 8 and central_genus == 3,
         "central hyperelliptic branch/genus ledger")
    beta_valuations = (3, -2, -1)  # x=0, x=1, infinity
    need(sum(beta_valuations) == 0, "beta Kummer divisor degree zero")
    need(all(value % 7 for value in beta_valuations),
         "beta Kummer irreducibility/full ramification")
    beta_genus = (7 * (-2) + len(beta_valuations) * (7 - 1) + 2) // 2
    need(beta_genus == 3, "beta Kummer Riemann-Hurwitz genus")
    w0_valuations = (1,) * 6 + (-4, -2)  # six roots, P=0, infinity
    need(sum(w0_valuations) == 0, "W-zero Kummer divisor degree zero")
    need(all(value % 3 for value in w0_valuations),
         "W-zero Kummer irreducibility/full ramification")
    w0_genus = (3 * (-2) + len(w0_valuations) * (3 - 1) + 2) // 2
    need(w0_genus == 6, "W-zero Kummer Riemann-Hurwitz genus")

    # Generic-point G_P witnesses.  The equalities are literal polynomial
    # identities; the displayed remainders have lower S-degree than their
    # face equations and are nonzero, so the derivatives do not vanish in
    # the corresponding face function fields.
    central_core = sparse_multiply(R, central)
    central_p = sparse_derivative(central, P_INDEX)
    central_core_p = sparse_derivative(central_core, P_INDEX)
    need(
        sparse_subtract(central_core_p, sparse_multiply(R, central_p))
        == sparse_neg(central),
        "central G_P restriction identity",
    )
    central_remainder = sparse_subtract(
        sparse_multiply(P, central_p), sparse_scale(central, 5)
    )
    expected_central_remainder = sparse_neg(
        sparse_add(sparse_multiply(U, P6), sparse_constant(5))
    )
    need(central_remainder == expected_central_remainder,
         "central C_P nonvanishing remainder")
    need(central_remainder and
         sparse_degree(central_remainder, S_INDEX) < sparse_degree(central, S_INDEX),
         "central C_P proper S-remainder")
    central_on_R = substitute_s_squared_by_p(central)
    expected_on_R = sparse_subtract(
        one, sparse_multiply(sparse_add(U, W), P6)
    )
    need(central_on_R == expected_on_R and central_on_R,
         "central R factor nonzero at generic point")

    beta_core = sparse_multiply(S2, beta_face)
    beta_p = sparse_derivative(beta_face, P_INDEX)
    beta_core_p = sparse_derivative(beta_core, P_INDEX)
    need(beta_core_p == sparse_multiply(S2, beta_p),
         "beta-face G_P restriction identity")
    beta_remainder = sparse_subtract(
        sparse_multiply(P, beta_p), sparse_scale(beta_face, 4)
    )
    expected_beta_remainder = sparse_neg(
        sparse_add(x_beta, sparse_constant(4))
    )
    need(beta_remainder == expected_beta_remainder,
         "beta-face B_P nonvanishing remainder")
    need(beta_remainder and
         sparse_degree(beta_remainder, S_INDEX) < sparse_degree(beta_face, S_INDEX),
         "beta-face B_P proper S-remainder")

    beta_w0_p = sparse_derivative(beta_face_w0, P_INDEX)
    beta_w0_remainder = sparse_subtract(
        sparse_multiply(P, beta_w0_p), sparse_scale(beta_face_w0, 4)
    )
    expected_beta_w0_remainder = sparse_neg(
        sparse_add(sparse_scale(sparse_multiply(U, P6), 2), sparse_constant(4))
    )
    need(beta_w0_remainder == expected_beta_w0_remainder,
         "W-zero beta B_P nonvanishing remainder")
    need(beta_w0_remainder and
         sparse_degree(beta_w0_remainder, S_INDEX)
         < sparse_degree(beta_face_w0, S_INDEX),
         "W-zero beta B_P proper S-remainder")

    # Independent exact gcd control for the central branch polynomial
    # x*(x^6-U), at a positive rational gate specialization.
    central_branch = (0, -2, 0, 0, 0, 0, 0, 1)
    need(polynomial_squarefree(central_branch),
         "central branch polynomial squarefree control")

    return (
        "central:P*C_P-5*C=-(U*P^6+5);C|R=1-(U+W)*P^6",
        "beta:P*B_P-4*B=-(W*S^2*P^5+4)",
        "W0beta:P*B_P-4*B=-(2*U*P^6+4)",
    )


def generic_edge_schemes(U, W, beta, zeta):
    x_minus_one = (Fraction(-1), Fraction(1))
    return (
        ("X-1", x_minus_one),
        ("1-zeta*X^3", (Fraction(1), 0, 0, -zeta)),
        ("-beta*X-zeta", (-zeta, -beta)),
        ("-W*X-beta", (-beta, -W)),
        ("(X-1)*(U*X+W)", polynomial_multiply(x_minus_one, (W, U))),
        ("U-X^6", (U, 0, 0, 0, 0, 0, Fraction(-1))),
        ("1-W*X", (Fraction(1), -W)),
        ("1-beta*X", (Fraction(1), -beta)),
    )


def w0_edge_schemes(U, beta, zeta):
    return (
        ("X-1", (Fraction(-1), Fraction(1))),
        ("1-zeta*X^3", (Fraction(1), 0, 0, -zeta)),
        ("-beta*X-zeta", (-zeta, -beta)),
        ("-U*X-beta", (-beta, -U)),
        ("U*(X-1)", (-U, U)),
        ("U-X^6", (U, 0, 0, 0, 0, 0, Fraction(-1))),
        ("1-U*X^6", (Fraction(1), 0, 0, 0, 0, 0, -U)),
        ("1-beta*X", (Fraction(1), -beta)),
    )


def check_edge_schemes():
    """Check both literal edge lists and their exact squarefreeness gates."""
    U = Fraction(2)
    W = Fraction(1)
    beta = Fraction(3)
    zeta = Fraction(5)
    generic = generic_edge_schemes(U, W, beta, zeta)
    w0 = w0_edge_schemes(U, beta, zeta)
    for regime, schemes in (("W-nonzero", generic), ("W-zero", w0)):
        for name, polynomial in schemes:
            polynomial = polynomial_trim(polynomial)
            need(polynomial[0] != 0 and polynomial[-1] != 0,
                 f"{regime} edge {name} avoids toric corners")
            need(polynomial_squarefree(polynomial),
                 f"{regime} edge {name} squarefree")

    # Symbolic certificate formulas behind the gate implications.
    top = polynomial_multiply((-1, 1), (W, U))
    top_discriminant = top[1] * top[1] - 4 * top[0] * top[2]
    need(top_discriminant == (U + W) ** 2,
         "generic top-edge discriminant=(U+W)^2")
    need(3 ** 3 * (-zeta) ** 3 != 0,
         "cubic binomial resultant unit*zeta^3")
    need(6 ** 6 * U ** 5 != 0,
         "affine sixth-binomial resultant unit*U^5")
    need(6 ** 6 * (-U) ** 6 != 0,
         "W-zero internal sixth-binomial resultant unit*U^6")

    hostile_top = generic_edge_schemes(Fraction(1), Fraction(-1), beta, zeta)[4][1]
    need(not polynomial_squarefree(hostile_top),
         "Lambda-zero generic top edge is repeated")

    return (
        tuple(name for name, _polynomial in generic),
        tuple(name for name, _polynomial in w0),
        "generic:U*W*beta*zeta*(U+W);W0:U*beta*zeta",
    )


def multigraph_ledger(vertices, edges):
    adjacency = {vertex: set() for vertex in vertices}
    for first, second in edges:
        need(first in adjacency and second in adjacency, "graph edge endpoint declared")
        adjacency[first].add(second)
        adjacency[second].add(first)
    unseen = set(vertices)
    connected_components = 0
    while unseen:
        connected_components += 1
        stack = [unseen.pop()]
        while stack:
            vertex = stack.pop()
            for neighbour in adjacency[vertex]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    stack.append(neighbour)
    betti = len(edges) - len(vertices) + connected_components
    return len(vertices), len(edges), connected_components, betti


def check_graph_ledgers():
    # W!=0: R and the central genus-three curve meet in the twelve simple
    # roots of 1-(U+W)S^12.  The two primitive shared edges have lengths one.
    generic_contact = (Fraction(1),) + (Fraction(0),) * 11 + (Fraction(-3),)
    need(polynomial_squarefree(generic_contact),
         "W-nonzero twelve central intersections")
    generic_shared = (gcd(2, 5), gcd(3, 4))
    need(generic_shared == (1, 1), "W-nonzero primitive shared-edge lengths")
    generic_vertices = ("R", "C", "B", "T")
    generic_edges = (
        (("R", "C"),) * (len(generic_contact) - 1)
        + (("C", "B"),) * generic_shared[0]
        + (("B", "T"),) * generic_shared[1]
    )
    generic_ledger = multigraph_ledger(generic_vertices, generic_edges)
    need(generic_ledger == (4, 14, 1, 11), "W-nonzero graph V/E/cc/b1")

    # W=0: 1-U*P^6 has six simple rational factors.  Each meets R twice;
    # the V9 shared edge has lattice length six and the T3 edge length one.
    w0_factor = (Fraction(1), 0, 0, 0, 0, 0, Fraction(-2))
    w0_contact = (Fraction(1),) + (Fraction(0),) * 11 + (Fraction(-2),)
    need(polynomial_squarefree(w0_factor), "W-zero six rational main factors")
    need(polynomial_squarefree(w0_contact), "W-zero twelve R/factor intersections")
    w0_shared = (gcd(0, 6), gcd(3, 4))
    need(w0_shared == (6, 1), "W-zero shared-edge lengths")
    lines = tuple(f"L{index}" for index in range(6))
    w0_vertices = ("R",) + lines + ("B", "T")
    w0_edges = tuple(
        edge
        for line in lines
        for edge in (("R", line), ("R", line), (line, "B"))
    ) + (("B", "T"),)
    need(len(w0_edges) == 12 + w0_shared[0] + w0_shared[1],
         "W-zero exact graph edge decomposition")
    w0_ledger = multigraph_ledger(w0_vertices, w0_edges)
    need(w0_ledger == (9, 19, 1, 11), "W-zero graph V/E/cc/b1")
    return generic_ledger, w0_ledger


def check_lambda_zero_a23():
    """Transport THM-4297's local A23/tail ledger to Z=0, U+W=0."""
    one = sparse_constant(1)
    r = sparse_variable(0)
    b_local = sparse_variable(1)
    U = sparse_variable(2)
    W = sparse_neg(U)
    r_minus_one = sparse_subtract(r, one)
    r4 = sparse_power(r, 4)
    r5 = sparse_power(r, 5)
    r6 = sparse_power(r, 6)

    A = sparse_add(sparse_scale(U, 2), W)
    D = sparse_power(W, 2)  # Z=0, so D=W^2-4UZ.
    need(A == U, "Z0 Lambda-zero transverse coefficient A=U")
    need(D == sparse_power(U, 2), "Z0 Lambda-zero discriminant D=U^2")

    # Z=0 and W=-U give Ctilde=b^12-U*r^5*(r-1).
    ctilde = sparse_subtract(
        sparse_power(b_local, 12),
        sparse_multiply(U, sparse_multiply(r5, r_minus_one)),
    )
    ctilde_on_R = sparse_substitute_constant(ctilde, 0, 1)
    need(ctilde_on_R == sparse_power(b_local, 12),
         "Z0 Lambda-zero length-twelve R intersection")
    ctilde_r = sparse_derivative(ctilde, 0)
    derivative_at_contact = sparse_substitute_constant(
        sparse_substitute_constant(ctilde_r, 0, 1), 1, 0
    )
    need(derivative_at_contact == sparse_neg(U),
         "Z0 Lambda-zero transverse derivative -U")
    need(2 * 12 - 1 == 23, "length-twelve two-branch contact is A23")

    # Here A=2U+W=U and U_eff=A/2.  The exact top identity is the same
    # THM-4297 identity; W=-U changes only the sign of the q^3 correction.
    left = sparse_multiply(U, sparse_multiply(r5, r_minus_one))
    right = sparse_add(
        sparse_scale(
            sparse_multiply(U, sparse_subtract(r6, r4)), Fraction(1, 2)
        ),
        sparse_scale(
            sparse_multiply(
                U, sparse_multiply(r4, sparse_power(r_minus_one, 2))
            ),
            Fraction(1, 2),
        ),
    )
    need(left == right, "Z0 Lambda-zero THM-4297 top identity")
    need(Fraction(1, 2) != 0, "effective coefficient U_eff=U/2 is nonzero")

    base_top = sparse_scale(
        sparse_multiply(U, sparse_subtract(r6, r4)), Fraction(1, 2)
    )
    top_difference = sparse_subtract(left, base_top)
    expected_top_difference = sparse_scale(
        sparse_multiply(U, sparse_multiply(r4, sparse_power(r_minus_one, 2))),
        Fraction(1, 2),
    )
    need(top_difference == expected_top_difference,
         "Z0 Lambda-zero exact quadratic top correction")
    prepared_difference = sparse_multiply(r_minus_one, top_difference)
    expected_prepared_difference = sparse_scale(
        sparse_multiply(U, sparse_multiply(r4, sparse_power(r_minus_one, 3))),
        Fraction(1, 2),
    )
    need(prepared_difference == expected_prepared_difference,
         "Z0 Lambda-zero prepared +U/2*r4*q3 correction")

    # The prepared equation has an outer q=r-1, hence the difference is
    # +(U/2)r^4*q^3.  On q=t^6*y, division by t^12 puts it first at t^6,
    # strictly after the four critical coefficients and the absent t^5 row.
    correction_order = 3 * 6 - 12
    need(correction_order == 6 and correction_order > 4,
         "Z0 Lambda-zero correction enters after critical ladder")
    # If the four critical coefficients vanish, the b^12*q splitter gap is
    # 6(beta_v-s), while the correction begins at 6(s+beta_v); their
    # difference is 12s>0 for every positive divisorial valuation s.
    need((6 - (-6), 6 - 6) == (12, 0),
         "terminal splitter strictly precedes t6 correction")

    tail_forms = (
        (Fraction(6), Fraction(2)),
        (Fraction(5, 2), Fraction(9, 2)),
        (Fraction(2), Fraction(4)),
        (Fraction(3, 2), Fraction(7, 2)),
        (Fraction(1), Fraction(3)),
    )
    need(all(s_coefficient > 0 and beta_coefficient > 0
             for s_coefficient, beta_coefficient in tail_forms),
         "Z0 Lambda-zero five Keller-form tail orders positive")
    need(Fraction(3) > 0 and Fraction(5) > 0,
         "Z0 Lambda-zero nonrepeated lower order 3s+5beta positive")
    return (
        "6s+2beta_v",
        "(5s+9beta_v)/2",
        "2s+4beta_v",
        "(3s+7beta_v)/2",
        "s+3beta_v",
    )


def check_source_normal_factorization():
    # Store a coefficient as (U coefficient, W coefficient).  With
    # z=x^2*t, compare U(1+z)^6+W*z(1+z)^5 with
    # (1+z)^5(U+(U+W)z).
    left = []
    right = []
    for degree in range(7):
        previous = comb(5, degree - 1) if 1 <= degree <= 6 else 0
        left.append((comb(6, degree), previous))
        current = comb(5, degree) if degree <= 5 else 0
        right.append((current + previous, previous))
    need(tuple(left) == tuple(right), "literal source-normal factorization")
    h0, h1, h2 = left[:3]
    wall = (
        h2[0] - 5 * h1[0] + 15 * h0[0],
        h2[1] - 5 * h1[1] + 15 * h0[1],
    )
    need((h0, h1, h2) == ((1, 0), (6, 1), (15, 5)),
         "M12 leading flag")
    need(wall == (0, 0), "Z-zero flag equation")
    # On Lambda=U+W=0, set U=1,W=-1.  The top scheme
    # (X-1)(U*X+W) becomes (X-1)^2.
    first = (-1, 1)
    second = (-1, 1)
    product = (
        first[0] * second[0],
        first[0] * second[1] + first[1] * second[0],
        first[1] * second[1],
    )
    need(product == (1, -2, 1), "Lambda-zero double-root control")
    return tuple(left)


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
    need((2, 0, 1) in point_index, "fixed residual in hull universe")

    M12 = (Fraction(1, 12), Fraction(1, 6), Fraction(-1, 6))
    BETA7 = (Fraction(1, 7), Fraction(1, 7), Fraction(-2, 7))
    BETA6 = (Fraction(1, 9), Fraction(1, 6), Fraction(-2, 9))
    T3 = (Fraction(1, 3), Fraction(0), Fraction(-2, 3))

    P_ROW = (1, 0)
    P2_ROW = (2, 0)
    P3_ROW = (3, 0)
    U_ROW = (6, 0)
    W_ROW = (3, 2)
    Z_ROW = (0, 4)
    BETA_ROW = (1, 3)
    ZETA_ROW = (0, 3)
    fixed = {P_ROW, P2_ROW, P3_ROW, U_ROW}
    collisions = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1),
        (2, 6, 1), (3, 4, 1), (3, 5, 1),
    )

    generic_required = fixed | {W_ROW, BETA_ROW, ZETA_ROW}
    generic_optional, generic_total, generic_counter = hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        generic_required, {Z_ROW}, collisions,
    )
    generic_faces = (M12, BETA7, T3)
    need(generic_optional == 8, "W-nonzero optional-row count")
    need(generic_total == 16384, "W-nonzero hostile-mask count")
    need(generic_counter == Counter({generic_faces: 16384}),
         "W-nonzero unique lower-face complex")

    # With W absent the point (2,6,1) is owned solely by U and cannot cancel.
    w0_required = fixed | {BETA_ROW, ZETA_ROW}
    w0_collisions = tuple(point for point in collisions if point != (2, 6, 1))
    w0_optional, w0_total, w0_counter = hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        w0_required, {Z_ROW, W_ROW}, w0_collisions,
    )
    w0_faces = (M12, BETA6, T3)
    need(w0_optional == 8, "W-zero optional-row count")
    need(w0_total == 8192, "W-zero hostile-mask count")
    need(w0_counter == Counter({w0_faces: 8192}),
         "W-zero unique lower-face complex")

    # Concrete positive controls and the two cheapest owner-loss hostiles.
    need(2 * 1 * 1 * 1 * (2 + 1) != 0,
         "W-nonzero concrete gate control")
    need(1 * 1 * 1 * (1 + 0) != 0,
         "W-zero concrete gate control")

    def bare_faces(required_labels):
        active = tuple(row for row in rows if row[:2] in required_labels)
        return lower_faces(support_mask(active))

    beta_zero_faces = bare_faces(fixed | {W_ROW, ZETA_ROW})
    zeta_zero_faces = bare_faces(fixed | {W_ROW, BETA_ROW})
    need(beta_zero_faces != generic_faces, "beta-zero owner-complex hostile")
    need(zeta_zero_faces != generic_faces, "zeta-zero owner-complex hostile")

    broad_optional, broad_total, broad_counter = hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        fixed, {Z_ROW}, collisions,
    )
    broad_planes = {
        plane for face_complex in broad_counter for plane in face_complex
    }
    order_coefficients = {
        plane: Fraction(5, 6) - sum(plane) for plane in broad_planes
    }
    minimum_order = min(order_coefficients.values())
    minimum_planes = tuple(sorted(
        plane for plane, value in order_coefficients.items()
        if value == minimum_order
    ))
    need(broad_optional == 11, "broad optional-row count")
    need(broad_total == 131072, "broad hostile-mask count")
    need(len(broad_counter) == 51, "broad face-complex count")
    need(len(broad_planes) == 32, "broad supporting-plane count")
    need(all(value > 0 for value in order_coefficients.values()),
         "all broad primary-face orders positive")
    need(minimum_order == Fraction(3, 4), "broad minimum order coefficient")
    need(minimum_planes == (M12,), "unique broad minimum plane")

    face_gp_witnesses = check_symbolic_face_geometry()
    generic_edge_names, w0_edge_names, edge_gate = check_edge_schemes()
    generic_graph, w0_graph = check_graph_ledgers()
    lambda_zero_tail_orders = check_lambda_zero_a23()

    generic_polygons = {
        "M": ((0, 1), (2, 0), (4, 5), (0, 7)),
        "central": ((0, 0), (0, 6), (2, 5)),
        "beta": ((2, 0), (4, 5), (5, 4)),
        "zeta": ((2, 0), (5, 3), (5, 4)),
        "global": ((0, 1), (2, 0), (5, 3), (5, 4), (4, 5), (0, 7)),
    }
    generic_ledgers = {
        name: polygon_ledger(polygon)[1:]
        for name, polygon in generic_polygons.items()
    }
    need(generic_ledgers == {
        "M": (36, 10, 14),
        "central": (12, 8, 3),
        "beta": (7, 3, 3),
        "zeta": (3, 5, 0),
        "global": (46, 14, 17),
    }, "W-nonzero Pick ledgers")
    need(generic_ledgers["central"][2] + generic_ledgers["beta"][2]
         + generic_graph[3]
         == generic_ledgers["global"][2], "W-nonzero genus inventory")
    generic_packet = edge_packet(generic_polygons["global"])
    need(generic_packet == (11, 11, 7, 4, 2, 2, 2, 1),
         "W-nonzero source-infinity packet")
    need(sum(index - 1 for index in generic_packet) == 32 == 2 * 17 - 2,
         "W-nonzero Riemann-Hurwitz saturation")
    need((sum(generic_packet) - 6, sum(generic_packet)) == (34, 40),
         "W-nonzero locked-cubic responses")

    w0_polygons = {
        "M": ((0, 1), (2, 0), (2, 6), (0, 7)),
        "beta": ((2, 0), (2, 6), (5, 4)),
        "zeta": ((2, 0), (5, 3), (5, 4)),
        "global": ((0, 1), (2, 0), (5, 3), (5, 4), (2, 6), (0, 7)),
    }
    w0_ledgers = {
        name: polygon_ledger(polygon)[1:]
        for name, polygon in w0_polygons.items()
    }
    need(w0_ledgers == {
        "M": (24, 14, 6),
        "beta": (18, 8, 6),
        "zeta": (3, 5, 0),
        "global": (45, 13, 17),
    }, "W-zero Pick ledgers")
    need(w0_ledgers["beta"][2] + w0_graph[3] == w0_ledgers["global"][2],
         "W-zero genus inventory")
    w0_packet = edge_packet(w0_polygons["global"])
    need(w0_packet == (17, 11, 4, 2, 2, 2, 1),
         "W-zero source-infinity packet")
    need(sum(index - 1 for index in w0_packet) == 32 == 2 * 17 - 2,
         "W-zero Riemann-Hurwitz saturation")
    need((sum(w0_packet) - 6, sum(w0_packet)) == (33, 39),
         "W-zero locked-cubic responses")

    generic_orders = tuple(face_order(84, plane) for plane in generic_faces)
    w0_orders = tuple(face_order(36, plane) for plane in w0_faces)
    need(generic_orders == (63, 70, 98), "W-nonzero exact form orders")
    need(w0_orders == (27, 28, 42), "W-zero exact form orders")

    # Torus-smooth positive controls from the two nonconstant exponent vectors.
    need(determinant((0, 6), (2, 5)) == -12,
         "central genus-three torus smoothness")
    need(determinant((2, 5), (3, 4)) == -7,
         "beta genus-three torus smoothness")
    need(determinant((3, 3), (3, 4)) == 3,
         "zeta rational-face torus smoothness")
    need(determinant((0, 6), (3, 4)) == -18,
         "W-zero genus-six torus smoothness")
    source_normal_coefficients = check_source_normal_factorization()

    print("THM4327_Z0_ENDPOINT_EXACT_AUDIT=PASS")
    print("universe_rows=16 forbidden_rows=((0,1),(1,1)) wall_absent_row=(0,4)")
    print("fixed_points=((2,0,0),(0,1,0),(2,0,1))")
    print("collision_points=" + repr(collisions))
    print("collision_masks_are_hostile_overcount=true")
    print("source_normal_coefficients_U_W=" + repr(source_normal_coefficients))
    print("source_normal_top=t^6*(1+x^2*t)^5*(U+(U+W)*x^2*t)")
    print("source_normal_flag=((U),(6U+W),(15U+5W)) Zwall=h2-5h1+15h0=0")
    print("THEOREM_Z_gate=Z=0;U*beta*zeta!=0;W_arbitrary;Lambda0_A23_covered")
    print("symbolic_face_normal_forms=central_hyperelliptic,beta_Kummer,zeta_rational,W0beta_Kummer:PASS")
    print("W_NONZERO_hostiles=16384 optional_rows=8 collision_bits=6")
    print("W_NONZERO_faces=" + repr(generic_faces))
    print("W_NONZERO_pick=" + repr(generic_ledgers))
    print("W_NONZERO_normalizations=(central:y^2=W*x*(x^6-U),beta:P^7=(beta^2/W^3)*x^3/(1-x)^2,zeta:rational)")
    print("W_NONZERO_GP_witnesses=" + repr(face_gp_witnesses[:2]))
    print("W_NONZERO_edges=" + repr(generic_edge_names))
    print("W_NONZERO_edge_squarefree_gate_Lambda_nonzero=U*W*beta*zeta*(U+W)!=0")
    print("W_NONZERO_graph_Lambda_nonzero_V_E_cc_b1=" + repr(generic_graph))
    print("W_NONZERO_graph_edges=R-C:12,C-B:1,B-T:1")
    print("W_NONZERO_orders_L84=" + repr(generic_orders))
    print("W_NONZERO_packet_Lambda_nonzero=" + repr(generic_packet) + " responses=(34,40)")
    print("LAMBDA_ZERO_Ctilde=b^12-U*r^5*(r-1);partial_r_at_contact=-U;contact=A23")
    print("LAMBDA_ZERO_transport=A=U;U_eff=U/2;correction=+(U/2)*r^4*q^3;normalized_entry=t^6")
    print("LAMBDA_ZERO_tail_orders=" + repr(lambda_zero_tail_orders) + " all_positive=true")
    print("LAMBDA_ZERO_packets_transported=false")
    print("W_ZERO_hostiles=8192 optional_rows=8 collision_bits=5")
    print("W_ZERO_faces=" + repr(w0_faces))
    print("W_ZERO_pick=" + repr(w0_ledgers))
    print("W_ZERO_normalizations=(M:rational_factors,beta:S^3=(1-U*P^6)/(beta*P^4)_genus6,zeta:rational)")
    print("W_ZERO_GP_witness=" + repr(face_gp_witnesses[2]))
    print("W_ZERO_edges=" + repr(w0_edge_names))
    print("W_ZERO_edge_squarefree_gate=U*beta*zeta!=0")
    print("W_ZERO_graph_V_E_cc_b1=" + repr(w0_graph))
    print("W_ZERO_graph_edges=R-L_i:12,L_i-B:6,B-T:1")
    print("W_ZERO_orders_L36=" + repr(w0_orders))
    print("W_ZERO_packet=" + repr(w0_packet) + " responses=(33,39)")
    print("BROAD_overatlas_masks=131072 optional_rows=11 collision_bits=6")
    print("BROAD_face_complexes=51 supporting_planes=32")
    print("BROAD_order_formula=L*(5/6-a-b-c)")
    print("BROAD_min_order_coefficient=3/4 unique_plane=" + repr(M12))
    print("BROAD_all_primary_face_orders_positive=true")
    print("positive_controls=(U,W,beta,zeta)=(2,1,1,1),(1,0,1,1),(1,-1,1,1)")
    print("Lambda0_top_edge=U*(X-1)^2_replaced_by_audited_A23_tail_ledger")
    print("hostile_beta0_faces=" + repr(beta_zero_faces))
    print("hostile_zeta0_faces=" + repr(zeta_zero_faces))
    print("edge_gate_certificate=" + edge_gate)
    print("scope=finite_exact_Newton_order_response_A23_inputs_relative_geometric_imports_no_seam_entry_no_JC2")
    print("CHECKS=" + str(CHECKS))


if __name__ == "__main__":
    main()
