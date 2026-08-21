#!/usr/bin/env python3
"""Exact companion for THM-3642.

This file verifies the two universal fourth-stable endpoint identities,
their actual local target-ring controls, and the deterministic Q6 actual
J0/J1/J2 pivot certificate.  The frozen theorem package passed independent
hostile audit before promotion.
"""

import ast
import hashlib
from math import factorial
from pathlib import Path

import sympy as sp
from flint import fmpq, fmpq_mat, nmod_mat, nmod_poly


CHECKS = 0
MODULUS = 1000003


def require(label, condition):
    """Record one exact gate which remains active under python -O."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def as_fmpq(value):
    """Convert a SymPy rational to python-flint fmpq exactly."""
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def rational_pair(value):
    """Return an integer numerator/denominator pair."""
    if isinstance(value, fmpq):
        return int(value.numerator), int(value.denominator)
    value = sp.Rational(value)
    return int(value.p), int(value.q)


def rational_vector_hash(values, separator=";"):
    """SHA-256 of a complete ordered rational coefficient vector."""
    payload = separator.join(
        f"{numerator}/{denominator}"
        for numerator, denominator in map(rational_pair, values)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def integer_vector_hash(values, separator=";"):
    """SHA-256 of a complete ordered integer index vector."""
    payload = separator.join(str(int(value)) for value in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def target_coefficient_hash(values, packet):
    """Bind every nonzero rational coefficient to its target monomial."""
    entries = []
    for value, (metadata, _restriction, _delta) in zip(values, packet):
        if not value:
            continue
        numerator, denominator = rational_pair(value)
        entries.append(
            f"{metadata[0]},{metadata[1]},{metadata[2]}:"
            f"{numerator}/{denominator}"
        )
    return hashlib.sha256(";".join(entries).encode("ascii")).hexdigest()


def polynomial_hash(polynomial):
    """Hash all descending rational coefficients, including internal zeros."""
    polynomial = sp.Poly(polynomial, x, domain=sp.QQ)
    return rational_vector_hash(polynomial.all_coeffs())


def modular_rational(value):
    """Reduce a rational modulo the fixed audit prime."""
    numerator, denominator = rational_pair(value)
    return numerator * pow(denominator, -1, MODULUS) % MODULUS


def to_nmod(polynomial):
    """Convert a rational univariate polynomial to nmod_poly."""
    polynomial = sp.Poly(polynomial, x, domain=sp.QQ)
    if polynomial.is_zero:
        return nmod_poly([], MODULUS)
    coefficients = [0] * (polynomial.degree() + 1)
    for (degree,), coefficient in polynomial.terms():
        coefficients[degree] = modular_rational(coefficient)
    return nmod_poly(coefficients, MODULUS)


def truncated(expr, point, x, t, source_order=2, stable_order=4):
    """Bivariate Taylor dictionary at (x,t)=(point,0)."""
    answer = {}
    for source_degree in range(source_order + 1):
        for stable_degree in range(stable_order + 1):
            value = (
                sp.diff(expr, x, source_degree, t, stable_degree)
                .subs({x: point, t: 0})
                / factorial(source_degree)
                / factorial(stable_degree)
            )
            if value:
                answer[(source_degree, stable_degree)] = sp.factor(value)
    return answer


def truncated_product(left, right, source_order=2, stable_order=4):
    """Multiply two truncated bivariate Taylor dictionaries."""
    answer = {}
    for (left_source, left_stable), left_value in left.items():
        for (right_source, right_stable), right_value in right.items():
            source_degree = left_source + right_source
            stable_degree = left_stable + right_stable
            if source_degree > source_order or stable_degree > stable_order:
                continue
            key = (source_degree, stable_degree)
            answer[key] = answer.get(key, 0) + left_value * right_value
    return {key: sp.factor(value) for key, value in answer.items() if value}


def truncated_power(value, exponent):
    """Power in the fixed source-order-two/stable-order-four jet ring."""
    answer = {(0, 0): sp.S.One}
    for _ in range(exponent):
        answer = truncated_product(answer, value)
    return answer


print("THM-3642 exact companion -- zero-debt fourth closure")
print("status=PROVED VERIFIED-EXACT INDEPENDENTLY HOSTILE-AUDITED")


print("SECTION compiler and zero-second-debt folds")
x, q, t = sp.symbols("x q t")
y, z, w = sp.symbols("y z w")
points = (-1, 0, 1)

Q1 = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
Q6 = sp.expand(Q1 - sp.Rational(259, 36) * x**2 * (x**2 - 1) ** 2)
Q7 = (
    -x**7
    - sp.Rational(27, 4) * x**6
    + 3 * x**5
    + 18 * x**4
    - 3 * x**3
    - sp.Rational(27, 2) * x**2
    + x
    - sp.Rational(3, 4)
)

D_general = 1 + x**2 * q
b_general = sp.expand((D_general - 1) * (D_general + 2) ** 2)
c_general = sp.expand(x * D_general * (D_general + 2))
e_general = sp.expand(q * (D_general + 3))
y_general = c_general / 3
z_general = e_general + 3

require(
    "compiler surface relation",
    sp.expand(c_general**2 * e_general - b_general * (b_general + 4)) == 0,
)


def fold_packet(label, polynomial, expected_curvature):
    """Verify the common ordinary triple and the Q-curvature row."""
    values = tuple(sp.factor(polynomial.subs(x, point)) for point in points)
    slopes = tuple(sp.factor(sp.diff(polynomial, x).subs(x, point)) for point in points)
    curvatures = tuple(
        sp.factor(sp.diff(polynomial, x, 2).subs(x, point)) for point in points
    )
    require(f"{label} values", values == (-3, sp.Rational(-3, 4), -3))
    require(
        f"{label} slopes",
        slopes == (sp.Rational(-9, 2), 1, sp.Rational(9, 2)),
    )
    require(f"{label} curvatures", curvatures == expected_curvature)
    require(
        f"{label} zero second debt row",
        5 * curvatures[0] + 13 * curvatures[2] + 243 == 0,
    )

    restricted_y = sp.expand(y_general.subs(q, polynomial))
    restricted_z = sp.expand(z_general.subs(q, polynomial))
    target_values = tuple(
        (
            sp.factor(restricted_y.subs(x, point)),
            sp.factor(restricted_z.subs(x, point)),
        )
        for point in points
    )
    tangent_y = sp.Matrix(
        [sp.diff(restricted_y, x).subs(x, point) for point in points]
    )
    tangent_z = sp.Matrix(
        [sp.diff(restricted_z, x).subs(x, point) for point in points]
    )
    require(f"{label} target triple", target_values == ((0, 0),) * 3)
    require(f"{label} tangent y", tangent_y == sp.ones(3, 1))
    require(f"{label} tangent z", tangent_z == sp.Matrix([-9, 4, 9]))
    require(
        f"{label} ordinary triple",
        tuple(
            sp.det(
                sp.Matrix(
                    [
                        [tangent_y[i], tangent_z[i]],
                        [tangent_y[j], tangent_z[j]],
                    ]
                )
            )
            for i, j in ((0, 1), (1, 2), (0, 2))
        )
        == (13, 5, 18),
    )
    return restricted_y, restricted_z


fold_packet(
    "Q6",
    Q6,
    (sp.Rational(-451, 18), sp.Rational(-251, 9), sp.Rational(-163, 18)),
)
fold_packet(
    "Q7",
    Q7,
    (sp.Rational(-27, 2), -27, sp.Rational(-27, 2)),
)
print("PASS Q6_Q7=same_ordinary_triple_and_zero_second_debt")


print("SECTION universal fourth-stable two-form identities")
lambda_coefficients = (
    sp.Rational(5, 18),
    -1,
    sp.Rational(13, 18),
)


Q6_IDENTITY = {
    (0, -1, 0): sp.Rational(16246280, 531441),
    (0, -1, 1): -sp.Rational(4489, 6561),
    (0, -1, 2): -sp.Rational(10, 81),
    (0, 0, 1): -sp.Rational(64, 81),
    (0, 1, 0): sp.Rational(13390648, 531441),
    (0, 1, 1): -sp.Rational(6559, 6561),
    (0, 1, 2): -sp.Rational(26, 81),
    (2, -1, 0): sp.Rational(2012, 2187),
    (2, -1, 1): sp.Rational(5, 27),
    (2, 1, 0): -sp.Rational(2012, 2187),
    (2, 1, 1): -sp.Rational(13, 27),
}

Q7_IDENTITY = {
    (0, -1, 0): sp.Rational(2300, 81),
    (0, -1, 1): -sp.Rational(1, 9),
    (0, -1, 2): -sp.Rational(10, 81),
    (0, 1, 0): sp.Rational(3140, 81),
    (0, 1, 1): -sp.Rational(7, 9),
    (0, 1, 2): -sp.Rational(26, 81),
    (2, -1, 0): sp.Rational(4, 9),
    (2, -1, 1): sp.Rational(5, 27),
    (2, 1, 0): -sp.Rational(4, 9),
    (2, 1, 1): -sp.Rational(13, 27),
}


def verify_two_form_identity(label, polynomial, identity, expected_debt):
    """Verify the identity on every relevant target two-form monomial.

    Write omega=A dy^dz+B dy^dw+C dz^dw.  No closedness or decomposability
    equation is injected: the gate therefore covers arbitrary target
    two-forms and, a fortiori, every dF^dG from an actual target-ring pair.
    """
    pulled_y = sp.expand(y_general.subs(q, polynomial + t**2))
    pulled_z = sp.expand(z_general.subs(q, polynomial + t**2))
    pulled_y_x = sp.diff(pulled_y, x)
    pulled_y_t = sp.diff(pulled_y, t)
    pulled_z_x = sp.diff(pulled_z, x)
    pulled_z_t = sp.diff(pulled_z, t)
    area = sp.expand(pulled_y_x * pulled_z_t - pulled_y_t * pulled_z_x)

    packets = {}
    for point in points:
        y_jet = truncated(pulled_y, point, x, t)
        z_jet = truncated(pulled_z, point, x, t)
        bases = (
            truncated(area, point, x, t),
            truncated(pulled_y_x, point, x, t),
            truncated(pulled_z_x, point, x, t),
        )
        packets[point] = (
            [truncated_power(y_jet, exponent) for exponent in range(5)],
            [truncated_power(z_jet, exponent) for exponent in range(5)],
            bases,
        )

    def monomial_pullback(kind, y_degree, z_degree, w_degree, point):
        y_powers, z_powers, bases = packets[point]
        value = truncated_product(y_powers[y_degree], z_powers[z_degree])
        value = truncated_product(value, {(0, w_degree): sp.S.One})
        return truncated_product(value, bases[kind])

    def residual(kind, y_degree, z_degree, w_degree, coefficients):
        values = {
            point: monomial_pullback(
                kind, y_degree, z_degree, w_degree, point
            )
            for point in points
        }
        left = sum(
            lambda_coefficients[index] * values[point].get((0, 4), 0)
            for index, point in enumerate(points)
        )
        right = sum(
            coefficient * values[point].get((source_degree, stable_degree), 0)
            for (stable_degree, point, source_degree), coefficient in coefficients.items()
        )
        return sp.factor(left - right)

    monomial_count = 0
    mutation_detected = False
    mutated = dict(identity)
    mutated[(0, -1, 0)] += 1
    for kind in range(3):
        for y_degree in range(5):
            for z_degree in range(5 - y_degree):
                for w_degree in range(5 - y_degree - z_degree):
                    require(
                        f"{label} two-form monomial {kind,y_degree,z_degree,w_degree}",
                        residual(kind, y_degree, z_degree, w_degree, identity) == 0,
                    )
                    if residual(kind, y_degree, z_degree, w_degree, mutated) != 0:
                        mutation_detected = True
                    monomial_count += 1
    require(f"{label} active coefficient mutation", mutation_detected)

    forced_debt = sp.factor(
        sum(
            coefficient
            for (stable_degree, _point, source_degree), coefficient in identity.items()
            if stable_degree == 0 and source_degree == 0
        )
    )
    require(f"{label} forced debt", forced_debt == expected_debt)
    print(
        f"PASS {label}_two_form_monomials={monomial_count} "
        f"forced_lambda_J4={forced_debt}"
    )


verify_two_form_identity("Q6", Q6, Q6_IDENTITY, sp.Rational(365888, 6561))
verify_two_form_identity("Q7", Q7, Q7_IDENTITY, sp.Rational(5440, 81))
print("PASS identities_require_neither_J1_nor_J3_and_no_rank_division")


print("SECTION actual local target-ring controls")
F0_Q6 = (
    y
    + sp.Rational(81, 65) * y**2
    - sp.Rational(371, 972) * y * z
    + sp.Rational(277, 15795) * z**2
    - sp.Rational(561271753, 118272960) * y**3
    - sp.Rational(1067, 9477) * y**2 * z
    - sp.Rational(10544747, 23654592) * y * z**2
)
F2_Q6 = (
    -sp.Rational(13, 243) * y
    + sp.Rational(4024, 15795) * z
    + sp.Rational(227269286, 3080025) * y**2
    - sp.Rational(2564065739, 399171240) * y * z
    - sp.Rational(180336166, 249482025) * z**2
)
F0_Q7 = (
    y
    + sp.Rational(81, 65) * y**2
    - sp.Rational(1, 4) * y * z
    - sp.Rational(1, 65) * z**2
    - sp.Rational(84419, 20800) * y**3
    - sp.Rational(2741, 1755) * y**2 * z
    - sp.Rational(302927, 561600) * y * z**2
)
F2_Q7 = (
    y
    + sp.Rational(8, 65) * z
    + sp.Rational(7406968, 114075) * y**2
    - sp.Rational(469273, 70200) * y * z
    - sp.Rational(11632, 12675) * z**2
)


def verify_local_control(label, polynomial, F0, F2, expected_j4, expected_debt):
    """Check F=F0+w^2 F2, G=w in the actual local target ring."""
    general_F0 = sp.expand(F0.subs({y: y_general, z: z_general}))
    general_F2 = sp.expand(F2.subs({y: y_general, z: z_general}))
    U = sp.expand(general_F0.subs(q, polynomial))
    delta_F0 = sp.expand(sp.diff(general_F0, q).subs(q, polynomial))
    delta2_F0 = sp.expand(sp.diff(general_F0, q, 2).subs(q, polynomial))
    restricted_F2 = sp.expand(general_F2.subs(q, polynomial))
    delta_F2 = sp.expand(sp.diff(general_F2, q).subs(q, polynomial))
    C = sp.expand(delta_F0 + restricted_F2)
    K = sp.expand(sp.Rational(1, 2) * delta2_F0 + delta_F2)

    require(
        f"{label} J0 retained two-jet",
        tuple(
            tuple(sp.factor(sp.diff(U, x, order).subs(x, point)) for order in (1, 2, 3))
            for point in points
        )
        == ((1, 0, 0),) * 3,
    )
    require(
        f"{label} J2 retained one-jet",
        tuple(
            tuple(sp.factor(sp.diff(C, x, order).subs(x, point)) for order in (1, 2))
            for point in points
        )
        == ((0, 0),) * 3,
    )
    j4_values = tuple(sp.factor(sp.diff(K, x).subs(x, point)) for point in points)
    require(f"{label} J4 vector", j4_values == expected_j4)
    debt = sp.factor(
        sum(
            lambda_coefficients[index] * j4_values[index]
            for index in range(3)
        )
    )
    require(f"{label} local debt", debt == expected_debt)

    # The pullback of F0+w^2F2 is even in t and G=w, so the Jacobian is
    # d_x(F0+w^2F2), also even.  Verify the stable coefficients directly.
    pulled_f = sp.expand(
        general_F0.subs(q, polynomial + t**2)
        + t**2 * general_F2.subs(q, polynomial + t**2)
    )
    pulled_jacobian = sp.Poly(sp.diff(pulled_f, x), t)
    require(f"{label} J1 parity", pulled_jacobian.nth(1) == 0)
    require(f"{label} J3 parity", pulled_jacobian.nth(3) == 0)
    require(f"{label} J0 formula", pulled_jacobian.nth(0) == sp.diff(U, x))
    require(f"{label} J2 formula", pulled_jacobian.nth(2) == sp.diff(C, x))
    require(f"{label} J4 formula", pulled_jacobian.nth(4) == sp.diff(K, x))
    print(f"PASS {label}_actual_local_control_J1=J3=0 lambda_J4={debt}")


verify_local_control(
    "Q6",
    Q6,
    F0_Q6,
    F2_Q6,
    (
        sp.Rational(97158947953, 997928100),
        -sp.Rational(55841390863, 997928100),
        -sp.Rational(37631662223, 997928100),
    ),
    sp.Rational(365888, 6561),
)
verify_local_control(
    "Q7",
    Q7,
    F0_Q7,
    F2_Q7,
    (
        sp.Rational(7993879, 91260),
        -sp.Rational(29539213, 456300),
        -sp.Rational(13841293, 456300),
    ),
    sp.Rational(5440, 81),
)
print("PASS local_controls_are_finite_polynomials_in_y=c/3_z=e+3_w")


print("SECTION Q6 ACTUAL LIFT")

restricted_b = sp.Poly(sp.expand(b_general.subs(q, Q6)), x, domain=sp.QQ)
restricted_c = sp.Poly(sp.expand(c_general.subs(q, Q6)), x, domain=sp.QQ)
restricted_e = sp.Poly(sp.expand(e_general.subs(q, Q6)), x, domain=sp.QQ)
delta_b = sp.Poly(sp.diff(b_general, q).subs(q, Q6), x, domain=sp.QQ)
delta_c = sp.Poly(sp.diff(c_general, q).subs(q, Q6), x, domain=sp.QQ)
delta_e = sp.Poly(sp.diff(e_general, q).subs(q, Q6), x, domain=sp.QQ)
restriction_generators = (restricted_b, restricted_c, restricted_e)
delta_generators = (delta_b, delta_c, delta_e)
restriction_degrees = tuple(poly.degree() for poly in restriction_generators)
require("Q6 restriction degrees", restriction_degrees == (24, 17, 14))


def exact_monomials(cutoff):
    """Raw target monomials ordered by nested b,c,e exponents.

    The cutoff is the weighted restriction degree using (24,17,14), not
    target total degree.  Return metadata, restriction, and the restriction
    of the chosen representative's q derivative.
    """
    powers = []
    for generator, degree in zip(restriction_generators, restriction_degrees):
        row = [sp.Poly(1, x, domain=sp.QQ)]
        for _ in range(cutoff // degree):
            row.append(row[-1] * generator)
        powers.append(row)

    answer = []
    for b_exponent in range(cutoff // restriction_degrees[0] + 1):
        for c_exponent in range(cutoff // restriction_degrees[1] + 1):
            for e_exponent in range(cutoff // restriction_degrees[2] + 1):
                metadata = (b_exponent, c_exponent, e_exponent)
                weighted_degree = sum(
                    exponent * degree
                    for exponent, degree in zip(metadata, restriction_degrees)
                )
                if weighted_degree > cutoff:
                    continue
                monomial = (
                    powers[0][b_exponent]
                    * powers[1][c_exponent]
                    * powers[2][e_exponent]
                )
                delta_monomial = sp.Poly(0, x, domain=sp.QQ)
                if b_exponent:
                    delta_monomial += (
                        b_exponent
                        * powers[0][b_exponent - 1]
                        * powers[1][c_exponent]
                        * powers[2][e_exponent]
                        * delta_b
                    )
                if c_exponent:
                    delta_monomial += (
                        c_exponent
                        * powers[0][b_exponent]
                        * powers[1][c_exponent - 1]
                        * powers[2][e_exponent]
                        * delta_c
                    )
                if e_exponent:
                    delta_monomial += (
                        e_exponent
                        * powers[0][b_exponent]
                        * powers[1][c_exponent]
                        * powers[2][e_exponent - 1]
                        * delta_e
                    )
                answer.append((metadata, monomial, delta_monomial))
    return answer


def modular_monomials(cutoff):
    """Modular copy of the exact raw-monomial packet."""
    return [
        (metadata, to_nmod(polynomial), to_nmod(delta_polynomial))
        for metadata, polynomial, delta_polynomial in exact_monomials(cutoff)
    ]


def exact_rref_solution(columns, right_hand_side):
    """Set free RREF variables to zero in an exact rational system."""
    row_count = max(
        [column.degree() for column in columns] + [right_hand_side.degree()]
    ) + 1
    data = []
    for row in range(row_count):
        data.extend(as_fmpq(column.nth(row)) for column in columns)
        data.append(as_fmpq(right_hand_side.nth(row)))
    matrix = fmpq_mat(row_count, len(columns) + 1, data)
    reduced, rank = matrix.rref()
    solution = [fmpq(0)] * len(columns)
    for row in range(rank):
        pivot = next(
            (
                column
                for column in range(len(columns) + 1)
                if reduced[row, column]
            ),
            None,
        )
        if pivot == len(columns):
            return None, (row_count, len(columns), rank)
        if pivot is not None:
            solution[pivot] = reduced[row, len(columns)]
    return solution, (row_count, len(columns), rank)


base_packet = exact_monomials(126)
base_polynomials = [item[1] for item in base_packet]
require("Q6 cutoff126 raw monomial count", len(base_packet) == 105)
j0_columns = [restricted_c.diff() * polynomial for polynomial in base_polynomials]
j0_columns += [-restricted_e.diff() * polynomial for polynomial in base_polynomials]
j0_solution, j0_shape = exact_rref_solution(
    j0_columns, sp.Poly(1, x, domain=sp.QQ)
)
require("Q6 J0 exact solution exists", j0_solution is not None)
require("Q6 J0 exact shape rank", j0_shape == (143, 210, 142))

F1_restriction = sp.Poly(0, x, domain=sp.QQ)
G1_restriction = sp.Poly(0, x, domain=sp.QQ)
delta_F1 = sp.Poly(0, x, domain=sp.QQ)
delta_G1 = sp.Poly(0, x, domain=sp.QQ)
G1_target_vector = list(j0_solution[: len(base_packet)])
F1_target_vector = list(j0_solution[len(base_packet) :])
for index, (_metadata, monomial, delta_monomial) in enumerate(base_packet):
    g1_coefficient = G1_target_vector[index]
    if g1_coefficient:
        numerator, denominator = rational_pair(g1_coefficient)
        coefficient = sp.Rational(numerator, denominator)
        G1_restriction += coefficient * monomial
        delta_G1 += coefficient * delta_monomial
    f1_coefficient = F1_target_vector[index]
    if f1_coefficient:
        numerator, denominator = rational_pair(f1_coefficient)
        coefficient = sp.Rational(numerator, denominator)
        F1_restriction += coefficient * monomial
        delta_F1 += coefficient * delta_monomial

require(
    "Q6 exact J0 identity",
    restricted_c.diff() * G1_restriction
    - F1_restriction * restricted_e.diff()
    == sp.Poly(1, x, domain=sp.QQ),
)
require(
    "Q6 F1 retained values",
    tuple(F1_restriction.eval(point) for point in points) == (0, 0, 0),
)
require(
    "Q6 G1 retained values",
    tuple(G1_restriction.eval(point) for point in points)
    == (sp.Rational(1, 3),) * 3,
)

F1_target_hash = target_coefficient_hash(F1_target_vector, base_packet)
G1_target_hash = target_coefficient_hash(G1_target_vector, base_packet)
require(
    "Q6 F1 target representative hash",
    F1_target_hash
    == "88c69cb90ca318c7f84e58e83bc3e66104c9c56963818f7cc56f870b3ae4da77",
)
require(
    "Q6 G1 target representative hash",
    G1_target_hash
    == "a21b37f603b738ac89c70d0fee4ce1c73c117f6c1c3081fe16b0031302857267",
)
require("Q6 F1 target support", sum(bool(value) for value in F1_target_vector) == 71)
require("Q6 G1 target support", sum(bool(value) for value in G1_target_vector) == 68)
require(
    "Q6 F1 restriction hash",
    polynomial_hash(F1_restriction)
    == "648b6d09c9ada7160f4cf6d8637d3ce6cbd9d27a26f87f894da3b8a5a8bbbd75",
)
require(
    "Q6 G1 restriction hash",
    polynomial_hash(G1_restriction)
    == "e9363ac4075038f87f82daa5b4b83e2f9e4da57d4650177ee86b83a2e1934587",
)
require(
    "Q6 delta F1 representative hash",
    polynomial_hash(delta_F1)
    == "8e44affbd6dab1d95644d15f86437f3b1ebef0a6ce33b0cae9942ba96b01eb22",
)
require(
    "Q6 delta G1 representative hash",
    polynomial_hash(delta_G1)
    == "2fe85085fea62bee81ac10a0eb88908060d6d33a6495d152469d84d77d405dc6",
)
print(
    "PASS Q6_J0_exact_shape=143x210_rank142 "
    f"F1_target_support={sum(bool(value) for value in F1_target_vector)} "
    f"F1_target_hash={F1_target_hash} "
    f"G1_target_support={sum(bool(value) for value in G1_target_vector)} "
    f"G1_target_hash={G1_target_hash}"
)
print(
    "PASS Q6_J0_restrictions "
    f"F1_degree={F1_restriction.degree()} F1_hash={polynomial_hash(F1_restriction)} "
    f"G1_degree={G1_restriction.degree()} G1_hash={polynomial_hash(G1_restriction)}"
)
print(
    "PASS Q6_J0_representative_deltas "
    f"deltaF1_degree={delta_F1.degree()} deltaF1_hash={polynomial_hash(delta_F1)} "
    f"deltaG1_degree={delta_G1.degree()} deltaG1_hash={polynomial_hash(delta_G1)}"
)

# A full modular solve chooses a deterministic small characteristic-zero
# certificate.  The modular ranks are not asserted as rational ranks.
modular_base_packet = modular_monomials(126)
modular_F1 = to_nmod(F1_restriction)
modular_G1 = to_nmod(G1_restriction)
modular_delta_F1 = to_nmod(delta_F1)
modular_delta_G1 = to_nmod(delta_G1)
modular_c = to_nmod(restricted_c)
modular_e = to_nmod(restricted_e)
modular_delta_c = to_nmod(delta_c)
modular_delta_e = to_nmod(delta_e)

lift_packet = modular_monomials(240)
lift_polynomials = [item[1] for item in lift_packet]
require("Q6 cutoff240 raw monomial count", len(lift_packet) == 561)
coupled_columns = []
for polynomial in lift_polynomials:
    coupled_columns.append(
        (
            -2 * polynomial * modular_e.derivative(),
            polynomial.derivative() * modular_G1
            - 2 * polynomial * modular_G1.derivative(),
        )
    )
for polynomial in lift_polynomials:
    coupled_columns.append(
        (
            2 * modular_c.derivative() * polynomial,
            2 * modular_F1.derivative() * polynomial
            - modular_F1 * polynomial.derivative(),
        )
    )
for polynomial in lift_polynomials:
    coupled_columns.append(
        (nmod_poly([], MODULUS), -3 * polynomial * modular_e.derivative())
    )
for polynomial in lift_polynomials:
    coupled_columns.append(
        (nmod_poly([], MODULUS), 3 * modular_c.derivative() * polynomial)
    )

known_j1_modular = (
    2 * modular_c.derivative() * modular_delta_e
    + modular_F1.derivative() * modular_G1
    - modular_F1 * modular_G1.derivative()
    - 2 * modular_delta_c * modular_e.derivative()
)
known_j2_modular = (
    3 * modular_c.derivative() * modular_delta_G1
    + 2 * modular_F1.derivative() * modular_delta_e
    + modular_delta_c.derivative() * modular_G1
    - modular_F1 * modular_delta_e.derivative()
    - 2 * modular_delta_c * modular_G1.derivative()
    - 3 * modular_delta_F1 * modular_e.derivative()
)
rhs_j1_modular = -known_j1_modular
rhs_j2_modular = -known_j2_modular
row_count_j1 = max(
    [known_j1_modular.degree()] + [pair[0].degree() for pair in coupled_columns]
) + 1
row_count_j2 = max(
    [known_j2_modular.degree()] + [pair[1].degree() for pair in coupled_columns]
) + 1
require("Q6 modular J1 row count", row_count_j1 == 257)
require("Q6 modular J2 row count", row_count_j2 == 366)

operator_data = []
augmented_data = []
for row in range(row_count_j1):
    row_values = [int(pair[0][row]) for pair in coupled_columns]
    operator_data.extend(row_values)
    augmented_data.extend(row_values + [int(rhs_j1_modular[row])])
for row in range(row_count_j2):
    row_values = [int(pair[1][row]) for pair in coupled_columns]
    operator_data.extend(row_values)
    augmented_data.extend(row_values + [int(rhs_j2_modular[row])])

operator_modular = nmod_mat(
    row_count_j1 + row_count_j2,
    len(coupled_columns),
    operator_data,
    MODULUS,
)
augmented_modular = nmod_mat(
    row_count_j1 + row_count_j2,
    len(coupled_columns) + 1,
    augmented_data,
    MODULUS,
)
reduced_augmented, augmented_rank = augmented_modular.rref()
operator_rank = operator_modular.rank()
require("Q6 modular operator rank", operator_rank == 618)
require("Q6 modular augmented rank", augmented_rank == 618)

pivot_columns = []
for row in range(augmented_rank):
    pivot = next(
        (
            column
            for column in range(len(coupled_columns) + 1)
            if reduced_augmented[row, column]
        ),
        None,
    )
    require(f"Q6 modular pivot is operator column row={row}", pivot != len(coupled_columns))
    pivot_columns.append(pivot)
require("Q6 modular pivot column count", len(pivot_columns) == 618)

selected_modular_data = []
for row in range(row_count_j1 + row_count_j2):
    selected_modular_data.extend(
        int(operator_modular[row, column]) for column in pivot_columns
    )
selected_modular = nmod_mat(
    row_count_j1 + row_count_j2,
    len(pivot_columns),
    selected_modular_data,
    MODULUS,
)
transpose_reduced, transpose_rank = selected_modular.transpose().rref()
require("Q6 selected modular column independence", transpose_rank == 618)
pivot_rows = []
for row in range(transpose_rank):
    pivot_rows.append(
        next(
            column
            for column in range(row_count_j1 + row_count_j2)
            if transpose_reduced[row, column]
        )
    )
require("Q6 modular pivot row count", len(pivot_rows) == 618)

pivot_column_hash = integer_vector_hash(pivot_columns, separator=",")
pivot_row_hash = integer_vector_hash(pivot_rows, separator=",")
require(
    "Q6 modular pivot column hash",
    pivot_column_hash
    == "38236094240c4954f5bb640c2bc12b942ea3a3e3c524e39c11223c9223d91b1b",
)
require(
    "Q6 modular pivot row hash",
    pivot_row_hash
    == "46fcca02583f1d60d52a2a9f615d28a1c33a5f0fe0d8f658b4fae24c1b04492d",
)
print(
    f"PASS Q6_mod_{MODULUS}_operator=623x2244_rank618_augmented_rank618 "
    f"pivot_column_hash={pivot_column_hash} pivot_row_hash={pivot_row_hash}"
)

# Rebuild only the selected pivot skeleton over Q.  This exact square solve,
# followed by the full Q[x] residual identities below, is the characteristic-
# zero proof.  We make no claim about the rank of the full rational operator.
exact_lift_packet = exact_monomials(240)
exact_lift_polynomials = [item[1] for item in exact_lift_packet]


def exact_coupled_column(column):
    """Return the selected rational (J1,J2) operator column."""
    block_size = len(exact_lift_polynomials)
    block, index = divmod(column, block_size)
    polynomial = exact_lift_polynomials[index]
    zero = sp.Poly(0, x, domain=sp.QQ)
    if block == 0:
        return (
            -2 * polynomial * restricted_e.diff(),
            polynomial.diff() * G1_restriction
            - 2 * polynomial * G1_restriction.diff(),
        )
    if block == 1:
        return (
            2 * restricted_c.diff() * polynomial,
            2 * F1_restriction.diff() * polynomial
            - F1_restriction * polynomial.diff(),
        )
    if block == 2:
        return zero, -3 * polynomial * restricted_e.diff()
    if block == 3:
        return zero, 3 * restricted_c.diff() * polynomial
    raise RuntimeError(f"bad coupled block {block}")


selected_exact_columns = [exact_coupled_column(column) for column in pivot_columns]
known_j1_exact = (
    2 * restricted_c.diff() * delta_e
    + F1_restriction.diff() * G1_restriction
    - F1_restriction * G1_restriction.diff()
    - 2 * delta_c * restricted_e.diff()
)
known_j2_exact = (
    3 * restricted_c.diff() * delta_G1
    + 2 * F1_restriction.diff() * delta_e
    + delta_c.diff() * G1_restriction
    - F1_restriction * delta_e.diff()
    - 2 * delta_c * G1_restriction.diff()
    - 3 * delta_F1 * restricted_e.diff()
)

exact_rhs = []
for row in pivot_rows:
    if row < row_count_j1:
        exact_rhs.append(as_fmpq(-known_j1_exact.nth(row)))
    else:
        exact_rhs.append(as_fmpq(-known_j2_exact.nth(row - row_count_j1)))

square_data = []
for row in pivot_rows:
    if row < row_count_j1:
        square_data.extend(
            as_fmpq(pair[0].nth(row)) for pair in selected_exact_columns
        )
    else:
        degree = row - row_count_j1
        square_data.extend(
            as_fmpq(pair[1].nth(degree)) for pair in selected_exact_columns
        )

exact_square = fmpq_mat(618, 618, square_data)
exact_rhs_matrix = fmpq_mat(618, 1, exact_rhs)
exact_solution = exact_square.solve(exact_rhs_matrix)
require(
    "Q6 exact selected square residual",
    exact_square * exact_solution == exact_rhs_matrix,
)

target_vectors = [[fmpq(0)] * len(exact_lift_polynomials) for _ in range(4)]
restrictions = [sp.Poly(0, x, domain=sp.QQ) for _ in range(4)]
for local_index, column in enumerate(pivot_columns):
    coefficient = exact_solution[local_index, 0]
    if not coefficient:
        continue
    block, index = divmod(column, len(exact_lift_polynomials))
    target_vectors[block][index] = coefficient
    numerator, denominator = rational_pair(coefficient)
    restrictions[block] += (
        sp.Rational(numerator, denominator) * exact_lift_polynomials[index]
    )

F2_restriction, G2_restriction, F3_restriction, G3_restriction = restrictions
F2_target_vector, G2_target_vector, F3_target_vector, G3_target_vector = target_vectors

j1_exact = (
    known_j1_exact
    - 2 * F2_restriction * restricted_e.diff()
    + 2 * restricted_c.diff() * G2_restriction
)
j2_exact = (
    known_j2_exact
    + F2_restriction.diff() * G1_restriction
    - 2 * F2_restriction * G1_restriction.diff()
    + 2 * F1_restriction.diff() * G2_restriction
    - F1_restriction * G2_restriction.diff()
    - 3 * F3_restriction * restricted_e.diff()
    + 3 * restricted_c.diff() * G3_restriction
)
require("Q6 full exact J1 identity", j1_exact.is_zero)
require("Q6 full exact J2 identity", j2_exact.is_zero)

expected_lift_data = {
    "F2": (
        185,
        "efa9ce6798214abbb07584ef133184b85b999f7cf5752113952b81023ee0c953",
        235,
        "b9166bc8740ca5e0263f8aa72c422f26b9aeb411571249ba7f24cec8711e7d09",
    ),
    "G2": (
        185,
        "4e69d81eee145f85ca5ef0b8d6e7ce13ea808089ab0732d964f6b5c52662b449",
        232,
        "c85105a341edfacfe7d36adf342a0c0ffce66e8ae2dbf41cc42e6bccb948f347",
    ),
    "F3": (
        185,
        "1a2ee587687ebccd1e6fccaa05e159235b16f020a0a3c353c3a541e3a193b756",
        240,
        "f8ee4cfd4b8bbc0d5a918aecff353b972379bba53b6e37348a6dab17d49c495c",
    ),
    "G3": (
        63,
        "6751af82e78165b96cf26f00b6f70cd5654b403137268892eabaaad01ca8a267",
        238,
        "64a5648703793adc0ef2a29d1eb8630ee442750220f4d26042bd4483e8b3fb11",
    ),
}
for label, vector, restriction in (
    ("F2", F2_target_vector, F2_restriction),
    ("G2", G2_target_vector, G2_restriction),
    ("F3", F3_target_vector, F3_restriction),
    ("G3", G3_target_vector, G3_restriction),
):
    support = sum(bool(value) for value in vector)
    target_hash = target_coefficient_hash(vector, exact_lift_packet)
    restriction_hash = polynomial_hash(restriction)
    expected_support, expected_target_hash, expected_degree, expected_hash = (
        expected_lift_data[label]
    )
    require(f"Q6 {label} target support", support == expected_support)
    require(f"Q6 {label} target hash", target_hash == expected_target_hash)
    require(f"Q6 {label} restriction degree", restriction.degree() == expected_degree)
    require(f"Q6 {label} restriction hash", restriction_hash == expected_hash)
    print(
        f"PASS Q6_{label}_target_support={support} "
        f"target_hash={target_hash} "
        f"restriction_degree={restriction.degree()} "
        f"restriction_hash={restriction_hash}"
    )
print("PASS Q6_exact_selected_square=618x618 full_QQ_J1=J2=0")


print("SECTION source AST gate")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")
print("PASS scope=two_fixed_folds_actual_local_controls_no_global_pair_JC2_OPEN")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
