#!/usr/bin/env python3
"""Exact algebraic referee for THM-2065.

THM-2051 supplies the analytic bounded-relation alternative.  This script
audits the subsequent integer-linear pullback: internal template circuits,
primitive orthogonal rays, fixed-N uniqueness, and finite circuit-ray
intersection with a hereditary-primitivity filter.

Tournament Analysis is not useful here.  Relation patterns are symmetric
hyperedges and their signed coefficients are load-bearing.  An orientation
would forget the linear form whose projective kernel is the theorem.
"""

from itertools import combinations, product
from math import comb, gcd
from random import Random


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def gcd_many(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def dot(left, right):
    return left[0] * right[0] + left[1] * right[1]


def relation_row(rows, coefficients):
    return (
        sum(k * row[0] for k, row in zip(coefficients, rows)),
        sum(k * row[1] for k, row in zip(coefficients, rows)),
    )


def canonical_direction(vector):
    require(vector != (0, 0), "zero vector has no projective direction")
    divisor = gcd(abs(vector[0]), abs(vector[1]))
    answer = (vector[0] // divisor, vector[1] // divisor)
    if answer[0] < 0 or (answer[0] == 0 and answer[1] < 0):
        answer = (-answer[0], -answer[1])
    return answer


def orthogonal_direction(row):
    require(row != (0, 0), "zero relation row is internal")
    return canonical_direction((row[1], -row[0]))


def coefficient_patterns(n, height):
    values = tuple(range(-height, 0)) + tuple(range(1, height + 1))
    answer = []
    for size in range(3, min(5, n) + 1):
        for support in combinations(range(n), size):
            for nonzero in product(values, repeat=size):
                if nonzero[0] < 0 or gcd_many(nonzero) != 1:
                    continue
                coefficients = [0] * n
                for index, value in zip(support, nonzero):
                    coefficients[index] = value
                answer.append(tuple(coefficients))
    return answer


def fixed_n_candidate(row, n_value):
    """Return the unique integral M solving A*N+B*M=0, if it exists."""
    a_value, b_value = row
    require(row != (0, 0) and n_value > 0, "invalid fixed-N problem")
    if b_value == 0 or (-a_value * n_value) % b_value:
        return None
    return (-a_value * n_value) // b_value


def determinant(left, right):
    return left[0] * right[1] - left[1] * right[0]


def hereditary_primitive(rows, parameter):
    return all(
        gcd_many(
            dot(row, parameter)
            for position, row in enumerate(rows)
            if position != deleted
        )
        == 1
        for deleted in range(len(rows))
    )


# Every nonzero integer covector has exactly one primitive projective kernel.
kernel_checks = 0
primitive_box = [
    (n_value, m_value)
    for n_value in range(-10, 11)
    for m_value in range(-10, 11)
    if (n_value, m_value) != (0, 0) and gcd(abs(n_value), abs(m_value)) == 1
]
for a_value in range(-20, 21):
    for b_value in range(-20, 21):
        row = (a_value, b_value)
        if row == (0, 0):
            continue
        expected = orthogonal_direction(row)
        for parameter in primitive_box:
            require(
                (dot(row, parameter) == 0)
                == (canonical_direction(parameter) == expected),
                ("projective kernel mismatch", row, parameter, expected),
            )
            kernel_checks += 1


# Fixed N cuts a nonzero relation row at at most one longitudinal integer.
fixed_n_checks = 0
for a_value in range(-12, 13):
    for b_value in range(-12, 13):
        row = (a_value, b_value)
        if row == (0, 0):
            continue
        for n_value in range(1, 9):
            candidate = fixed_n_candidate(row, n_value)
            direct = [
                m_value
                for m_value in range(-100, 101)
                if dot(row, (n_value, m_value)) == 0
            ]
            if candidate is None:
                require(not direct, ("missed fixed-N solution", row, n_value))
            elif -100 <= candidate <= 100:
                require(direct == [candidate], ("fixed-N uniqueness", row, n_value))
            else:
                require(not direct, ("out-of-box fixed-N solution", row, n_value))
            fixed_n_checks += 1


# Exhaust the complete height-two packet on deterministic five-row templates.
patterns = coefficient_patterns(5, 2)
rng = Random(2065)
template_checks = 0
internal_templates = 0
circuit_free_templates = 0
wheel_intersections = 0
for _ in range(50):
    rows = []
    while len(rows) < 5:
        row = (rng.randint(-60, 60), rng.randint(-60, 60))
        if row != (0, 0):
            rows.append(row)

    pulled = {relation_row(rows, coefficients) for coefficients in patterns}
    has_internal = (0, 0) in pulled
    if has_internal:
        internal_templates += 1
    else:
        circuit_free_templates += 1
    rays = {orthogonal_direction(row) for row in pulled if row != (0, 0)}

    for parameter in primitive_box:
        scalar_relation = any(dot(row, parameter) == 0 for row in pulled)
        predicted = has_internal or canonical_direction(parameter) in rays
        require(
            scalar_relation == predicted,
            ("finite packet pullback mismatch", rows, parameter),
        )
        template_checks += 1

    # A direct hereditary filter can be intersected after the ray packet with
    # no loss.  This is the set-theoretic content of the THM-2062 composition.
    for n_value in range(1, 9):
        direct = {
            m_value
            for m_value in range(-30, 31)
            if gcd(n_value, abs(m_value)) == 1
            and hereditary_primitive(rows, (n_value, m_value))
            and any(dot(row, (n_value, m_value)) == 0 for row in pulled)
        }
        candidates = {
            candidate
            for row in pulled
            if row != (0, 0)
            for candidate in [fixed_n_candidate(row, n_value)]
            if candidate is not None
            and -30 <= candidate <= 30
            and gcd(n_value, abs(candidate)) == 1
            and hereditary_primitive(rows, (n_value, candidate))
        }
        if has_internal:
            candidates.update(
                m_value
                for m_value in range(-30, 31)
                if gcd(n_value, abs(m_value)) == 1
                and hereditary_primitive(rows, (n_value, m_value))
            )
        require(direct == candidates, ("wheel intersection mismatch", rows, n_value))
        wheel_intersections += 1


# The internal-circuit guardrail really leaves a full two-dimensional plane.
guardrail_rows = [(1, 0), (0, 1), (1, 1), (2, 3), (5, 8)]
guardrail_coefficients = (1, 1, -1, 0, 0)
require(
    relation_row(guardrail_rows, guardrail_coefficients) == (0, 0),
    "internal circuit guardrail failed",
)
guardrail_points = 0
for parameter in primitive_box:
    speeds = [dot(row, parameter) for row in guardrail_rows]
    require(
        sum(k * speed for k, speed in zip(guardrail_coefficients, speeds)) == 0,
        ("internal relation did not persist", parameter),
    )
    guardrail_points += 1


# The three-row circuit coefficients are the signed 2x2 minors.
minor_checks = 0
for _ in range(10000):
    triple = [
        (rng.randint(-100, 100), rng.randint(-100, 100)) for _ in range(3)
    ]
    coefficients = (
        determinant(triple[1], triple[2]),
        determinant(triple[2], triple[0]),
        determinant(triple[0], triple[1]),
    )
    require(
        relation_row(triple, coefficients) == (0, 0),
        ("signed-minor circuit failed", triple, coefficients),
    )
    minor_checks += 1


height = 2**20
crude_pattern_bound = sum(comb(13, size) * (2 * height) ** size for size in range(3, 6))
canonical_ray_bound = crude_pattern_bound // 2

print("THM-2065 exact two-anchor circuit-ray referee")
print(f"THM-2051 height: {height}")
print(f"crude 13-row relation-pattern bound: {crude_pattern_bound}")
print(f"canonical projective-ray bound: {canonical_ray_bound}")
print(f"projective kernel checks: {kernel_checks}")
print(f"fixed-N uniqueness checks: {fixed_n_checks}")
print(f"height-two five-row patterns: {len(patterns)}")
print(
    "five-row templates: "
    f"{internal_templates} internal, {circuit_free_templates} circuit-free"
)
print(f"finite packet/parameter checks: {template_checks}")
print(f"hereditary-wheel intersections: {wheel_intersections}")
print(f"internal-circuit guardrail points: {guardrail_points}")
print(f"signed-minor circuit checks: {minor_checks}")
print("PASS")
