#!/usr/bin/env python3
"""Exact finite controls for the THM-3213 transcendental-radius corollary.

The analytic obstruction is proved on paper.  This script checks only its
finite combinatorial hypotheses, normalization maps, and required hostiles.
"""

from ast import Assert, parse, walk
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from math import comb, factorial


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def ordered_stirling_profile(n):
    row = [1]
    for m in range(1, n + 1):
        nxt = [0] * (m + 1)
        for k in range(1, m + 1):
            left = row[k] if k < len(row) else 0
            right = row[k - 1] if k - 1 < len(row) else 0
            nxt[k] = k * (left + right)
        row = nxt
    return row


@lru_cache(None)
def numerator_power(d):
    base = {
        (1, 0, 0): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (0, 1, 1): 1,
        (1, 1, 1): 3,
    }
    poly = {(0, 0, 0): 1}
    for _ in range(d):
        out = {}
        for a, va in poly.items():
            for b, vb in base.items():
                exponent = tuple(a[i] + b[i] for i in range(3))
                out[exponent] = out.get(exponent, 0) + va * vb
        poly = out
    return poly


def quotient_kernel(d, c1, c2, c3):
    answer = 0
    for q in range(min(c1, c2, c3) + 1):
        coeff = numerator_power(d).get((c1 - q, c2 - q, c3 - q), 0)
        answer += coeff * comb(d + q - 1, d - 1)
    return answer


def output_coordinate(r, d):
    profile = ordered_stirling_profile(r)
    answer = 0
    for c1 in range(1, r + 1):
        for c2 in range(1, r + 1):
            for c3 in range(1, r + 1):
                weight = quotient_kernel(d, c1, c2, c3)
                answer += weight * profile[c1] * profile[c2] * profile[c3]
    return answer


def permutation_sign(perm):
    inversions = sum(
        perm[i] > perm[j]
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
    )
    return -1 if inversions % 2 else 1


def det_i_minus_ax(adjacency):
    """Exact polynomial det(I-A diag(X)) as exponent-tuple -> coefficient."""
    q = len(adjacency)
    answer = {}
    for perm in permutations(range(q)):
        coeff = permutation_sign(perm)
        exponent = [0] * q
        for i, j in enumerate(perm):
            if i == j:
                continue
            if not adjacency[i][j]:
                coeff = 0
                break
            coeff *= -1
            exponent[j] += 1
        if coeff:
            key = tuple(exponent)
            answer[key] = answer.get(key, 0) + coeff
    return {key: value for key, value in answer.items() if value}


def uniform_polynomial(poly):
    answer = {}
    for exponent, coeff in poly.items():
        degree = sum(exponent)
        answer[degree] = answer.get(degree, 0) + coeff
    return {degree: coeff for degree, coeff in answer.items() if coeff}


def log_derivative_uniform(poly, variable):
    answer = {}
    for exponent, coeff in poly.items():
        if exponent[variable]:
            degree = sum(exponent)
            answer[degree] = answer.get(degree, 0) + coeff * exponent[variable]
    return {degree: coeff for degree, coeff in answer.items() if coeff}


def polynomial_value(poly, value):
    return sum(Fraction(coeff) * value**degree for degree, coeff in poly.items())


def polynomial_derivative_value(poly, value):
    return sum(
        Fraction(degree * coeff) * value ** (degree - 1)
        for degree, coeff in poly.items()
        if degree
    )


def require_tournament_and_strong(adjacency, name):
    q = len(adjacency)
    require(all(adjacency[i][i] == 0 for i in range(q)),
            ("loop", name, adjacency))
    require(all(adjacency[i][j] + adjacency[j][i] == 1
                for i in range(q) for j in range(i + 1, q)),
            ("not tournament", name, adjacency))
    reach = [[bool(adjacency[i][j]) or i == j for j in range(q)] for i in range(q)]
    for k in range(q):
        for i in range(q):
            for j in range(q):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    require(all(reach[i][j] for i in range(q) for j in range(q)),
            ("not strong", name, adjacency))


records = []
positive = zero = 0
for r in range(1, 13):
    scale = factorial(r) ** 3
    for d in range(1, 13):
        raw = output_coordinate(r, d)
        expected_positive = d <= 3 * r
        require((raw > 0) == expected_positive,
                ("Hamilton-cut positivity", r, d, raw))
        positive += raw > 0
        zero += raw == 0
        normalized = Fraction(raw, scale)
        require(normalized * scale == raw,
                ("factorial normalization", r, d, normalized, raw))
        if r >= d:
            twice = normalized / comb(r - 1, d - 1)
            require(twice * comb(r - 1, d - 1) == normalized,
                    ("binomial normalization", r, d, twice, normalized))
        records.append((r, d, raw, normalized))


def c3_hamilton(profile):
    profile = profile + (0,)
    return 3 * sum(
        profile[c] ** 3
        + profile[c + 1] * profile[c] ** 2
        + profile[c + 1] ** 2 * profile[c]
        for c in range(1, len(profile) - 1)
    )


same_h_profiles = (
    (0, 15, 78, 198, 240, 120),
    (0, 15, 90, 210, 240, 120),
)
require(same_h_profiles[0][-1] == same_h_profiles[1][-1] == 120,
        ("same-H input", same_h_profiles))
hostile_outputs = tuple(c3_hamilton(row) for row in same_h_profiles)
require(hostile_outputs == (178036299, 193215375),
        ("different-jet hostile", hostile_outputs))


# Regular C5: exact Perron balance and the W_Q resolvent orientation.
q = 5
c5 = tuple(
    tuple(int((j - i) % q in (1, 2)) for j in range(q))
    for i in range(q)
)
require_tournament_and_strong(c5, "C5")
c5_rows = tuple(sum(row) for row in c5)
c5_columns = tuple(sum(c5[i][j] for i in range(q)) for j in range(q))
require(c5_rows == (2,) * q and c5_columns == (2,) * q,
        ("C5 Perron balance", c5_rows, c5_columns))
c5_det = uniform_polynomial(det_i_minus_ax(c5))
c5_alpha = Fraction(1, 2)
require(polynomial_value(c5_det, c5_alpha) == 0,
        ("C5 Perron pole", c5_det, c5_alpha))
require(polynomial_derivative_value(c5_det, c5_alpha) != 0,
        ("C5 nonsimple Perron pole", c5_det, c5_alpha))
for t in (Fraction(1, 10), Fraction(1, 3)):
    lhs = tuple(1 - t * sum(row) for row in c5)
    require(lhs == (1 - 2 * t,) * q,
            ("C5 resolvent orientation", t, lhs))
    require(Fraction(q) * t / (1 - 2 * t) > 0,
            ("C5 resolvent numerator", t))


# Strong Q4 hostile: its equal positive pole is not diagonal-critical.
q4 = (
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (1, 0, 0, 1),
    (1, 1, 0, 0),
)
require_tournament_and_strong(q4, "Q4")
q4_det = det_i_minus_ax(q4)
expected_q4_det = {
    (0, 0, 0, 0): 1,
    (1, 1, 1, 0): -1,
    (0, 1, 1, 1): -1,
    (1, 1, 1, 1): -1,
}
require(q4_det == expected_q4_det, ("Q4 determinant", q4_det))
q4_uniform = uniform_polynomial(q4_det)
q4_log_gradients = tuple(log_derivative_uniform(q4_det, i) for i in range(4))
expected_q4_gradients = (
    {3: -1, 4: -1},
    {3: -2, 4: -1},
    {3: -2, 4: -1},
    {3: -1, 4: -1},
)
require(q4_uniform == {0: 1, 3: -2, 4: -1},
        ("Q4 uniform determinant", q4_uniform))
require(q4_log_gradients == expected_q4_gradients,
        ("Q4 unequal critical gradients", q4_log_gradients))


primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
for n in range(0, 50):
    left = (n + 1) * Fraction(1, factorial(n + 1))
    right = Fraction(1, factorial(n))
    require(left == right, ("P-recursive hostile recurrence", n, left, right))
for p in primes:
    denominator = Fraction(1, factorial(p)).denominator
    require(denominator % p == 0,
            ("fresh denominator hostile", p, denominator))


tree = parse(open(__file__, encoding="utf-8").read(), filename=__file__)
require(not any(isinstance(node, Assert) for node in walk(tree)),
        ("assert node", __file__))

digest = sha256("\n".join(map(str, records)).encode()).hexdigest()
print("TOURNAMENT P-RECURSIVE RADIUS EXACT CONTROLS")
print(f"grid=r_1_12,d_1_12;positive={positive};zero={zero};record_sha256={digest}")
print(f"same_H=120;different_C3_lifts={hostile_outputs}")
print(f"normalizations=factorial_cube_and_fixed_d_binomial;checked_records={len(records)}")
print(f"prime_denominator_P_recursive_hostile=1/n!;primes={primes}")
print(f"regular_C5=row_col_sum_2;alpha=1/2;uniform_det={c5_det};W_equal=5t/(1-2t)")
print(f"strong_Q4_uniform_det={q4_uniform};log_gradients={q4_log_gradients}")
print("analytic_dependency=THM3213_asymptotic_and_transcendental_radius_lemma")
print("FAILED_CHECKS=NONE")
