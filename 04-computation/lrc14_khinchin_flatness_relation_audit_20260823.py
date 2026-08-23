#!/usr/bin/env python3
"""Exact audit of the Khinchin-flatness relation reduction for LRC(14).

The cited geometric input is the lonely-runner zonotope

    Z(n) = [1/(k+1), k/(k+1)]^k | n^perp

with projected lattice Lambda = Z^k | n^perp.  A counterexample makes Z(n)
Lambda-free.  Khinchin's Flatness Theorem then supplies a nonzero dual-lattice
direction of width at most Flt(k-1).

This script verifies the exact quotient-lattice algebra, the width formula,
and the explicit dimension-12 numerical consequence.  It does not prove the
cited flatness theorem and does not promote reserved THM-3743.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb, gcd, isqrt


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def dot(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> Fraction:
    require(len(left) == len(right), "dot-product dimensions")
    return sum((a * b for a, b in zip(left, right)), Fraction(0, 1))


def primitive_pair_relation(
    speeds: tuple[int, ...], first: int, second: int
) -> tuple[int, ...]:
    common_divisor = gcd(speeds[first], speeds[second])
    relation = [0] * len(speeds)
    relation[first] = speeds[second] // common_divisor
    relation[second] = -speeds[first] // common_divisor
    return tuple(relation)


def triple_relation(
    speeds: tuple[int, ...], first: int, second: int, third: int
) -> tuple[int, ...]:
    first_speed = speeds[first]
    second_speed = speeds[second]
    third_speed = speeds[third]
    relation = [0] * len(speeds)
    relation[first] = second_speed - third_speed
    relation[second] = third_speed - first_speed
    relation[third] = first_speed - second_speed
    common_divisor = 0
    for coefficient in relation:
        common_divisor = gcd(common_divisor, abs(coefficient))
    require(common_divisor > 0, "distinct speeds give a nonzero triple relation")
    return tuple(coefficient // common_divisor for coefficient in relation)


def projected_vertex(
    vertex: tuple[Fraction, ...], speeds: tuple[int, ...]
) -> tuple[Fraction, ...]:
    rational_speeds = tuple(Fraction(speed, 1) for speed in speeds)
    numerator = dot(vertex, rational_speeds)
    denominator = dot(rational_speeds, rational_speeds)
    return tuple(
        coordinate - numerator * speed / denominator
        for coordinate, speed in zip(vertex, rational_speeds)
    )


def exact_projected_width(
    speeds: tuple[int, ...], relation: tuple[int, ...]
) -> tuple[Fraction, int]:
    runner_count = len(speeds)
    lower = Fraction(1, runner_count + 1)
    upper = Fraction(runner_count, runner_count + 1)
    rational_relation = tuple(Fraction(coefficient, 1) for coefficient in relation)
    rational_speeds = tuple(Fraction(speed, 1) for speed in speeds)
    require(dot(rational_relation, rational_speeds) == 0, "dual vector is a speed relation")

    values: list[Fraction] = []
    projection_checks = 0
    for endpoint_choices in product((lower, upper), repeat=runner_count):
        vertex = tuple(endpoint_choices)
        projected = projected_vertex(vertex, speeds)
        require(dot(rational_relation, projected) == dot(rational_relation, vertex), "projection preserves dual pairing")
        values.append(dot(rational_relation, projected))
        projection_checks += 1

    measured_width = max(values) - min(values)
    predicted_width = (upper - lower) * sum(abs(coefficient) for coefficient in relation)
    require(measured_width == predicted_width, "projected-box width formula")
    return measured_width, projection_checks


checks = 0

print("=== Quotient-lattice dual and width audit ===")
speed_families: dict[int, tuple[tuple[int, ...], ...]] = {}
for runner_count in range(2, 7):
    speed_families[runner_count] = (
        tuple(range(1, runner_count + 1)),
        tuple(2 * index + 1 for index in range(runner_count)),
        tuple(range(1, runner_count)) + (runner_count + 2,),
    )

relation_count = 0
projection_count = 0
for runner_count, families in speed_families.items():
    for speeds in families:
        require(gcd(*speeds) == 1, "test speed vector is primitive")
        relations: list[tuple[int, ...]] = []
        for first in range(runner_count):
            for second in range(first + 1, runner_count):
                relations.append(primitive_pair_relation(speeds, first, second))
                for third in range(second + 1, runner_count):
                    relations.append(triple_relation(speeds, first, second, third))
        for relation in relations:
            width, local_projection_count = exact_projected_width(speeds, relation)
            side_length = Fraction(runner_count - 1, runner_count + 1)
            require(width == side_length * sum(abs(coefficient) for coefficient in relation), "exact width replay")
            relation_count += 1
            projection_count += local_projection_count
            checks += local_projection_count + 3

print(f"Primitive speed vectors checked: {sum(len(families) for families in speed_families.values())}.")
print(f"Pair/triple dual relations checked: {relation_count}.")
print(f"Projected cube vertices checked: {projection_count}.")
print("For every audited a in Z^k intersect n^perp:")
print("  <a, projection(x)> = <a, x>,")
print("  width_a(Z(n)) = ((k-1)/(k+1)) ||a||_1.")


print("\n=== Explicit LRC(14) consequence ===")
dimension = 12
flatness_bound_squared = (
    Fraction((dimension + 1) * (2 * dimension + 1), 6) * dimension**3
)
require(flatness_bound_squared == 93600, "dimension-12 explicit flatness square")
relation_bound_squared = Fraction(7, 6) ** 2 * flatness_bound_squared
require(relation_bound_squared == 127400, "dimension-12 relation square")
integer_cap = isqrt(relation_bound_squared.numerator // relation_bound_squared.denominator)
require(integer_cap == 356, "integer relation cap")
require(integer_cap**2 < relation_bound_squared < (integer_cap + 1) ** 2, "sharp integer rounding")
checks += 4

print("Using Flt(d) <= sqrt((d+1)(2d+1)/6) d^(3/2):")
print("  Flt(12) <= sqrt(93600) = 60 sqrt(26).")
print("The LRC(14) zonotope side length is 12/14 = 6/7, so a counterexample would force")
print("  (6/7)||a||_1 <= 60 sqrt(26), hence ||a||_1 <= 70 sqrt(26) < 357.")
print("Exact integer consequence: some nonzero a in Z^13 with a dot n = 0 has ||a||_1 <= 356.")


print("\n=== Minimal-width Graver split ===")
reduced_pair_ratio_count = sum(
    1
    for first in range(1, integer_cap)
    for second in range(first + 1, integer_cap + 1 - first)
    if gcd(first, second) == 1
)
require(reduced_pair_ratio_count == 19314, "bounded reduced pair-ratio atlas")

triple_speeds = (3, 4, 5)
triple_minimum = None
triple_minimizers: list[tuple[int, int, int]] = []
for relation in product(range(-4, 5), repeat=3):
    if relation == (0, 0, 0) or sum(
        coefficient * speed for coefficient, speed in zip(relation, triple_speeds)
    ):
        continue
    relation_norm = sum(abs(coefficient) for coefficient in relation)
    if triple_minimum is None or relation_norm < triple_minimum:
        triple_minimum = relation_norm
        triple_minimizers = [relation]
    elif relation_norm == triple_minimum:
        triple_minimizers.append(relation)
    checks += 1
require(triple_minimum == 4, "genuine triple shortest relation")
require((1, -2, 1) in triple_minimizers, "genuine triple witness")
require(
    all(sum(coefficient != 0 for coefficient in relation) == 3 for relation in triple_minimizers),
    "all shortest 3,4,5 relations have genuine support three",
)
checks += 3

print("Choose an l1-minimal nonzero speed relation; it is conformally indecomposable (a Graver element).")
print("If it has support two, its primitive coefficients are the reduced speed pair and their sum is <=356.")
print(f"Unordered distinct coprime pair ratios with numerator+denominator<=356: {reduced_pair_ratio_count}.")
print("The other branch is genuinely higher: for speeds (3,4,5), the shortest relation is")
print("  (1,-2,1), of l1 norm 4, while every pair relation is longer.")


print("\n=== Join with the THM-2052 rank-eleven atlas ===")
relation_candidate_count = sum(
    2**support * comb(13, support) * comb(integer_cap, support)
    for support in range(1, 14)
)
require(
    relation_candidate_count == 1978967793896659449022201064,
    "short ambient relation count",
)
relation_height = 91**6
mixed_cofactor_square_bound = integer_cap**2 * 3**11 * relation_height**22
mixed_cofactor_cap = isqrt(mixed_cofactor_square_bound)
require(
    mixed_cofactor_cap
    == 296721347184071259951513500572385227832063299530012809732642802018052529730741866571256016315918699080945038248181476620474893833141605,
    "mixed Hadamard cofactor cap",
)
require(
    mixed_cofactor_cap**2 <= mixed_cofactor_square_bound < (mixed_cofactor_cap + 1) ** 2,
    "mixed cofactor exact floor",
)
checks += 3

print(f"All nonzero ambient vectors with l1 norm <=356: {relation_candidate_count} (a finite upper universe).")
print("Let W be THM-2052's support-at-most-three, height-91^6 relation span.")
print("If dim(W)=11 and the flatness relation a is outside W, then rank becomes 12.")
print("Hadamard on eleven sparse rows plus a gives the primitive speed cap")
print(f"  max(n_i) <= {mixed_cofactor_cap}.")
print("If a lies in W, the unresolved rank-eleven code itself contains a nonzero l1-at-most-356 vector.")


print("\n=== Safe hostile control: arithmetic progression ===")
ap_speeds = tuple(range(1, 14))
ap_relation = (1, -2, 1) + (0,) * 10
ap_width, ap_projection_checks = exact_projected_width(ap_speeds, ap_relation)
require(ap_width == Fraction(24, 7), "AP short-relation width")
lonely_time = Fraction(1, 14)
minimum_distance = min(
    min((lonely_time * speed) % 1, 1 - ((lonely_time * speed) % 1))
    for speed in ap_speeds
)
require(minimum_distance == Fraction(1, 14), "AP boundary lonely witness")
checks += ap_projection_checks + 2

print("Speeds (1,2,...,13) have relation (1,-2,1,0,...,0), ||a||_1=4, width=24/7.")
print("They are nevertheless a valid boundary instance at t=1/14.")
print("Therefore a short flatness relation is necessary for a counterexample, not sufficient.")


print("\n=== Scope verdict ===")
print("This exact audit supports a conditional/cited reduction only:")
print("  counterexample => lattice-free projected zonotope => one relation with l1 norm <=356.")
print("It does not preserve support, sign pattern, endpoint owner, or semantic arrival.")
print("THM-3743 remains RESERVED / UNPROVED and is not used here.")
print(f"CHECKS={checks}")
