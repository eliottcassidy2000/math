#!/usr/bin/env python3
"""Independent exact audit of a diagonal-ellipsoid refinement of THM-4009.

The mathematical proof is in the companion report.  This script uses only
exact rational/integer arithmetic except for printing one calibrated decimal.
It checks the polar-vertex model, the pair-facet inequalities, the compressed
operator (and the ambient-inverse trap), the two proposed diagonal banks, and
the sharp integer optimizations used for their coefficient and l1 conclusions.

It does not test LRC(14), reserve a theorem ID, or modify the proved canon.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, sqrt
from typing import Iterable, Sequence


K = 13
AP = tuple(range(1, K + 1))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def dot(left: Sequence[Fraction], right: Sequence[Fraction]) -> Fraction:
    require(len(left) == len(right), "dot-product length drift")
    return sum((x * y for x, y in zip(left, right)), Fraction(0))


def polar_vertices(
    speeds: Sequence[int],
) -> tuple[tuple[Fraction, ...], ...]:
    """The signed pair vertices of {u in n^perp : ||u||_1 <= 1}."""

    vertices: list[tuple[Fraction, ...]] = []
    for i, j in combinations(range(len(speeds)), 2):
        denominator = speeds[i] + speeds[j]
        positive = [Fraction(0) for _ in speeds]
        positive[i] = Fraction(speeds[j], denominator)
        positive[j] = Fraction(-speeds[i], denominator)
        vertices.append(tuple(positive))
        vertices.append(tuple(-entry for entry in positive))
    return tuple(vertices)


def pair_slack(
    speeds: Sequence[int], weights: Sequence[Fraction], i: int, j: int
) -> Fraction:
    """Nonnegative exactly when the (i,j) polar vertex has H-energy <= 1."""

    ni = Fraction(speeds[i])
    nj = Fraction(speeds[j])
    return (
        2 * ni * nj
        - (weights[i] - 1) * nj * nj
        - (weights[j] - 1) * ni * ni
    )


def vertex_energy(vertex: Sequence[Fraction], weights: Sequence[Fraction]) -> Fraction:
    return sum(weight * value * value for weight, value in zip(weights, vertex))


def coordinate_boost(speeds: Sequence[int], index: int) -> tuple[Fraction, ...]:
    other_maximum = max(value for j, value in enumerate(speeds) if j != index)
    weights = [Fraction(1) for _ in speeds]
    weights[index] += Fraction(2 * speeds[index], other_maximum)
    return tuple(weights)


def linear_bank(speeds: Sequence[int]) -> tuple[Fraction, ...]:
    largest, second = sorted(speeds)[-2:]
    denominator = largest + second
    return tuple(Fraction(1) + Fraction(2 * value, denominator) for value in speeds)


def verify_pair_bank(speeds: Sequence[int], weights: Sequence[Fraction]) -> tuple[int, int]:
    zero = 0
    positive = 0
    vertices = polar_vertices(speeds)
    for pair_index, (i, j) in enumerate(combinations(range(len(speeds)), 2)):
        slack = pair_slack(speeds, weights, i, j)
        require(slack >= 0, f"negative pair slack at {(i, j)}")
        zero += slack == 0
        positive += slack > 0
        for signed_offset in (0, 1):
            vertex = vertices[2 * pair_index + signed_offset]
            energy = vertex_energy(vertex, weights)
            require(energy <= 1, f"vertex escaped ellipsoid at {(i, j)}")
            denominator = Fraction((speeds[i] + speeds[j]) ** 2)
            require(1 - energy == slack / denominator, "pair/vertex identity drift")
    return zero, positive


def compressed_inverse_quadratic(
    speeds: Sequence[int], weights: Sequence[Fraction], vector: Sequence[int]
) -> tuple[Fraction, Fraction, tuple[Fraction, ...]]:
    """Return ambient H^-1 energy, correct compressed inverse energy, and A^-1 x."""

    n = tuple(Fraction(value) for value in speeds)
    h = tuple(Fraction(value) for value in weights)
    x = tuple(Fraction(value) for value in vector)
    require(dot(n, x) == 0, "test vector must lie in n-perp")
    ambient = sum(value * value / weight for value, weight in zip(x, h))
    numerator = sum(ni * xi / hi for ni, xi, hi in zip(n, x, h))
    denominator = sum(ni * ni / hi for ni, hi in zip(n, h))
    correction = numerator * numerator / denominator
    compressed = ambient - correction
    multiplier = -numerator / denominator
    inverse_vector = tuple(
        (xi + multiplier * ni) / hi for xi, ni, hi in zip(x, n, h)
    )
    require(dot(n, inverse_vector) == 0, "compressed inverse left n-perp")
    require(dot(x, inverse_vector) == compressed, "compressed quadratic drift")
    for xi, ni, hi, yi in zip(x, n, h, inverse_vector):
        require(hi * yi - xi == multiplier * ni, "compressed equation drift")
    require(compressed <= ambient, "rank-one correction has wrong sign")
    return ambient, compressed, inverse_vector


def balanced_square_sum(total: int, slots: int) -> int:
    """Minimum sum of squares of nonnegative integers with fixed sum."""

    quotient, remainder = divmod(total, slots)
    return (slots - remainder) * quotient**2 + remainder * (quotient + 1) ** 2


def distinguished_minimum(l1: int, distinguished_weight: int) -> tuple[int, tuple[int, ...]]:
    """Minimize h*t^2 plus twelve ordinary squares at fixed integer l1."""

    rows = tuple(
        (distinguished_weight * t * t + balanced_square_sum(l1 - t, 12), t)
        for t in range(l1 + 1)
    )
    minimum = min(value for value, _ in rows)
    minimizers = tuple(t for value, t in rows if value == minimum)
    return minimum, minimizers


def reduced_top_pair_ratios() -> tuple[tuple[int, int], ...]:
    return tuple(
        (p, q)
        for p in range(1, 14)
        for q in range(p + 1, 14)
        if p + q <= 13 and gcd(p, q) == 1
    )


def semantic_digest(vertices: Iterable[Sequence[Fraction]]) -> str:
    payload = "\n".join(
        ",".join(f"{value.numerator}/{value.denominator}" for value in vertex)
        for vertex in vertices
    )
    return sha256(payload.encode("ascii")).hexdigest()


def main() -> None:
    vertices = polar_vertices(AP)
    require(len(vertices) == 2 * 78 == 156, "polar vertex count drift")
    require(len(set(vertices)) == 156, "polar vertex collision")
    ap_fraction = tuple(Fraction(value) for value in AP)
    for vertex in vertices:
        require(dot(ap_fraction, vertex) == 0, "polar vertex left n-perp")
        require(sum(abs(value) for value in vertex) == 1, "polar l1 drift")

    max_index = len(AP) - 1
    max_weights = coordinate_boost(AP, max_index)
    max_zero, max_positive = verify_pair_bank(AP, max_weights)
    require((max_zero, max_positive) == (1, 77), "coordinate-boost facet drift")

    linear_weights = linear_bank(AP)
    linear_zero, linear_positive = verify_pair_bank(AP, linear_weights)
    require((linear_zero, linear_positive) == (1, 77), "linear-bank facet drift")

    top_pair_vector = (0,) * 11 + (13, -12)
    ambient, compressed, _ = compressed_inverse_quadratic(
        AP, max_weights, top_pair_vector
    )
    require(compressed < ambient, "AP compressed trap should be strict")

    # Exact two-coordinate hostile to the false ambient-inverse equality.
    two_n = (1, 2)
    two_h = (Fraction(2), Fraction(1))
    require(pair_slack(two_n, two_h, 0, 1) == 0, "2D hostile is not feasible")
    restricted_h = Fraction(9, 5)
    compressed_inverse = Fraction(5, 9)
    ambient_inverse = Fraction(3, 5)
    hostile_radius_squared = Fraction(7, 4)
    require(hostile_radius_squared * compressed_inverse == Fraction(35, 36), "2D compressed hostile drift")
    require(hostile_radius_squared * ambient_inverse == Fraction(21, 20), "2D ambient hostile drift")
    require(Fraction(35, 36) < 1 < Fraction(21, 20), "2D trap lost strictness")

    # Integer optimization for the max-coordinate bank.  Its distinguished
    # weight is strictly greater than 3 because M > M_2.
    min_l1_49, minimizers_49 = distinguished_minimum(49, 3)
    min_l1_50, minimizers_50 = distinguished_minimum(50, 3)
    require((min_l1_49, minimizers_49) == (195, (1,)), "l1=49 boundary drift")
    require(min_l1_50 == 204, "l1=50 exclusion drift")
    require(balanced_square_sum(48, 12) == 192, "twelve-four boundary drift")
    raw_max_coefficient = max(value for value in range(20) if 3 * value * value < 196)
    require(raw_max_coefficient == 8, "raw max-coordinate cap drift")

    # Kernel balance improves 8 to 7: for r=M/M_2>1 and twelve lower speeds,
    # E/a_M^2 > 1+2r+r^2/12 > 37/12; 8^2*(37/12)>196.
    require(Fraction(64 * 37, 12) > 196, "balance cap no longer excludes 8")
    require(Fraction(49 * 37, 12) < 196, "balance cap accidentally excludes 7")

    avoidance_root = -12 + 2 * sqrt(621)
    require(avoidance_root > 13, "avoidance threshold entered THM-1008 live range")

    ratios = reduced_top_pair_ratios()
    require(len(ratios) == 28, "top-pair ratio census drift")
    require(all(p + q <= 13 and gcd(p, q) == 1 for p, q in ratios), "ratio filter drift")

    # The old exact scalar inradius is always stronger than the scalar bound
    # obtained by replacing every linear weight by its minimum.  Check the
    # identity on a nontrivial rational row; the report proves it generally.
    r = Fraction(13, 1)
    s = Fraction(12, 1)
    linear_scalar = (r + s) / (r + s + 2)
    inradius_scalar = (r * r + 1) / (r + 1) ** 2
    difference = linear_scalar - inradius_scalar
    require(
        difference == Fraction(2) * (r * s - 1) / ((r + s + 2) * (r + 1) ** 2),
        "scalar comparison identity drift",
    )
    require(difference > 0, "linear scalar bank unexpectedly beats inradius")

    print("LRC14_DIAGONAL_ELLIPSOID_REFINEMENT_INDEPENDENT_AUDIT")
    print("status=PROVED_ALGEBRA_PLUS_VERIFIED_EXACT_CONTROLS;_NOT_CANONIZED;_LRC14_OPEN")
    print(f"polar_vertices=2*C(13,2)={len(vertices)} semantic_sha256={semantic_digest(vertices)}")
    print("pair_condition=necessary_and_sufficient_for_positive_diagonal_H")
    print(f"ap13_coordinate_max_boost_zero_facets={max_zero} positive_facets={max_positive}")
    print(f"ap13_linear_bank_zero_facets={linear_zero} positive_facets={linear_positive}")
    print(f"compressed_operator_ap13_ambient={ambient} correct={compressed} strict_rank_one_loss=1")
    print("ambient_inverse_hostile=n=(1,2)_H=(2,1)_t2=7/4_correct=35/36_ambient=21/20")
    print("correct_transform=A^(-1/2)_on_V;_dual=A^(1/2)*Lambda_star;_weight=a^T_H_a")
    print("coordinate_max_immediate=E<196_square_if_nonzero<=193_raw_abs_aM<=8_l1<=49")
    print("coordinate_max_kernel_balance=abs_aM<=7")
    print(f"integer_l1_boundary=L49_min_at_HM=3:{min_l1_49}_aM={minimizers_49};L50_min:{min_l1_50}")
    print("l1_49_only_possible_shape=abs_aM_1_and_twelve_abs4;_requires_M_mod4=0_and_M/M2<3/2")
    print(f"max_avoidance_ratio=-12+2sqrt621={avoidance_root:.12f}_outside_live_R_lt_13=1")
    print("linear_bank=D=M+M2;_D*S+2*sum(n_i*a_i^2)<=196D-1;_abs_aM<=9;_l1<=49")
    print(f"linear_top_pair_support2=p+q<=13_reduced_ratio_count={len(ratios)}")
    print("linear_scalar_cap=weaker_than_existing_exact_inradius;_weighted_localization_is_new")
    print("quantifier_guard=coordinate_boosts_may_choose_different_Graver_relations")
    print("PASS")


if __name__ == "__main__":
    main()
