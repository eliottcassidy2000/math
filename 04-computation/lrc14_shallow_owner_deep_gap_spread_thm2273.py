#!/usr/bin/env python3
"""Exact folded and transport ledger for THM-2273.

The analytic proof uses THM-2080's guard/danger fold, THM-2263's
valuation-gap danger-pair cap, and the root-occupancy inequalities already
proved in THM-2255.  This companion checks:

* the parity-sharp guard/deep debit and both equality families;
* all 150 strict first-depth-one profiles and the unique worst row;
* compatibility of the two load-bearing pair equalities;
* the common-time two-shallow transport optimization;
* the uniform and boundary-profile image floors; and
* the exact deepest-successor safe-gap counts.

Every check remains active under ``python -O``.
"""

from fractions import Fraction as Q
from math import gcd


P = 13
DELTA5 = Q(961, 6930)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def thm2080_fold(x: Q, y: Q) -> Q:
    return min(x, y) + max(Q(0), x + y - 1) - 2 * x * y


def guard_normalized_debit(parity: int, a: int, k: int) -> Q:
    """p times the correction above 5/49 for C_a intersect D_(p*k)."""
    require(parity in (0, 1), "parity must be zero or one")
    require(a > 0 and a % 2 == 1, "reduced guard factor must be odd")
    require(k > 0 and gcd(a, k) == 1, "guard/deep row must be reduced")
    require(a % P and k % P, "normalized factors must be thirteen-units")
    x = Q(a % 14, 14)
    signed_k = k if parity == 0 else -k
    y = Q(signed_k % 7, 7)
    return -2 * thm2080_fold(x, y) / (a * k)


def guard_deep_cap(depth: int) -> Q:
    power = P**depth
    correction = Q(5, 49 * power) if depth % 2 else Q(5, 294 * power)
    return Q(5, 49) + correction


def guard_deep_overlap(guard: int, deep: int) -> Q:
    """Exact mu(C_guard intersect D_deep), using THM-2080's fold."""
    common = gcd(guard, deep)
    a = guard // common
    b = deep // common
    require(a % 2 == 1, "reduced guard is not odd")
    x = Q(a % 14, 14)
    y = Q(b % 7, 7)
    eligibility_overlap = Q(2, 49) + Q(2, a * b) * thm2080_fold(x, y)
    return Q(1, 7) - eligibility_overlap


def fold14(value: int) -> int:
    residue = value % 14
    return residue * (14 - residue)


def danger_overlap(left: int, right: int) -> Q:
    common = gcd(left, right)
    a = left // common
    b = right // common
    if a > b:
        a, b = b, a
    return Q(
        4 * a * b + fold14(a + b) - fold14(b - a),
        196 * a * b,
    )


def danger_gap_upper(depth: int) -> Q:
    power = P**depth
    correction = Q(6, 49 * power) if depth % 2 == 0 else Q(5, 588 * power)
    return Q(1, 49) + correction


def shallow_flow_floor(middle: int, deepest: int) -> Q:
    return (
        DELTA5
        - guard_deep_cap(deepest)
        - danger_gap_upper(middle - 1)
    )


def ceil_fraction(value: Q) -> int:
    return -(-value.numerator // value.denominator)


def safe_gap_count_floor(image_mass: Q, successor_depth: int) -> int:
    require(successor_depth >= 0, "successor depth must be nonnegative")
    return ceil_fraction(Q(7 * P**successor_depth, 6) * image_mass)


def main() -> None:
    # The exact normalized guard/deep endpoint.  The universal
    # |F|<=1/8 tail leaves two odd-depth cells and twenty-two even-depth
    # cells.
    odd_bank = [
        (guard_normalized_debit(1, a, k), a, k)
        for a in range(1, 3, 2)
        for k in range(1, 3 // a)
        if a * k <= 2 and gcd(a, k) == 1 and a % P and k % P
    ]
    # The range above produces k=1,2.
    require(len(odd_bank) == 2, "odd guard/deep bank size drift")
    odd_max = max(row[0] for row in odd_bank)
    odd_equal = {(a, k) for value, a, k in odd_bank if value == odd_max}
    require(odd_max == Q(5, 49), "odd guard/deep endpoint drift")
    require(odd_equal == {(1, 1)}, "odd guard/deep equality drift")
    require(Q(1, 12) < odd_max, "odd analytic tail does not separate")

    even_bank = [
        (guard_normalized_debit(0, a, k), a, k)
        for a in range(1, 15, 2)
        for k in range(1, 15 // a + 1)
        if a * k <= 14 and gcd(a, k) == 1 and a % P and k % P
    ]
    require(len(even_bank) == 22, "even guard/deep bank size drift")
    even_max = max(row[0] for row in even_bank)
    even_equal = {(a, k) for value, a, k in even_bank if value == even_max}
    require(even_max == Q(5, 294), "even guard/deep endpoint drift")
    require(even_equal == {(1, 6)}, "even guard/deep equality drift")
    require(Q(1, 60) < even_max, "even analytic tail does not separate")

    for depth in range(1, 20):
        if depth % 2:
            witness = (1, P**depth)
        else:
            witness = (1, 6 * P**depth)
        require(
            guard_deep_overlap(*witness) == guard_deep_cap(depth),
            f"guard/deep equality family drift at depth {depth}",
        )

    # All strict profiles, including both geometric boundary classes.
    rows = []
    for deepest in range(5, 20):
        for middle in range(2, deepest):
            rows.append(
                (
                    shallow_flow_floor(middle, deepest),
                    (1, middle, deepest),
                )
            )
    require(len(rows) == 150, "strict profile census drift")
    uniform_floor = min(value for value, _ in rows)
    worst_profiles = {
        profile for value, profile in rows if value == uniform_floor
    }
    require(worst_profiles == {(1, 3, 5)}, "uniform worst profile drift")
    require(
        uniform_floor == Q(5696989, 367580070),
        "uniform shallow-flow floor drift",
    )

    proposed_coarse = Q(558290567, 36022846860)
    require(
        uniform_floor - proposed_coarse == Q(29, 72773428),
        "sharp guard/deep gain drift",
    )

    # The two pair maxima coexist on one scalar ladder.
    guard = 1
    c1, c2, c3 = P, P**3, P**5
    require(
        danger_overlap(c1, c2) == danger_gap_upper(2),
        "shallow-pair equality compatibility drift",
    )
    require(
        guard_deep_overlap(guard, c3) == guard_deep_cap(5),
        "guard/deep equality compatibility drift",
    )

    # Common time b+1.  A is the c1-owner expansion through time two,
    # C is the c2-owner expansion through time two, and B is the latter's
    # additional expiration expansion.
    a_factor = Q(169, 20)
    c_factor = Q(169, 120)
    b_factor = Q(2197, 240)
    split_e1 = Q(11, 23)
    split_e2 = Q(12, 23)
    common_factor = Q(2197, 460)
    require(split_e1 + split_e2 == 1, "transport split is not normalized")
    require(
        a_factor * split_e1 + c_factor * split_e2 == common_factor,
        "noncontracted-union branch mismatch",
    )
    require(
        b_factor * split_e2 == common_factor,
        "second-owner expiration branch mismatch",
    )

    uniform_image = common_factor * uniform_floor
    require(
        uniform_image == Q(5696989, 76962600),
        "uniform common-time image floor drift",
    )
    require(
        uniform_image - Q(1, 14) == Q(1397623, 538738200),
        "one-fourteenth image margin drift",
    )
    require(
        uniform_image - Q(6, 91) == Q(4357723, 538738200),
        "one depth-positive safe-gap margin drift",
    )
    require(uniform_image < Q(1, 13), "uniform image unexpectedly crosses 1/13")

    # Boundary b=2: the minimum is c=5 and the successor depth is at least
    # two, forcing at least fifteen safe gaps.
    b2_rows = [
        (
            shallow_flow_floor(2, deepest),
            (1, 2, deepest),
        )
        for deepest in range(5, 20)
    ]
    b2_floor = min(value for value, _ in b2_rows)
    b2_worst = {profile for value, profile in b2_rows if value == b2_floor}
    require(b2_worst == {(1, 2, 5)}, "b=2 worst profile drift")
    require(b2_floor == Q(80120351, 5146120980), "b=2 floor drift")
    b2_image = common_factor * b2_floor
    require(b2_image == Q(80120351, 1077476400), "b=2 image drift")
    require(
        safe_gap_count_floor(b2_image, 2) == 15,
        "b=2 safe-gap count drift",
    )

    # Adjacent-deep boundary c=b+1: there is no forced 13 in the successor
    # speed.  Its sharp uniform row is (1,5,6), but the image crosses 1/13.
    adjacent_rows = [
        (
            shallow_flow_floor(middle, middle + 1),
            (1, middle, middle + 1),
        )
        for middle in range(4, 19)
    ]
    adjacent_floor = min(value for value, _ in adjacent_rows)
    adjacent_worst = {
        profile for value, profile in adjacent_rows if value == adjacent_floor
    }
    require(adjacent_worst == {(1, 5, 6)}, "adjacent worst profile drift")
    require(
        adjacent_floor == Q(271263857, 16724893185),
        "adjacent shallow-flow floor drift",
    )
    adjacent_image = common_factor * adjacent_floor
    require(
        adjacent_image == Q(271263857, 3501798300),
        "adjacent common-time image drift",
    )
    require(
        adjacent_image - Q(1, 13) == Q(1894757, 3501798300),
        "adjacent one-thirteenth margin drift",
    )

    # Full profile-gap count supplied by the uniform image floor.
    gap_counts = {
        gap: safe_gap_count_floor(uniform_image, gap - 1)
        for gap in range(1, 18)
    }
    expected_prefix = {1: 1, 2: 2, 3: 15, 4: 190, 5: 2467}
    require(
        {gap: gap_counts[gap] for gap in expected_prefix} == expected_prefix,
        "safe-gap count prefix drift",
    )

    print("THM-2273 SHALLOW-OWNER / DEEP-SUCCESSOR GAP SPREAD")
    print(f"guard_deep_normalized_odd_max={odd_max} equality={sorted(odd_equal)}")
    print(f"guard_deep_normalized_even_max={even_max} equality={sorted(even_equal)}")
    print(f"strict_profiles={len(rows)} worst={sorted(worst_profiles)}")
    print(f"uniform_shallow_exact_one_floor={uniform_floor}")
    print(f"coarse_floor_gain={uniform_floor - proposed_coarse}")
    print(
        "transport_factors="
        f"A:{a_factor} C:{c_factor} B:{b_factor} optimized:{common_factor}"
    )
    print(f"uniform_common_time_image={uniform_image}")
    print(f"uniform_image_minus_1/14={uniform_image - Q(1, 14)}")
    print(f"uniform_image_minus_6/91={uniform_image - Q(6, 91)}")
    print(
        f"b2_floor/image/gaps={b2_floor} {b2_image} "
        f"{safe_gap_count_floor(b2_image, 2)}"
    )
    print(
        "adjacent_floor/image/minus_1/13="
        f"{adjacent_floor} {adjacent_image} {adjacent_image-Q(1, 13)}"
    )
    print("uniform_successor_gap_counts:")
    for gap in range(1, 18):
        print(f"  c-b={gap}: {gap_counts[gap]}")
    print("status=PROVED_VERIFIED_EXACT")


if __name__ == "__main__":
    main()
