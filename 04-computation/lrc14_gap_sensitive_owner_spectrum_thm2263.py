#!/usr/bin/env python3
"""Exact folded-formula and profile-ledger audit for THM-2263."""

from fractions import Fraction
from math import gcd


P = 13
DELTA5 = Fraction(961, 6930)
OLD_STRICT = Fraction(88159, 1171170)
OLD_REPEATED = Fraction(14627, 585585)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fold14(value):
    residue = value % 14
    return residue * (14 - residue)


def defect(a, k):
    """The even-gap folded numerator F(a+k)-F(k-a)."""
    return fold14(a + k) - fold14(k - a)


def quotient(a, k):
    return Fraction(defect(a, k), a * k)


def rho(a, b):
    common = gcd(a, b)
    a //= common
    b //= common
    return Fraction(
        4 * a * b + fold14(a + b) - fold14(b - a),
        196 * a * b,
    )


def gap_upper(depth):
    power = P**depth
    if depth % 2 == 0:
        return Fraction(1, 49) + Fraction(6, 49 * power)
    return Fraction(1, 49) + Fraction(5, 588 * power)


def gap_lower(depth):
    power = P**depth
    if depth % 2 == 0:
        return Fraction(1, 49) - Fraction(5, 588 * power)
    return Fraction(1, 49) - Fraction(6, 49 * power)


def cross_score_even(g, u):
    """196*13^d times the two cross-pair correction for shallow ratio 1:2."""
    require(gcd(g, u) == 1, "cross-score inputs must be coprime")
    if u % 2:
        return quotient(g, u) + quotient(2 * g, u)
    require(g % 2 == 1, "coprime even-u input has odd g")
    return quotient(g, u) + quotient(g, u // 2)


def cross_upper(depth):
    power = P**depth
    if depth % 2 == 0:
        return Fraction(34, 196 * power)
    return Fraction(32, 11 * 196 * power)


def main():
    # The universal folded quotient interval. The analytic tails are
    # |defect|<=49; only these small banks can meet an endpoint.
    lower_bank = [
        (quotient(a, k), a, k)
        for a in range(1, 30)
        for k in range(1, 30)
        if a * k <= 29 and gcd(a, k) == 1 and a % P and k % P
    ]
    require(len(lower_bank) == 75, "lower endpoint bank drift")
    lower_value = min(row[0] for row in lower_bank)
    lower_equal = {
        (a, k) for value, a, k in lower_bank if value == lower_value
    }
    require(lower_value == Fraction(-5, 3), "folded lower endpoint drift")
    require(
        lower_equal == {(1, 12), (12, 1)},
        "folded lower equality locus drift",
    )

    upper_bank = [
        (quotient(a, k), a, k)
        for a in range(1, 3)
        for k in range(1, 3)
        if a * k < 3 and gcd(a, k) == 1
    ]
    upper_value = max(row[0] for row in upper_bank)
    upper_equal = {
        (a, k) for value, a, k in upper_bank if value == upper_value
    }
    require(upper_value == 24, "folded upper endpoint drift")
    require(upper_equal == {(1, 1)}, "folded upper equality locus drift")
    require(Fraction(-49, 30) > lower_value, "lower analytic tail failed")
    require(Fraction(49, 3) < upper_value, "upper analytic tail failed")

    # Check every gap used by the 165-row ledger and both equality families.
    for depth in range(1, 19):
        power = P**depth
        require(
            rho(1, power) == (
                gap_upper(depth) if depth % 2 == 0 else gap_lower(depth)
            ),
            f"1:p equality drift at gap {depth}",
        )
        require(
            rho(12, power) == (
                gap_lower(depth) if depth % 2 == 0 else gap_upper(depth)
            ),
            f"12:p equality drift at gap {depth}",
        )
        require(
            rho(1, 12 * power) == rho(12, power),
            f"reciprocal equality family drift at gap {depth}",
        )

    # Strict profiles: exact gap-sensitive pair ledger.
    strict_rows = []
    for deepest in range(5, 20):
        for middle in range(2, deepest):
            cap_sum = (
                gap_upper(middle - 1)
                + gap_upper(deepest - middle)
                + gap_upper(deepest - 1)
            )
            exclusive = DELTA5 - cap_sum
            strict_rows.append((exclusive, (1, middle, deepest), cap_sum))
    require(len(strict_rows) == 150, "strict profile census drift")
    strict_floor = min(row[0] for row in strict_rows)
    strict_worst = {
        profile for value, profile, _ in strict_rows if value == strict_floor
    }
    strict_cap = max(row[2] for row in strict_rows)
    require(strict_worst == {(1, 3, 5)}, "strict worst profile drift")
    require(strict_cap == Fraction(12531, 199927), "strict cap sum drift")
    require(
        strict_floor == Fraction(15041431, 197927730),
        "strict exclusive floor drift",
    )
    require(
        strict_floor - OLD_STRICT == Fraction(144, 199927),
        "strict improvement drift",
    )
    strict_owner = strict_floor / 3
    strict_expiration = Fraction(169, 20) * strict_owner
    require(
        strict_owner == Fraction(15041431, 593783190),
        "strict owner floor drift",
    )
    require(
        strict_expiration == Fraction(15041431, 70270200),
        "strict expiration floor drift",
    )
    require(strict_expiration > Fraction(1, 7), "strict threshold not crossed")

    # The strict equality carrier is compatible on the common-core ladder.
    strict_speeds = (1, P**2, P**4)
    strict_pair_sum = sum(
        (
            rho(strict_speeds[i], strict_speeds[j])
            for i in range(3)
            for j in range(i + 1, 3)
        ),
        Fraction(),
    )
    require(strict_pair_sum == strict_cap, "strict equality carrier drift")

    # Universal second maximum for a distinct same-valuation pair.
    same_bank = [
        (rho(a, b), a, b)
        for a in range(1, 10)
        for b in range(a + 1, 10)
        if a * b <= 9 and gcd(a, b) == 1
    ]
    same_nonmax = [row for row in same_bank if row[1:] != (1, 2)]
    require(
        max(row[0] for row in same_nonmax) == Fraction(1, 21),
        "same-valuation second maximum drift",
    )
    require(
        Fraction(1, 49) + Fraction(1, 40) < Fraction(1, 21),
        "same-valuation analytic tail failed",
    )

    # If the shallow pair is 1:2, the two deep cross-pair corrections
    # have their own exact parity spectrum.
    cross_bank = [
        (cross_score_even(g, u), g, u)
        for g in range(1, 51)
        for u in range(1, 51)
        if g * u <= 50 and gcd(g, u) == 1 and g % P and u % P
    ]
    cross_max = max(row[0] for row in cross_bank)
    cross_min = min(row[0] for row in cross_bank)
    require(cross_max == 34, "even cross-score maximum drift")
    require(
        {(g, u) for value, g, u in cross_bank if value == cross_max}
        == {(1, 1), (1, 2)},
        "even cross-score equality locus drift",
    )
    require(cross_min == Fraction(-32, 11), "odd cross-score maximum drift")
    require(
        {(g, u) for value, g, u in cross_bank if value == cross_min}
        == {(1, 11), (11, 2)},
        "odd cross-score equality locus drift",
    )

    # Repeated-first profiles. The 1:2 branch beats the universal
    # non-1:2 bound, and gap four is the unique maximum.
    repeated_rows = []
    for deepest in range(5, 20):
        depth = deepest - 1
        ratio_two_cap = Fraction(1, 14) + Fraction(2, 49) + cross_upper(depth)
        non_ratio_two_cap = Fraction(1, 21) + 2 * gap_upper(depth)
        require(
            ratio_two_cap > non_ratio_two_cap,
            f"same-pair branch comparison drift at {deepest}",
        )
        repeated_rows.append(
            (DELTA5 - ratio_two_cap, (1, 1, deepest), ratio_two_cap)
        )
    require(len(repeated_rows) == 15, "repeated profile census drift")
    repeated_floor = min(row[0] for row in repeated_rows)
    repeated_worst = {
        profile for value, profile, _ in repeated_rows if value == repeated_floor
    }
    repeated_cap = max(row[2] for row in repeated_rows)
    require(repeated_worst == {(1, 1, 5)}, "repeated worst profile drift")
    require(repeated_cap == Fraction(3206, 28561), "repeated cap sum drift")
    require(
        repeated_floor == Fraction(5229541, 197927730),
        "repeated exclusive floor drift",
    )
    require(
        repeated_floor - OLD_REPEATED == Fraction(577, 399854),
        "repeated improvement drift",
    )
    repeated_owner = repeated_floor / 3
    repeated_expiration = Fraction(169, 20) * repeated_owner
    require(
        repeated_expiration == Fraction(5229541, 70270200),
        "repeated expiration floor drift",
    )
    require(
        repeated_expiration - Fraction(1, 14)
        == Fraction(210241, 70270200),
        "repeated half-comb margin drift",
    )
    require(
        repeated_expiration < Fraction(1, 7),
        "repeated one-comb threshold unexpectedly crossed",
    )

    repeated_speeds = (1, 2, P**4)
    repeated_pair_sum = sum(
        (
            rho(repeated_speeds[i], repeated_speeds[j])
            for i in range(3)
            for j in range(i + 1, 3)
        ),
        Fraction(),
    )
    require(repeated_pair_sum == repeated_cap, "repeated equality carrier drift")

    print("THM-2263 GAP-SENSITIVE OWNER SPECTRUM -- exact audit")
    print(f"folded quotient interval: [{lower_value}, {upper_value}]")
    print(f"lower endpoint bank/equalities: {len(lower_bank)} {sorted(lower_equal)}")
    print(f"strict profiles/worst: {len(strict_rows)} {sorted(strict_worst)}")
    print(f"strict pair cap sum: {strict_cap}")
    print(f"strict exclusive/owner: {strict_floor} {strict_owner}")
    print(f"strict expiration image: {strict_expiration}")
    print(f"repeated profiles/worst: {len(repeated_rows)} {sorted(repeated_worst)}")
    print(f"repeated pair cap sum: {repeated_cap}")
    print(f"repeated exclusive/expiration: {repeated_floor} {repeated_expiration}")
    print(f"cross-score interval: [{cross_min}, {cross_max}]")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
