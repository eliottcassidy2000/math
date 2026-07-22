#!/usr/bin/env python3
"""Exact referee for THM-2061's strict LRC(14) dyadic seam."""

import argparse
from fractions import Fraction
from itertools import combinations
from math import comb, gcd


DELTA = Fraction(1, 14)
ELIGIBILITY = Fraction(1, 7)
DIAMOND = Fraction(6, 7)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def norm_fraction(value):
    residue = value.numerator % value.denominator
    return Fraction(min(residue, value.denominator - residue), value.denominator)


def nearest_integer(value):
    quotient, residue = divmod(value.numerator, value.denominator)
    if 2 * residue == value.denominator:
        return None
    return quotient if 2 * residue < value.denominator else quotient + 1


def folded_clearance(theta, x, y):
    lifts = (theta / 2, (theta + 1) / 2)
    return max(
        min(norm_fraction(x * lift), norm_fraction(y * lift)) for lift in lifts
    )


def folded_measure_formula(alpha, beta):
    r = alpha - beta
    s = alpha + beta
    n = (2 * alpha + 7) // 14
    total = sum(
        min(Fraction(2 * r, 7), Fraction(2 * alpha, 7) - (2 * j - 1))
        for j in range(1, n + 1)
    )
    return Fraction(2, r * s) * total


def diamond_measure_piecewise(alpha, beta):
    breakpoints = sorted(
        {Fraction(k, 2 * alpha) for k in range(2 * alpha + 1)}
        | {Fraction(k, 2 * beta) for k in range(2 * beta + 1)}
    )
    answer = Fraction(0)
    for left, right in zip(breakpoints, breakpoints[1:]):
        value_left = norm_fraction(alpha * left) + norm_fraction(beta * left)
        value_right = norm_fraction(alpha * right) + norm_fraction(beta * right)
        if value_left > DIAMOND and value_right > DIAMOND:
            answer += right - left
        elif value_left <= DIAMOND and value_right <= DIAMOND:
            continue
        else:
            crossing = left + (DIAMOND - value_left) * (right - left) / (
                value_right - value_left
            )
            answer += right - crossing if value_right > DIAMOND else crossing - left
    return answer


def merge_closed(intervals):
    answer = []
    for left, right in sorted(intervals):
        if answer and left <= answer[-1][1]:
            answer[-1] = (answer[-1][0], max(answer[-1][1], right))
        else:
            answer.append((left, right))
    return answer


def intersect_closed(left, right):
    answer = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if hi >= lo:
            answer.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return merge_closed(answer)


def weak_safe_components(speeds):
    current = [(Fraction(0), Fraction(1))]
    for speed in speeds:
        safe = [
            ((Fraction(k) + DELTA) / speed,
             (Fraction(k + 1) - DELTA) / speed)
            for k in range(speed)
        ]
        current = intersect_closed(current, safe)
        if not current:
            break
    return tuple(current)


def universal_open_signature(components, speed):
    word = []
    minimum_slack = None
    for left, right in components:
        centre = (left + right) / 2
        half_width = (right - left) / 2
        owner = nearest_integer(speed * centre)
        if owner is None:
            return None
        slack = ELIGIBILITY - abs(speed * centre - owner) - speed * half_width
        if slack <= 0:
            return None
        word.append(owner & 1)
        minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
    return tuple(word), minimum_slack


def all_unit_eligible_odd_residues(modulus):
    units = [r for r in range(modulus) if gcd(r, modulus) == 1]
    return [
        z
        for z in range(1, 2 * modulus, 2)
        if all(7 * min((z * r) % modulus, modulus - (z * r) % modulus) < modulus
               for r in units)
    ]


def finite_core_slice(bound):
    scanned = 0
    primitive = 0
    pinned = 0
    universal_tails = 0
    complementary_pairs = 0
    max_cap = 0
    hardest = None
    for core in combinations(range(1, bound + 1), 11):
        scanned += 1
        if gcd(*core) != 1:
            continue
        primitive += 1
        if not all(any(c % modulus == 0 for c in core) for modulus in range(2, 15)):
            continue
        pinned += 1
        components = weak_safe_components(core)
        require(components, ("empty weak core", core))
        widest = max(right - left for left, right in components)
        require(widest > 0, ("no bulk component", core))
        strict_ratio = Fraction(2, 7) / widest
        cap = (strict_ratio.numerator - 1) // strict_ratio.denominator
        max_cap = max(max_cap, cap)
        if hardest is None or (cap, -widest, core) > (hardest[0], -hardest[1], hardest[2]):
            hardest = (cap, widest, core, len(components))

        by_word = {}
        for speed in range(1, cap + 1, 2):
            signature = universal_open_signature(components, speed)
            if signature is None:
                continue
            word, _ = signature
            universal_tails += 1
            by_word.setdefault(word, []).append(speed)
        for word, speeds in by_word.items():
            complement = tuple(1 - bit for bit in word)
            if complement not in by_word or word >= complement:
                continue
            complementary_pairs += len(speeds) * len(by_word[complement])

    return {
        "scanned": scanned,
        "expected": comb(bound, 11),
        "primitive": primitive,
        "pinned": pinned,
        "universal_tails": universal_tails,
        "complementary_pairs": complementary_pairs,
        "max_cap": max_cap,
        "hardest": hardest,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--finite-N", type=int, default=19)
    args = parser.parse_args()

    fold_rows = 0
    for denominator in range(2, 61):
        for numerator in range(denominator):
            theta = Fraction(numerator, denominator)
            for x in range(3, 24, 2):
                for y in range(1, x, 2):
                    a = (x + y) // 2
                    b = (x - y) // 2
                    folded = folded_clearance(theta, x, y)
                    diamond_value = norm_fraction(a * theta) + norm_fraction(b * theta)
                    require(folded == (1 - diamond_value) / 2,
                            ("fold identity", theta, x, y))
                    parity_obstruction = False
                    if norm_fraction(x * theta) < ELIGIBILITY and norm_fraction(y * theta) < ELIGIBILITY:
                        n_x = nearest_integer(x * theta)
                        n_y = nearest_integer(y * theta)
                        require(n_x is not None and n_y is not None,
                                ("eligible half-integer", theta, x, y))
                        parity_obstruction = (n_x - n_y) % 2 == 1
                        if parity_obstruction:
                            determinant = abs(n_x * y - n_y * x)
                            require(determinant > 0 and determinant % 2 == 1,
                                    ("odd determinant", theta, x, y))
                            require(determinant % gcd(x, y) == 0,
                                    ("gcd determinant", theta, x, y))
                    require((folded < DELTA) == parity_obstruction,
                            ("strict parity iff", theta, x, y))
                    require((folded < DELTA) == (diamond_value > DIAMOND),
                            ("strict diamond iff", theta, x, y))
                    fold_rows += 1

    residue_table = {}
    for modulus in range(2, 15):
        eligible = all_unit_eligible_odd_residues(modulus)
        expected = [] if modulus % 2 == 0 else [modulus]
        require(eligible == expected, ("unit residue table", modulus, eligible))
        residue_table[modulus] = eligible

    measure_rows = 0
    cap_rows = 0
    maximizers = []
    cap = Fraction(4, 63)
    for alpha in range(2, 1001):
        for beta in range(1, alpha):
            if gcd(alpha, beta) != 1 or (alpha - beta) % 2 == 0:
                continue
            formula = folded_measure_formula(alpha, beta)
            require(formula <= cap, ("measure cap", alpha, beta, formula))
            if formula == cap:
                maximizers.append((alpha, beta))
            cap_rows += 1
            if alpha <= 160:
                direct = diamond_measure_piecewise(alpha, beta)
                require(formula == direct,
                        ("measure formula", alpha, beta, formula, direct))
                measure_rows += 1
    require(maximizers == [(5, 4)], ("unique reduced maximizer", maximizers))

    finite = finite_core_slice(args.finite_N)
    require(finite["scanned"] == finite["expected"], "finite census count")
    require(finite["complementary_pairs"] == 0,
            ("finite seam survivor", finite))

    print("THM-2061 STRICT DYADIC TWO-TAIL SEAM AUDIT")
    print("fold/parity/diamond rational rows checked:", fold_rows)
    print("all-unit eligible odd residues mod 2N:")
    print(" even N=2,4,6,8,10,12,14: none")
    print(" odd N=3,5,7,9,11,13: [N]")
    print("independent exact measure rows checked:", measure_rows)
    print("reduced measure-cap rows checked:", cap_rows)
    print("sharp reduced maximizer:", maximizers)
    print(
        "finite core slice N={}: all={} primitive={} divisor_pinned={}".format(
            args.finite_N, finite["scanned"], finite["primitive"], finite["pinned"]
        )
    )
    print(
        "finite universal odd tails={} complementary pairs={}".format(
            finite["universal_tails"], finite["complementary_pairs"]
        )
    )
    print("finite max intrinsic odd cap / hardest core:",
          finite["max_cap"], finite["hardest"])
    print("endpoint polarity: closed weak core / strict open tail teeth")
    print("carrier: ordered core components + parity words + rational tooth slacks")
    print("tournament: binary ownership only; widths and determinant remain sidecars")
    print("PASS")


if __name__ == "__main__":
    main()
