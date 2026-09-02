#!/usr/bin/env python3
"""Exact r_7=7 minority-anchor endpoint-pierce method hostile.

The displayed row is safe by the c=7 specialization of THM-771.  This
artifact does not prove LRC(14).  It verifies that three currently used
sufficient banks all fail on the row, then records the exact THM-771
capacity-slack/owner-overlap calculation at one core-safe endpoint.

All checks use integer or fractions.Fraction arithmetic only.
"""

from collections import Counter
from fractions import Fraction
from math import gcd, isqrt


H = 420
ANCHOR = 2 * H
HALF_TURN_MODULUS = 28 * H
NONSEVEN_TAILS = tuple(11 + 1680 * k for k in range(7))
SEVEN_MULTIPLES = (525, 945, 1365, 1575, 9009)
W = NONSEVEN_TAILS + SEVEN_MULTIPLES
S = (ANCHOR,) + W

ENDPOINT_OWNER = 11
ENDPOINT_NUMERATOR = 1
ENDPOINT_DENOMINATOR = 22
OUTER_SCALE = 7
LIFT_DENOMINATOR = OUTER_SCALE * ENDPOINT_DENOMINATOR


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def circular_numerator(v: int, numerator: int, denominator: int) -> int:
    residue = (v * numerator) % denominator
    return min(residue, denominator - residue)


def factorization(n: int) -> tuple[tuple[int, int], ...]:
    factors = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            factors.append((p, exponent))
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


def product_divisors(n: int) -> set[int]:
    result = {1}
    for p, exponent in factorization(n):
        result = {d * p**e for d in result for e in range(exponent + 1)}
    return result


def trial_divisors(n: int) -> set[int]:
    result = {1, n}
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            result.add(d)
            result.add(n // d)
    return result


def mc7_capacity(a: int) -> Fraction:
    total = Fraction()
    for v in S:
        if v % a:
            q = a // gcd(a, v)
            total += Fraction((q + 6) // 7, q)
    return total


def danger_sheets(v: int) -> tuple[int, ...]:
    return tuple(
        k
        for k in range(OUTER_SCALE)
        if 14
        * circular_numerator(
            v,
            ENDPOINT_NUMERATOR + ENDPOINT_DENOMINATOR * k,
            LIFT_DENOMINATOR,
        )
        < LIFT_DENOMINATOR
    )


def tight_sheets(v: int) -> tuple[int, ...]:
    return tuple(
        k
        for k in range(OUTER_SCALE)
        if 14
        * circular_numerator(
            v,
            ENDPOINT_NUMERATOR + ENDPOINT_DENOMINATOR * k,
            LIFT_DENOMINATOR,
        )
        == LIFT_DENOMINATOR
    )


def main() -> None:
    require(H % 420 == 0, "the row left the 420 wall")
    require(len(S) == 13 and len(W) == len(set(W)) == 12, "row cardinality")
    require(all(w > 0 and w % 2 for w in W), "W is not positive odd")
    require(gcd(*S) == 1, "row is not primitive")

    core = tuple(v // OUTER_SCALE for v in S if v % OUTER_SCALE == 0)
    tails = tuple(v for v in S if v % OUTER_SCALE)
    require(core == (120, 75, 135, 195, 225, 1287), "c=7 core changed")
    require(tails == NONSEVEN_TAILS and len(tails) == 7, "r_7 is not seven")
    missing_gates = tuple(
        m for m in range(2, 15) if not any(v % m == 0 for v in S)
    )
    require(not missing_gates, "a THM-366 denominator gate is missing")

    # Complete integer half-turn grid.  The first block partitions j with
    # 7 not dividing j; the second partitions j=7s with s odd.
    modulus = HALF_TURN_MODULUS
    d_partition = Counter(
        sum(5040 < (j * d) % modulus < 6720 for d in NONSEVEN_TAILS)
        for j in range(modulus)
        if j % 7
    )
    c_block = SEVEN_MULTIPLES[:4]
    c_partition = Counter(
        sum(5040 < (j * c) % modulus < 6720 for c in c_block)
        for j in range(modulus)
        if j % 7 == 0 and j % 14
    )
    require(d_partition == Counter({1: 10080}), "D block is not a partition")
    require(c_partition == Counter({1: 840}), "C block is not a partition")

    grid_survivors = []
    anchor_blocked = 0
    odd_blocked = 0
    for j in range(modulus):
        numerator = modulus // 2 + j
        distances = tuple(circular_numerator(v, numerator, modulus) for v in S)
        if 14 * min(distances) >= modulus:
            grid_survivors.append(j)
        if j % 14 == 0:
            require(distances[0] == 0, "anchor failed to block j=0 mod 14")
            anchor_blocked += 1
        else:
            require(14 * min(distances[1:]) < modulus, "odd grid blocker missing")
            odd_blocked += 1
    require(not grid_survivors, "integer half-turn grid has a witness")

    # Direct unit-bank replay at every numerator, with the projective sign
    # classes independently checked when the anchor does not kill the bank.
    bank_rows = []
    for p in range(8, 15):
        bank_modulus = 2 * p
        units = tuple(a for a in range(1, bank_modulus) if gcd(a, bank_modulus) == 1)
        survivors = tuple(
            a
            for a in units
            if all(
                14 * circular_numerator(v, a, bank_modulus) >= bank_modulus
                for v in S
            )
        )
        require(not survivors, f"p={p} unit bank has a witness")
        represented = {
            min(w % bank_modulus, bank_modulus - (w % bank_modulus))
            for w in W
            if gcd(w, bank_modulus) == 1
        }
        target = {
            min(a, bank_modulus - a)
            for a in units
        }
        anchor_killed = H % p == 0
        if not anchor_killed:
            require(represented == target, f"p={p} sign classes are incomplete")
        bank_rows.append((p, len(units), len(represented), anchor_killed))

    # Every represented adaptive divisor scale, using independent divisor
    # constructors.  Strict inequality is required by MC7, so equality fails.
    divisors_trial = {
        d for v in S for d in trial_divisors(v) if d >= 2
    }
    divisors_product = {
        d for v in S for d in product_divisors(v) if d >= 2
    }
    require(divisors_trial == divisors_product, "divisor constructors disagree")
    capacities = {a: mc7_capacity(a) for a in sorted(divisors_trial)}
    capacity_minimum = min(capacities.values())
    capacity_minimizers = tuple(
        a for a, value in capacities.items() if value == capacity_minimum
    )
    capacity_closures = tuple(a for a, value in capacities.items() if value < 1)
    require(len(capacities) == 78, "represented-divisor count changed")
    require(not capacity_closures, "adaptive capacity closes the row")
    require(capacity_minimum == 1, "adaptive capacity minimum changed")
    require(capacity_minimizers == (7, 21), "capacity minimizers changed")

    # THM-771 endpoint event at theta=1/22.  The six core distances are all
    # closed-safe, while owner 11 is on the half-integer boundary.
    core_numerators = tuple(
        circular_numerator(c, ENDPOINT_NUMERATOR, ENDPOINT_DENOMINATOR)
        for c in core
    )
    require(core_numerators == (10, 9, 3, 3, 5, 11), "core endpoint vector")
    require(
        all(14 * d >= ENDPOINT_DENOMINATOR for d in core_numerators),
        "endpoint is not core-safe",
    )
    require(
        ENDPOINT_OWNER * ENDPOINT_NUMERATOR * 2 == ENDPOINT_DENOMINATOR,
        "selected owner is not at a half-integer endpoint",
    )

    danger = {w: danger_sheets(w) for w in tails}
    require(danger[ENDPOINT_OWNER] == (), "event owner still has a bad sheet")
    require(
        all(len(danger[w]) == 1 for w in tails if w != ENDPOINT_OWNER),
        "a non-event owner does not have one bad sheet",
    )
    expected_danger = {
        11: (),
        1691: (0,),
        3371: (2,),
        5051: (2,),
        6731: (4,),
        8411: (6,),
        10091: (6,),
    }
    require(danger == expected_danger, "owner-location multiset changed")
    require(tight_sheets(ENDPOINT_OWNER) == (0, 5), "two-sheet endpoint pair")

    owner_multiplicity = Counter(k for sheets in danger.values() for k in sheets)
    capacity_slack = sum(1 - len(sheets) for sheets in danger.values())
    overlap_debt = sum(max(multiplicity - 1, 0) for multiplicity in owner_multiplicity.values())
    free_sheets = tuple(k for k in range(OUTER_SCALE) if owner_multiplicity[k] == 0)
    require(capacity_slack == 1, "THM-771 capacity slack changed")
    require(overlap_debt == 2, "THM-771 overlap debt changed")
    require(free_sheets == (1, 3, 5), "free-sheet set changed")
    require(
        len(free_sheets) == capacity_slack + overlap_debt,
        "THM-771 F=Q+Omega identity failed",
    )

    safe_lifts = []
    lift_clearances = []
    for k in range(OUTER_SCALE):
        numerator = ENDPOINT_NUMERATOR + ENDPOINT_DENOMINATOR * k
        clearance = min(
            circular_numerator(v, numerator, LIFT_DENOMINATOR) for v in S
        )
        if 14 * clearance >= LIFT_DENOMINATOR:
            safe_lifts.append(k)
            lift_clearances.append((k, clearance))
    require(tuple(safe_lifts) == free_sheets, "direct safe lifts differ from owner deck")
    require(tuple(lift_clearances) == ((1, 15), (3, 21), (5, 11)), "lift margins")

    # The coarse metric corollary is also available, but it is inherited
    # verbatim from THM-771 rather than a new theorem.
    core_radius = max(core)
    reduced_winding = max(tails)
    require(reduced_winding >= 7 * core_radius, "THM-771 scale cutoff failed")

    print("LRC14 r7=7 MINORITY ENDPOINT PIERCE HOSTILE -- FINITE-EXACT")
    print(f"h={H} anchor={ANCHOR} primitive={gcd(*S)} missing_thm366={missing_gates}")
    print("W=" + ",".join(map(str, W)))
    print("core=" + ",".join(map(str, core)) + " tails=" + ",".join(map(str, tails)))
    print(
        f"thm771_metric=max_tail:{reduced_winding} seven_max_core:{7 * core_radius} "
        f"passes={reduced_winding >= 7 * core_radius}"
    )
    print(
        f"half_turn_grid={modulus} anchor_blocked={anchor_blocked} "
        f"odd_blocked={odd_blocked} survivors={len(grid_survivors)}"
    )
    print(f"D_partition={dict(d_partition)} C_partition={dict(c_partition)}")
    for p, clock_count, class_count, anchor_killed in bank_rows:
        print(
            f"unit_bank_p={p} clocks={clock_count} represented_sign_classes={class_count} "
            f"anchor_killed={anchor_killed} survivors=0"
        )
    print(
        f"mc7_scales={len(capacities)} closures={len(capacity_closures)} "
        f"minimum={capacity_minimum} minimizers={','.join(map(str, capacity_minimizers))}"
    )
    print(
        f"endpoint_theta={ENDPOINT_NUMERATOR}/{ENDPOINT_DENOMINATOR} "
        f"core_numerators={','.join(map(str, core_numerators))} "
        f"core_min={min(core_numerators)}/{ENDPOINT_DENOMINATOR}"
    )
    print(
        "owner_danger="
        + ";".join(
            f"{w}:{','.join(map(str, danger[w])) if danger[w] else '-'}" for w in tails
        )
    )
    print(
        f"endpoint_owner={ENDPOINT_OWNER} tight_sheets="
        f"{','.join(map(str, tight_sheets(ENDPOINT_OWNER)))}"
    )
    print(
        f"thm771_defect=F:{len(free_sheets)} Q:{capacity_slack} "
        f"Omega:{overlap_debt} sigma:0 free_sheets={','.join(map(str, free_sheets))}"
    )
    print(
        "safe_lifts="
        + ",".join(
            f"k{k}:t={(ENDPOINT_NUMERATOR + ENDPOINT_DENOMINATOR * k)}/{LIFT_DENOMINATOR}:"
            f"clearance={clearance}/{LIFT_DENOMINATOR}"
            for k, clearance in lift_clearances
        )
    )
    print("STATUS FINITE-EXACT METHOD HOSTILE; SAFETY INHERITED FROM THM-771")
    print("RESULT PASS")


if __name__ == "__main__":
    main()
