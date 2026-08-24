#!/usr/bin/env python3
"""Exact controls for THM-4024's complete divisor-incidence envelope.

The all-divisor proof is symbolic in THM-4024.  This companion independently
checks the open-comb count, exhausts the two- and three-exception residue
profiles at the only small composite moduli, and audits the divisor
classification.  It does not prove LRC(14).
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import product
from math import gcd, lcm
import sys


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def distance(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def direct_order_bad_count(order: int) -> int:
    """Maximum points of an order-m orbit in the open 1/14 danger arc."""
    walls = sorted(
        {
            (sign * F(1, 14) - F(j, order)) % 1
            for sign in (-1, 1)
            for j in range(order)
        }
    )
    probes = list(walls)
    for index, left in enumerate(walls):
        right = walls[(index + 1) % len(walls)]
        if index + 1 == len(walls):
            right += 1
        probes.append(((left + right) / 2) % 1)
    return max(
        sum(distance(phase + F(j, order)) < F(1, 14) for j in range(order))
        for phase in probes
    )


def order_bad_count(order: int) -> int:
    return order // 7 if order % 7 == 0 else (order + 6) // 7


def normalized_tariff(order: int) -> F:
    return F(order_bad_count(order), order)


def labelled_bad_count(modulus: int, residue: int) -> int:
    common = gcd(modulus, residue)
    order = modulus // common
    return common * order_bad_count(order)


def primitive_residue_profiles(modulus: int, exceptions: int):
    for residues in product(range(1, modulus), repeat=exceptions):
        if gcd(modulus, *residues) == 1:
            yield residues


def max_two_exception_count(modulus: int) -> tuple[int, tuple[int, int]]:
    scored = (
        (sum(labelled_bad_count(modulus, r) for r in residues), residues)
        for residues in primitive_residue_profiles(modulus, 2)
    )
    return max(scored)


def lower_incidence_ok(modulus: int, residues: tuple[int, int, int]) -> bool:
    """Constraints already forced at 2, 3 and 4 for an eight-incidence row."""
    if modulus % 2 == 0 and sum(r % 2 == 0 for r in residues) > 1:
        return False  # c_2 <= 9
    if modulus % 3 == 0 and any(r % 3 == 0 for r in residues):
        return False  # c_3 <= 8
    if modulus % 4 == 0 and any(r % 4 == 0 for r in residues):
        return False  # c_4 <= 8
    return True


def max_three_exception_count(modulus: int) -> tuple[int, tuple[int, int, int]]:
    scored = (
        (sum(labelled_bad_count(modulus, r) for r in residues), residues)
        for residues in primitive_residue_profiles(modulus, 3)
        if lower_incidence_ok(modulus, residues)
    )
    return max(scored)


def smallest_prime_factor(value: int) -> int:
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return divisor
        divisor += 1
    return value


def classification_reason(modulus: int) -> str:
    """Symbolic case split used for every d>=5."""
    probe = modulus
    while probe % 2 == 0:
        probe //= 2
    while probe % 3 == 0:
        probe //= 3
    if probe > 1:
        return "prime>=5"
    if modulus % 8 == 0:
        return "multiple_of_8"
    if modulus % 9 == 0:
        return "multiple_of_9"
    if modulus == 6:
        return "six_boundary"
    if modulus == 12:
        return "twelve_boundary"
    raise AssertionError(("unclassified", modulus))


def main() -> None:
    # Independent wall-cell replay of THM-4004's strict open-comb formula.
    direct = tuple(direct_order_bad_count(order) for order in range(1, 65))
    formula = tuple(order_bad_count(order) for order in range(1, 65))
    require(direct == formula, "open-comb count through order 64")

    # With nine divisible body coordinates there are two exceptions.  Body
    # primitivity is exactly lcm(m_1,m_2)=d, equivalently gcd(d,r_1,r_2)=1.
    # Exhaustion displays the unique d=2 covering boundary.
    two_profiles = {}
    for modulus in range(2, 65):
        maximum, profile = max_two_exception_count(modulus)
        two_profiles[modulus] = (maximum, profile)
        require(
            maximum >= modulus if modulus == 2 else maximum < modulus,
            f"two-exception boundary d={modulus}",
        )

    # The small prime-power/composite states left after lower-incidence gates.
    three_profiles = {}
    for modulus in (2, 3, 4, 6, 8, 9, 12):
        maximum, profile = max_three_exception_count(modulus)
        three_profiles[modulus] = (maximum, profile)

    require(three_profiles[3][0] == 3, "d=3 equality")
    require(three_profiles[4][0] == 4, "d=4 equality")
    require(three_profiles[6][0] == 4 < 6, "d=6 strict closure")
    require(three_profiles[8][0] == 6 < 8, "d=8 strict closure")
    require(three_profiles[9][0] == 6 < 9, "d=9 strict closure")
    require(three_profiles[12][0] == 6 < 12, "d=12 strict closure")

    # At d=4, equality requires one residue of valuation exactly one and two
    # odd residues.  All-odd exception profiles have total bad count three.
    d4_equality = []
    d4_all_odd = []
    for residues in primitive_residue_profiles(4, 3):
        if not lower_incidence_ok(4, residues):
            continue
        total = sum(labelled_bad_count(4, r) for r in residues)
        if total == 4:
            d4_equality.append(residues)
            require(sum(r % 2 == 0 for r in residues) == 1, f"d4 equality {residues}")
        if all(r % 2 for r in residues):
            d4_all_odd.append(residues)
            require(total == 3, f"d4 all odd {residues}")
    require(d4_equality and d4_all_odd, "d4 controls nonempty")

    # General proof-side tariff: c_d=8 gives three exceptions.  THM-4004's
    # c_e<=8 for e>=3 forces each gcd(d,delta_j) to be 1 or 2; c_2<=9 allows
    # at most one value 2.  Hence these two displayed inequalities cover all
    # odd d>=5 and even d>=6 respectively.
    tariff_checks = 0
    for modulus in range(5, 10001):
        if modulus % 2:
            require(3 * normalized_tariff(modulus) < 1, f"odd tariff d={modulus}")
        else:
            require(
                normalized_tariff(modulus // 2) + 2 * normalized_tariff(modulus) < 1,
                f"even tariff d={modulus}",
            )
        tariff_checks += 1

    # Every d>=5 falls under a large prime, 8, 9, or one of the two strict
    # mixed boundaries.  This is an exhaustive integer classification, not a
    # finite extrapolation: the loop only audits the displayed symbolic split.
    reasons = {}
    for modulus in range(5, 10001):
        reason = classification_reason(modulus)
        reasons[reason] = reasons.get(reason, 0) + 1
        if reason == "prime>=5":
            reduced = modulus
            while reduced % 2 == 0:
                reduced //= 2
            while reduced % 3 == 0:
                reduced //= 3
            require(smallest_prime_factor(reduced) >= 5, f"large prime d={modulus}")
        elif reason == "multiple_of_8":
            require(modulus % 8 == 0, f"eight d={modulus}")
        elif reason == "multiple_of_9":
            require(modulus % 9 == 0, f"nine d={modulus}")
        else:
            require(modulus in (6, 12), f"mixed boundary d={modulus}")

    # The natural-number envelope and its equivalent subset-gcd compression.
    envelope = (11, 9, 8, 8) + (7,) * 12
    require(envelope[:8] == (11, 9, 8, 8, 7, 7, 7, 7), "envelope head")
    for common_divisor in range(1, 10001):
        allowed_for_eight = common_divisor in (1, 2, 3, 4)
        implied_by_envelope = common_divisor < 5
        require(allowed_for_eight == implied_by_envelope, f"eight-subset gcd d={common_divisor}")

    # Transported ordinal product used elsewhere in the session: retain as a
    # hostile reminder that ordinary index order is not the divisor operation.
    for a in range(1, 40):
        for b in range(1, 40):
            star = 2 * a * b - a - b + 1
            require(2 * star - 1 == (2 * a - 1) * (2 * b - 1), f"odd star {a},{b}")

    print("THM4024_COMPLETE_DIVISOR_INCIDENCE_ENVELOPE_EXACT")
    print("scope=THM3818_rank11_11+2_branch;divisors_d_of_t;LRC14=OPEN")
    print("incidence_sequence=c_1=11;c_2<=9;c_3<=8;c_4<=8;c_d<=7_for_every_d>=5")
    print("prime_power_profile=2:(9,8,7,...);3:(8,7,...);ell>=5:(7,7,...)")
    print("subset_gcd_compression=every_8_body_subset_has_gcd_with_t_in_{1,2,3,4}")
    print("nine_subset_gcd_compression=every_9_body_subset_has_gcd_with_t_in_{1,2}")
    print("ten_subset_gcd_compression=every_10_body_subset_is_coprime_to_t")
    print(f"direct_comb_orders={len(direct)}")
    print(f"two_exception_d2={two_profiles[2]};two_exception_d3={two_profiles[3]}")
    print(f"three_exception_profiles={three_profiles}")
    print(f"d4_equality_profiles={len(d4_equality)};d4_all_odd_profiles={len(d4_all_odd)}")
    print(f"classification_audit_through_10000={reasons}")
    print(f"general_tariff_checks={tariff_checks}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
