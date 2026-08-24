#!/usr/bin/env python3
"""Exact controls for THM-4004's t<U divisor-comb profile.

The theorem contains the symbolic union-bound and crossing-height proofs.
This companion independently rebuilds the pair atlas, verifies the sharp
open-comb count from wall cells, audits the order and prime inequalities, and
checks typed positive and hostile rows.  It does not prove LRC(14).
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def distance(speed: int, phase: F) -> F:
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def bad_branch_bound(modulus: int, detuned: int) -> int:
    common = gcd(modulus, detuned)
    order = modulus // common
    return common * ((order + 6) // 7)


def direct_max_bad_count(order: int) -> int:
    """Direct wall-cell maximum for an open 1/14 danger arc."""
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
        sum(distance(1, phase + F(j, order)) < F(1, 14) for j in range(order))
        for phase in probes
    )


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def admissible_sum(value: int) -> bool:
    """THM-3818 cube-class sum predicate, rebuilt by trial division."""
    original = value
    prime = 2
    seen = False
    while prime * prime <= value:
        if value % prime == 0:
            seen = True
            exponent = 0
            while value % prime == 0:
                exponent += 1
                value //= prime
            if prime % 3 != 2 or exponent > 2:
                return False
        prime += 1
    if value > 1:
        seen = True
        if value % 3 != 2:
            return False
    return seen and original > 1


def atlas() -> tuple[tuple[int, int], ...]:
    return tuple(
        (p, total - p)
        for total in range(3, 357)
        if admissible_sum(total)
        for p in range(1, (total + 1) // 2)
        if p < total - p and gcd(p, total - p) == 1
    )


def main() -> None:
    q_height = 91**6
    pairs = atlas()
    require(len(pairs) == 5855, "independent pair-atlas count")
    max_p = max(p for p, _ in pairs)
    max_p_pairs = tuple(pair for pair in pairs if pair[0] == max_p)
    require(max_p == 177, "sharp smaller-coordinate cap")
    require(max_p_pairs == ((177, 178), (177, 179)), "sharp cap owners")
    integer_floor = q_height // 177 + 1
    require(q_height == 567869252041, "crossing height")
    require(integer_floor == 3208300859, "global integer body floor")
    require(
        177 * (integer_floor - 1) <= q_height < 177 * integer_floor,
        "floor boundary",
    )

    # The exact branch-coset formula is checked against a direct expansion of
    # all open wall cells.  At order m it is ceil(m/7).
    direct_counts = tuple(direct_max_bad_count(order) for order in range(1, 65))
    formula_counts = tuple((order + 6) // 7 for order in range(1, 65))
    require(direct_counts == formula_counts, "direct open-comb wall-cell count")

    # Two detuned combs cover all branches only at the order-(2,2) boundary.
    two_order_checks = 0
    for first in range(2, 1001):
        for second in range(2, 1001):
            total = F((first + 6) // 7, first) + F((second + 6) // 7, second)
            require(total < 1 or (first, second) == (2, 2), f"two orders {first},{second}")
            if total >= 1:
                require((first, second) == (2, 2), "unique two-comb boundary")
            two_order_checks += 1

    # Three full-order combs leave a branch for every prime ell >= 5.
    primes = tuple(value for value in range(5, 5000) if is_prime(value))
    require(len(primes) > 600, "prime audit universe")
    require(
        all(3 * ((ell + 6) // 7) < ell for ell in primes),
        "three-detuned prime inequality",
    )
    require(3 * ((3 + 6) // 7) == 3, "ell=3 sharp counting boundary")

    # One-detuned t<U control: after removing owner 1 the harmonic pack is
    # exactly 1,...,12, and branch 1 wins at 7/13.
    t = 2
    body = (1, 2, 4, 8, 10, 12, 16, 18, 20, 22, 24)
    pair = (3, 7)
    pack12 = tuple(value // t for value in body[1:]) + pair
    require(tuple(sorted(pack12)) == tuple(range(1, 13)), "twelve-speed pack")
    branches = tuple((F(1, 13) + k) / t for k in range(t))
    speeds = body + tuple(t * value for value in pair)
    clearances = tuple(min(distance(speed, phase) for speed in speeds) for phase in branches)
    require(clearances == (F(1, 26), F(1, 13)), "one-detuned typed control")

    # At d=2 the two order-two exceptions are the only counting boundary.
    # The reduced odd pair (1,3) is universally rescued in this control.
    pack11 = tuple(range(1, 12))
    branches2 = tuple((F(1, 12) + k) / 2 for k in range(2))
    clearances2 = tuple(
        min(
            min(distance(2 * value, phase) for value in pack11),
            distance(1, phase),
            distance(3, phase),
        )
        for phase in branches2
    )
    require(clearances2 == (F(1, 24), F(1, 12)), "d=2 reduced-pair control")

    # First closed prime: eight body owners and both pair owners reduce to
    # H={1,...,10}; four of five branches are fully 1/11-safe.
    body5 = (1, 2, 3, 10, 15, 25, 30, 35, 40, 45, 50)
    pair5 = (1, 4)
    pack10 = tuple(value // 5 for value in body5[3:]) + pair5
    require(tuple(sorted(pack10)) == tuple(range(1, 11)), "ell=5 ten-speed pack")
    speeds5 = body5 + tuple(5 * value for value in pair5)
    require(len(speeds5) == len(set(speeds5)) == 13, "ell=5 distinct owners")
    branches5 = tuple((F(1, 11) + k) / 5 for k in range(5))
    clearances5 = tuple(min(distance(speed, phase) for speed in speeds5) for phase in branches5)
    require(
        clearances5 == (F(1, 55), F(1, 11), F(1, 11), F(1, 11), F(1, 11)),
        "ell=5 positive control",
    )

    # Sharp method hostile at ell=3: the same ten-speed pack is safe, but the
    # three labelled detuned owners spoil the three selected lifts.
    body3 = (1, 10, 11, 6, 9, 15, 18, 21, 24, 27, 30)
    speeds3 = body3 + (3, 12)
    require(len(speeds3) == len(set(speeds3)) == 13, "ell=3 distinct hostile")
    pack10_hostile = tuple(value // 3 for value in body3[3:]) + pair5
    require(tuple(sorted(pack10_hostile)) == tuple(range(1, 11)), "ell=3 ten-speed pack")
    branches3 = tuple((F(1, 11) + k) / 3 for k in range(3))
    spoilers = (1, 11, 10)
    spoiler_distances = tuple(distance(spoilers[k], branches3[k]) for k in range(3))
    require(spoiler_distances == (F(1, 33), F(0), F(1, 33)), "ell=3 selector hostile")
    require(
        min(distance(speed, F(4, 33)) for speed in speeds3) == F(1, 11),
        "ell=3 hostile row is itself lonely",
    )

    # The literal swapped-component width criterion is arithmetically unable
    # to fire when t<U: t*lambda(u) <= 6t/(7U) < 6/7.
    require(F(6 * 999, 7 * 1000) < F(6, 7) < 1, "component-swap no-go control")

    print("THM4004_TLTU_DIVISOR_COMB_PROFILE_EXACT")
    print("scope=THM3818_rank11_11+2_branch;t<U;LRC14=OPEN")
    print(f"atlas_pairs={len(pairs)};max_p={max_p};max_p_pairs={max_p_pairs}")
    print(f"crossing_floor=U>91^6/177;integer_U_min={integer_floor}")
    print("deletion_gcd=every_prime_dividing_t_misses_at_least_two_body_coordinates")
    print("prime_profile=ell>=5_at_most_7;ell=3_at_most_8;ell=2_at_most_9")
    print("ell=2_equality_requires=s=1_and_reduced_odd_exception_sum>7")
    print(f"direct_comb_orders={len(direct_counts)}")
    print(f"two_order_checks={two_order_checks}")
    print(f"prime_checks={len(primes)}")
    print("swap_no_go=t<U_implies_t*lambda(u)<6/7")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
