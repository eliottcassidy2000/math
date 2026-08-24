#!/usr/bin/env python3
"""Independent hostile audit of THM-4024.

This path works with gcd/order profiles rather than owner residues.  It also
constructs a primitive distinct d=4 selector hostile and a deliberately
nonprimitive d=6 equality control.  It does not import the primary companion.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def divisors(value: int) -> tuple[int, ...]:
    return tuple(d for d in range(1, value + 1) if value % d == 0)


def comb_upper_fraction(order: int) -> F:
    """Normalized open-comb count, with a safe ceiling at multiples of seven."""
    return F((order + 6) // 7, order)


def inherited_profile_ok(modulus: int, gs: tuple[int, int, int]) -> bool:
    """Apply THM-4004's c_2<=9 and c_e<=8 (e>=3) at proper divisors."""
    for divisor in divisors(modulus):
        if divisor in (1, modulus):
            continue
        incidence = 8 + sum(g % divisor == 0 for g in gs)
        cap = 9 if divisor == 2 else 8
        if incidence > cap:
            return False
    return True


def distance(speed: int, phase: F) -> F:
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def branch_bad_sets(
    modulus: int, y: F, exceptions: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            label
            for label in range(modulus)
            if distance(delta, (y + label) / modulus) < F(1, 14)
        )
        for delta in exceptions
    )


def full_clearance(speeds: tuple[int, ...], phase: F) -> F:
    return min(distance(speed, phase) for speed in speeds)


def main() -> None:
    require(
        all(comb_upper_fraction(order) <= F(1, 3) for order in range(3, 10001)),
        "tariff <=1/3 for order >=3",
    )
    require(
        all(comb_upper_fraction(order) <= F(1, 4) for order in range(5, 10001)),
        "tariff <=1/4 for order >=5",
    )

    unresolved: dict[int, list[tuple[tuple[int, int, int], tuple[int, int, int], F]]] = {}
    for modulus in range(3, 501):
        proper = tuple(g for g in divisors(modulus) if g < modulus)
        for gs in combinations_with_replacement(proper, 3):
            if not inherited_profile_ok(modulus, gs):
                continue
            orders = tuple(modulus // g for g in gs)
            tariff = sum((comb_upper_fraction(order) for order in orders), F(0))
            if tariff >= 1:
                unresolved.setdefault(modulus, []).append((gs, orders, tariff))

    require(tuple(unresolved) == (3, 4), "only d=3,4 escape the strict tariff")
    require(
        unresolved[3] == [((1, 1, 1), (3, 3, 3), F(1))],
        "d=3 equality profile",
    )
    require(
        unresolved[4] == [((1, 1, 2), (4, 4, 2), F(1))],
        "d=4 equality profile",
    )

    # Primitive, distinct, syntactically typed d=4 control.  The divided pack
    # is H={1,...,10}; its y=1/11 lifts are all spoiled, but the row itself is
    # 1/11-safe elsewhere.  Hence this is a selector hostile, not an LRC
    # counterexample.
    d4_y = F(1, 11)
    d4_exceptions = (2, 9, 11)
    d4_sets = branch_bad_sets(4, d4_y, d4_exceptions)
    require(tuple(gcd(4, value) for value in d4_exceptions) == (2, 1, 1), "d4 gcds")
    require(d4_sets == ((0, 2), (3,), (1,)), "d4 branch partition")
    require(set().union(*map(set, d4_sets)) == set(range(4)), "d4 covers labels")
    d4_body = (8, 12, 20, 24, 28, 32, 36, 40) + d4_exceptions
    d4_pair = (4, 16)
    d4_speeds = d4_body + d4_pair
    require(gcd(*d4_body) == 1, "d4 primitive body")
    require(len(set(d4_speeds)) == 13, "d4 distinct row")
    require(full_clearance(d4_speeds, F(21, 22)) == F(1, 11), "d4 safe control")

    # The tempting d=6 order-(3,3,3) equality is geometrically real, but its
    # three gcd-two exceptions make the whole body even.  This is the exact
    # lower-incidence obstruction used by the proof.
    d6_y = F(1, 11)
    d6_exceptions = (2, 22, 46)
    d6_sets = branch_bad_sets(6, d6_y, d6_exceptions)
    require(tuple(gcd(6, value) for value in d6_exceptions) == (2, 2, 2), "d6 gcds")
    require(d6_sets == ((0, 3), (1, 4), (2, 5)), "d6 branch partition")
    d6_body = (12, 18, 30, 36, 42, 48, 54, 60) + d6_exceptions
    require(gcd(*d6_body) == 2, "d6 equality nonprimitive")

    # Exact prime-power consequence, audited well beyond the theorem's first
    # relevant levels without treating this finite loop as its proof.
    prime_power_checks = 0
    for exponent in range(3, 13):
        require(2**exponent not in unresolved, f"2^{exponent}")
        prime_power_checks += 1
    for exponent in range(2, 9):
        require(3**exponent not in unresolved, f"3^{exponent}")
        prime_power_checks += 1
    for prime in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        for exponent in range(1, 5):
            if prime**exponent <= 500:
                require(prime**exponent not in unresolved, f"{prime}^{exponent}")
                prime_power_checks += 1

    print("THM4024_COMPLETE_DIVISOR_INCIDENCE_INDEPENDENT_AUDIT")
    print("abstract_unresolved_moduli_through_500=3,4")
    print("d3_profile=(3,3,3);tariff=1")
    print("d4_profile=(2,4,4);tariff=1;primitive_hostile=(2,9,11)")
    print("d4_branch_sets=((0,2),(3),(1));safe_phase=21/22;clearance=1/11")
    print("d6_naive_profile=(3,3,3);branch_sets=((0,3),(1,4),(2,5));body_gcd=2")
    print("proved_profile=c2<=9;c3<=8;c4<=8;c_d<=7_for_every_d>=5_dividing_t")
    print(f"prime_power_checks={prime_power_checks}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
