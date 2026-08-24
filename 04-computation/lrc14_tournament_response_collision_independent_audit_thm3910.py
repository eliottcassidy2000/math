#!/usr/bin/env python3
"""Independent THM-3910 wall-cell audit of the AP11 cubic-response collision."""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import sys


sys.stdout.reconfigure(newline="\n")


def boundaries(n: int) -> set[Fraction]:
    modulus = 14 * n
    return {
        Fraction((14 * k + sign) % modulus, modulus)
        for k in range(n)
        for sign in (-1, 1)
    }


def is_danger(n: int, x: Fraction) -> bool:
    phase = (n * x) % 1
    return min(phase, 1 - phase) < Fraction(1, 14)


def ledger(pair_frequencies: tuple[int, int]) -> dict[tuple[int, int, int], Fraction]:
    speeds = tuple(range(1, 12)) + pair_frequencies
    walls = {Fraction(0), Fraction(1)}
    for speed in speeds:
        walls.update(boundaries(speed))
    ordered = sorted(walls)
    answer: dict[tuple[int, int, int], Fraction] = defaultdict(Fraction)
    for left, right in zip(ordered, ordered[1:]):
        if left == right:
            continue
        middle = (left + right) / 2
        g = int(all(not is_danger(speed, middle) for speed in range(1, 12)))
        a = int(is_danger(pair_frequencies[0], middle))
        b = int(is_danger(pair_frequencies[1], middle))
        answer[(g, a, b)] += right - left
    if sum(answer.values(), Fraction(0)) != 1:
        raise RuntimeError("atom ledger does not have total mass one")
    return dict(answer)


def moments(atom: dict[tuple[int, int, int], Fraction]) -> tuple[Fraction, ...]:
    m0 = sum(mass * bits[0] for bits, mass in atom.items())
    ma = sum(mass * bits[1] for bits, mass in atom.items())
    mb = sum(mass * bits[2] for bits, mass in atom.items())
    m0a = sum(mass * bits[0] * bits[1] for bits, mass in atom.items())
    m0b = sum(mass * bits[0] * bits[2] for bits, mass in atom.items())
    mab = sum(mass * bits[1] * bits[2] for bits, mass in atom.items())
    m0ab = sum(mass * bits[0] * bits[1] * bits[2] for bits, mass in atom.items())
    target = atom.get((1, 0, 0), Fraction(0))
    return m0, ma, mb, m0a, m0b, mab, m0ab, target


def fmt(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    first = moments(ledger((72, 492)))
    second = moments(ledger((216, 1476)))
    expected_common = (
        Fraction(10931, 194040),
        Fraction(1, 7),
        Fraction(1, 7),
        Fraction(11, 1176),
        Fraction(1, 123),
        Fraction(5, 246),
    )
    if first[:6] != expected_common or second[:6] != expected_common:
        raise RuntimeError("one/two response mismatch")
    if first[6:] != (Fraction(1, 861), Fraction(79579, 1988910)):
        raise RuntimeError("first cubic/target response mismatch")
    if second[6:] != (Fraction(89, 61992), Fraction(961493, 23866920)):
        raise RuntimeError("second cubic/target response mismatch")
    if second[6] - first[6] != Fraction(17, 61992):
        raise RuntimeError("cubic split mismatch")
    if second[7] - first[7] != Fraction(17, 61992):
        raise RuntimeError("target split mismatch")

    print("LRC14_TOURNAMENT_RESPONSE_COLLISION_INDEPENDENT_AUDIT_20260823")
    print("method=sorted_rational_wall_cells+direct_midpoint_bits;no_interval_union_or_intersection")
    print("common_one_two=" + ",".join(fmt(x) for x in expected_common))
    print(f"cubic={fmt(first[6])},{fmt(second[6])};difference={fmt(second[6]-first[6])}")
    print(f"target={fmt(first[7])},{fmt(second[7])};difference={fmt(second[7]-first[7])}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
