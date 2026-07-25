#!/usr/bin/env python3
"""Exact arithmetic ledger for the septimal/slope addendum to THM-2193.

Reproduction:
  python3 04-computation/lrc14_septimal_torsion_slope_addendum_thm2193.py
  python3 -O 04-computation/lrc14_septimal_torsion_slope_addendum_thm2193.py
"""

from fractions import Fraction
from math import comb, lcm, prod


PROFILE = (29, 105, 178, 258, 416, 718)
PROFILE_PRODUCT = prod(PROFILE)
MINOR_BOUND = 216 * PROFILE_PRODUCT
TORSION_THRESHOLD = 1372 * MINOR_BOUND * MINOR_BOUND


def danger_residues(q: int) -> tuple[int, ...]:
    """Strict radius-1/14 residues in Z/qZ."""
    return tuple(r for r in range(q) if 14 * min(r, q - r) < q)


def check_septimal_cosets(k: int) -> None:
    q = 7**k
    danger = danger_residues(q)
    assert len(danger) == q // 7
    for s in range(k):
        counts = [0] * (7**s)
        for r in danger:
            counts[r % (7**s)] += 1
        assert len(set(counts)) == 1
        assert counts[0] == q // (7 * 7**s)


def check_lucas_shadow(k: int) -> None:
    q = 7**k
    for j in range(q):
        assert comb(q - 1, j) % 7 == (1 if j % 2 == 0 else 6)
    for j in range(1, q):
        assert comb(q, j) % 7 == 0


def main() -> None:
    assert PROFILE_PRODUCT == 41_768_105_783_040
    assert MINOR_BOUND == 9_021_910_849_136_640
    assert TORSION_THRESHOLD == (
        111_673_769_007_323_628_596_227_411_751_731_200
    )
    assert 1372 == 196 * 7

    for k in range(1, 7):
        check_septimal_cosets(k)
    for k in range(1, 5):
        check_lucas_shadow(k)

    # THM-2188 ambient product versus primitive runner slope.
    core_mass = Fraction(4319, 51480)
    ambient_mass = Fraction(6, 7) * core_mass
    endpoint_numerator = 6 * core_mass
    assert ambient_mass == Fraction(617, 8580)
    assert endpoint_numerator == Fraction(4319, 8580)
    assert ambient_mass - endpoint_numerator / 7 == 0

    bank_height = 4
    finite_bank = 25
    period = 720720
    for q in range(1, finite_bank + 1):
        period = lcm(period, q)
    period = lcm(period, 7**bank_height)
    pumped_speed = 7 + period
    pumped_mass = ambient_mass - endpoint_numerator / pumped_speed
    assert pumped_mass == Fraction(
        4319 * (pumped_speed - 7), 60060 * pumped_speed
    )
    assert pumped_mass > 0
    for q in range(1, finite_bank + 1):
        assert pumped_speed % q == 7 % q
    assert pumped_speed % (7**bank_height) == 7

    print("THM-2193 septimal/slope addendum exact referee")
    print(f"profile={PROFILE}")
    print(f"profile_product={PROFILE_PRODUCT}")
    print(f"minor_bound={MINOR_BOUND}")
    print(f"torsion_threshold={TORSION_THRESHOLD}")
    print("septimal_coset_depths=1..6 PASS")
    print("lucas_depths=1..4 PASS")
    print(f"hostile_ambient_mass={ambient_mass}")
    print(f"hostile_endpoint_numerator={endpoint_numerator}")
    print(f"hostile_test_speed={pumped_speed}")
    print(f"hostile_pumped_mass={pumped_mass}")
    print("ALL EXACT CHECKS PASS")


if __name__ == "__main__":
    main()
