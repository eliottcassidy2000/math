#!/usr/bin/env python3
"""Exact hostile controls for THM-2167.

The theorem is symbolic.  This referee checks its mixed-height prime bound,
sign-sharp state constants, full binary rank-two digit map, simultaneous
carries, and a pump which preserves two relations and positivity while
destroying distinctness.
"""

from collections import Counter
from itertools import product
from math import prod


PRIMES = (2, 3, 5, 7, 11, 13)
ANCHOR_HEIGHT = 29
COMPANION_HEIGHT = 105


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def dot(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x * y for x, y in zip(a, b))


def independent_mod_q(
    r: tuple[int, ...], s: tuple[int, ...], q: int
) -> bool:
    return any(
        (r[i] * s[j] - r[j] * s[i]) % q
        for i in range(len(r))
        for j in range(i + 1, len(r))
    )


def carry_path(
    row: tuple[int, ...], relation: tuple[int, ...], q: int
) -> tuple[int, ...]:
    path = []
    power = 1
    while True:
        quotient = tuple(v // power for v in row)
        path.append(-dot(relation, quotient))
        if not any(quotient):
            return tuple(path)
        power *= q


def delete_digit_block(
    row: tuple[int, ...], q: int, j: int, k: int
) -> tuple[int, ...]:
    """Delete base-q positions j,...,k-1 and concatenate."""
    low_modulus = q**j
    return tuple(v % low_modulus + low_modulus * (v // (q**k)) for v in row)


def owner(row: tuple[int, ...], q: int, j: int) -> tuple[int, ...]:
    return tuple(i for i, v in enumerate(row) if v >= q**j)


def main() -> None:
    prime_product = prod(PRIMES)
    minor_bound = 2 * ANCHOR_HEIGHT * COMPANION_HEIGHT
    require(minor_bound < prime_product, "minor bound does not force a good prime")

    sharp_r = (22, 0)
    sharp_s = (0, 105)
    sharp_det = 22 * 105
    sharp_profile = tuple(
        (q, independent_mod_q(sharp_r, sharp_s, q)) for q in PRIMES
    )
    require(sharp_det == 2310, "wrong endpoint determinant")
    require(
        sharp_profile
        == (
            (2, False),
            (3, False),
            (5, False),
            (7, False),
            (11, False),
            (13, True),
        ),
        "endpoint prime profile changed",
    )

    anchor_l1 = 12 * 29 + 28
    height43_l1 = 12 * 43 + 42
    height105_l1 = 12 * 105 + 104
    support_ge_three_pairs = (anchor_l1 - 1) * (height43_l1 - 1)
    support_two_pairs = (29 + 28 - 1) * (height105_l1 - 1)
    owner_states = 14 * support_ge_three_pairs
    require(support_ge_three_pairs == 208875, "support>=3 state count changed")
    require(support_two_pairs == 76328, "support-two state count changed")
    require(owner_states == 2924250, "carry-owner state count changed")

    ap = tuple(range(1, 14))
    r = (1, 1, -1) + (0,) * 10
    s = (1, 0, 1, -1) + (0,) * 9
    require(dot(r, ap) == dot(s, ap) == 0, "AP controls are not relations")
    require(independent_mod_q(r, s, 2), "AP controls lost binary rank")

    fibres = Counter()
    for digit in product(range(2), repeat=13):
        fibres[(dot(r, digit) % 2, dot(s, digit) % 2)] += 1
    require(
        fibres == Counter({key: 2**11 for key in product(range(2), repeat=2)}),
        "binary affine fibre count changed",
    )

    carry_r = carry_path(ap, r, 2)
    carry_s = carry_path(ap, s, 2)
    require(carry_r[0] == carry_s[0] == 0, "initial carries are nonzero")
    require(carry_r[-1] == carry_s[-1] == 0, "terminal carries are nonzero")
    for relation, path in ((r, carry_r), (s, carry_s)):
        power = 1
        for j in range(len(path) - 1):
            digit = tuple((v // power) % 2 for v in ap)
            require(
                2 * path[j + 1] == path[j] + dot(relation, digit),
                "binary carry recurrence failed",
            )
            power *= 2

    pump_row = (13, 16, 1, 2, 9, 10, 11, 12, 14, 15, 17, 18, 19)
    pump_r = (-2, 1, 0, 0, 0, 1) + (0,) * 7
    pump_s = (0, 0, 0, 0, 1, 0, 0, -2, 0, 1) + (0,) * 3
    require(len(set(pump_row)) == 13, "pump source is not distinct")
    require(independent_mod_q(pump_r, pump_s, 3), "pump relations lost ternary rank")
    require(
        dot(pump_r, pump_row) == dot(pump_s, pump_row) == 0,
        "pump source relations fail",
    )

    pumped = delete_digit_block(pump_row, 3, 1, 2)
    require(min(pumped) > 0, "pump destroyed positivity")
    require(
        owner(pump_row, 3, 1) == owner(pump_row, 3, 2),
        "pump endpoints have different owner masks",
    )
    pump_carry_r = carry_path(pump_row, pump_r, 3)
    pump_carry_s = carry_path(pump_row, pump_s, 3)
    require(pump_carry_r[1] == pump_carry_r[2], "first pump carry does not close")
    require(pump_carry_s[1] == pump_carry_s[2], "second pump carry does not close")
    require(
        dot(pump_r, pumped) == dot(pump_s, pumped) == 0,
        "pumped row lost a relation",
    )
    require(len(set(pumped)) < 13, "pump did not destroy distinctness")
    require(pumped[0] == pumped[1], "designated coordinates did not merge")

    print("THM-2167 RANK-TWO ADAPTIVE CARRY REFEREE")
    print(f"primes={PRIMES}")
    print(f"prime_product={prime_product}")
    print(f"mixed_heights={(ANCHOR_HEIGHT, COMPANION_HEIGHT)}")
    print(f"minor_bound={minor_bound}")
    print(f"abstract_endpoint_determinant={sharp_det}")
    print(f"abstract_endpoint_independence={sharp_profile}")
    print(f"support_ge_three_carry_pairs={support_ge_three_pairs}")
    print(f"support_two_carry_pairs={support_two_pairs}")
    print(f"sorted_carry_owner_states={owner_states}")
    print("binary_rank_two_digit_fibres=" + str(sorted(fibres.items())))
    print(f"binary_fibre_size={2**11}")
    print(f"ap_carry_r={carry_r}")
    print(f"ap_carry_s={carry_s}")
    print(f"pump_owner={owner(pump_row, 3, 1)}")
    print(f"pump_carry_r={pump_carry_r}")
    print(f"pump_carry_s={pump_carry_s}")
    print(f"pump_row={pump_row}")
    print(f"pumped_row={pumped}")
    print("pump_preserves=positivity,two_relations")
    print("pump_destroys=distinctness")
    print("status=PASS")


if __name__ == "__main__":
    main()
