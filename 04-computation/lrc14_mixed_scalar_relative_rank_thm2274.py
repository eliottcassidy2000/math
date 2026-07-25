"""Exact certificate for THM-2274.

The script verifies two independent branches.

1. A primitive mixed-scalar relation of support at least three is retained
   as a whole Fourier face.  The signed Selberg tensor is positive on that
   face at degree 35 for every support type.
2. A support-two relation is separated by a positive crossing polynomial on
   the adaptive cut consisting of its two coordinates.

All load-bearing arithmetic uses Fraction and explicit exceptions so that
``python`` and ``python -O`` execute the same validity checks.
"""

from fractions import Fraction as F
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative_packet_bound(height: int, support: int, guard_in_support: int) -> F:
    """Coefficient-free lower bound on the known-relation face."""

    require(height >= 1, "height must be positive")
    require(3 <= support <= 9, "support outside theorem range")
    require(guard_in_support in (0, 1), "bad guard flag")
    require(support - guard_in_support <= 8, "too many ordinary coordinates")

    epsilon = F(1, height + 1)
    defect = 2 * epsilon
    u_guard = F(5, 7) + epsilon
    u_ordinary = F(6, 7) + epsilon

    zero_mode = (
        u_guard
        * u_ordinary**8
        * (1 - defect / u_guard - 8 * defect / u_ordinary)
    )

    zero_coordinate_product = (
        u_guard ** (1 - guard_in_support)
        * u_ordinary ** (8 - support + guard_in_support)
    )

    tail = F(0)
    for ell in range(1, height + 1):
        active = F(5, 16 * ell) + epsilon
        coefficient = active**support * zero_coordinate_product
        coefficient += (
            support
            * defect
            * active ** (support - 1)
            * zero_coordinate_product
        )
        coefficient += (
            defect
            * active**support
            * zero_coordinate_product
            * (
                F(1 - guard_in_support, 1) / u_guard
                + F(8 - support + guard_in_support, 1) / u_ordinary
            )
        )
        tail += coefficient

    return zero_mode - 2 * tail


def packet_spectrum(height: int) -> list[tuple[F, int, int]]:
    rows = []
    for support in range(3, 10):
        for guard_in_support in (0, 1):
            if support - guard_in_support <= 8:
                rows.append(
                    (
                        relative_packet_bound(
                            height, support, guard_in_support
                        ),
                        support,
                        guard_in_support,
                    )
                )
    return sorted(rows)


def crossing_margin(N: int, floor_a: F, floor_b: F) -> F:
    """Two-factor product margin with 2+7 Jackson telescope errors."""

    error_a = F(3, N)
    error_b = F(21, 2 * N)
    return (
        (floor_a - error_a) * (floor_b - error_b)
        - error_a
        - error_b
    )


rows_34 = packet_spectrum(34)
rows_35 = packet_spectrum(35)

require(len(rows_34) == 13 and len(rows_35) == 13, "support-type count")
require(rows_34[0][1:] == (3, 1), "degree-34 worst type")
require(rows_35[0][1:] == (3, 1), "degree-35 worst type")
require(rows_34[0][0] < 0, "degree 34 must fail this majorant")
require(all(value > 0 for value, _, _ in rows_35), "degree 35 positivity")

expected_34 = F(
    -2764922237396438494134094765517911935722912519,
    3044225083073044036025020827143413039104000000000,
)
expected_35 = F(
    10912708836373489079295440740131626642463500500311420793,
    1804923419549521964583407381958492581322740462946091008000,
)
require(rows_34[0][0] == expected_34, "degree-34 exact fraction")
require(rows_35[0][0] == expected_35, "degree-35 exact fraction")
require(expected_35 > F(1, 166), "advertised degree-35 floor")

guard_pair_floor = F(4, 7)
seven_ordinary_floor = F(15, 154)
ordinary_pair_floor = F(66, 91)
guard_six_floor = F(191, 6930)

guard_pair_fail = crossing_margin(
    354, guard_pair_floor, seven_ordinary_floor
)
guard_pair_pass = crossing_margin(
    355, guard_pair_floor, seven_ordinary_floor
)
ordinary_pair_fail = crossing_margin(
    1058, ordinary_pair_floor, guard_six_floor
)
ordinary_pair_pass = crossing_margin(
    1059, ordinary_pair_floor, guard_six_floor
)

require(guard_pair_fail == F(-3, 15010072), "N=354 exact margin")
require(
    guard_pair_pass == F(21177, 135854950), "N=355 exact margin"
)
require(
    ordinary_pair_fail == F(-861505, 47060301288),
    "N=1058 exact margin",
)
require(
    ordinary_pair_pass == F(44021, 78582173670),
    "N=1059 exact margin",
)
require(guard_pair_fail < 0 < guard_pair_pass, "guard-pair boundary")
require(
    ordinary_pair_fail < 0 < ordinary_pair_pass,
    "ordinary-pair boundary",
)

guard_pair_types = 8
ordinary_pair_types = comb(8, 2)
require(
    guard_pair_types + ordinary_pair_types == comb(9, 2) == 36,
    "complete scalar pair split",
)

height_guard_pair = 2 * 355 - 2
height_ordinary_pair = 2 * 1059 - 2
scalar_rank_height = max(35, height_guard_pair, height_ordinary_pair)
original_rank_height = 2 * scalar_rank_height

require(height_guard_pair == 708, "guard-pair crossing height")
require(height_ordinary_pair == 2116, "ordinary-pair crossing height")
require(scalar_rank_height == 2116, "uniform scalar rank height")
require(original_rank_height == 4232, "fixed-section lift height")

print("THM-2274 MIXED SCALAR RELATIVE-RANK AUDIT")
print(
    "relative_packet_types=13 "
    f"H34_min={expected_34} H35_min={expected_35}"
)
print(
    "relative_packet_worst=(support=3,guard_in_support=1) "
    "H35_floor_gt=1/166"
)
print(
    "guard_pair_cut: types=8 "
    f"N354={guard_pair_fail} N355={guard_pair_pass} height=708"
)
print(
    "ordinary_pair_cut: types=28 "
    f"N1058={ordinary_pair_fail} N1059={ordinary_pair_pass} height=2116"
)
print(
    "pair_partition=36/36 scalar_rank_height=2116 "
    "fixed_section_original_height=4232"
)
print("ALL EXACT CHECKS PASSED")
