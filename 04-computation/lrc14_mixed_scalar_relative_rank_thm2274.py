"""Exact certificate for THM-2274.

The script verifies two independent branches.

1. A primitive mixed-scalar relation of support at least three is retained
   as a whole Fourier face.  The signed Selberg tensor is positive on that
   face at degree 35 for every support type.
2. A support-two relation is separated by a positive crossing polynomial on
   the adaptive cut consisting of its two coordinates.  Here the exact
   Jackson first-moment formula from THM-2145 replaces the preliminary
   coarse 3/(2N) error bound.

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


def jackson_coefficient(N: int, k: int) -> int:
    """Integral Fourier numerator C_k for the normalized Jackson kernel."""

    require(N >= 2, "Jackson order must be at least two")
    require(0 <= k <= 2 * N - 2, "Jackson frequency outside support")
    if k <= N:
        numerator = (
            4 * N**3
            - 6 * N * k**2
            + 2 * N
            + 3 * k**3
            - 3 * k
        )
    else:
        numerator = (2 * N - k) ** 3 - (2 * N - k)
    require(numerator % 6 == 0, "Jackson coefficient lost integrality")
    return numerator // 6


def jackson_eta_upper(N: int) -> F:
    """Rigorous eta_N upper bound using pi < 355/113."""

    c_zero_numerator = N * (2 * N**2 + 1)
    require(c_zero_numerator % 3 == 0, "Jackson zero mode not integral")
    c_zero = c_zero_numerator // 3
    odd_sum = sum(
        F(jackson_coefficient(N, k), k**2)
        for k in range(1, 2 * N - 2, 2)
    )
    # eta_N = 1/2 - 4*odd_sum/(pi^2*c_zero).
    # The strict upper bound pi < 355/113 gives the displayed rational
    # upper bound.
    eta_upper = F(1, 2) - 4 * F(113, 355) ** 2 * odd_sum / c_zero
    require(0 < eta_upper < F(1, 2), "bad Jackson error bound")
    return eta_upper


def crossing_margin(eta: F, floor_a: F, floor_b: F) -> F:
    """Two-factor product margin with 2+7 Jackson telescope errors."""

    error_a = 2 * eta
    error_b = 7 * eta
    return (floor_a - error_a) * (floor_b - error_b) - 9 * eta


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

eta_99 = jackson_eta_upper(99)
eta_100 = jackson_eta_upper(100)
eta_297 = jackson_eta_upper(297)
eta_298 = jackson_eta_upper(298)

guard_pair_fail = crossing_margin(
    eta_99, guard_pair_floor, seven_ordinary_floor
)
guard_pair_pass = crossing_margin(
    eta_100, guard_pair_floor, seven_ordinary_floor
)
ordinary_pair_fail = crossing_margin(
    eta_297, ordinary_pair_floor, guard_six_floor
)
ordinary_pair_pass = crossing_margin(
    eta_298, ordinary_pair_floor, guard_six_floor
)

require(guard_pair_fail < -F(1, 100000), "N=99 ledger must fail")
require(guard_pair_pass > F(1, 100000), "N=100 ledger must pass")
require(
    ordinary_pair_fail < -F(1, 100000),
    "N=297 ledger must fail",
)
require(
    ordinary_pair_pass > F(1, 100000),
    "N=298 ledger must pass",
)
require(guard_pair_fail < 0 < guard_pair_pass, "guard-pair boundary")
require(
    ordinary_pair_fail < 0 < ordinary_pair_pass,
    "ordinary-pair boundary",
)
require(
    guard_pair_floor - 2 * eta_100 > 0
    and seven_ordinary_floor - 7 * eta_100 > 0,
    "guard-pair crossing factors",
)
require(
    ordinary_pair_floor - 2 * eta_298 > 0
    and guard_six_floor - 7 * eta_298 > 0,
    "ordinary-pair crossing factors",
)

guard_pair_types = 8
ordinary_pair_types = comb(8, 2)
require(
    guard_pair_types + ordinary_pair_types == comb(9, 2) == 36,
    "complete scalar pair split",
)

height_guard_pair = 2 * 100 - 2
height_ordinary_pair = 2 * 298 - 2
scalar_rank_height = max(35, height_guard_pair, height_ordinary_pair)
original_rank_height = 2 * scalar_rank_height

require(height_guard_pair == 198, "guard-pair crossing height")
require(height_ordinary_pair == 594, "ordinary-pair crossing height")
require(scalar_rank_height == 594, "uniform scalar rank height")
require(original_rank_height == 1188, "fixed-section lift height")

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
    "N99_margin_lt=-1/100000 N100_margin_gt=1/100000 height=198"
)
print(
    "ordinary_pair_cut: types=28 "
    "N297_margin_lt=-1/100000 N298_margin_gt=1/100000 height=594"
)
print(
    "pair_partition=36/36 scalar_rank_height=594 "
    "fixed_section_original_height=1188"
)
print("ALL EXACT CHECKS PASSED")
