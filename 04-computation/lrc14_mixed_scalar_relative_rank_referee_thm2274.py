#!/usr/bin/env python3
"""Independent exact referee and sharpening for THM-2274.

The current theorem uses the coarse Jackson error 3/(2N) in its
support-two branch.  This referee reconstructs THM-2145's exact Jackson
first-moment bound and proves the sharper bandwidths N=100 and N=298.

Every check uses explicit exceptions, so optimized Python executes the same
audit.
"""

from fractions import Fraction as F
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def arctan_partial(reciprocal: int, last_index: int) -> F:
    """Alternating partial sum for arctan(1/reciprocal)."""

    total = F(0)
    for index in range(last_index + 1):
        total += F(
            (-1) ** index,
            (2 * index + 1) * reciprocal ** (2 * index + 1),
        )
    return total


def tangent_double(value: F) -> F:
    return 2 * value / (1 - value * value)


# Self-contained rational pi brackets from Machin's identity
# pi/4=4 arctan(1/5)-arctan(1/239).
tan_one_fifth = F(1, 5)
tan_four_fifths = tangent_double(tangent_double(tan_one_fifth))
tan_difference = (
    tan_four_fifths - F(1, 239)
) / (
    1 + tan_four_fifths * F(1, 239)
)
require(tan_difference == 1, "Machin tangent identity")

# Odd truncation is a lower bound and even truncation is an upper bound.
pi_lower_certificate = 4 * (
    4 * arctan_partial(5, 7)
    - arctan_partial(239, 2)
)
pi_upper_certificate = 4 * (
    4 * arctan_partial(5, 4)
    - arctan_partial(239, 1)
)

pi_lower = F(103993, 33102)
pi_upper = F(355, 113)
require(pi_lower_certificate > pi_lower, "lower pi certificate")
require(pi_upper_certificate < pi_upper, "upper pi certificate")
require(pi_lower < pi_upper, "ordered pi brackets")


def jackson_coefficient(N: int, frequency: int) -> int:
    """Integer C_k from THM-2145."""

    k = abs(frequency)
    require(0 <= k <= 2 * N - 2, "frequency outside Jackson support")
    if k <= N:
        numerator = (
            4 * N**3
            - 6 * N * k**2
            + 2 * N
            + 3 * k**3
            - 3 * k
        )
        require(numerator % 6 == 0, "first Jackson branch nonintegral")
        return numerator // 6
    numerator = (2 * N - k) ** 3 - (2 * N - k)
    require(numerator % 6 == 0, "second Jackson branch nonintegral")
    return numerator // 6


def reconstruct_coefficient(N: int, frequency: int) -> int:
    """Direct convolution of the two triangular Fejer coefficient rows."""

    k = abs(frequency)
    total = 0
    for left in range(-N + 1, N):
        right = k - left
        if abs(right) < N:
            total += (N - abs(left)) * (N - abs(right))
    return total


def jackson_eta_bound(N: int, pi_bound: F) -> F:
    """Substitute a rational pi bound in the exact first-moment formula."""

    c_zero = F(N * (2 * N**2 + 1), 3)
    odd_sum = sum(
        F(jackson_coefficient(N, k), k * k)
        for k in range(1, 2 * N - 2, 2)
    )
    return F(1, 2) - 4 * odd_sum / (pi_bound * pi_bound * c_zero)


def crossing_margin(error: F, floor_a: F, floor_b: F) -> F:
    """Exact 2+7 block crossing margin."""

    return (
        (floor_a - 2 * error) * (floor_b - 7 * error)
        - 9 * error
    )


def margin_is_decreasing_below_centipercent(
    floor_a: F, floor_b: F
) -> bool:
    # derivative = 28e-(7a+2b+9), evaluated at e<=1/100
    return F(28, 100) - (7 * floor_a + 2 * floor_b + 9) < 0


# Independently reconstruct every Jackson coefficient used at the four
# boundary bandwidths.
coefficient_checks = 0
for bandwidth in (99, 100, 297, 298):
    for frequency in range(0, 2 * bandwidth - 1):
        require(
            jackson_coefficient(bandwidth, frequency)
            == reconstruct_coefficient(bandwidth, frequency),
            f"Jackson convolution mismatch N={bandwidth} k={frequency}",
        )
        coefficient_checks += 1

guard_pair_floor = F(4, 7)
seven_ordinary_floor = F(15, 154)
ordinary_pair_floor = F(66, 91)
guard_six_floor = F(191, 6930)

require(
    F(5, 7) + F(6, 7) - 1 == guard_pair_floor,
    "guard+ordinary union floor",
)
require(
    1 - F(2, 7) + F(1, 91) == ordinary_pair_floor,
    "ordinary pair overlap floor",
)
require(
    margin_is_decreasing_below_centipercent(
        guard_pair_floor, seven_ordinary_floor
    ),
    "guard-pair margin monotonicity",
)
require(
    margin_is_decreasing_below_centipercent(
        ordinary_pair_floor, guard_six_floor
    ),
    "ordinary-pair margin monotonicity",
)


def eta_upper(N: int) -> F:
    # pi<355/113 makes 1/pi^2 larger and eta smaller only after checking
    # the sign: eta=1/2-positive/pi^2.  Substituting the upper pi therefore
    # gives an upper bound on eta.
    return jackson_eta_bound(N, pi_upper)


def eta_lower(N: int) -> F:
    # Substituting the lower pi gives a strict lower bound on eta.
    return jackson_eta_bound(N, pi_lower)


def certified_crossing(
    N: int, floor_a: F, floor_b: F
) -> bool:
    error = eta_upper(N)
    return (
        floor_a - 2 * error > 0
        and floor_b - 7 * error > 0
        and crossing_margin(error, floor_a, floor_b) > 0
    )


guard_fail = crossing_margin(
    eta_lower(99), guard_pair_floor, seven_ordinary_floor
)
guard_pass = crossing_margin(
    eta_upper(100), guard_pair_floor, seven_ordinary_floor
)
ordinary_fail = crossing_margin(
    eta_lower(297), ordinary_pair_floor, guard_six_floor
)
ordinary_pass = crossing_margin(
    eta_upper(298), ordinary_pair_floor, guard_six_floor
)

require(guard_fail < 0 < guard_pass, "guard pair 99/100 boundary")
require(
    ordinary_fail < 0 < ordinary_pass,
    "ordinary pair 297/298 boundary",
)
require(eta_upper(100) < F(1, 100), "guard pass error range")
require(eta_upper(298) < F(1, 100), "ordinary pass error range")

first_guard_pass = next(
    N
    for N in range(2, 101)
    if certified_crossing(N, guard_pair_floor, seven_ordinary_floor)
)
first_ordinary_pass = next(
    N
    for N in range(2, 299)
    if certified_crossing(N, ordinary_pair_floor, guard_six_floor)
)
require(first_guard_pass == 100, "first guard-pair certified bandwidth")
require(
    first_ordinary_pass == 298,
    "first ordinary-pair certified bandwidth",
)

guard_pair_types = 8
ordinary_pair_types = comb(8, 2)
require(
    guard_pair_types + ordinary_pair_types == comb(9, 2) == 36,
    "complete adaptive pair split",
)

guard_height = 2 * first_guard_pass - 2
ordinary_height = 2 * first_ordinary_pass - 2
scalar_rank_height = max(35, guard_height, ordinary_height)
original_rank_height = 2 * scalar_rank_height

require(guard_height == 198, "guard-pair height")
require(ordinary_height == 594, "ordinary-pair height")
require(scalar_rank_height == 594, "sharpened scalar rank height")
require(original_rank_height == 1188, "fixed-section lift height")

print("THM-2274 INDEPENDENT EXACT JACKSON REFEREE")
print("pi_brackets_certified=103993/33102<pi<355/113")
print(f"jackson_coefficients_reconstructed={coefficient_checks}")
print(
    "guard_pair_boundary="
    "N99:true_margin_negative,"
    "N100:certified_positive,"
    "height:198"
)
print(
    "ordinary_pair_boundary="
    "N297:true_margin_negative,"
    "N298:certified_positive,"
    "height:594"
)
print("interval_error=2*first_Jackson_moment_independent_of_interval_length")
print("adaptive_pair_partition=guard:8,ordinary:28,total:36")
print("scalar_rank_height=594")
print("fixed_section_original_height=1188")
print("ALL_EXACT_CHECKS_PASSED")
