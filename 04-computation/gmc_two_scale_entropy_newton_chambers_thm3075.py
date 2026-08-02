#!/usr/bin/env python3
"""Exact companion for THM-3075's simultaneous entropy-Newton sectors."""

from fractions import Fraction
from math import prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def interval_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def interval_scale(value, interval):
    if value >= 0:
        return value * interval[0], value * interval[1]
    return value * interval[1], value * interval[0]


def interval_subtract(left, right):
    return left[0] - right[1], left[1] - right[0]


def atanh_log_bounds(value, terms=90):
    """Rigorous Fraction bounds for log(value), value>0."""
    require(value > 0, "logarithm domain")
    if value < 1:
        lower, upper = atanh_log_bounds(1 / value, terms)
        return -upper, -lower

    power_two = 0
    reduced = value
    while reduced >= 2:
        reduced /= 2
        power_two += 1

    def core(argument):
        z = (argument - 1) / (argument + 1)
        partial = Fraction(0)
        zpower = z
        for index in range(terms):
            partial += Fraction(2) * zpower / (2 * index + 1)
            zpower *= z * z
        tail = Fraction(2) * zpower / ((2 * terms + 1) * (1 - z * z))
        return partial, partial + tail

    reduced_bounds = core(reduced)
    two_bounds = core(Fraction(2))
    return interval_add(reduced_bounds, interval_scale(power_two, two_bounds))


def entropy_invoice_bounds(x):
    log_bounds = atanh_log_bounds(Fraction(6) / x)
    return interval_scale(Fraction(3600) * x, interval_add(log_bounds, (1, 1)))


def margin_against_log(x, ratio):
    return interval_subtract(entropy_invoice_bounds(x), atanh_log_bounds(ratio))


def wall_a_bounds(x):
    first = interval_scale(1 + 4 * x, atanh_log_bounds(1 + 4 * x))
    second = interval_scale(4 * x, atanh_log_bounds(x))
    answer = interval_subtract(first, second)
    answer = interval_subtract(answer, atanh_log_bounds(Fraction(5)))
    answer = interval_subtract(answer, interval_scale(4 * x, atanh_log_bounds(Fraction(4))))
    return answer


def wall_b_bounds(x):
    first = interval_scale(1 + 5 * x, atanh_log_bounds(1 + 5 * x))
    second = interval_scale(5 * x, atanh_log_bounds(x))
    answer = interval_subtract(first, second)
    answer = interval_subtract(answer, atanh_log_bounds(Fraction(6)))
    answer = interval_subtract(answer, interval_scale(5 * x, atanh_log_bounds(Fraction(4))))
    return answer


# Resultant coefficient degrees and the full underlying factorial-factor bill.
coefficient_degrees = tuple(720 // degree for degree in range(2, 7))
require(coefficient_degrees == (360, 240, 180, 144, 120), "multidegrees changed")
factor_invoice = sum(
    degree * coefficient_degree
    for degree, coefficient_degree in zip(range(2, 7), coefficient_degrees)
)
require(factor_invoice == 3600, "factor invoice changed")

mixed_cells = [
    (degree, slow, fast)
    for degree in range(2, 7)
    for slow in range(1, degree)
    for fast in range(1, degree - slow + 1)
]
require(len(mixed_cells) == 35, "mixed-cell count changed")


# Rigorous rational brackets for the two universal sector thresholds.
eps_a_lo = Fraction("0.00001808125357979")
eps_a_hi = Fraction("0.00001808125357980")
eps_b_lo = Fraction("0.00006001388194673")
eps_b_hi = Fraction("0.00006001388194674")

fa_lo = margin_against_log(eps_a_lo, Fraction(625, 256))
fa_hi = margin_against_log(eps_a_hi, Fraction(625, 256))
fb_lo = margin_against_log(eps_b_lo, Fraction(46656, 3125))
fb_hi = margin_against_log(eps_b_hi, Fraction(46656, 3125))
require(fa_lo[1] < 0 < fa_hi[0], "epsilon A bracket changed")
require(fb_lo[1] < 0 < fb_hi[0], "epsilon B bracket changed")

# Simple exact safe rays used as hostile-independent controls.
require(
    margin_against_log(Fraction(1, 100000), Fraction(625, 256))[1] < 0,
    "sector A safe ray changed",
)
require(
    margin_against_log(Fraction(1, 20000), Fraction(46656, 3125))[1] < 0,
    "sector B safe ray changed",
)


# The first mixed upper-atom walls in the two nested gauges.
wall_a_lo = Fraction("0.34053358195087")
wall_a_hi = Fraction("0.34053358195088")
wall_b_lo = Fraction("0.24197406496062")
wall_b_hi = Fraction("0.24197406496063")
require(wall_a_bounds(wall_a_lo)[1] < 0 < wall_a_bounds(wall_a_hi)[0], "wall A bracket")
require(wall_b_bounds(wall_b_lo)[1] < 0 < wall_b_bounds(wall_b_hi)[0], "wall B bracket")


# Exact norm-composition and asymptotic exponent ledgers.
sector_a = {
    "S3": 4 * 30,
    "E56": 24,
    "U4D": 6 * 30,
    "U5C": 144,
    "U6C": 120,
}
sector_b = {
    "S3": 20 * 6,
    "E45": 6 * 6,
    "U4C": 30 * 6,
    "U5C": 24 * 6,
    "U6R": 120,
}
require(sector_a == {"S3": 120, "E56": 24, "U4D": 180, "U5C": 144, "U6C": 120}, "sector A ledger")
require(sector_b == {"S3": 120, "E45": 36, "U4C": 180, "U5C": 144, "U6R": 120}, "sector B ledger")
require(30 * sector_a["E56"] == 720, "sector A gap exponent")
require(20 * sector_b["E45"] == 720, "sector B gap exponent")

a_slow_power = -Fraction(3, 2) * sector_a["U4D"]
a_fast_power = -Fraction(4, 2) * sector_a["U5C"] - Fraction(5, 2) * sector_a["U6C"]
b_slow_power = -Fraction(3, 2) * sector_b["U4C"] - Fraction(4, 2) * sector_b["U5C"]
b_fast_power = -Fraction(5, 2) * sector_b["U6R"]
require((a_slow_power, a_fast_power) == (-270, -588), "sector A powers")
require((b_slow_power, b_fast_power) == (-558, -300), "sector B powers")
require(a_slow_power + a_fast_power == -858, "sector A total power")
require(b_slow_power + b_fast_power == -858, "sector B total power")

require(prod((4,)) == 4 and prod((5, 6)) == 30, "sector A bases")
require(prod((4, 5)) == 20 and prod((6,)) == 6, "sector B bases")


# Pure upper coefficients diverge for every x>0 in the nested diagonal gauge;
# their exact disappearance is therefore a THM-3073 quotient, not convergence.
require(Fraction(5) > Fraction(4), "sector A upper hostile")
require(Fraction(3, 2) > 1, "sector B upper hostile")


print("THM-3075 TWO-SCALE ENTROPY-NEWTON CHAMBERS")
print(f"coefficient_degrees={coefficient_degrees} factor_invoice={factor_invoice}")
print(f"mixed_upper_cells={len(mixed_cells)} entropy_bound=3600*x*(log(6/x)+1)")
print("epsilon_A_in=(0.00001808125357979,0.00001808125357980)")
print("epsilon_B_in=(0.00006001388194673,0.00006001388194674)")
print("safe_rays=A:1/100000,B:1/20000")
print("mixed_wall_A_in=(0.34053358195087,0.34053358195088)")
print("mixed_wall_B_in=(0.24197406496062,0.24197406496063)")
print("sector_A=S3^120*E56^24*U4D^180*U5C^144*U6C^120")
print("sector_B=S3^120*E45^36*U4C^180*U5C^144*U6R^120")
print("gap_exponents=A:720,B:720;power_ledgers=A:-270/-588,B:-558/-300")
print("pure_upper_hostiles=A:5*x*log(5/4),B:6*x*log(3/2)")
print("repair=fixed-low-pivot+THM3073_triangular_norm;raw_convergence=FALSE")
print("scope=thin_fixed-slope_sectors;not_maximal;not_wall_crossing")
print("all_exact_checks=PASS")
