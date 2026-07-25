#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2265.

The script checks the squared-Fejer normalization and moment bound, the
six- and seven-speed safe-mass inputs, the equal-bandwidth cutoff, and several
anisotropic bandwidth controls.  The Fourier factorization itself is proved
symbolically in the theorem.
"""

from fractions import Fraction as F


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fejer_square_mass(n):
    # Parseval for coefficients 1-|k|/n, |k|<n.
    return F(1) + 2 * sum((F(n - k, n) ** 2 for k in range(1, n)), F(0))


def centered_first_moment_bound(n):
    # Split at 1/(2n), using F_n<=n on the center and
    # F_n<=1/(4n|t|^2) on the tail.
    numerator_bound = F(1, 2) - F(1, 4 * n * n)
    return numerator_bound / fejer_square_mass(n)


def coupling_margin(n_f, n_e):
    safe_f = F(18, 91)
    safe_e = F(15, 154)
    err_f = F(9, n_f)
    err_e = F(21, 2 * n_e)
    return (safe_f - err_f) * (safe_e - err_e) - err_f - err_e


for n in range(2, 201):
    mass = fejer_square_mass(n)
    require(mass == F(2 * n * n + 1, 3 * n), "Fejer mass mismatch")
    require(
        centered_first_moment_bound(n) < F(3, 4 * n),
        "first-moment estimate failed",
    )

six_safe = 1 - F(6, 7) + 5 * F(1, 91)
seven_safe = F(15, 154)
require(six_safe == F(18, 91), "six-speed Hunter floor mismatch")

margin1161 = coupling_margin(1161, 1161)
margin1162 = coupling_margin(1162, 1162)
require(margin1161 == -F(2293, 699620922), "N=1161 margin mismatch")
require(margin1162 == F(465, 35106344), "N=1162 margin mismatch")
require(margin1161 < 0 < margin1162, "equal-bandwidth cutoff mismatch")

first_equal = None
for n in range(108, 1163):
    if coupling_margin(n, n) > 0:
        first_equal = n
        break
require(first_equal == 1162, "first positive equal bandwidth mismatch")

anisotropic_rows = [
    (1161, 1162, F(449, 75023949)),
    (1130, 1187, F(4569, 4699279585)),
    (1019, 1304, F(93, 423215156)),
    (600, 4427, F(9, 326526200)),
]
for n_f, n_e, expected in anisotropic_rows:
    margin = coupling_margin(n_f, n_e)
    require(margin == expected and margin > 0, "anisotropic margin mismatch")

coefficient_height = 2 * 1162 - 2
largest_seven_core_sum = sum(range(7, 14))
core_carry_bound = coefficient_height * largest_seven_core_sum
require(coefficient_height == 2322, "coefficient height mismatch")
require(largest_seven_core_sum == 70, "seven-core sum mismatch")
require(core_carry_bound == 162540, "core carry mismatch")

print("THM-2265 exact referee")
print("fejer_square_mass_and_moment_bound=PASS N=2..200")
print(f"six_speed_safe_floor={six_safe}")
print(f"seven_speed_safe_floor={seven_safe}")
print(f"equal_N_1161_margin={margin1161}")
print(f"equal_N_1162_margin={margin1162}")
print(f"first_positive_equal_bandwidth={first_equal}")
for n_f, n_e, expected in anisotropic_rows:
    print(
        f"anisotropic_bandwidth=({n_f},{n_e}) "
        f"heights=({2*n_f-2},{2*n_e-2}) margin={expected}"
    )
print(f"balanced_relation_height={coefficient_height}")
print(f"largest_seven_core_sum={largest_seven_core_sum}")
print(f"defect_six_core_carry_bound={core_carry_bound}")
print("ALL_CHECKS_PASS")
