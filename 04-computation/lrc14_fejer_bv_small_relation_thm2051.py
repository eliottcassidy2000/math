#!/usr/bin/env python3
"""Exact arithmetic audit for THM-2051's centered BONF5 budget."""

from fractions import Fraction as F
from math import comb, factorial


H = 2**20
p = F(1, 7)


def centered_coefficient(m: int) -> F:
    return sum(
        F((-1) ** j * comb(13 - m, j - m)) * p ** (j - m)
        for j in range(m, 6)
    )


baseline = sum(F((-1) ** j * comb(13, j)) * p**j for j in range(6))
coefficients = {m: centered_coefficient(m) for m in range(2, 6)}
assert baseline == F(2052, 16807)
assert coefficients == {2: F(24, 343), 3: F(-24, 49), 4: F(-2, 7), 5: F(-1)}

# Exact lower bound e^14 > sum_{n=0}^{18}14^n/n! > H+1, hence log(H+1)<14.
exp14_lower = sum(F(14**n, factorial(n)) for n in range(19))
assert exp14_lower > H + 1
epsilon_bound = F(29, 2 * (H + 1))

telescoping_factor = sum(
    F(comb(13, m)) * abs(coefficients[m]) * m * F(6, 7) ** (m - 1)
    for m in range(2, 6)
)
error_bound = epsilon_bound * telescoping_factor
margin = baseline - error_bound
assert telescoping_factor == F(1477008, 343)
assert error_bound == F(21416616, 359661911)
assert margin == F(1102265820, 17623433639) > 0

print("THM-2051 exact centered BONF5 audit")
print(f"H={H}")
print(f"B5_equilibrium={baseline}")
for m in range(2, 6):
    print(f"centered coefficient c_{m}={coefficients[m]}")
print(f"finite lower bound for exp(14)={exp14_lower}>{H + 1}")
print(f"Fejer BV epsilon<{epsilon_bound}")
print(f"telescoping factor={telescoping_factor}")
print(f"total centered error<{error_bound}")
print(f"strict BONF5 margin>{margin}")
print("Conclusion: no support-2..5 relation with 0<|k_i|<=2^20 implies a strict lonely interval.")
