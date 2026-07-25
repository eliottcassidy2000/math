#!/usr/bin/env python3
"""Exact arithmetic companion for THM-2279.

This checks the uniform shallow-owner/target floors, the first ideal
Zakharov word-length threshold, and the delayed-tail hostile family on all
150 strict first-depth-one profiles.  The BV transfer estimate itself is
proved symbolically in the theorem.
"""

from fractions import Fraction
from math import gcd


def circle_norm(x: Fraction) -> Fraction:
    y = x - (x.numerator // x.denominator)
    return min(y, 1 - y)


def T(x: Fraction, k: int = 1) -> Fraction:
    y = x * (13**k)
    return y - (y.numerator // y.denominator)


def zakharov_constant(n: int) -> Fraction:
    return Fraction(1, n) * Fraction(n, n + 1) ** (n + 1)


delta5 = Fraction(961, 6930)
L0 = Fraction(5696989, 367580070)
e0 = L0 / 2
G1 = Fraction(10, 91)
eta0 = delta5 - G1
handoff_floor = e0 * eta0 / 2
expiration_support_floor = Fraction(169, 20) * e0
word_product = expiration_support_floor * eta0

assert eta0 == Fraction(2593, 90090)
assert handoff_floor == Fraction(14772292477, 132461154025200)
assert expiration_support_floor == Fraction(5696989, 87001200)
assert word_product == Fraction(14772292477, 7837938108000)

# The parity-sharp guard/deep cap used in THM-2273 is largest at depth one.
def G(d: int) -> Fraction:
    if d % 2:
        return Fraction(5, 49) + Fraction(5, 49 * 13**d)
    return Fraction(5, 49) + Fraction(5, 294 * 13**d)


depth_caps = [(d, G(d)) for d in range(1, 20)]
assert max(depth_caps, key=lambda row: row[1]) == (1, G1)

# If the post-expiration support and opposite residual were honest,
# tail-independent word families, n=195 would be the first Zakharov
# threshold before cylinder-boundary losses.
first_word_threshold = next(
    n for n in range(1, 10000) if zakharov_constant(n) < word_product
)
assert first_word_threshold == 195
assert zakharov_constant(194) >= word_product

# Exact pointwise hostile family.  The packet is a local membership control,
# not a global scalar cover.
strict_profiles = [(b, c) for c in range(5, 20) for b in range(2, c)]
assert len(strict_profiles) == 150

x_infinity = Fraction(79, 338)
y_infinity = Fraction(1, 2)
assert T(x_infinity, 2) == y_infinity
assert int(13 * x_infinity) == 3
assert int(13 * T(x_infinity)) == 0

q = [1 + 338 * i for i in range(1, 6)]
assert all(gcd(a, 13) == 1 for a in q)

checks = 0
for b, c in strict_profiles:
    for delay in (c, c + 1, c + 7):
        y = Fraction(1, 2) + Fraction(1, 2 * 13**delay)
        x = (39 + y) / 169

        # Same two inverse-root digits as x_infinity.
        assert int(13 * x) == 3
        assert int(13 * T(x)) == 0
        assert T(x, 2) == y

        # Strict c_1-exclusive source membership.
        assert circle_norm(x) > Fraction(1, 7)
        assert all(circle_norm(a * x) >= Fraction(1, 14) for a in q)
        assert circle_norm(13 * x) < Fraction(1, 14)
        assert circle_norm(13**b * x) >= Fraction(1, 14)
        assert circle_norm(13**c * x) >= Fraction(1, 14)

        # Normalized c_1-core future: 0^delay 1^infinity.
        for k in range(delay):
            assert circle_norm(T(y, k)) >= Fraction(1, 14)
        assert T(y, delay) == 0
        for k in range(delay, delay + 4):
            assert circle_norm(T(y, k)) < Fraction(1, 14)

        # The limiting tail is 0^infinity.
        for k in range(delay + 4):
            assert T(y_infinity, k) == y_infinity
            assert circle_norm(T(y_infinity, k)) >= Fraction(1, 14)
        checks += 1

print("THM-2279 exact companion")
print(f"strict_profiles={len(strict_profiles)}")
print(f"tail_family_checks={checks}")
print(f"L0={L0}")
print(f"e0=L0/2={e0} ~= {float(e0):.15f}")
print(f"G1={G1}")
print(f"eta0=delta5-G1={eta0} ~= {float(eta0):.15f}")
print(
    "handoff_floor=e0*eta0/2="
    f"{handoff_floor} ~= {float(handoff_floor):.15f}"
)
print(
    "BV_horizon_condition: 13^k >= "
    f"(4/e0) S = {Fraction(4, 1) / e0} S"
)
print(f"expiration_support_floor={expiration_support_floor}")
print(f"ideal_word_product={word_product}")
print(f"first_ideal_Zakharov_threshold={first_word_threshold}")
print(f"C_194 ~= {float(zakharov_constant(194)):.15f} >= product")
print(f"C_195 ~= {float(zakharov_constant(195)):.15f} < product")
print("hostile_root_word=(3,0)")
print("hostile_futures=0^infinity versus 0^R 1^infinity for arbitrary R>=c")
print("scope=no fixed-expiration return, no profile exclusion, no LRC14 proof")
