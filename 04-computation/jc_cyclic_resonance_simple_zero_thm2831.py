#!/usr/bin/env python3
"""Exact referee for THM-2831.

The theorem is symbolic.  This companion independently checks the
zero-source specialization of the exact THM-2827 Faber formula, the
simple-zero valuation dichotomy, the resonance exponent arithmetic, and
the sharp formal c=0 boundary.  It uses only exact integer/Fraction
arithmetic and explicit failures; Python optimization removes no check.
"""

from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def generalized_binomial(a, n):
    value = Fraction(1)
    for j in range(n):
        value *= a - j
        value /= j + 1
    return value


def active_row(j):
    """Return (odd T exponent, coefficient) at d=s=0, or None."""
    if j < 2 or (j - 2) % 3:
        return None
    h = (j - 2) // 3
    coefficient = 4 * generalized_binomial(Fraction(2 * j - 1, 2), j + 1 + h)
    require(coefficient != 0, f"active row {j} lost its coefficient")
    return 2 * h + 1, coefficient


print("THM-2831 CYCLIC RESONANCE SIMPLE-ZERO AUDIT")
print()
print("zero-source active Faber rows")

for k in range(2, 11):
    R = 3 * k + 2
    D = 4 * k + 3
    e = 2 * k + 1
    rows = [(j, active_row(j)) for j in range(1, R + 1)]
    live = [(j, data[0], data[1]) for j, data in rows if data is not None]
    require(live[-1][0] == R, f"top row missing at R={R}")
    require(live[-1][1] == e, f"wrong top exponent at R={R}")
    require(len({exponent for _, exponent, _ in live}) == len(live),
            f"nonunique T exponents at R={R}")
    for j, data in rows:
        if data is None:
            require(j % 3 != 2, f"missed residue-two row j={j}")
    if k <= 6:
        print(f"R={R:2d} D={D:2d} live={len(live):2d} "
              f"top=T^{e} coeff={live[-1][2]}")

require(active_row(8)[1] == Fraction(-195, 131072),
        "R=8 hostile coefficient changed")

print()
print("simple-zero valuation dichotomy")
for k in range(2, 11):
    e = 2 * k + 1
    # a=0: the unique highest exponent controls.
    require(-e != 1, f"polar order accidentally equals one for k={k}")
    # a>=1: every possible least live odd exponent gives order at least two.
    for a in range(1, 13):
        for ell in range(1, e + 1, 2):
            order = a + ell * (2 * a - 1)
            require(order >= 2 and order != 1,
                    f"regular order one at k={k}, a={a}, ell={ell}")
print("a=0 gives -e; a>=1 gives a+ell(2a-1)>=2: PASS")

print()
print("top-only UFD resonance arithmetic")
for k in range(2, 11):
    D = 4 * k + 3
    e = (D - 1) // 2
    require(2 * e + 1 == D, f"wrong odd degree at k={k}")
    require(0 < e + 1 < D and (e + 1) % D != 0,
            f"S exponent became a D-multiple at k={k}")
    for delta in range(1, 8):
        N = D * delta
        y_exponent = 1 + e * (N + 2)
        formal_c_zero_exponent = (e + 1) * N + y_exponent
        require(formal_c_zero_exponent == D * (N + 1),
                f"formal c=0 boundary failed at k={k}, delta={delta}")
print("S exponent=(D+1)/2 is nondivisible; c=0 merges to D(N+1): PASS")

print()
print("RESULT: PASS")
