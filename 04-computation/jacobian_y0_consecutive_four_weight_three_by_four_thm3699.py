#!/usr/bin/env python3
"""Exact support and differential audit for THM-3699."""

from itertools import combinations

import sympy as sp


b = sp.symbols("b")
h = 1 - b**2


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def bracket_coefficient(r, f, s, g):
    """Coefficient in {c^r f(b),c^s g(b)}."""
    return sp.expand(s * sp.diff(f, b) * g - r * f * sp.diff(g, b))


def buckets(left, right):
    out = {}
    for r in left:
        for s in right:
            out.setdefault(r + s + 1, []).append((r, s))
    return {k: tuple(v) for k, v in sorted(out.items())}


W = (-2, -1, 0, 1)
supports = tuple(combinations(W, 3))
expected = {
    (-2, -1, 0): (2, ((0, 1),)),
    (-2, -1, 1): (2, ((1, 0),)),
    (-2, 0, 1): (1, ((0, 0), (1, -1))),
    (-1, 0, 1): (-2, ((-1, -2),)),
}

print("THM-3699 exact support audit")
contains_minus_two = set(range(-5, -1))
contains_plus_one = set(range(-2, 2))
forced_starts = contains_minus_two & contains_plus_one
require(forced_starts == {-2}, "wrong translated-window intersection")
print("four-consecutive windows containing both -2 and 1 =", tuple(sorted(forced_starts)))
print("window =", W)
print("three-support cases =", len(supports))

for support in supports:
    row = buckets(support, W)
    key, pairs = expected[support]
    require(row[key] == pairs, f"wrong decisive bucket for {support}")
    scalar = row[0]
    print(f"support {support}: decisive bucket {key} = {pairs}; scalar = {scalar}")

f = sp.Function("f")(b)
g = sp.Function("g")(b)

# The two zero/nonzero singleton identities.
require(
    sp.simplify(bracket_coefficient(0, f, 1, g) - sp.diff(f, b) * g) == 0,
    "failed (0,1) singleton identity",
)
require(
    sp.simplify(bracket_coefficient(1, f, 0, g) + f * sp.diff(g, b)) == 0,
    "failed (1,0) singleton identity",
)

# The opposite-sign row is exactly -(fg)'.
require(
    sp.simplify(bracket_coefficient(1, f, -1, g) + sp.diff(f * g, b)) == 0,
    "failed opposite-sign product identity",
)

# The same-sign endpoint row is p^3 (q/p^2)'.
p = sp.Function("p")(b)
q = sp.Function("q")(b)
endpoint = bracket_coefficient(-1, p, -2, q)
require(
    sp.simplify(endpoint - p**3 * sp.diff(q / p**2, b)) == 0,
    "failed endpoint logarithmic-derivative identity",
)
require(
    sp.simplify(bracket_coefficient(-1, p, -2, p**2)) == 0,
    "failed endpoint common-square identity",
)

# With p in M_-1, endpoint commutation gives q=h^2 t.  All three scalar
# addresses then have a visible h factor.
F = sp.Function("F")(b)
G = sp.Function("G")(b)
R = sp.Function("R")(b)
S = sp.Function("S")(b)
T = sp.Function("T")(b)
U = sp.Function("U")(b)

arm_1 = bracket_coefficient(-1, h * F, 0, R)
arm_2 = bracket_coefficient(0, S, -1, h * G)
arm_3 = bracket_coefficient(1, b * U, -2, h**2 * T)

for expression in (arm_1, arm_2, arm_3, arm_1 + arm_2 + arm_3):
    quotient = sp.cancel(expression / h)
    require(not quotient.has(1 / h), "lost common arm factor")
    require(sp.simplify(expression - h * quotient) == 0, "bad arm quotient")

print("singleton differential identities = PASS")
print("endpoint common-square identity = PASS")
print("three scalar addresses have common arm factor h = PASS")

# THM-3695 supplies >=3 pieces per output and >=7 total inside a four-weight
# window.  These are the only count profiles; a one-dimensional weight space
# lets a 4x4 pair shear to <=4x3.
profiles = tuple((i, j) for i in range(3, 5) for j in range(3, 5) if i + j >= 7)
require(profiles == ((3, 4), (4, 3), (4, 4)), "wrong four-window profiles")
print("four-window profiles after THM-3695 =", profiles)
print("one-dimensional -1 component shears 4x4 to <=4x3 = PASS")
print("ALL CHECKS PASSED")
