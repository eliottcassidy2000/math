#!/usr/bin/env python3
"""Exact support-placement and differential audit for THM-3707."""

from collections import defaultdict

import sympy as sp


b = sp.symbols("b")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def bracket_coefficient(r, f, s, g):
    """Coefficient in {c^r f(b),c^s g(b)}."""
    return sp.expand(s * sp.diff(f, b) * g - r * f * sp.diff(g, b))


def buckets(left, right):
    out = defaultdict(list)
    for r in left:
        for s in right:
            out[r + s + 1].append((r, s))
    return {weight: tuple(out[weight]) for weight in sorted(out)}


print("THM-3707 exact aligned-AP support audit")

# At a simple arm only a complementary (-2,1) or (1,-2) address can
# contribute.  Solve both oriented placement equations exactly.  The index
# distance times the common step is 3, so only the positive divisors of 3
# occur; the loop below records orientations before quotienting by support.
oriented = []
for i in range(3):
    for j in range(4):
        if j > i and 3 % (j - i) == 0:
            step = 3 // (j - i)
            start = -2 - i * step
            if start + j * step == 1:
                oriented.append((step, start, "P:-2,Q:+1", i, j))
        if i > j and 3 % (i - j) == 0:
            step = 3 // (i - j)
            start = 1 - i * step
            if start + j * step == -2:
                oriented.append((step, start, "P:+1,Q:-2", i, j))

expected_oriented = {
    (1, -2, "P:-2,Q:+1", 0, 3),
    (3, -2, "P:-2,Q:+1", 0, 1),
    (3, -2, "P:+1,Q:-2", 1, 0),
    (3, -5, "P:-2,Q:+1", 1, 2),
    (3, -5, "P:+1,Q:-2", 2, 1),
    (3, -8, "P:-2,Q:+1", 2, 3),
}
require(
    len(oriented) == len(expected_oriented) and set(oriented) == expected_oriented,
    "oriented scalar-address enumeration changed",
)

placements = defaultdict(list)
for step, start, orientation, i, j in oriented:
    placements[(step, start)].append((orientation, i, j))

expected_placements = {
    (1, -2): [("P:-2,Q:+1", 0, 3)],
    (3, -2): [("P:-2,Q:+1", 0, 1), ("P:+1,Q:-2", 1, 0)],
    (3, -5): [("P:-2,Q:+1", 1, 2), ("P:+1,Q:-2", 2, 1)],
    (3, -8): [("P:-2,Q:+1", 2, 3)],
}
require(dict(placements) == expected_placements, "support-placement quotient changed")

print("oriented arm placements =")
for row in sorted(oriented, key=lambda item: (item[0], item[1], item[3], item[4])):
    print(" ", row)
print("support placements =", tuple(sorted(placements)))

expected_scalar = {
    (1, -2): ((-2, 1), (-1, 0), (0, -1)),
    (3, -2): ((-2, 1), (1, -2)),
    (3, -5): ((-5, 4), (-2, 1), (1, -2)),
    (3, -8): ((-2, 1),),
}

for step, start in sorted(placements):
    left = tuple(start + k * step for k in range(3))
    right = tuple(start + k * step for k in range(4))
    table = buckets(left, right)
    require(table[0] == expected_scalar[(step, start)], "wrong scalar fibre")
    sizes = tuple(len(table[w]) for w in table)
    print(
        f"(d,a)=({step},{start}): P={left}; Q={right}; "
        f"fibre_sizes={sizes}; scalar={table[0]}"
    )

# d=1,a=-2: output weight 2 is the singleton (0,1).  Its vanishing
# coefficient is p_0' q_1, forcing the retained weight-zero piece to be a
# scalar because q_1 is active.
p0 = sp.Function("p0")(b)
q1 = sp.Function("q1")(b)
require(
    sp.simplify(bracket_coefficient(0, p0, 1, q1) - sp.diff(p0, b) * q1) == 0,
    "failed d=1 singleton identity",
)

# d=3,a=-2: the lowest singleton (-2,-2) forces g=lambda*f.  The whole
# scalar double then compresses to one homogeneous complementary bracket.
f = sp.Function("f")(b)
p = sp.Function("p")(b)
q = sp.Function("q")(b)
lam = sp.symbols("lambda", nonzero=True)
require(bracket_coefficient(-2, f, -2, lam * f) == 0, "failed low commutation")
scalar_double = bracket_coefficient(-2, f, 1, q) + bracket_coefficient(
    1, p, -2, lam * f
)
compressed = bracket_coefficient(-2, f, 1, q - lam * p)
require(sp.simplify(scalar_double - compressed) == 0, "failed scalar compression")

# d=3,a=-8 has exactly one scalar address.  The d=3,a=-5 placement is
# literally the universal THM-3702 support word.
require(expected_scalar[(3, -8)] == ((-2, 1),), "lost top singleton")
require(
    tuple(-5 + 3 * k for k in range(3)) == (-5, -2, 1)
    and tuple(-5 + 3 * k for k in range(4)) == (-5, -2, 1, 4),
    "THM-3702 support identification changed",
)

print("d=1 retained-zero singleton identity = PASS")
print("d=3,a=-2 endpoint compression to one homogeneous bracket = PASS")
print("d=3,a=-5 identification with THM-3702 = PASS")
print("d=3,a=-8 scalar singleton = PASS")
print("ALL CHECKS PASSED")
