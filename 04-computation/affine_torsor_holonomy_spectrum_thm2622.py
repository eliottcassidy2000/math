#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2622.

The proof in the theorem is general.  This program exhausts the two live
finite specializations: Aff(F_13) for the LRC deck and
AGL(2,2)=V_4 semidirect S_3 for the quartic resolvent deck.
"""

from collections import Counter


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fixed_affine_line(a, c, p):
    return tuple(x for x in range(p) if (a * x + c - x) % p == 0)


def add2(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])


def mat2(a, x):
    return (
        (a[0][0] * x[0] + a[0][1] * x[1]) % 2,
        (a[1][0] * x[0] + a[1][1] * x[1]) % 2,
    )


def mul2(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(2)) % 2
              for j in range(2))
        for i in range(2)
    )


I2 = ((1, 0), (0, 1))
V4 = ((0, 0), (0, 1), (1, 0), (1, 1))
GL2 = tuple(
    ((a, b), (c, d))
    for a in range(2) for b in range(2)
    for c in range(2) for d in range(2)
    if (a * d - b * c) % 2 == 1
)


def affine2(g, x):
    a, c = g
    return add2(mat2(a, x), c)


def compose2(g, h):
    """Return g after h."""
    a, c = g
    b, d = h
    return mul2(a, b), add2(mat2(a, d), c)


def inverse2(g):
    for h in AFF2:
        if compose2(g, h) == (I2, (0, 0)) and compose2(h, g) == (I2, (0, 0)):
            return h
    raise RuntimeError("affine inverse missing")


def fixed2(g):
    return tuple(x for x in V4 if affine2(g, x) == x)


def matrix_order(a):
    power = I2
    for order in range(1, 7):
        power = mul2(power, a)
        if power == I2:
            return order
    raise RuntimeError("GL2 matrix order exceeded six")


def cycle_type(g):
    unseen = set(V4)
    lengths = []
    while unseen:
        start = min(unseen)
        x = start
        length = 0
        while True:
            unseen.remove(x)
            length += 1
            x = affine2(g, x)
            if x == start:
                break
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


AFF2 = tuple((a, c) for a in GL2 for c in V4)


print("THM-2622 exact affine-torsor holonomy controls")

c13_hist = Counter()
c13_fixed_action = Counter()
c13_twisted_unique = 0
for a in range(1, 13):
    for c in range(13):
        count = len(fixed_affine_line(a, c, 13))
        c13_hist[count] += 1
        if a == 1:
            c13_fixed_action[count] += 1
        else:
            require(count == 1, "nonidentity C13 linear part lost uniqueness")
            c13_twisted_unique += 1
require(c13_hist == Counter({1: 143, 0: 12, 13: 1}),
        "C13 affine fixed-section spectrum changed")
require(c13_fixed_action == Counter({0: 12, 13: 1}),
        "C13 fixed-action spectrum changed")
require(c13_twisted_unique == 143, "C13 twisted atlas count changed")
print(f"C13_affine_maps=156 fixed_section_hist={sorted(c13_hist.items())}")
print(f"C13_fixed_action_hist={sorted(c13_fixed_action.items())} "
      f"generator_twisted_unique={c13_twisted_unique}")

require(len(GL2) == 6 and len(AFF2) == 24, "V4 affine group size changed")
v4_fixed_hist = Counter(len(fixed2(g)) for g in AFF2)
v4_cycle_hist = Counter(cycle_type(g) for g in AFF2)
require(v4_fixed_hist == Counter({0: 9, 1: 8, 2: 6, 4: 1}),
        "V4 affine fixed-section spectrum changed")
require(v4_cycle_hist == Counter({
    (4,): 6, (3, 1): 8, (2, 2): 3, (2, 1, 1): 6,
    (1, 1, 1, 1): 1,
}), "AGL(2,2) cycle atlas changed")
print(f"V4_GL_size={len(GL2)} affine_size={len(AFF2)} "
      f"fixed_section_hist={sorted(v4_fixed_hist.items())}")
print(f"V4_cycle_type_hist={sorted(v4_cycle_hist.items())}")

linear_fixed = Counter()
for a in GL2:
    order = matrix_order(a)
    for c in V4:
        linear_fixed[order, len(fixed2((a, c)))] += 1
require(linear_fixed == Counter({
    (1, 4): 1, (1, 0): 3,
    (2, 2): 6, (2, 0): 6,
    (3, 1): 8,
}), "linear-order/fixed-section atlas changed")
print(f"V4_linear_order_fixed_hist={sorted(linear_fixed.items())}")

# One order-two linear part and all translations form a Sylow-2 D4 subgroup.
order_two = next(a for a in GL2 if matrix_order(a) == 2)
d4 = tuple((a, c) for a in (I2, order_two) for c in V4)
require(all(compose2(g, h) in d4 for g in d4 for h in d4),
        "selected D4 atlas is not a subgroup")
d4_fixed_hist = Counter(len(fixed2(g)) for g in d4)
d4_cycle_hist = Counter(cycle_type(g) for g in d4)
require(d4_fixed_hist == Counter({0: 5, 2: 2, 4: 1}),
        "D4 fixed-section spectrum changed")
require(d4_cycle_hist == Counter({
    (4,): 2, (2, 2): 3, (2, 1, 1): 2, (1, 1, 1, 1): 1,
}), "D4 cycle atlas changed")
print(f"D4_subgroup_size={len(d4)} fixed_section_hist={sorted(d4_fixed_hist.items())} "
      f"cycle_hist={sorted(d4_cycle_hist.items())}")

# Changing local origins conjugates the affine holonomy and preserves sections.
conjugacy_checks = 0
for gauge in AFF2:
    gauge_inverse = inverse2(gauge)
    for holonomy in AFF2:
        conjugate = compose2(compose2(gauge, holonomy), gauge_inverse)
        require(len(fixed2(conjugate)) == len(fixed2(holonomy)),
                "affine gauge conjugacy changed the section count")
        conjugacy_checks += 1
require(conjugacy_checks == 576, "affine conjugacy check count changed")

# On a two-edge cycle, compatible sections are exactly fixed points of the
# composite holonomy.  This is the finite model of the general proof.
cycle_checks = 0
for first in AFF2:
    for second in AFF2:
        compatible = tuple(
            x for x in V4
            if affine2(second, affine2(first, x)) == x
        )
        require(compatible == fixed2(compose2(second, first)),
                "two-edge sections differ from holonomy fixed points")
        cycle_checks += 1
require(cycle_checks == 576, "two-edge cycle check count changed")
print(f"affine_gauge_conjugacy_checks={conjugacy_checks} "
      f"two_edge_cycle_checks={cycle_checks}")
print("verdict=PASS: cyclic sections are affine-holonomy fixed points; "
      "C13 and quartic V4 have the stated spectra")
print("SCOPE: abstract deck torsors and exact monodromy dictionary; no physical "
      "LRC connector, Keller exclusion, or Jacobian closure")
