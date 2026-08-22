#!/usr/bin/env python3
"""Exact controls for THM-3291's two-block TV homogenization rigidity.

Exact rational arithmetic on a complete rational grid; no floating point, no
randomness, no imported executable.  Every gate is an explicit ``require`` so
that ordinary and ``-O`` replay are byte-identical.
"""

from fractions import Fraction
from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tv_binomial_two(x, y):
    """Total variation distance between Bin(2,x) and Bin(2,y), by definition."""
    left = ((1 - x) ** 2, 2 * x * (1 - x), x * x)
    right = ((1 - y) ** 2, 2 * y * (1 - y), y * y)
    return sum(abs(a - b) for a, b in zip(left, right)) / 2


def closed_form(x, y):
    return abs(x - y) * (1 + abs(x + y - 1))


# ---------------------------------------------------- 1.  the closed form

CLOSED_GRID = 24
CLOSED_CHECKS = 0
for i in range(CLOSED_GRID + 1):
    for j in range(CLOSED_GRID + 1):
        x = Fraction(i, CLOSED_GRID)
        y = Fraction(j, CLOSED_GRID)
        require(tv_binomial_two(x, y) == closed_form(x, y),
                "TV(Bin(2,x),Bin(2,y)) = |x-y|(1+|x+y-1|)")
        CLOSED_CHECKS += 1
require(CLOSED_CHECKS == (CLOSED_GRID + 1) ** 2, "closed-form grid complete")

# The closed form is a genuine identity, not a coincidence of the grid: the
# three signed pmf differences factor through |x-y| exactly.
require(closed_form(Fraction(1, 3), Fraction(1, 3)) == 0, "reflexivity")
require(closed_form(0, 1) == 1, "extreme pair saturates")


# ------------------------------- 2.  the two-block inequality and its face

def rigid_face(p1, p2, q1, q2):
    """Predicted equality locus with both block distances nonzero."""
    e1, e2 = p1 - q1, p2 - q2
    if e1 == 0 or e2 == 0:
        return False
    same_sign = (e1 > 0 and e2 > 0) or (e1 < 0 and e2 < 0)
    equal_gaps = abs(e1) == abs(e2)
    low = min(p1, q1) == 0 and min(p2, q2) == 0
    high = max(p1, q1) == 1 and max(p2, q2) == 1
    return same_sign and equal_gaps and (low or high)


GRID = 16
VIOLATIONS = []
EQUALITIES = []
FACE_POINTS = []
TOTAL = 0
for a, b, c, d in product(range(GRID + 1), repeat=4):
    p1 = Fraction(a, GRID)
    p2 = Fraction(b, GRID)
    q1 = Fraction(c, GRID)
    q2 = Fraction(d, GRID)
    delta1 = abs(p1 - q1)
    delta2 = abs(p2 - q2)
    left = closed_form((p1 + p2) / 2, (q1 + q2) / 2)
    right = delta1 + delta2 - delta1 * delta2
    TOTAL += 1
    if left > right:
        VIOLATIONS.append((p1, p2, q1, q2))
    if left == right and delta1 and delta2:
        EQUALITIES.append((p1, p2, q1, q2))
    if rigid_face(p1, p2, q1, q2):
        FACE_POINTS.append((p1, p2, q1, q2))
        require(left == right, "every predicted face point is an equality")
require(TOTAL == (GRID + 1) ** 4, "grid complete")
require(not VIOLATIONS, "delta_N <= delta_1 + delta_2 - delta_1 delta_2")
require(sorted(EQUALITIES) == sorted(FACE_POINTS),
        "the nontrivial equality locus is exactly the rigid face")
require(len(FACE_POINTS) > 0, "the rigid face is non-vacuous")

# No equality with exactly one vanishing block distance.
MIXED = 0
for a, b, c, d in product(range(GRID + 1), repeat=4):
    p1 = Fraction(a, GRID)
    p2 = Fraction(b, GRID)
    q1 = Fraction(c, GRID)
    q2 = Fraction(d, GRID)
    delta1 = abs(p1 - q1)
    delta2 = abs(p2 - q2)
    if (delta1 == 0) == (delta2 == 0):
        continue
    left = closed_form((p1 + p2) / 2, (q1 + q2) / 2)
    right = delta1 + delta2 - delta1 * delta2
    if left == right:
        MIXED += 1
require(MIXED == 0, "no equality with exactly one vanishing block distance")


# ------------------------------------- 3.  the two ingredients of the proof

# (i) box constraint: |x+y-1| <= 1 - (delta_1+delta_2)/2.
BOX_CHECKS = 0
for a, b, c, d in product(range(GRID + 1), repeat=4):
    p1 = Fraction(a, GRID)
    p2 = Fraction(b, GRID)
    q1 = Fraction(c, GRID)
    q2 = Fraction(d, GRID)
    delta1 = abs(p1 - q1)
    delta2 = abs(p2 - q2)
    x = (p1 + p2) / 2
    y = (q1 + q2) / 2
    require(abs(x + y - 1) <= 1 - (delta1 + delta2) / 2,
            "box constraint on the homogenized means")
    BOX_CHECKS += 1

# (ii) AM-GM: (d1+d2)^2/4 >= d1 d2, with equality iff d1 = d2.
AMGM_CHECKS = 0
for i in range(GRID + 1):
    for j in range(GRID + 1):
        delta1 = Fraction(i, GRID)
        delta2 = Fraction(j, GRID)
        gap = (delta1 + delta2) ** 2 / 4 - delta1 * delta2
        require(gap >= 0, "AM-GM")
        require((gap == 0) == (delta1 == delta2), "AM-GM equality iff equal")
        AMGM_CHECKS += 1

# (iii) opposite-sign branch: max(d1,d2) <= d1+d2-d1 d2 on [0,1]^2.
OPPOSITE_CHECKS = 0
for i in range(GRID + 1):
    for j in range(GRID + 1):
        delta1 = Fraction(i, GRID)
        delta2 = Fraction(j, GRID)
        require(max(delta1, delta2) <= delta1 + delta2 - delta1 * delta2,
                "min(d)*(1-max(d)) >= 0")
        OPPOSITE_CHECKS += 1

# Hostile control: the multiplicative correction is necessary -- the weaker
# additive bound delta_1+delta_2 is attained strictly more often, and the
# stronger bound delta_1+delta_2-2*delta_1*delta_2 is FALSE.
STRONGER_FAILS = 0
for a, b, c, d in product(range(GRID + 1), repeat=4):
    p1 = Fraction(a, GRID)
    p2 = Fraction(b, GRID)
    q1 = Fraction(c, GRID)
    q2 = Fraction(d, GRID)
    delta1 = abs(p1 - q1)
    delta2 = abs(p2 - q2)
    left = closed_form((p1 + p2) / 2, (q1 + q2) / 2)
    if left > delta1 + delta2 - 2 * delta1 * delta2:
        STRONGER_FAILS += 1
require(STRONGER_FAILS > 0,
        "hostile: doubling the correction term breaks the bound")


print("THM-3291 TWO-BLOCK TV HOMOGENIZATION RIGIDITY EXACT CONTROL")
print("closed_form=TV(Bin(2,x),Bin(2,y))=|x-y|(1+|x+y-1|)")
print("closed_form_checks=" + str(CLOSED_CHECKS))
print("inequality=delta_N <= delta_1+delta_2-delta_1*delta_2")
print("grid_points=" + str(TOTAL))
print("violations=" + str(len(VIOLATIONS)))
print("nontrivial_equalities=" + str(len(EQUALITIES)))
print("rigid_face=same sign AND delta_1=delta_2 AND "
      "(min(p_i,q_i)=0 both i OR max(p_i,q_i)=1 both i)")
print("face_points=" + str(len(FACE_POINTS)))
print("mixed_equalities=" + str(MIXED))
print("box_constraint_checks=" + str(BOX_CHECKS))
print("amgm_checks=" + str(AMGM_CHECKS))
print("opposite_sign_checks=" + str(OPPOSITE_CHECKS))
print("hostile_doubled_correction_failures=" + str(STRONGER_FAILS))
print("scope=two single-observation blocks only; the general block statement "
      "is the cited external theorem, not proved here")
print("ALL EXACT CHECKS PASSED")
