#!/usr/bin/env python3
"""Exact verifier for the five-centre equal-step AP functional in THM-1132.

This script verifies the circle arrangement independently of the hand-sorted
centre-gap formula.  It does not claim that arbitrary killer offsets form an
AP, or that a frozen sigma supplies the missing core-landing/drift lemma.
"""

from fractions import Fraction as F

ONE = F(1)
RADIUS = F(1, 14)
ARC_WIDTH = F(1, 7)
THRESHOLD = F(1, 7)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def arrangement_G(sigma):
    """Largest complementary gap after five open radius-1/14 arcs."""
    sigma %= ONE
    arcs = []
    for m in range(5):
        centre = (m * sigma) % ONE
        lo, hi = centre - RADIUS, centre + RADIUS
        if lo < 0:
            arcs.extend(((lo + ONE, ONE), (F(0), hi)))
        elif hi > ONE:
            arcs.extend(((lo, ONE), (F(0), hi - ONE)))
        else:
            arcs.append((lo, hi))

    merged = []
    for lo, hi in sorted(arcs):
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))

    gaps = []
    for i, (_, hi) in enumerate(merged):
        next_lo = merged[(i + 1) % len(merged)][0]
        if i + 1 == len(merged):
            next_lo += ONE
        gaps.append(next_lo - hi)
    return max(gaps, default=F(0))


def centre_gap_formula(sigma):
    """Hand-sorted largest gap H among {0,sigma,...,4 sigma}."""
    sigma %= ONE
    u = min(sigma, ONE - sigma)
    if u <= F(1, 4):
        return max(u, ONE - 4 * u)
    if u <= F(1, 3):
        return u
    return max(3 * u - ONE, ONE - 2 * u)


def formula_G(sigma):
    return centre_gap_formula(sigma) - ARC_WIDTH


def in_strict_good_bands(sigma):
    sigma %= ONE
    return (
        F(0) <= sigma < F(5, 28)
        or F(2, 7) < sigma < F(5, 14)
        or F(3, 7) < sigma < F(4, 7)
        or F(9, 14) < sigma < F(5, 7)
        or F(23, 28) < sigma < ONE
    )


# Exact regression on many rational arrangements, including all formula walls.
points = {F(j, q) for q in range(1, 85) for j in range(q + 1)}
walls = [
    F(0), F(5, 28), F(2, 7), F(5, 14), F(3, 7), F(4, 7),
    F(9, 14), F(5, 7), F(23, 28), F(1),
]
points.update(walls)
for sigma in sorted(points):
    observed = arrangement_G(sigma)
    predicted = formula_G(sigma)
    require(observed == predicted, f"formula mismatch at sigma={sigma}")
    require((observed > THRESHOLD) == in_strict_good_bands(sigma),
            f"band mismatch at sigma={sigma}")

require(sum((F(5, 28), F(1, 14), F(1, 7), F(1, 14), F(5, 28)), F(0)) == F(9, 14),
        "strict-good measure mismatch")
for sigma in walls[1:-1]:
    require(arrangement_G(sigma) == THRESHOLD,
            f"threshold endpoint mismatch at sigma={sigma}")

minimisers = (F(1, 5), F(2, 5), F(3, 5), F(4, 5))
for sigma in minimisers:
    require(arrangement_G(sigma) == F(2, 35),
            f"minimum mismatch at sigma={sigma}")
require(arrangement_G(F(1, 3)) == F(4, 21), "G(1/3) mismatch")
require(arrangement_G(F(1, 2)) == F(5, 14), "G(1/2) mismatch")

# One exact step-two landing probe; this is an example, not a uniform bridge.
core = (1, 2, 4, 7, 9, 11, 12)
t = F(7, 13)
sigma = (2 * t) % ONE


def norm(x):
    x %= ONE
    return min(x, ONE - x)


core_min = min(norm(p * t) for p in core)
require(sigma == F(1, 13), "step-two sigma probe mismatch")
require(core_min == F(1, 13), "core-safe probe mismatch")
require(arrangement_G(sigma) == F(50, 91), "probe G mismatch")

print("EXACT equal-step five-centre AP functional")
print("H(u)=max(u,1-4u) for 0<=u<=1/4")
print("H(u)=u for 1/4<=u<=1/3")
print("H(u)=max(3u-1,1-2u) for 1/3<=u<=1/2")
print("G(sigma)=H(min(sigma,1-sigma))-1/7")
print("strict G(sigma)>1/7 bands:")
print("  [0,5/28) U (2/7,5/14) U (3/7,4/7)")
print("  U (9/14,5/7) U (23/28,1]")
print("strict-good circle measure=9/14")
print("G=1/7 at sigma=5/28,2/7,5/14,3/7,4/7,9/14,5/7,23/28")
print("global min G=2/35 at sigma=1/5,2/5,3/5,4/5")
print("G(1/3)=4/21")
print("G(1/2)=5/14")
print(f"exact rational arrangement regression points={len(points)}: PASS")
print("step-two probe t=7/13: sigma=1/13, G=50/91, core_min=1/13")
print("scope=equal-step five-centre AP functional only; landing/drift and arbitrary offsets OPEN")
print("tournament=not used: an orientation forgets the metric cyclic gaps determining G")
print("certificate=PASS")
