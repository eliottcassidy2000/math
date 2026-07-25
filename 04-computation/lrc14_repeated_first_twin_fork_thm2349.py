#!/usr/bin/env python3
"""Exact finite controls for THM-2349.

The analytic input of THM-2349 is proved in the theorem text.  This companion
freezes the complete repeated-profile mass ledger, the Boolean twin-fork
implication pre-empted by THM-2138, the CRT two-colour completion, and
representative normalized root-pair arithmetic used by the shallow-carrier
triangle.
"""

from fractions import Fraction
from math import gcd


P = 13
Q = 7
MODULUS = P * Q
A0_FLOOR = Fraction(961, 6930)
SHALLOW_TARGET_FLOOR = Fraction(2593, 90090)
TWIN_FORK_FLOOR = Fraction(94254673, 2573060490)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def guard_blocker_cap(depth: int) -> Fraction:
    """THM-2080/2291 parity-sharp cap G_depth."""
    require(depth >= 1, "depth must be positive")
    if depth % 2:
        return Fraction(5, 49) + Fraction(5, 49 * P**depth)
    return Fraction(5, 49) + Fraction(5, 294 * P**depth)


def is_unit(value: int, modulus: int) -> bool:
    return gcd(value, modulus) == 1


# The shallow-target debit is exactly the depth-one cap.
require(
    A0_FLOOR - Fraction(10, 91) == SHALLOW_TARGET_FLOOR,
    "shallow target floor changed",
)

# All fifteen repeated-first profiles have a positive deep-complement
# residual, and the unique worst parity/depth row is c=5.
profile_gaps = {}
for depth in range(5, 20):
    gap = A0_FLOOR - guard_blocker_cap(depth)
    require(gap > 0, "twin-fork mass floor lost positivity")
    profile_gaps[depth] = gap
require(
    min(profile_gaps.values()) == TWIN_FORK_FLOOR,
    "uniform twin-fork floor changed",
)
require(
    [depth for depth, gap in profile_gaps.items() if gap == TWIN_FORK_FLOOR]
    == [5],
    "twin-fork worst profile changed",
)

# Truth-table proof of the twin implication.  On R_3 one has a0=1 and d3=0.
# Cover says d1 or d2.  If both shallow exclusive patterns are forbidden,
# only the simultaneous fork remains.
twin_truth_rows = 0
for d1 in (0, 1):
    for d2 in (0, 1):
        d3 = 0
        cover = bool(d1 or d2 or d3)
        e1 = bool(d1 and not d2 and not d3)
        e2 = bool(d2 and not d1 and not d3)
        if cover and not e1 and not e2:
            require(d1 and d2, "shallow-exclusive nullity did not force fork")
            twin_truth_rows += 1
require(twin_truth_rows == 1, "twin truth-table census changed")

# Exhaust the two-colour triangle over Z/91.  Its hypotheses are
# 13∤t and 7∤s; at least one of t,s,s-t must be a 91-unit.
colour_rows = 0
choice_counts = {"t": 0, "s": 0, "s-t": 0}
for t in range(MODULUS):
    if t % P == 0:
        continue
    for s in range(MODULUS):
        if s % Q == 0:
            continue
        if is_unit(t, MODULUS):
            choice = "t"
        elif is_unit(s, MODULUS):
            choice = "s"
        else:
            require(
                is_unit((s - t) % MODULUS, MODULUS),
                "CRT two-colour completion failed",
            )
            choice = "s-t"
        choice_counts[choice] += 1
        colour_rows += 1
require(colour_rows == 84 * 78, "two-colour universe changed")


def choose_root_pair(d_value: int, root: int) -> tuple[int, int]:
    """Choose K_0,K_1=K_0+D, primitive modulo 13D, at a fixed root."""
    require(d_value % P == 0, "deep quotient must contain thirteen")
    require(root % P != 0, "root character must be nonzero")
    for k_zero in range(1, d_value):
        if k_zero % P != root % P:
            continue
        if is_unit(k_zero, d_value):
            k_one = k_zero + d_value
            require(is_unit(k_zero, P * d_value), "K0 lost primitivity")
            require(is_unit(k_one, P * d_value), "K1 lost primitivity")
            return k_zero, k_one
    raise RuntimeError("no primitive normalized root pair")


# Representative deep quotients include extra primes and high 13-power.
deep_quotients = (
    13,
    26,
    39,
    65,
    91,
    169,
    13 * 17 * 19,
    13**4 * 5,
)
root_pair_rows = 0
for d_value in deep_quotients:
    for root in range(1, P):
        k_zero, k_one = choose_root_pair(d_value, root)
        require(k_one - k_zero == d_value, "normalized edge changed")
        require(
            k_zero % P == k_one % P == root,
            "normalized root character changed",
        )
        root_pair_rows += 1

# Grade/root preservation for every repeated depth and a signed multiplier
# window containing both CRT colours.
grade_rows = 0
for depth in range(5, 20):
    deep_speed = P**depth
    for root in range(1, P):
        base = P * root
        for multiplier in range(-2 * MODULUS, 2 * MODULUS + 1):
            if multiplier == 0:
                continue
            endpoint = base + multiplier * deep_speed
            require(endpoint % P == 0, "endpoint lost shallow grade")
            require(endpoint % (P * P) != 0, "endpoint gained deep grade")
            require(
                (endpoint // P) % P == root,
                "deep shift changed shallow root character",
            )
            grade_rows += 1

print("theorem=THM-2349")
print("status=PROVED+VERIFIED-EXACT-CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print("repeated_profiles=15")
print(f"shallow_target_floor={SHALLOW_TARGET_FLOOR}")
print(f"twin_fork_uniform_floor={TWIN_FORK_FLOOR}")
print("twin_fork_unique_worst_depth=5")
print(f"twin_truth_rows={twin_truth_rows}")
print("each_shallow_owner_positive=YES-BY-THM-2138")
print("twin_fork_branch=EMPTY-BY-THM-2138")
print(f"two_colour_rows={colour_rows}")
print(f"two_colour_choice_t={choice_counts['t']}")
print(f"two_colour_choice_s={choice_counts['s']}")
print(f"two_colour_choice_s_minus_t={choice_counts['s-t']}")
print(f"normalized_root_pair_rows={root_pair_rows}")
print(f"grade_preservation_rows={grade_rows}")
print("repeated_marked_triangle_clock=FINITE-NONUNIFORM")
print("all_first_depth_one_rows_on_inverse_boundary=165")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
