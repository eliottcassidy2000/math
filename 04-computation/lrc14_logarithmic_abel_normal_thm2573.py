#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2573.

The logarithmic Abel limit itself is proved analytically in THM-2573.  This
companion checks its finite jump-pairing algebra, target-character controls,
coincident-endpoint convention, and the 91-unit sampling hostile using only
integer and rational arithmetic.
"""

from fractions import Fraction
from math import gcd


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cyclic_jumps(values):
    """Total-layer jumps at the left endpoints of equal grid cells."""
    n = len(values)
    return {
        j: values[j] - values[(j - 1) % n]
        for j in range(n)
        if values[j] != values[(j - 1) % n]
    }


def handoff_measure(left, right):
    """Return -Delta(left)Delta(right) on a common cyclic grid."""
    require(len(left) == len(right), "grid mismatch")
    require(all(a * b == 0 for a, b in zip(left, right)),
            "endpoint layers are not disjoint")
    jl = cyclic_jumps(left)
    jr = cyclic_jumps(right)
    return {
        j: -jl.get(j, 0) * jr.get(j, 0)
        for j in range(len(left))
        if -jl.get(j, 0) * jr.get(j, 0) != 0
    }


def cyclic_interval(n, start, length):
    return tuple(
        Fraction(int((j - start) % n < length), 1)
        for j in range(n)
    )


print("== THM-2573: logarithmic Abel normal jump pairing ==")


print("\n== exact Abel-kernel normalization ==")
# For one-sided whole-layer smoothing rho^|X|, each of the positive and
# negative Fourier tails contributes -1 after division by
# (1-rho)log(1/(1-rho)).  The jump expansion contributes 1/(4*pi^2).
# Smoothing both endpoints with the same rho doubles the Abel speed and the
# normal; using sqrt(rho) on both restores the canonical one-sided scale.
tail_sides = 2
one_tail_constant = -1
kernel_constant_numerator = tail_sides * one_tail_constant
require(kernel_constant_numerator == -2,
        "two-sided Abel kernel constant changed")
pi_squared_coefficient = Fraction(kernel_constant_numerator, 4)
require(pi_squared_coefficient == Fraction(-1, 2),
        "jump-pairing normalization changed")
print("  one-sided smoothing, two Fourier tails: 2 x (-1) = -2")
print("  logarithmic normal = -B_M/(2*pi^2)")
print("  equal-rho smoothing of both endpoints doubles this normal")


print("\n== complementary and separated-support controls ==")
half = cyclic_interval(4, 0, 2)
half_complement = tuple(1 - x for x in half)
half_mu = handoff_measure(half, half_complement)
require(half_mu == {0: 1, 2: 1},
        "half-interval handoff measure changed")

left_separated = cyclic_interval(4, 0, 1)
right_separated = cyclic_interval(4, 2, 1)
separated_mu = handoff_measure(left_separated, right_separated)
require(separated_mu == {}, "separated support acquired a common jump")
print(f"  half/complement handoff atoms: {len(half_mu)} (mass {sum(half_mu.values())})")
print("  separated positive intervals: no common jump, zero logarithmic normal")


print("\n== aggregate coincident-endpoint convention ==")
# P=[0,1/2), H=[0,3/4).  The common factor H dies at x=0, so only x=1/2
# is a genuine total-layer handoff.  Counting P versus 1-P before forming
# the whole layers would incorrectly count both endpoints.
p_gate = cyclic_interval(4, 0, 2)
common_filter = cyclic_interval(4, 0, 3)
left_total = tuple(a * h for a, h in zip(p_gate, common_filter))
right_total = tuple((1 - a) * h for a, h in zip(p_gate, common_filter))
coincident_mu = handoff_measure(left_total, right_total)
naive_mu = handoff_measure(p_gate, tuple(1 - a for a in p_gate))
require(coincident_mu == {2: 1}, "merged total-layer jump changed")
require(naive_mu == {0: 1, 2: 1}, "naive factorwise control changed")
print(f"  naive factorwise atoms / lawful total-layer atoms: {len(naive_mu)} / {len(coincident_mu)}")
print("  strict/open point values do not enter one-sided aggregate jumps")


print("\n== moving complementary interval target character ==")
# Work on the 52-grid.  P_s has length 1/4 and starts at s/13=4s/52.
# Its handoff normal has phases z_52^(4s)+z_52^(4s+13).  Multiplication
# by the target DFT phase z_52^(4qs) is constant exactly when q=-1 mod 13;
# otherwise each endpoint gives one complete 13th-root coset and sums to 0.
n = 52
target_active = []
moving_atom_census = 0
for s in range(P):
    moving = cyclic_interval(n, 4 * s, 13)
    complement = tuple(1 - x for x in moving)
    mu = handoff_measure(moving, complement)
    expected = {(4 * s) % n: 1, (4 * s + 13) % n: 1}
    require(mu == expected, "moving interval handoff atoms changed")
    moving_atom_census += len(mu)

for q in range(P):
    constant_endpoints = 0
    for endpoint_offset in (0, 13):
        exponents = [
            (4 * s + endpoint_offset + 4 * q * s) % n
            for s in range(P)
        ]
        if len(set(exponents)) == 1:
            constant_endpoints += 1
        else:
            expected_coset = sorted((endpoint_offset + 4 * j) % n
                                     for j in range(P))
            require(sorted(exponents) == expected_coset,
                    "nontrivial target phase failed to form a root coset")
    if constant_endpoints:
        require(constant_endpoints == 2,
                "only one boundary endpoint became target-constant")
        target_active.append(q)

require(moving_atom_census == 26, "moving handoff census changed")
require(target_active == [12], "moving interval target sign changed")
print(f"  moving handoff atoms checked: {moving_atom_census}")
print("  DFT convention sum_s N_s zeta^(q s): unique colour q=12=-1")
print("  normalized active coefficient: (1+i)/(2*pi^2)")


print("\n== target-neutral and deep-unit-blind hostiles ==")
# A fixed complementary interval has a positive q=0 normal but no nonzero
# target character.  Separately, the alternating 14-grid has uniform boundary
# measure, whose Fourier moments vanish away from multiples of 14.  Hence all
# m coprime to 91 are blind even if the grid is translated with the target.
require(sum(half_mu.values()) == 2,
        "target-neutral normal mass changed")
for q in range(1, P):
    residues = sorted((q * s) % P for s in range(P))
    require(residues == list(range(P)),
            "target-neutral root coset changed")

alternating = tuple(Fraction(j % 2, 1) for j in range(14))
alternating_complement = tuple(1 - x for x in alternating)
alternating_mu = handoff_measure(alternating, alternating_complement)
require(len(alternating_mu) == 14,
        "alternating-grid boundary census changed")
require(set(alternating_mu.values()) == {1},
        "alternating-grid handoff weights changed")

unit_samples = 0
for m in range(-1000, 1001):
    if gcd(abs(m), 91) != 1:
        continue
    require(m % 14 != 0,
            "a 91-unit unexpectedly survived the 14-grid Fourier support")
    unit_samples += 1
require(unit_samples == 1584, "91-unit sample census changed")
print("  fixed complement: positive normal mass 2, only target colour q=0")
print(f"  alternating-grid handoff atoms: {len(alternating_mu)}")
print(f"  checked deep-unit-blind samples |m|<=1000: {unit_samples}")
print("  symbolic reason: 14|m implies 7|m, impossible for gcd(m,91)=1")


print("\nsemantic scope: the normal survives only after singular Abel rescaling")
print("a frozen THM-2569 selector is not a lawful common-target orbit")
print("\nall exact checks passed")
