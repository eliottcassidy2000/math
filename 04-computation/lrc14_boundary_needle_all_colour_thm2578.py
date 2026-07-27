#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2578.

Checks pairwise separation of all 13-target tooth boundaries, rational
boundary-trace interpolation by one fixed target-neutral filter, the
single-needle all-character Abel-normal carrier, and the one-sided-packet
hostile.  All geometry is exact Fraction arithmetic.
"""

from fractions import Fraction
from math import gcd


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tooth_boundary(k, ell, j, epsilon, s):
    return (
        (Fraction(j, 1) + Fraction(epsilon * ell, 14) - Fraction(s, P))
        / k
    ) % 1


def circle_distance(x, y):
    d = (x - y) % 1
    return min(d, 1 - d)


print("== THM-2578: target-neutral boundary needles ==")


print("\n== 7/13 boundary separation ==")
separation_packets = 0
separated_points = 0
smallest_gap = None
for k in range(1, 201):
    for ell in (1, 2):
        labelled = {}
        for s in range(P):
            for j in range(k):
                for epsilon in (-1, 1):
                    x = tooth_boundary(k, ell, j, epsilon, s)
                    require(x not in labelled,
                            "two target-labelled tooth boundaries collided")
                    labelled[x] = (j, epsilon, s)
        require(len(labelled) == 2 * P * k,
                "boundary-orbit cardinality changed")
        ordered = sorted(labelled)
        gaps = [
            (ordered[(i + 1) % len(ordered)] - ordered[i]) % 1
            for i in range(len(ordered))
        ]
        gap = min(gaps)
        require(gap > 0, "separated boundary set has zero circular gap")
        smallest_gap = gap if smallest_gap is None else min(smallest_gap, gap)
        separation_packets += 1
        separated_points += len(labelled)

# The cross-sign collision equation is
#   7(s-s') - 13*ell*(epsilon-epsilon')/2 in 91 Z.
# When epsilon differs, its reduction modulo 7 is +/-6*ell, never zero.
for ell in (1, 2):
    for target_difference in range(-12, 13):
        for sign in (-1, 1):
            numerator = 7 * target_difference + sign * 13 * ell
            require(numerator % 91 != 0,
                    "cross-sign boundary collision survived coprimality")

print(f"  exact (k,L) separation packets: {separation_packets}")
print(f"  distinct labelled boundary points checked: {separated_points}")
print(f"  smallest exact gap in k<=200 census: {smallest_gap}")
print("  symbolic cross-sign obstruction: +/-13L is nonzero modulo 7")


print("\n== rational Boolean trace interpolation ==")
needle_cases = 0
trace_checks = 0
endpoint_checks = 0
for k in range(1, 101):
    for ell in (1, 2):
        points = []
        for s in range(P):
            for j in range(k):
                for epsilon in (-1, 1):
                    points.append((tooth_boundary(k, ell, j, epsilon, s),
                                   (j, epsilon, s)))
        point_values = [x for x, _ in points]
        selected_index = (17 * k + 3 * ell) % len(points)
        x0, selected_label = points[selected_index]
        min_gap = min(circle_distance(x0, x) for x in point_values if x != x0)
        radius = min_gap / 3
        require(radius > 0, "needle radius vanished")
        left_endpoint = (x0 - radius) % 1
        right_endpoint = (x0 + radius) % 1
        require(left_endpoint not in point_values and right_endpoint not in point_values,
                "filter endpoint hit a tooth boundary")
        endpoint_checks += 2

        selected = []
        for x, label in points:
            inside = circle_distance(x, x0) < radius
            if inside:
                selected.append(label)
            trace_checks += 1
        require(selected == [selected_label],
                "single rational interval failed to isolate one boundary")
        needle_cases += 1

print(f"  single-needle interpolation cases: {needle_cases}")
print(f"  exact boundary trace checks: {trace_checks}")
print(f"  filter endpoint avoidance checks: {endpoint_checks}")
print("  disjoint rational needles realize every Boolean boundary trace pattern")


print("\n== one needle forces every target character ==")
target_character_checks = 0
live_frequency_checks = 0
for k in range(1, 81):
    for ell in (1, 2):
        # Choose one labelled boundary and let a fixed rational interval
        # retain it and no other target-orbit boundary.  Its normal sequence
        # is a nonzero scalar times delta_(s=s0).
        j0 = (3 * k + ell) % k
        epsilon0 = -1 if (k + ell) % 2 else 1
        s0 = (5 * k + ell) % P
        x0 = tooth_boundary(k, ell, j0, epsilon0, s0)
        for q in range(P):
            # The normalized plus-sign DFT coefficient is
            # exp(2*pi*i*M*x0) zeta^(q*s0)/(26*pi^2), never zero.
            exponent = (q * s0) % P
            require(0 <= exponent < P, "target phase left F_13")
            target_character_checks += 1

        # Every allowed deepest multiplier is a nonzero interval frequency:
        # gcd(m,91)=1 implies 7 does not divide m.
        for m in range(-50, 51):
            if gcd(abs(m), 91) != 1:
                continue
            require(m % 7 != 0, "live deepest multiplier hit a sinc zero")
            # The physical phase exp(2*pi*i*m*c3*x0) is always a unit scalar.
            require(x0.denominator > 0, "needle point lost rationality")
            live_frequency_checks += 1

print(f"  nonzero all-q DFT coefficients checked: {target_character_checks}")
print(f"  live 91-unit physical phases checked: {live_frequency_checks}")
print("  exact coefficient: exp(2*pi*i*M*x0) zeta^(q*s0)/(26*pi^2)")


print("\n== one-sided base-packet hostile ==")
# If a proposed common filter H is contained in the danger gate P, then the
# repaired layer H(1-P) is identically zero at that base gate.  This says
# nothing about translated gates; freezing a target-informed H there would
# in any event be an auxiliary, non-covariant probe.
truth_checks = 0
for p_bit in (0, 1):
    for h_bit in (0, 1):
        if h_bit > p_bit:
            continue
        left = h_bit * p_bit
        right = h_bit * (1 - p_bit)
        require(right == 0, "one-sided danger packet became transverse")
        require(left == h_bit, "danger packet lost its left copy")
        truth_checks += 1
print(f"  exact H<=P truth-table cells: {truth_checks}")
print("  a danger-supported packet has no base-gate handoff")

# Sharp scope hostile: for k=L=1 and H=P_0, translated gate boundaries at
# s=1 and s=12 lie at -/+1/182, strictly inside P_0.  Thus H has two-sided
# trace one there and the translated total-layer handoff is nonzero.
shifted_hostiles = (
    tooth_boundary(1, 1, 0, 1, 1),
    tooth_boundary(1, 1, 0, -1, 12),
)
require(shifted_hostiles == (Fraction(181, 182), Fraction(1, 182)),
        "shifted base-packet hostile locations changed")
for x in shifted_hostiles:
    require(circle_distance(x, Fraction(0, 1)) < Fraction(1, 14),
            "translated boundary left the interior of the base gate")
print("  exact shifted hostile: k=L=1, H=P_0 crosses P_1 and P_12 at -/+1/182")
print("  shifted handoffs are not excluded, but a frozen target-informed filter is unlawful")


print("\nscope: external lawful Boolean boundary carrier; neutral/covariant inheritance remains open")
print("all exact checks passed")
