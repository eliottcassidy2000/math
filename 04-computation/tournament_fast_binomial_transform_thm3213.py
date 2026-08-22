#!/usr/bin/env python3
"""Exact fast-binomial-transform audit for THM-3213."""

from fractions import Fraction
from math import comb, factorial


def direct_profile(q_value, seed_order, r):
    n = seed_order * r
    row = [q_value(j) ** r for j in range(n + 1)]
    answer = [0] * (n + 1)
    for c in range(1, n + 1):
        row = [row[j + 1] - row[j] for j in range(len(row) - 1)]
        answer[c] = row[0]
    return answer


def convolution_profile(q_value, seed_order, r):
    n = seed_order * r
    left = [Fraction(q_value(j) ** r, factorial(j)) for j in range(n + 1)]
    right = [Fraction((-1) ** j, factorial(j)) for j in range(n + 1)]
    product = [Fraction(0) for _ in range(2 * n + 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j <= n:
                product[i + j] += a * b
    answer = [0] * (n + 1)
    for c in range(1, n + 1):
        scaled = product[c] * factorial(c)
        if scaled.denominator != 1:
            raise RuntimeError((seed_order, r, c, scaled))
        answer[c] = scaled.numerator
    return answer


controls = (
    ("K1", 1, lambda j: j, range(1, 13)),
    ("C3", 3, lambda j: j**3 + 2 * j, range(1, 8)),
)

checks = 0
for name, order, q_value, levels in controls:
    for r in levels:
        direct = direct_profile(q_value, order, r)
        convolution = convolution_profile(q_value, order, r)
        if direct != convolution:
            raise RuntimeError((name, r))
        checks += len(direct) - 1
    print(f"seed={name};levels={levels.start}..{levels.stop - 1}")


def q_from_ordered_profile(profile):
    return lambda j: sum(profile[c] * comb(j, c) for c in range(1, len(profile)))


hostiles = (
    ("mask40", (0, 15, 78, 198, 240, 120)),
    ("mask76", (0, 15, 90, 210, 240, 120)),
)
hostile_rows = {}
for name, profile in hostiles:
    q_value = q_from_ordered_profile(profile)
    for r in range(1, 5):
        direct = direct_profile(q_value, 5, r)
        convolution = convolution_profile(q_value, 5, r)
        if direct != convolution:
            raise RuntimeError((name, r))
        checks += len(direct) - 1
        if r == 1:
            hostile_rows[name] = convolution
    print(f"seed={name};levels=1..4")


def c3_hamilton(profile):
    profile = profile + [0]
    return 3 * sum(
        profile[c] ** 3
        + profile[c + 1] * profile[c] ** 2
        + profile[c + 1] ** 2 * profile[c]
        for c in range(1, len(profile) - 1)
    )


hostile_h = tuple(c3_hamilton(hostile_rows[name]) for name, _ in hostiles)
if hostile_h != (178036299, 193215375):
    raise RuntimeError(hostile_h)

print(f"exact_profile_coordinates_checked={checks}")
print(f"same_H_hostile_C3_lift_H={hostile_h}")
print("identity=F_r(c)/c!=[z^c](sum_j Q_A(j)^r z^j/j!)(sum_l (-z)^l/l!)")
print("algorithm=one truncated polynomial product plus endpoint scaling")

print('scope=complete_factor_moving_jet_and_fixed_output_depths_not_all_growing_cyclic_depths')
print('FAILED_CHECKS=NONE')
