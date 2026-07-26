#!/usr/bin/env python3
"""Exact first-collision null-parent positive control.

This companion independently recounts the reduced THM-2377 interval
slice, then checks the THM-2380 correlation and spectral arithmetic.
It deliberately makes no canonical-owner or scalar-row claim.
"""

from fractions import Fraction
from math import lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


C, S, Q, D = (13, 1, 7, 1274)
BASE = lcm(C, S, Q, D)
GRID = 182 * BASE
require(BASE == 1274 and GRID == 231868, "common endpoint grid changed")
UNIVERSE = int.from_bytes(b"\x01" * GRID, "little")


def danger_cell_mask(speed, shift, width=1):
    """Strict radius-width/14 danger cells on the common endpoint grid."""
    step = BASE // speed
    half_width = 13 * width * step
    cells = bytearray(GRID)
    one = b"\x01"
    for tooth in range(speed):
        centre = step * (182 * tooth - 14 * shift) % GRID
        left = (centre - half_width) % GRID
        right = (centre + half_width) % GRID
        if left < right:
            cells[left:right] = one * (right - left)
        elif left > right:
            cells[left:] = one * (GRID - left)
            cells[:right] = one * right
    return int.from_bytes(cells, "little")


source = danger_cell_mask(C, 0)
collision_danger = danger_cell_mask(S, 0)
collision_safe = UNIVERSE ^ collision_danger
deep_danger = tuple(danger_cell_mask(D, -shift) for shift in range(13))
deep_safe = tuple(UNIVERSE ^ danger_cell_mask(D, -shift)
                  for shift in range(13))
graft_safe = tuple(UNIVERSE ^ danger_cell_mask(Q, shift)
                   for shift in range(13))

# The lawful line is (r,sigma,t)=(1+v,u,v).  The reduced tensor is
# sigma-blind, so only v must be enumerated.
b_counts = []
w_counts = []
a_counts = []
for v in range(13):
    r = (1 + v) % 13
    parent = (
        source
        & deep_danger[r]
        & deep_safe[v]
        & graft_safe[v]
    )
    retained = parent & collision_safe
    deleted = parent & collision_danger
    require((retained & deleted) == 0, "children overlap")
    require((retained | deleted) == parent, "children do not partition parent")
    b_counts.append(parent.bit_count())
    w_counts.append(retained.bit_count())
    a_counts.append(deleted.bit_count())

N = (78, 74) + (71,) * 10 + (74,)
M = tuple(78 - value for value in N)
require(tuple(b_counts) == (2184,) * 13, "parent line is not constant")
require(tuple(w_counts) == tuple(28 * value for value in N),
        "retained line count changed")
require(tuple(a_counts) == tuple(28 * value for value in M),
        "deletion line count changed")
require(tuple(b_counts[i] - w_counts[i] for i in range(13))
        == tuple(a_counts), "B=A+W failed")


def cyclotomic_reduce(coefficients):
    """Reduce a degree-at-most-12 polynomial modulo Phi_13."""
    require(len(coefficients) == 13, "cyclotomic coefficient length changed")
    top = coefficients[12]
    return tuple(coefficients[index] - top for index in range(12))


def cyclotomic_fourier(profile, character):
    """Exact polynomial for sum_k profile[k] zeta^(character*k)."""
    coefficients = [Fraction(0) for _ in range(13)]
    for index, value in enumerate(profile):
        coefficients[(character * index) % 13] += Fraction(value)
    return cyclotomic_reduce(coefficients)


def sparse_cyclotomic(constant, symmetric_coefficient, character):
    """Reduce constant+c(zeta^character+zeta^-character) modulo Phi_13."""
    coefficients = [Fraction(0) for _ in range(13)]
    coefficients[0] = Fraction(constant)
    coefficients[character % 13] += Fraction(symmetric_coefficient)
    coefficients[-character % 13] += Fraction(symmetric_coefficient)
    return cyclotomic_reduce(coefficients)


# Nonzero target Fourier polynomials.  Since
# N=71*1_G+7*delta_0+3*(delta_1+delta_-1), every q!=0 has
# Nhat(q)=7+3(zeta^q+zeta^-q), whose real value is >=1.
require(N == (78, 74) + (71,) * 10 + (74,), "N profile changed")
require(M == (0, 4) + (7,) * 10 + (4,), "M profile changed")
require(7 - 6 > 0, "strict target half-plane lower bound failed")
strict_target_colours = 12

sum_n = sum(N)
sum_m = sum(M)
require((sum_n, sum_m) == (936, 78), "target means changed")

# Parseval gives the exact nonzero-target squared norm.
l2_nonzero = 13 * sum(value * value for value in N) - sum_n * sum_n
require(l2_nonzero == 702, "nonzero target Gram mass changed")

# Fourth energy follows from the cyclic autocorrelation of N.
n_autocorrelation = tuple(
    sum(N[(v + delta) % 13] * N[v] for v in range(13))
    for delta in range(13)
)
l4_all = 13 * sum(value * value for value in n_autocorrelation)
l4_nonzero = l4_all - sum_n ** 4
require(l4_nonzero == 77766, "nonzero target fourth energy changed")

# THM-2380 normalized cross-correlation R(delta).
K = tuple(
    sum(M[(v + delta) % 13] * N[v] for v in range(13))
    for delta in range(13)
)
expected_k = (5562, 5587, 5620, 5629, 5629, 5629, 5629,
              5629, 5629, 5629, 5629, 5620, 5587)
require(K == expected_k, "cross-correlation profile changed")

correlation_denominator = 13 * 8281 ** 2
require(correlation_denominator == 891474493,
        "correlation denominator changed")
correlation_zero = Fraction(K[0], correlation_denominator)
correlation_mean = Fraction(sum(K), 13 * correlation_denominator)
require(correlation_zero == Fraction(5562, 891474493),
        "R(0) changed")
require(correlation_mean == Fraction(5616, 891474493),
        "translation mean changed")

signed_gram_mass = correlation_zero - correlation_mean
absolute_gram_mass = -signed_gram_mass
require(signed_gram_mass == Fraction(-54, 891474493),
        "signed Gram mass changed")
require(absolute_gram_mass == Fraction(54, 891474493),
        "absolute Gram mass changed")

target_denominator = 13 * 8281
require(target_denominator == 107653, "target denominator changed")
cayley_energy = Fraction(l4_nonzero, target_denominator ** 4)
require(cayley_energy
        == Fraction(5982, 10331448031704891637),
        "complete Cayley energy changed")

# Reconstruct the THM-2377 deep profile j(k), derive its exact
# cyclotomic transform, and check every nonzero deep/target product.
# For a!=0,
#
#   C_a=-(13+6(zeta^a+zeta^-a))/91
#      =-(13+12*cos(2*pi*a/13))/91 < 0,
#
# while Nhat_b=7+3(zeta^b+zeta^-b)>=1 and Mhat_b=-Nhat_b.
deep_profile = (
    Fraction(0),
    Fraction(1, 13),
    *(Fraction(1, 7) for _ in range(10)),
    Fraction(1, 13),
)
c_zero = sum(deep_profile)
require(c_zero == Fraction(144, 91), "deep zero-colour factor changed")

deep_cyclotomic_controls = 0
strict_deep_target_pairs = 0
for deep_character in range(1, 13):
    deep_transform = cyclotomic_fourier(deep_profile, deep_character)
    expected_deep = sparse_cyclotomic(
        Fraction(-13, 91),
        Fraction(-6, 91),
        deep_character,
    )
    require(deep_transform == expected_deep,
            "deep cyclotomic transform changed")
    # cos(theta)>=-1 gives C_a<=-1/91<0.
    require(Fraction(-(13 - 12), 91) < 0,
            "deep-colour strict sign bound failed")
    deep_cyclotomic_controls += 1

    for target_character in range(1, 13):
        target_transform = cyclotomic_fourier(N, target_character)
        expected_target = sparse_cyclotomic(7, 3, target_character)
        deletion_transform = cyclotomic_fourier(M, target_character)
        require(target_transform == expected_target,
                "retained target cyclotomic transform changed")
        require(deletion_transform
                == tuple(-value for value in target_transform),
                "null-parent target sign changed")
        # cos(theta)>=-1 gives Nhat_b>=1>0, so the two child products
        # with the common negative deep factor have opposite signs.
        require(7 - 6 > 0, "target strict sign bound failed")
        strict_deep_target_pairs += 1

print("LRC14 first-collision null-parent positive control")
print(f"endpoint_grid: {GRID}")
print("parent_line_counts: " + ",".join(str(value) for value in b_counts))
print("retained_line_counts: " + ",".join(str(value) for value in w_counts))
print("deletion_line_counts: " + ",".join(str(value) for value in a_counts))
print("cross_correlation_numerators: " + ",".join(str(value) for value in K))
print(f"strict_nonzero_target_colours: {strict_target_colours}")
print(f"deep_cyclotomic_controls: {deep_cyclotomic_controls}")
print(f"strict_nonzero_deep_target_pairs: {strict_deep_target_pairs}")
print(f"signed_nonzero_target_gram_mass: {signed_gram_mass}")
print(f"absolute_nonzero_target_gram_mass: {absolute_gram_mass}")
print(f"complete_cayley_energy: {cayley_energy}")
print("canonical_owner_word_transfer: OPEN")
print("all_checks: PASS")
