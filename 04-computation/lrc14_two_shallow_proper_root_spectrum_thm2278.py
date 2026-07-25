#!/usr/bin/env python3
"""Primary exact checks for THM-2278."""

from fractions import Fraction
from math import ceil


def deep_guard_cap(depth: int) -> Fraction:
    if depth % 2:
        return Fraction(5, 49) + Fraction(5, 49 * 13**depth)
    return Fraction(5, 49) + Fraction(5, 294 * 13**depth)


def danger_pair_cap(gap: int) -> Fraction:
    if gap % 2 == 0:
        return Fraction(1, 49) + Fraction(6, 49 * 13**gap)
    return Fraction(1, 49) + Fraction(5, 588 * 13**gap)


delta_5 = Fraction(961, 6930)

# A proper nonempty 13-root mask can be replaced by its complement, so its
# reduced size n lies in 1..6.  The six conjugate squared magnitudes have
# sum n(13-n)/2 and positive integral product.  The largest possible sum is
# 21.  AM-GM on the other five factors gives the universal floor below.
half_energy_by_size = {
    n: Fraction(n * (13 - n), 2)
    for n in range(1, 7)
}
assert half_energy_by_size == {
    1: 6,
    2: 11,
    3: 15,
    4: 18,
    5: 20,
    6: 21,
}
assert max(half_energy_by_size.values()) == 21

proper_masks = 2**13 - 2
complement_reduced_masks = sum(
    __import__("math").comb(13, n)
    for n in range(1, 7)
)
assert proper_masks == 8190
assert complement_reduced_masks == 4095

cyclotomic_mode_floor = Fraction(5, 21) ** 5
assert cyclotomic_mode_floor == Fraction(3125, 4084101)
assert cyclotomic_mode_floor * Fraction(21, 5) ** 5 == 1

# Recompute THM-2273's strict profile bank and common-time mass floor.
profiles = []
for c in range(5, 20):
    for b in range(2, c):
        shallow_mass = (
            delta_5
            - deep_guard_cap(c)
            - danger_pair_cap(b - 1)
        )
        profiles.append((b, c, shallow_mass))

assert len(profiles) == 150
minimum_profile = min(profiles, key=lambda row: row[2])
assert minimum_profile == (
    3,
    5,
    Fraction(5696989, 367580070),
)

shallow_mass_floor = minimum_profile[2]
transport_factor = Fraction(2197, 460)
common_time_image_floor = transport_factor * shallow_mass_floor
assert common_time_image_floor == Fraction(5696989, 76962600)
assert common_time_image_floor > Fraction(1, 14)

vector_mode_integral_floor = (
    cyclotomic_mode_floor * common_time_image_floor
)
every_residue_energy_floor = vector_mode_integral_floor / 169

assert vector_mode_integral_floor == Fraction(
    712123625, 12572921264904
)
assert every_residue_energy_floor == Fraction(
    712123625, 2124823693768776
)

# Gap-count consequences.  The actual successor speed is
# s=13^(c-b-1)u_3.  The smallest count at each depth offset has u_3=1.
def gap_count(successor_speed: int) -> int:
    return ceil(Fraction(7 * successor_speed, 6) * common_time_image_floor)


minimum_gap_counts = {
    depth_offset: gap_count(13**depth_offset)
    for depth_offset in range(4)
}
assert minimum_gap_counts == {0: 1, 1: 2, 2: 15, 3: 190}

adjacent = [(b, c) for b, c, _ in profiles if c == b + 1]
nonadjacent = [(b, c) for b, c, _ in profiles if c >= b + 2]
boundary_b2 = [(b, c) for b, c, _ in profiles if b == 2]
interior = [
    (b, c)
    for b, c, _ in profiles
    if b >= 3 and c >= b + 2
]

assert len(adjacent) == 15
assert len(nonadjacent) == 135
assert len(boundary_b2) == 15
assert len(interior) == 120
assert min(gap_count(13 ** (c - b - 1)) for b, c in nonadjacent) == 2
assert min(gap_count(13 ** (c - b - 1)) for b, c in boundary_b2) == 15

print("THM-2278 primary exact verification")
print(f"proper_nonempty_root_masks={proper_masks}")
print(f"complement_reduced_root_masks={complement_reduced_masks}")
print(
    "half_energy_by_size="
    f"{dict((size, int(value)) for size, value in half_energy_by_size.items())}"
)
print(f"cyclotomic_mode_floor={cyclotomic_mode_floor}")
print(f"strict_profiles_checked={len(profiles)}")
print(
    "unique_shallow_mass_minimizer="
    f"(1,{minimum_profile[0]},{minimum_profile[1]})"
)
print(f"shallow_mass_floor={shallow_mass_floor}")
print(f"common_time_image_floor={common_time_image_floor}")
print(f"vector_mode_integral_floor={vector_mode_integral_floor}")
print(f"every_residue_energy_floor={every_residue_energy_floor}")
print(f"minimum_gap_counts_by_depth_offset={minimum_gap_counts}")
print(
    "profile_split="
    f"adjacent:{len(adjacent)},"
    f"nonadjacent:{len(nonadjacent)},"
    f"b2:{len(boundary_b2)},"
    f"interior:{len(interior)}"
)
print("exact_frequency_landing=NOT_PROVED")
print("ALL_CHECKS_PASSED")
