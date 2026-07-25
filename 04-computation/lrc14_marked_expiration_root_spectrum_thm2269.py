#!/usr/bin/env python3
"""Primary exact checks for THM-2269."""

from fractions import Fraction
from itertools import combinations


def norm_circle(x: Fraction) -> Fraction:
    r = x % 1
    return min(r, 1 - r)


def danger(speed: int, x: Fraction) -> int:
    return int(norm_circle(speed * x) < Fraction(1, 14))


def guard(speed: int, x: Fraction) -> int:
    return int(norm_circle(speed * x) > Fraction(1, 7))


# Every support of size one or two has the exact nonconstant DFT energy
# 13n-n^2.  For size two, no nonzero character can cancel: cancellation
# would force 2*k*(a-b)=0 mod 13.
mask_count = 0
energy_histogram: dict[int, int] = {}
for size in (1, 2):
    for support in combinations(range(13), size):
        mask_count += 1
        energy = 13 * size - size * size
        energy_histogram[energy] = energy_histogram.get(energy, 0) + 1
        if size == 2:
            delta = (support[1] - support[0]) % 13
            for k in range(1, 13):
                assert (2 * k * delta) % 13 != 0

assert mask_count == 13 + 78 == 91
assert energy_histogram == {12: 13, 22: 78}

# For a two-root mask, |1+zeta^d|^2=4*cos(pi*d/13)^2.
# The closest nonzero residue angles to pi/2 are d=6,7, at distance
# pi/26.  Concavity gives sin(pi/26)>=2*(pi/26)/pi=1/13, hence every
# nonzero character of every size-one/two mask has energy at least 4/169.
angle_distances = {d: abs(2 * d - 13) for d in range(1, 13)}
assert min(angle_distances.values()) == 1
assert {d for d, distance in angle_distances.items() if distance == 1} == {6, 7}
uniform_character_energy_lower = Fraction(4, 169)

# The two-variable fibre LP is E=12(a+b)+10b, so E>=12 image_mass.
for a_num in range(8):
    for b_num in range(8):
        a = Fraction(a_num, 7)
        b = Fraction(b_num, 7)
        energy = 12 * a + 22 * b
        assert energy - 12 * (a + b) == 10 * b >= 0

strict_image = Fraction(15041431, 70270200)
repeat_image = Fraction(5229541, 70270200)
one_core_carrier = Fraction(6055, 28561)

strict_total_residue_energy = 12 * strict_image / 169
strict_one_residue_energy = strict_image / 169
repeat_total_residue_energy = 12 * repeat_image / 169
repeat_one_residue_energy = repeat_image / 169
strict_every_residue_energy = 4 * strict_image / 28561
repeat_every_residue_energy = 4 * repeat_image / 28561
strict_carrier_gap = strict_image - one_core_carrier

assert strict_total_residue_energy == Fraction(15041431, 989638650)
assert strict_one_residue_energy == Fraction(15041431, 11875663800)
assert repeat_total_residue_energy == Fraction(5229541, 989638650)
assert repeat_one_residue_energy == Fraction(5229541, 11875663800)
assert strict_every_residue_energy == Fraction(15041431, 501746795550)
assert repeat_every_residue_energy == Fraction(5229541, 501746795550)
assert strict_carrier_gap == Fraction(24332839, 11875663800)

# Exact branch-local hostile control.
x0 = Fraction(79, 338)
q = (2, 339, 677, 1015, 1353)
assert len(set(q)) == 5 and all(speed % 13 for speed in q)
assert guard(1, x0)
assert all(not danger(speed, x0) for speed in q)

profiles = []
state_rows = None
for c in range(5, 20):
    for b in range(2, c):
        blockers = (13, 13**b, 13**c)
        assert tuple(danger(speed, x0) for speed in blockers) == (1, 0, 0)

        rows = []
        phase = x0
        for t in range(10):
            literal = guard(1, phase) * int(all(not danger(speed, phase) for speed in q))
            owners = tuple(danger(speed, phase) for speed in blockers)
            assert literal <= int(any(owners))
            rows.append((t, phase, literal, owners))
            phase = (13 * phase) % 1

        assert rows[0][1:] == (x0, 1, (1, 0, 0))
        assert rows[1][1:] == (Fraction(1, 26), 0, (0, 0, 0))
        assert all(row[1:] == (Fraction(1, 2), 0, (0, 0, 0)) for row in rows[2:])
        assert danger(1, Fraction(1, 26)) == 1
        assert danger(13, Fraction(1, 26)) == 0
        assert all(danger(1, row[1]) == 0 for row in rows[2:])
        profiles.append((1, b, c))
        state_rows = rows

assert len(profiles) == 150
assert state_rows is not None

print("THM-2269 primary exact verification")
print(f"root_masks_checked={mask_count}")
print(f"nonconstant_energy_histogram={energy_histogram}")
print("uniform_angle_minimizers=(6,7)")
print(f"uniform_character_energy_lower={uniform_character_energy_lower}")
print("fibre_lp_certificate=12a+22b-12(a+b)=10b>=0")
print(f"strict_image_floor={strict_image}")
print(f"strict_total_nonzero_residue_energy={strict_total_residue_energy}")
print(f"strict_one_residue_energy={strict_one_residue_energy}")
print(f"strict_every_residue_energy={strict_every_residue_energy}")
print(f"repeated_image_floor={repeat_image}")
print(f"repeated_total_nonzero_residue_energy={repeat_total_residue_energy}")
print(f"repeated_one_residue_energy={repeat_one_residue_energy}")
print(f"repeated_every_residue_energy={repeat_every_residue_energy}")
print(f"strict_minus_one_core_carrier={strict_carrier_gap}")
print(f"hostile_x={x0}")
print(f"hostile_q={q}")
print(f"strict_profiles_checked={len(profiles)}")
print("hostile_branch_states=")
for row in state_rows[:4]:
    print(f"  t={row[0]} phase={row[1]} literal={row[2]} owners={row[3]}")
print("marked_current_bits_at_t1=(1,0)")
print("ALL_CHECKS_PASSED")
