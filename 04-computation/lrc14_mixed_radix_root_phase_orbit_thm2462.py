#!/usr/bin/env python3
"""Exact referee for THM-2462.

The script is dependency-free and uses only integer/Fraction arithmetic.
Every truth-bearing check uses ``require`` so optimized Python verifies the
same claims.
"""

from fractions import Fraction as F
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13


def floor_fraction(x):
    return x.numerator // x.denominator


def nearest_integer(x):
    return floor_fraction(x + F(1, 2))


def circle_distance(x):
    return abs(x - nearest_integer(x))


def root_mask_from_phase(rho, phase, radius):
    return frozenset(
        r for r in range(P)
        if circle_distance((phase + rho * r) / P) < radius
    )


def physical_root_mask(speed, y, radius):
    return frozenset(
        r for r in range(P)
        if circle_distance(F(speed) * (y + r) / P) < radius
    )


def ordinary_interval(rho, start):
    n = -rho * start
    return (F(n) - F(1, 2), F(n) - F(1, 14))


def guard_interval(start):
    n = -1 - 8 * start
    return (F(n) - F(1, 2), F(n) - F(1, 7))


def interval_intersection(left, right):
    lo = max(left[0], right[0])
    hi = min(left[1], right[1])
    return (lo, hi) if lo < hi else None


def backward_preimage(source, target, multiplier):
    """Choose the first k with source intersect (target+13k)/b nonempty."""
    sl, sh = source
    tl, th = target
    lower = floor_fraction((multiplier * sl - th) / P) - 2
    upper = floor_fraction((multiplier * sh - tl) / P) + 3
    for k in range(lower, upper + 1):
        pullback = ((tl + P * k) / multiplier, (th + P * k) / multiplier)
        hit = interval_intersection(source, pullback)
        if hit is not None:
            return hit, k
    raise RuntimeError("no mixed-radix preimage")


def affine_preimages(source, target, speed, base_speed=14):
    """All x0 cells with x0 in source and speed*(x0+13)/14 in target mod13."""
    sl, sh = source
    tl, th = target
    value_lo = F(speed) * (sl + P) / base_speed
    value_hi = F(speed) * (sh + P) / base_speed
    lower = floor_fraction((value_lo - th) / P) - 2
    upper = floor_fraction((value_hi - tl) / P) + 3
    result = []
    for k in range(lower, upper + 1):
        pullback = (
            F(base_speed, speed) * (tl + P * k) - P,
            F(base_speed, speed) * (th + P * k) - P,
        )
        hit = interval_intersection(source, pullback)
        if hit is not None:
            result.append((hit, k))
    return result


def midpoint(interval):
    return (interval[0] + interval[1]) / 2


def width(interval):
    return interval[1] - interval[0]


# ---------------------------------------------------------------------------
# 1. Nearest-integer phase laws, both orientations.
# ---------------------------------------------------------------------------

ordinary_phase_checks = 0
guard_phase_checks = 0
for rho in range(1, P):
    inv = pow(rho, -1, P)
    for start in range(P):
        # Negative ordinary branch: N=-rho*start and delta=-1/4.
        n_neg = -rho * start
        phase_neg = F(n_neg) - F(1, 4)
        expected_neg = frozenset((start, (start + inv) % P))
        require(
            root_mask_from_phase(rho, phase_neg, F(1, 14)) == expected_neg,
            f"negative ordinary phase law rho={rho}, start={start}",
        )

        # Positive ordinary branch: N=-rho*start and delta=+1/4.
        phase_pos = F(n_neg) + F(1, 4)
        expected_pos = frozenset((start, (start - inv) % P))
        require(
            root_mask_from_phase(rho, phase_pos, F(1, 14)) == expected_pos,
            f"positive ordinary phase law rho={rho}, start={start}",
        )
        ordinary_phase_checks += 2

        # For the guard, choose N so that the displayed start is fixed.
        n_guard_neg = -rho * start - 1
        guard_neg = F(n_guard_neg) - F(1, 4)
        expected_guard = frozenset((start + j * inv) % P for j in range(4))
        require(
            root_mask_from_phase(rho, guard_neg, F(1, 7)) == expected_guard,
            f"negative guard phase law rho={rho}, start={start}",
        )

        n_guard_pos = -rho * start - 2
        guard_pos = F(n_guard_pos) + F(1, 4)
        require(
            root_mask_from_phase(rho, guard_pos, F(1, 7)) == expected_guard,
            f"positive guard phase law rho={rho}, start={start}",
        )
        guard_phase_checks += 2


# ---------------------------------------------------------------------------
# 2. The explicit thirteen-word atlas.
# ---------------------------------------------------------------------------

ROWS = (
    (0, 9, 7, 3, 12),
    (1, 2, 9, 5, 7),
    (2, 3, 10, 6, 3),
    (3, 2, 10, 6, 7),
    (4, 10, 7, 2, 11),
    (5, 6, 3, 9, 11),
    (6, 2, 9, 2, 7),
    (7, 8, 2, 2, 6),
    (8, 2, 6, 9, 11),
    (9, 8, 3, 2, 7),
    (10, 9, 4, 3, 3),
    (11, 2, 6, 9, 5),
    (12, 8, 6, 2, 11),
)

MULTIPLIERS = (46, 42, 44, 48, 50)
SPEEDS = [14]
for multiplier in MULTIPLIERS:
    SPEEDS.append(SPEEDS[-1] * multiplier)
SPEEDS = tuple(SPEEDS)
RESIDUES = (1, 7, 8, 1, 9, 8)
require(SPEEDS == (14, 644, 27048, 1190112, 57125376, 2856268800),
        "explicit speed chain")
require(tuple(u % P for u in SPEEDS) == RESIDUES, "explicit speed residues")

EXPECTED_CENTERS = (
    F(3034865843, 3332313600),
    F(9058389143, 9996940800),
    F(17899453493, 19993881600),
    F(18124187051, 19993881600),
    F(8996707229, 9996940800),
    F(501517651, 555385600),
    F(362322469, 399877632),
    F(3604728023, 3998776320),
    F(823790063, 908812800),
    F(360620951, 399877632),
    F(18217325957, 19993881600),
    F(9062455079, 9996940800),
    F(3003819443, 3332313600),
)

parent_intervals = []
guard_counts = [0] * P
clean_cover_checks = 0
for row_index, (g, q, a, b, c) in enumerate(ROWS):
    target_intervals = (
        ordinary_interval(1, 0),
        ordinary_interval(7, q),
        guard_interval(g),
        ordinary_interval(1, a),
        ordinary_interval(9, b),
        ordinary_interval(8, c),
    )

    current = target_intervals[-1]
    chosen_k_reversed = []
    for coordinate in range(4, -1, -1):
        current, k = backward_preimage(
            target_intervals[coordinate], current, MULTIPLIERS[coordinate]
        )
        chosen_k_reversed.append(k)
    require(width(current) == F(1, 476044800),
            f"common x0 width row {row_index}")

    x0 = midpoint(current)
    y = (x0 + P) / SPEEDS[0]
    require(y == EXPECTED_CENTERS[row_index], f"midpoint row {row_index}")
    y_interval = ((current[0] + P) / SPEEDS[0],
                  (current[1] + P) / SPEEDS[0])
    require(width(y_interval) == F(1, 6664627200),
            f"parent width row {row_index}")
    parent_intervals.append(y_interval)

    expected_masks = (
        frozenset((0, 1)),
        frozenset((q, (q + 2) % P)),
        frozenset((g, (g + 5) % P, (g + 10) % P, (g + 15) % P)),
        frozenset((a, (a + 1) % P)),
        frozenset((b, (b + 3) % P)),
        frozenset((c, (c + 5) % P)),
    )
    radii = (F(1, 14), F(1, 14), F(1, 7),
             F(1, 14), F(1, 14), F(1, 14))
    actual_masks = tuple(
        physical_root_mask(speed, y, radius)
        for speed, radius in zip(SPEEDS, radii)
    )
    require(actual_masks == expected_masks, f"physical masks row {row_index}")

    # Quarter-point controls remain in the same open phase cell.
    for fraction in (F(1, 4), F(3, 4)):
        y_test = y_interval[0] + fraction * width(y_interval)
        require(
            tuple(
                physical_root_mask(speed, y_test, radius)
                for speed, radius in zip(SPEEDS, radii)
            ) == expected_masks,
            f"whole-cell control row {row_index}, fraction={fraction}",
        )

    gate, q_mask, guard, a_mask, b_mask, c_mask = expected_masks
    require(guard.isdisjoint(q_mask | a_mask | b_mask | c_mask),
            f"guard complete atom row {row_index}")
    union = gate | q_mask | guard | a_mask | b_mask | c_mask
    incidence = [
        sum(r in mask for mask in expected_masks)
        for r in range(P)
    ]
    require(union == frozenset(range(P)), f"clean cover row {row_index}")
    require(sorted(incidence) == [1] * 12 + [2],
            f"excess-one incidence row {row_index}")
    for r in guard:
        guard_counts[r] += 1
    clean_cover_checks += 1

require(guard_counts == [4] * P, "uniform guard density")

sorted_parent = sorted(parent_intervals)
minimum_gap = min(
    right[0] - left[1]
    for left, right in zip(sorted_parent, sorted_parent[1:])
)
center_minimum_gap = min(
    right - left for left, right in zip(
        sorted(EXPECTED_CENTERS), sorted(EXPECTED_CENTERS)[1:]
    )
)
require(center_minimum_gap == F(3341, 102009600), "minimum center gap")
require(minimum_gap > 0, "disjoint parent intervals")


# ---------------------------------------------------------------------------
# 3. Two appended root-constant labelled bits.
# ---------------------------------------------------------------------------


def blocker_bit_interval(danger):
    if danger:
        return (F(-1, 28), F(1, 28))
    return (F(1, 7), F(3, 14))


WORD_BITS = (0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0)
LABEL_MULTIPLIERS = MULTIPLIERS + (1305, 1302)
LABEL_UNIT_SPEEDS = [SPEEDS[0]]
for multiplier in LABEL_MULTIPLIERS:
    LABEL_UNIT_SPEEDS.append(LABEL_UNIT_SPEEDS[-1] * multiplier)
LABEL_UNIT_SPEEDS = tuple(LABEL_UNIT_SPEEDS)
require(tuple(u % P for u in LABEL_UNIT_SPEEDS) == RESIDUES + (1, 2),
        "labelled speed residues")

labelled_bit_checks = 0
labelled_width = None
for row_index, (g, q, a, b, c) in enumerate(ROWS):
    desired_bits = (g % 2 == 0, bool(WORD_BITS[row_index]))
    targets = (
        ordinary_interval(1, 0),
        ordinary_interval(7, q),
        guard_interval(g),
        ordinary_interval(1, a),
        ordinary_interval(9, b),
        ordinary_interval(8, c),
        blocker_bit_interval(desired_bits[0]),
        blocker_bit_interval(desired_bits[1]),
    )
    current = targets[-1]
    for coordinate in range(6, -1, -1):
        current, _ = backward_preimage(
            targets[coordinate], current, LABEL_MULTIPLIERS[coordinate]
        )
    if labelled_width is None:
        labelled_width = width(current)
    require(width(current) == labelled_width, f"labelled width row {row_index}")
    y = (midpoint(current) + P) / LABEL_UNIT_SPEEDS[0]

    # Physical blocker speeds are thirteen times the appended unit speeds.
    for desired, unit_speed in zip(desired_bits, LABEL_UNIT_SPEEDS[-2:]):
        blocker_speed = P * unit_speed
        values = {
            circle_distance(F(blocker_speed) * (y + r) / P) < F(1, 14)
            for r in range(P)
        }
        require(values == {desired}, f"root-constant bit row {row_index}")
        labelled_bit_checks += 1

require(labelled_width == F(1, 4853114880768000), "labelled common width")
require(labelled_bit_checks == 26, "labelled bit census")


# ---------------------------------------------------------------------------
# 4. Primitive affine speed control.
# ---------------------------------------------------------------------------

PRIMITIVE_SPEEDS = (14, 657, 27607, 1214721, 58306621, 2915331063)
require(tuple(u % P for u in PRIMITIVE_SPEEDS) == RESIDUES,
        "primitive speed residues")
require(gcd(*PRIMITIVE_SPEEDS) == 1, "primitive speed gcd")

primitive_final_cells = 0
primitive_max_cells = 0
for row_index, (g, q, a, b, c) in enumerate(ROWS):
    targets = (
        ordinary_interval(1, 0),
        ordinary_interval(7, q),
        guard_interval(g),
        ordinary_interval(1, a),
        ordinary_interval(9, b),
        ordinary_interval(8, c),
    )
    states = [(targets[0], ())]
    for target, speed in zip(targets[1:], PRIMITIVE_SPEEDS[1:]):
        next_states = []
        for interval, address in states:
            for pulled, k in affine_preimages(interval, target, speed):
                next_states.append((pulled, address + (k,)))
        states = next_states
        primitive_max_cells = max(primitive_max_cells, len(states))
        require(states, f"primitive phase cell row {row_index}")
    primitive_final_cells += len(states)

    states.sort(key=lambda item: (item[0][0], item[0][1], item[1]))
    interval, _ = states[0]
    y = (midpoint(interval) + P) / PRIMITIVE_SPEEDS[0]
    expected_masks = (
        frozenset((0, 1)),
        frozenset((q, (q + 2) % P)),
        frozenset((g, (g + 5) % P, (g + 10) % P, (g + 15) % P)),
        frozenset((a, (a + 1) % P)),
        frozenset((b, (b + 3) % P)),
        frozenset((c, (c + 5) % P)),
    )
    radii = (F(1, 14), F(1, 14), F(1, 7),
             F(1, 14), F(1, 14), F(1, 14))
    require(
        tuple(
            physical_root_mask(speed, y, radius)
            for speed, radius in zip(PRIMITIVE_SPEEDS, radii)
        ) == expected_masks,
        f"primitive direct masks row {row_index}",
    )

require(primitive_final_cells == 56, "primitive final-cell census")
require(primitive_max_cells == 16, "primitive maximum-cell census")


print("THM2462 MIXED-RADIX ROOT-PHASE AUDIT")
print(f"ordinary_phase_checks={ordinary_phase_checks}")
print(f"guard_phase_checks={guard_phase_checks}")
print("speed_chain=" + ",".join(map(str, SPEEDS)))
print("speed_residues=" + ",".join(map(str, RESIDUES)))
print(f"atlas_rows={len(ROWS)} clean_cover_checks={clean_cover_checks}")
print(f"x0_width={F(1,476044800)}")
print(f"parent_width={F(1,6664627200)}")
print(f"minimum_center_gap={center_minimum_gap}")
print("guard_counts=" + ",".join(map(str, guard_counts)))
print(f"normalized_guard_density={F(4,13)}")
print(f"labelled_bit_checks={labelled_bit_checks}")
print(f"labelled_x0_width={labelled_width}")
print("primitive_speeds=" + ",".join(map(str, PRIMITIVE_SPEEDS)))
print(f"primitive_gcd={gcd(*PRIMITIVE_SPEEDS)}")
print(f"primitive_final_cells={primitive_final_cells}")
print(f"primitive_max_cells={primitive_max_cells}")
print("mixed_radix_realization=PASS")
print("optimized_require_path=PASS")
