#!/usr/bin/env python3
"""Exact finite controls for THM-2367.

This dependency-free companion checks the root-average profiles, the
thirteen-window converse, two independent circulant hostiles, the full
165-profile valuation ledger, and the additive order-seven cover constraints.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac(value: Fraction) -> Fraction:
    return value - value.numerator // value.denominator


def circle_norm(value: Fraction) -> Fraction:
    residue = frac(value)
    return min(residue, 1 - residue)


def danger(value: Fraction, length_units: int = 1) -> bool:
    return circle_norm(value) < Fraction(length_units, 14)


def safe(value: Fraction, length_units: int = 1) -> bool:
    return not danger(value, length_units)


def interval_intersection(
    left: Fraction,
    right: Fraction,
    other_left: Fraction,
    other_right: Fraction,
) -> Fraction:
    return max(
        Fraction(0),
        min(right, other_right) - max(left, other_left),
    )


def circular_arc_mass(
    start: Fraction, length: Fraction, delta: Fraction
) -> Fraction:
    """Mass of [start,start+length] mod 1 inside ||x||<delta."""
    start = frac(start)
    end = start + length
    source = (
        [(start, end)]
        if end <= 1
        else [(start, Fraction(1)), (Fraction(0), end - 1)]
    )
    target = [(Fraction(0), delta), (1 - delta, Fraction(1))]
    return sum(
        interval_intersection(a, b, c, d)
        for a, b in source
        for c, d in target
    )


def valuation(value: int, prime: int) -> int:
    require(value > 0, "valuation requires a positive integer")
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


# ---------------------------------------------------------------------------
# Root-average step formula for unit-safe and guard-safe grafts.
# ---------------------------------------------------------------------------

root_average_cases = 0
nonconstant_root_cases: list[tuple[int, int, Fraction]] = []
for length_units in (1, 2):
    for c_value in range(1, 15):
        for cell in range(14):
            y_value = Fraction(2 * cell + 1, 28)
            actual = Fraction(
                sum(
                    safe(
                        Fraction(y_value + root, c_value),
                        length_units,
                    )
                    for root in range(c_value)
                ),
                c_value,
            )
            if c_value % 7 == 0:
                expected = Fraction(7 - length_units, 7)
            else:
                residue = length_units * c_value % 14
                fold = min(residue, 14 - residue)
                require(1 <= fold <= 6, "wrong folded root residue")
                delta = Fraction(fold, 14)
                epsilon = 1 if 1 <= residue <= 6 else -1
                centered_bit = int(circle_norm(y_value) < delta)
                expected = (
                    Fraction(7 - length_units, 7)
                    - Fraction(epsilon, c_value)
                    * (centered_bit - 2 * delta)
                )
                nonconstant_root_cases.append(
                    (length_units, c_value, delta)
                )
            require(actual == expected, "root-average step formula failed")
        root_average_cases += 1
require(root_average_cases == 28, "wrong root-average case count")


# ---------------------------------------------------------------------------
# Exact thirteen-window proof of noncirculancy.
# ---------------------------------------------------------------------------

window_controls = 0
for length_units in (1, 2):
    for c_residue in range(1, 15):
        if c_residue % 7 == 0:
            continue
        residue = length_units * c_residue % 14
        fold = min(residue, 14 - residue)
        delta = Fraction(fold, 14)
        require(0 < delta < Fraction(1, 2), "wrong window radius")
        for rho in range(1, P):
            x_zero = Fraction(rho, 14)
            cell_masses = tuple(
                circular_arc_mass(
                    x_zero + Fraction(index, P),
                    Fraction(1, P),
                    delta,
                )
                for index in range(P)
            )
            require(
                len(set(cell_masses)) > 1,
                "offset grid unexpectedly has equal cell masses",
            )
            window_sums = tuple(
                sum(
                    cell_masses[(index + offset) % P]
                    for offset in range(rho)
                )
                for index in range(P)
            )
            require(
                len(set(window_sums)) > 1,
                "thirteen-window drift witness disappeared",
            )
            window_controls += 1
require(window_controls == 288, "wrong converse-window control count")


def isolated_circulant(c_value: int, q_value: int) -> bool:
    return valuation(c_value, 7) > valuation(q_value, 7)


require(isolated_circulant(91, 1), "(91,1) hostile changed")
require(not isolated_circulant(91, 7), "(91,7) escape changed")
require(not isolated_circulant(13, 1), "(13,1) escape changed")


# ---------------------------------------------------------------------------
# Exact 90/91-mass Boolean cancellation hostile.
# ---------------------------------------------------------------------------

D = 16562
EXCLUDED = (
    (16555, 16562),
    (0, 13),
    (1625, 1651),
    (2457, 2463),
    (3263, 3289),
    (4907, 4927),
    (7449, 7455),
    (9087, 9113),
    (10725, 10751),
    (12363, 12389),
)
excluded_cells = sum(right - left for left, right in EXCLUDED)
require(excluded_cells == 182, "Boolean hostile deletion count changed")

segment_string = (
    "D=16562;excluded=[16555,16562);[0,13);[1625,1651);"
    "[2457,2463);[3263,3289);[4907,4927);[7449,7455);"
    "[9087,9113);[10725,10751);[12363,12389)"
)
segment_hash = sha256(segment_string.encode()).hexdigest()
require(
    segment_hash
    == "4e458401005911de96cd61a26638cee0c2c75b8aa033df0c0298a3438c0514eb",
    "Boolean hostile segment hash changed",
)


def grid_danger(
    cell: int,
    denominator: int,
    speed: int,
    shift: int,
    length_units: int = 1,
) -> bool:
    """Test ||speed*x+shift/13||<length_units/14 at a cell midpoint."""
    modulus = 2 * denominator * P
    numerator = (
        speed * (2 * cell + 1) * P
        + shift * 2 * denominator
    )
    residue = numerator % modulus
    distance = min(residue, (-numerator) % modulus)
    return 14 * distance < length_units * modulus


all_bits = (1 << D) - 1
mask_bits = 0
for cell in range(D):
    if not any(left <= cell < right for left, right in EXCLUDED):
        mask_bits |= 1 << cell
require(mask_bits.bit_count() == D - 182, "wrong Boolean hostile mass")

deep_danger_bits: list[int] = []
deep_safe_bits: list[int] = []
graft_safe_bits: list[int] = []
for shift in range(P):
    deep = 0
    graft = 0
    for cell in range(D):
        if grid_danger(cell, D, 13, -shift):
            deep |= 1 << cell
        if grid_danger(cell, D, 7, shift):
            graft |= 1 << cell
    deep_danger_bits.append(deep)
    deep_safe_bits.append(all_bits ^ deep)
    graft_safe_bits.append(all_bits ^ graft)

masked_counts: list[int] = []
unmasked_toothpick: list[int] = []
expected_profile = (0, 1078) + (2002,) * 10 + (1078,)
for t_shift in range(P):
    row = []
    for difference in range(P):
        r_shift = (t_shift + difference) % P
        row.append(
            (
                mask_bits
                & deep_danger_bits[r_shift]
                & deep_safe_bits[t_shift]
                & graft_safe_bits[t_shift]
            ).bit_count()
        )
    require(tuple(row) == expected_profile, "masked hostile lost circulancy")
    masked_counts.extend(row)
    r_shift = (t_shift + 1) % P
    unmasked_toothpick.append(
        (
            deep_danger_bits[r_shift]
            & deep_safe_bits[t_shift]
            & graft_safe_bits[t_shift]
        ).bit_count()
    )

require(
    unmasked_toothpick
    == [1098, 1084, 1104, 1078, 1104, 1078, 1104, 1078, 1104, 1078, 1104, 1084, 1098],
    "unmasked escaping pair became circulant",
)
matrix_string = ",".join(str(value) for value in masked_counts)
matrix_hash = sha256(matrix_string.encode()).hexdigest()
require(
    matrix_hash
    == "e79dc643a13e9f4aafecd3f6b007952885902449c7aa8e5825f7f26bae4c7825",
    "Boolean hostile matrix hash changed",
)


# ---------------------------------------------------------------------------
# Full nine-factor bare-owner hostile on every one of the 165 profiles.
# ---------------------------------------------------------------------------

f0_cells = 0
for cell in range(420):
    x_value = Fraction(2 * cell + 1, 840)
    if (
        safe(x_value, 2)
        and safe(2 * x_value)
        and safe(3 * x_value)
        and safe(5 * x_value)
    ):
        f0_cells += 1
require(f0_cells == 182, "low-frequency owner cell count changed")
f0_mass = Fraction(f0_cells, 420)
require(f0_mass == Fraction(13, 30), "low-frequency mass changed")

q_two = 2940
c_one = 13 * 7 * q_two
profiles = [
    (1, middle, deepest)
    for deepest in range(5, 20)
    for middle in range(1, deepest)
]
require(len(profiles) == 165, "valuation ledger changed")
profile_checks = 0
for _, middle, deepest in profiles:
    c_two = c_one * 7 * 13 ** (middle - 1)
    c_three = c_two * 7 * 13 ** (deepest - middle)
    require(
        (
            valuation(c_one, 13),
            valuation(c_two, 13),
            valuation(c_three, 13),
        )
        == (1, middle, deepest),
        "full hostile thirteen-profile changed",
    )
    require(q_two // 420 == 7, "q peeling ratio changed")
    require(c_one // q_two == 7 * 13, "c1 peeling ratio changed")
    require(c_two // c_one % 7 == 0, "c2 peeling lost seven")
    require(c_three // c_two % 7 == 0, "c3 peeling lost seven")
    profile_checks += 1

hostile_prefactor = f0_mass * Fraction(6, 7) ** 3 * Fraction(1, 7)
require(hostile_prefactor == Fraction(468, 12005), "hostile A changed")
owner_mass = hostile_prefactor * Fraction(6, 7)
require(owner_mass == Fraction(2808, 84035), "owner mass changed")
uncovered_mass = f0_mass * Fraction(6, 7) ** 5
require(
    uncovered_mass == Fraction(16848, 84035),
    "uncovered hostile mass changed",
)

j_profile = tuple(
    Fraction(0)
    if difference == 0
    else Fraction(1, 13)
    if difference in (1, P - 1)
    else Fraction(1, 7)
    for difference in range(P)
)
require(
    tuple(hostile_prefactor * value for value in j_profile)
    == (
        Fraction(0),
        Fraction(36, 12005),
        *(Fraction(468, 84035) for _ in range(10)),
        Fraction(36, 12005),
    ),
    "full hostile overlap profile changed",
)


# ---------------------------------------------------------------------------
# Additive Phi_7 top-layer and only-c3 role constraints.
# ---------------------------------------------------------------------------

top_weight_controls = 0
for guard_top in (0, 1):
    for top_units in range(6):
        weight = 2 * guard_top + top_units
        if weight == 0:
            continue
        if weight % 7 == 0:
            require(
                guard_top == 1 and top_units == 5,
                "top-layer constant mask has a new weight",
            )
        top_weight_controls += 1
require(top_weight_controls == 11, "wrong top-weight control count")

c3_role_patterns: list[tuple[int, int, int]] = []
for low_blockers in range(3):
    for top_units in range(1, 6):
        for top_low_blockers in range(low_blockers + 1):
            weight = top_units + top_low_blockers
            cover_permitted = weight <= low_blockers or weight >= 7
            if cover_permitted:
                c3_role_patterns.append(
                    (low_blockers, top_units, top_low_blockers)
                )
require(
    c3_role_patterns
    == [
        (1, 1, 0),
        (2, 1, 0),
        (2, 1, 1),
        (2, 2, 0),
        (2, 5, 2),
    ],
    "only-c3 role alternatives changed",
)


# ---------------------------------------------------------------------------
# Exact N=49 hard-role chamber.
# ---------------------------------------------------------------------------

orbit_size = 49
x_zero = Fraction(1, 686)
hard_speeds = {
    "H": 1,
    "q1": 7,
    "q2": 148,
    "q3": 171,
    "q4": 172,
    "q5": 243,
    "c1": 195,
    "c2": 16562,
    "c3": 215306,
}


def orbit_mask(speed: int, length_units: int = 1) -> set[int]:
    return {
        index
        for index in range(orbit_size)
        if danger(
            speed * (x_zero + Fraction(index, orbit_size)),
            length_units,
        )
    }


hard_masks = {
    name: orbit_mask(speed, 2 if name == "H" else 1)
    for name, speed in hard_speeds.items()
}
require(
    hard_masks["H"] == set(range(7)) | set(range(42, 49)),
    "hard guard mask changed",
)
require(hard_masks["c1"] == set(range(11, 18)), "hard c1 mask changed")
require(not hard_masks["c2"] and not hard_masks["c3"], "deep masks changed")
require(hard_masks["q1"] == {0, 7, 14, 21, 28, 35, 42}, "top q mask")
require(hard_masks["q2"] == set(range(35, 42)), "q2 mask changed")
require(hard_masks["q3"] == {18, 20, 22, 24, 26, 28, 30}, "q3 mask")
require(hard_masks["q4"] == {19, 21, 23, 25, 27, 29, 31}, "q4 mask")
require(hard_masks["q5"] == {7, 8, 9, 10, 32, 33, 34}, "q5 mask")

base_tile_names = ("H", "c1", "q2", "q3", "q4", "q5")
base_counts = [
    sum(index in hard_masks[name] for name in base_tile_names)
    for index in range(orbit_size)
]
require(base_counts == [1] * orbit_size, "hard base masks do not tile")
full_counts = [
    base_counts[index] + int(index in hard_masks["q1"])
    for index in range(orbit_size)
]
require(full_counts.count(1) == 42 and full_counts.count(2) == 7, "hard cover")

exclusive_c1 = (
    hard_masks["c1"]
    - hard_masks["H"]
    - hard_masks["q1"]
    - hard_masks["q2"]
    - hard_masks["q3"]
    - hard_masks["q4"]
    - hard_masks["q5"]
    - hard_masks["c2"]
    - hard_masks["c3"]
)
require(
    exclusive_c1 == {11, 12, 13, 15, 16, 17},
    "hard exclusive-owner sites changed",
)

boundary_radius: Fraction | None = None
for name, speed in hard_speeds.items():
    length_units = 2 if name == "H" else 1
    boundary = Fraction(length_units, 14)
    for index in range(orbit_size):
        phase = speed * (x_zero + Fraction(index, orbit_size))
        phase_distance = min(
            circle_norm(phase - boundary),
            circle_norm(phase + boundary),
        )
        x_distance = phase_distance / speed
        if boundary_radius is None or x_distance < boundary_radius:
            boundary_radius = x_distance
require(boundary_radius is not None, "missing hard chamber radius")
require(
    boundary_radius >= Fraction(1, 3014284),
    "hard chamber radius fell below the certified value",
)

require(
    (
        valuation(hard_speeds["c1"], 13),
        valuation(hard_speeds["c2"], 13),
        valuation(hard_speeds["c3"], 13),
    )
    == (1, 2, 3),
    "hard role thirteen-profile changed",
)


# ---------------------------------------------------------------------------
# Signed lower-event current and cross-chamber capacity screen.
# ---------------------------------------------------------------------------

hard_lower_roles = (
    ("H", hard_speeds["H"], 2),
    ("q2", hard_speeds["q2"], 1),
    ("q3", hard_speeds["q3"], 1),
    ("q4", hard_speeds["q4"], 1),
    ("q5", hard_speeds["q5"], 1),
    ("c1", hard_speeds["c1"], 1),
)
hard_absorbers = (
    hard_speeds["q1"],
    hard_speeds["c2"],
    hard_speeds["c3"],
)


def ordered_handoff_mass(
    lower_roles: tuple[tuple[str, int, int], ...],
) -> int:
    total = 0
    for name_i, speed_i, width_i in lower_roles:
        for name_j, speed_j, width_j in lower_roles:
            if name_i == name_j:
                continue
            common = gcd(speed_i, speed_j)
            if (
                width_i * speed_j + width_j * speed_i
            ) % (14 * common) == 0:
                total += common
    return total


def absorber_entry_capacity(
    lower_roles: tuple[tuple[str, int, int], ...],
    absorbers: tuple[int, ...],
) -> int:
    total = 0
    for _, speed, _ in lower_roles:
        for absorber in absorbers:
            common = gcd(speed, absorber)
            reduced = speed // common
            total += common * (reduced // 7 + 1)
    return total


hard_event_mass = sum(speed for _, speed, _ in hard_lower_roles)
hard_handoff_mass = ordered_handoff_mass(hard_lower_roles)
hard_absorber_capacity = absorber_entry_capacity(
    hard_lower_roles,
    hard_absorbers,
)
require(hard_event_mass == 930, "hard lower event mass changed")
require(hard_handoff_mass == 0, "hard chamber gained a wall handoff")
require(hard_absorber_capacity == 432, "hard absorber capacity changed")
require(
    hard_event_mass
    > hard_handoff_mass + hard_absorber_capacity,
    "hard chamber no longer fails the global event screen",
)

shield_c_two = 98 * 13 ** 2 * 60
shield_c_three = 13 * shield_c_two
shield_lower_roles = (
    ("H", 1, 2),
    ("q2", 2, 1),
    ("q3", 3, 1),
    ("q4", 4, 1),
    ("q5", 5, 1),
    ("c1", 13, 1),
)
shield_absorbers = (7, shield_c_two, shield_c_three)
shield_event_mass = sum(speed for _, speed, _ in shield_lower_roles)
shield_handoff_mass = ordered_handoff_mass(shield_lower_roles)
shield_absorber_capacity = absorber_entry_capacity(
    shield_lower_roles,
    shield_absorbers,
)
require(shield_event_mass == 28, "shield event mass changed")
require(shield_handoff_mass == 0, "shield control gained a handoff")
require(shield_absorber_capacity == 63, "shield capacity changed")
require(
    shield_event_mass
    <= shield_handoff_mass + shield_absorber_capacity,
    "lawful shield control unexpectedly fails event capacity",
)
require(
    (
        valuation(13, 13),
        valuation(shield_c_two, 13),
        valuation(shield_c_three, 13),
    )
    == (1, 2, 3),
    "shield control thirteen-profile changed",
)
require(
    valuation(7, 7) == 1
    and valuation(13, 7) == 0
    and valuation(shield_c_two, 7) > 1
    and valuation(shield_c_three, 7) > 1,
    "shield control septimal roles changed",
)

shield_witness = Fraction(319, 2000)
require(
    circle_norm(shield_witness) > Fraction(1, 7),
    "shield noncover witness left the guard complement",
)
require(
    all(
        not danger(speed * shield_witness)
        for speed in (
            7,
            2,
            3,
            4,
            5,
            13,
            shield_c_two,
            shield_c_three,
        )
    ),
    "shield control became a cover at the exact witness",
)

print("theorem=THM-2367")
print("status=PROVED+VERIFIED-EXACT")
print(f"root_average_cases={root_average_cases}")
print(f"converse_window_controls={window_controls}")
print("graft_examples=(91,1):circulant;(91,7):drift;(13,1):drift")
print(f"boolean_hostile_cells={D}")
print(f"boolean_hostile_deleted={excluded_cells}")
print("boolean_hostile_mass=90/91")
print(f"boolean_segment_sha256={segment_hash}")
print(f"boolean_matrix_sha256={matrix_hash}")
print(f"low_frequency_cells_selected={f0_cells}/420")
print(f"full_hostile_profiles={profile_checks}")
print(f"full_hostile_prefactor={hostile_prefactor}")
print(f"full_hostile_owner_mass={owner_mass}")
print(f"full_hostile_uncovered_mass={uncovered_mass}")
print(f"top_weight_controls={top_weight_controls}")
print(f"c3_role_permitted_patterns={len(c3_role_patterns)}")
print(f"hard_role_exclusive_owner_sites={len(exclusive_c1)}")
print(f"hard_role_chamber_radius={boundary_radius}")
print(
    "hard_event_screen="
    f"{hard_event_mass}>{hard_handoff_mass + hard_absorber_capacity}"
)
print(
    "shield_event_screen="
    f"{shield_event_mass}<={shield_handoff_mass + shield_absorber_capacity}"
)
print(f"shield_noncover_witness={shield_witness}")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
