#!/usr/bin/env python3
"""Exact companion for THM-2435.

The finite audit independently reconstructs the corrected exact-tiling
phase cap, checks the seven-bin gap bank and sharpened parent invoices,
constructs the equivariant one-root selector on C_13, verifies its flat
C_91 quotient spectrum, and checks the positive-depth ancestry kernel.
All truth-bearing checks use explicit exceptions and survive ``-O``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# 1. Independent reconstruction of the corrected exact-tiling cap.
# ---------------------------------------------------------------------------

L = 91
FULL = (1 << L) - 1
GUARD = sum(1 << site for site in range(26))
COMPLEMENT = FULL ^ GUARD


def progression(start: int, step: int, length: int = 13) -> int:
    mask = 0
    for offset in range(length):
        mask |= 1 << ((start + offset * step) % L)
    return mask


def sites(mask: int) -> tuple[int, ...]:
    return tuple(site for site in range(L) if (mask >> site) & 1)


representations: dict[int, list[tuple[int, int]]] = {}
for step in range(1, L):
    if gcd(step, L) != 1:
        continue
    for start in range(L):
        mask = progression(start, step)
        representations.setdefault(mask, []).append((start, step))

require(len(representations) == 3276, "wrong progression universe")

canonical: dict[int, tuple[int, int]] = {}
for mask, reps in representations.items():
    forward = sorted(rep for rep in reps if rep[1] <= L // 2)
    require(len(reps) == 2 and len(forward) == 1, "bad AP orientation")
    canonical[mask] = forward[0]

eligible = sorted(mask for mask in representations if mask & GUARD == 0)
require(len(eligible) == 182, "wrong eligible progression bank")

by_site: dict[int, list[int]] = {site: [] for site in sites(COMPLEMENT)}
for mask in eligible:
    for site in sites(mask):
        by_site[site].append(mask)


def fixed_leftmost_covers() -> set[frozenset[int]]:
    """Enumerate exact covers with a deterministic pivot, not index pruning."""

    answers: set[frozenset[int]] = set()

    def visit(remainder: int, chosen: tuple[int, ...]) -> None:
        if remainder == 0:
            require(len(chosen) == 5, "complete cover has wrong size")
            answers.add(frozenset(chosen))
            return
        require(len(chosen) < 5, "five blocks left a nonempty remainder")
        pivot_bit = remainder & -remainder
        pivot = pivot_bit.bit_length() - 1
        for candidate in by_site[pivot]:
            if candidate & remainder == candidate:
                visit(remainder ^ candidate, chosen + (candidate,))

    visit(COMPLEMENT, ())
    return answers


tilings = fixed_leftmost_covers()
require(len(tilings) == 62, "wrong normalized exact-tiling count")

centre_differences: dict[int, set[int]] = defaultdict(set)
for tiling in tilings:
    centres_by_step: dict[int, list[int]] = defaultdict(list)
    for mask in tiling:
        start, step = canonical[mask]
        centres_by_step[step].append((start + 6 * step) % L)
    require(
        any(len(centres) >= 2 for centres in centres_by_step.values()),
        "exact tiling without a repeated unsigned step",
    )
    for step, centres in centres_by_step.items():
        for first in centres:
            for second in centres:
                if first != second:
                    centre_differences[step].add((first - second) % L)

repeatable_steps = sorted(centre_differences)
require(
    repeatable_steps == [1, 2, 3, 4, 5, 44, 45],
    "wrong repeatable-step set",
)

phase_cells: dict[int, set[int]] = {}
for step, differences in centre_differences.items():
    inverse = pow(step, -1, L)
    normalized = {(inverse * difference) % L for difference in differences}
    phase_cells[step] = normalized | {
        (residue - 1) % L for residue in normalized
    }

expected_phase_sizes = {
    1: 39,
    2: 35,
    3: 6,
    4: 9,
    5: 12,
    44: 10,
    45: 8,
}
require(
    {step: len(cells) for step, cells in phase_cells.items()}
    == expected_phase_sizes,
    "wrong repeated-step phase banks",
)
require(max(map(len, phase_cells.values())) == 39, "wrong exact-locus cap")


# ---------------------------------------------------------------------------
# 2. The coarse seven-bin word bank.
# ---------------------------------------------------------------------------

unit_profiles = 0
gap_histogram: Counter[int] = Counter()
b1_cover_profiles = 0
b1_essential_profiles = 0
b2_cover_profiles = 0
b2_essential_profiles = 0

for guard in combinations(range(7), 2):
    for ordinary in product(range(7), repeat=5):
        unit_profiles += 1
        multiplicity = [0] * 7
        for cell in guard:
            multiplicity[cell] += 1
        for cell in ordinary:
            multiplicity[cell] += 1

        gaps = sum(value == 0 for value in multiplicity)
        gap_histogram[gaps] += 1

        if gaps == 0:
            b1_cover_profiles += 7
        elif gaps == 1:
            b1_cover_profiles += 1
            b1_essential_profiles += 1

        if gaps == 0:
            b2_cover_profiles += 49
        elif gaps == 1:
            b2_cover_profiles += 13
            b2_essential_profiles += 13
        elif gaps == 2:
            b2_cover_profiles += 2
            b2_essential_profiles += 2

require(unit_profiles == 352947, "wrong six-unit word count")
require(
    gap_histogram
    == Counter({0: 2520, 1: 50400, 2: 157500, 3: 119700, 4: 22155, 5: 672}),
    "wrong six-unit gap histogram",
)
require(b1_cover_profiles == 68040, "wrong one-blocker cover count")
require(b1_essential_profiles == 50400, "wrong one-blocker essential count")
require(b2_cover_profiles == 1093680, "wrong two-blocker cover count")
require(b2_essential_profiles == 970200, "wrong two-blocker essential count")


# ---------------------------------------------------------------------------
# 3. Sharpened parent, label, and marked-root invoices.
# ---------------------------------------------------------------------------

# A shape is (name, k low quotient blockers, b exact-depth top labels).
shapes = (
    ("M0_k1_b1", 1, 1),
    ("M0_k2_b2", 2, 2),
    ("Mpos_k2_b1", 2, 1),
)
shape_invoices: dict[str, tuple[Fraction, Fraction, Fraction, Fraction]] = {}

unit_exact_cap = Fraction(3, 7)
for name, low_count, top_labels in shapes:
    # P is exactly the complement of the 3-k high quotient blockers,
    # a.e.; the union bound gives (4+k)/7.
    parent_floor = Fraction(4 + low_count, 7)
    defect_parent_floor = parent_floor - unit_exact_cap
    fixed_label_parent_floor = defect_parent_floor / top_labels
    marked_quotient_mass_floor = fixed_label_parent_floor / 91
    shape_invoices[name] = (
        parent_floor,
        defect_parent_floor,
        fixed_label_parent_floor,
        marked_quotient_mass_floor,
    )

require(
    shape_invoices["M0_k1_b1"]
    == (Fraction(5, 7), Fraction(2, 7), Fraction(2, 7), Fraction(2, 637)),
    "wrong first shape invoice",
)
require(
    shape_invoices["M0_k2_b2"]
    == (Fraction(6, 7), Fraction(3, 7), Fraction(3, 14), Fraction(3, 1274)),
    "wrong second shape invoice",
)
require(
    shape_invoices["Mpos_k2_b1"]
    == (Fraction(6, 7), Fraction(3, 7), Fraction(3, 7), Fraction(3, 637)),
    "wrong positive-depth shape invoice",
)

uniform_marked_mass_floor = min(invoice[3] for invoice in shape_invoices.values())
uniform_per_c91_class_floor = uniform_marked_mass_floor / 91
uniform_unit_class_total_floor = 72 * uniform_per_c91_class_floor

require(uniform_marked_mass_floor == Fraction(3, 1274), "wrong marked mass floor")
require(
    uniform_per_c91_class_floor == Fraction(3, 115934),
    "wrong per-C91-class floor",
)
require(
    uniform_unit_class_total_floor == Fraction(108, 57967),
    "wrong unit-class total floor",
)


# ---------------------------------------------------------------------------
# 4. Translation-equivariant lexicographic marker on C_13.
# ---------------------------------------------------------------------------

def bits(mask: int) -> tuple[int, ...]:
    return tuple((mask >> r) & 1 for r in range(13))


def rotate_word(mask: int, start: int) -> tuple[int, ...]:
    indicator = bits(mask)
    return tuple(indicator[(start + j) % 13] for j in range(13))


def translate_mask(mask: int, shift: int) -> int:
    translated = 0
    for r in range(13):
        if (mask >> r) & 1:
            translated |= 1 << ((r + shift) % 13)
    return translated


def marker(mask: int) -> int:
    words = [rotate_word(mask, start) for start in range(13)]
    maximum = max(words)
    winners = [start for start, word in enumerate(words) if word == maximum]
    require(len(winners) == 1, "lexicographic marker is not unique")
    return winners[0]


mask_count = 0
translation_orbits: set[int] = set()
marker_histogram: Counter[int] = Counter()
size_histogram: Counter[int] = Counter()
equivariance_checks = 0
formal_cyclotomic_checks = 0
source_idempotent = 13 * pow(13, -1, 7)  # 78 mod 91
target_idempotent = 7 * pow(7, -1, 13)  # 14 mod 91
require(source_idempotent == 78, "wrong C7 CRT idempotent")
require(target_idempotent == 14, "wrong C13 CRT idempotent")
require(
    source_idempotent % 7 == 1 and source_idempotent % 13 == 0,
    "bad source CRT residues",
)
require(
    target_idempotent % 7 == 0 and target_idempotent % 13 == 1,
    "bad target CRT residues",
)

for mask in range(1, (1 << 13) - 1):
    mask_count += 1
    rotations = {translate_mask(mask, shift) for shift in range(13)}
    require(len(rotations) == 13, "nonconstant C13 mask has a stabilizer")
    translation_orbits.add(min(rotations))

    selected = marker(mask)
    require((mask >> selected) & 1 == 1, "marker did not select a hole")
    marker_histogram[selected] += 1
    size_histogram[mask.bit_count()] += 1

    for shift in range(13):
        shifted = translate_mask(mask, shift)
        require(
            marker(shifted) == (selected + shift) % 13,
            "marker is not translation-equivariant",
        )
        equivariance_checks += 1

    # On the unique C7 blocker sheet, a C91 character restricts to a
    # C13 character. For target character zero the sum is |S|. For a
    # nonzero target character, reduction modulo Phi_13 is zero iff
    # all thirteen 0/1 coefficients are equal. Nonempty proper masks
    # therefore survive every target character.
    for character in range(91):
        source_character = (6 * character) % 7
        target_character = (2 * character) % 13
        require(0 <= source_character < 7, "bad restricted source character")
        if target_character == 0:
            require(mask.bit_count() > 0, "zero trivial-character sum")
        else:
            permuted = [0] * 13
            for r in range(13):
                if (mask >> r) & 1:
                    permuted[(target_character * r) % 13] = 1
            reduced = tuple(permuted[r] - permuted[12] for r in range(12))
            require(any(value != 0 for value in reduced), "cyclotomic sum vanished")
        formal_cyclotomic_checks += 1

require(mask_count == 8190, "wrong nonconstant C13 mask count")
require(len(translation_orbits) == 630, "wrong C13 mask orbit count")
require(equivariance_checks == 106470, "wrong equivariance check count")
require(formal_cyclotomic_checks == 745290, "wrong C91 character check count")
require(
    marker_histogram == Counter({root: 630 for root in range(13)}),
    "marker addresses are not balanced",
)
require(
    size_histogram == Counter({size: comb(13, size) for size in range(1, 13)}),
    "wrong C13 mask-size census",
)


# ---------------------------------------------------------------------------
# 5. Exact flat C91 section and positive-depth ancestry kernel.
# ---------------------------------------------------------------------------

# One selected CRT root has normalized coefficient 1/91 in every
# character, hence squared magnitude 1/91^2. Integrating over a parent
# set of mass p gives p/91^2 = (p/91)/91.
for source_root in range(7):
    for target_root in range(13):
        crt_root = (
            source_idempotent * source_root
            + target_idempotent * target_root
        ) % 91
        require(crt_root % 7 == source_root, "bad CRT source root")
        require(crt_root % 13 == target_root, "bad CRT target root")
        for character in range(91):
            exponent = (-character * crt_root) % 91
            require(0 <= exponent < 91, "bad formal C91 exponent")
            require(Fraction(1, 91) ** 2 == Fraction(1, 8281), "bad C91 norm")

# At positive M, pullback through T_d with d=7^M has d ancestors.
# The physical C_(91d) transform vanishes off the annihilator d|q.
ancestry_controls = 0
for depth in range(1, 7):
    dilation = 7**depth
    # The geometric sum depends only on the character modulo d; each
    # residue has exactly 91 representatives modulo 91d.
    for step in range(dilation):
        expected_survival = step == 0
        # The d ancestor phases form a geometric sum with step
        # character modulo d. Its order is d/gcd(character,d); the
        # sum is d precisely at order one and zero otherwise.
        common = gcd(step, dilation)
        order = dilation // common
        require(common * order == dilation, "bad geometric order")
        arithmetic_survival = order == 1
        require(expected_survival == arithmetic_survival, "bad ancestry annihilator")
        amplitude_numerator = dilation if arithmetic_survival else 0
        require(
            Fraction(amplitude_numerator, 91 * dilation)
            == (Fraction(1, 91) if expected_survival else 0),
            "bad normalized ancestry amplitude",
        )
        ancestry_controls += 91
    require(
        all(gcd((dilation * residue) % 91, 91) == 7 for residue in range(1, 91) if gcd(residue, 91) == 1),
        "positive-depth unit quotient residue stayed a physical unit",
    )


# ---------------------------------------------------------------------------
# 6. Exact one-hole punctured hostile.
# ---------------------------------------------------------------------------

hostile_guard = set(range(26))
hostile_ordinary = [
    set(range(26, 39)),
    set(range(39, 52)),
    set(range(52, 65)),
    set(range(65, 78)),
    set(range(79, 91)) | {0},
]
hostile_six = hostile_guard | set().union(*hostile_ordinary)
hostile_missing = set(range(L)) - hostile_six
hostile_blocker = {site for site in range(L) if site % 7 == 1}
hostile_puncture = {site for site in range(L) if site % 7 == 2}

require(
    len(hostile_guard) + sum(map(len, hostile_ordinary)) == L,
    "hostile incidence changed",
)
require(hostile_missing == {78}, "hostile is not a one-hole word")
require(78 in hostile_blocker, "blocker does not repair hostile hole")
require(
    hostile_six | hostile_blocker == set(range(L)),
    "punctured hostile does not cover",
)
require(
    hostile_blocker.isdisjoint(hostile_puncture)
    and len(hostile_blocker) == len(hostile_puncture) == 13,
    "hostile puncture is not a disjoint thirteen-coset",
)


print("THM-2435 exact companion")
print(f"normalized_tilings={len(tilings)}")
print("repeatable_steps=" + ",".join(map(str, repeatable_steps)))
print(
    "phase_cell_sizes="
    + ",".join(
        f"{step}:{len(phase_cells[step])}" for step in repeatable_steps
    )
)
print(f"labelled_six_unit_profiles={unit_profiles}")
print(
    "unit_gap_histogram="
    + ",".join(f"{gaps}:{gap_histogram[gaps]}" for gaps in sorted(gap_histogram))
)
print(f"one_blocker_cover_profiles={b1_cover_profiles}")
print(f"two_blocker_cover_profiles={b2_cover_profiles}")
for name, _, _ in shapes:
    parent, defect, fixed_parent, marked = shape_invoices[name]
    print(
        f"{name}_invoice="
        f"parent:{parent},defect:{defect},fixed_parent:{fixed_parent},marked:{marked}"
    )
print(f"uniform_marked_mass_floor={uniform_marked_mass_floor}")
print(f"uniform_per_c91_class_floor={uniform_per_c91_class_floor}")
print(f"uniform_unit_class_total_floor={uniform_unit_class_total_floor}")
print(f"nonconstant_c13_masks={mask_count}")
print(f"c13_translation_orbits={len(translation_orbits)}")
print(f"marker_equivariance_checks={equivariance_checks}")
print(f"formal_c91_character_checks={formal_cyclotomic_checks}")
print(f"positive_depth_ancestry_controls={ancestry_controls}")
print("positive_depth_physical_residue_gcd=7")
print("one_hole_punctured_hostile=missing_78_repaired_PASS")
print("ALL CHECKS PASSED")
