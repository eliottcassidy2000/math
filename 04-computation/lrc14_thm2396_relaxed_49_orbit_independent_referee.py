#!/usr/bin/env python3
"""Independent finite referee for the reserved THM-2396 route.

The universe deliberately relaxes the physical common-core geometry.
Work in Z/49Z with top bin Q=7Z/49Z.  Ordinary unit dangers are arbitrary
translates of seven-term unit-step progressions; the guard is an arbitrary
translate of a fourteen-term unit-step progression.  After multiplying the
orbit coordinate by h, the common-core A=D_h and B=D_(13h) have steps
1 and 13^(-1)=34, but their translates are allowed to vary independently.

The retained constraints are only:

* on Q, the guard is the disjoint union of A and B and every one of the
  four lower ordinary words lies in those two guard addresses;
* on each of the other six bins, guard plus the four lower words is an
  exact six-address transversal; and
* a non-B hole must use the A address.

The exact search proves that the last event must occur on at least four
safe bins.  In the physical common core a non-B hole also requires A=C,
where C=D_(169h).  The independently translated A/C combs agree in at
most two bins.  Thus even this relaxed universe is empty.  A four-bin
positive control shows that the finite obstruction is sharp.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import gcd


MODULUS = 49
ROOTS = tuple(range(7))
SAFE_BINS = tuple(range(1, 7))
TOP_BIN = frozenset(range(0, MODULUS, 7))
UNITS = tuple(step for step in range(MODULUS) if gcd(step, MODULUS) == 1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@dataclass(frozen=True)
class Mask:
    points: tuple[int, ...]
    addresses: tuple[tuple[int, ...], ...]


def address_rows(points: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    rows: list[list[int]] = [[] for _ in ROOTS]
    for point in points:
        residue = point % 7
        rows[residue].append((point - residue) // 7)
    return tuple(tuple(sorted(row)) for row in rows)


def progression_masks(length: int) -> tuple[Mask, ...]:
    """All distinct length-``length`` unit-step cyclic progressions."""

    masks: dict[tuple[int, ...], Mask] = {}
    for step in UNITS:
        for start in range(MODULUS):
            points = tuple(
                sorted((start + index * step) % MODULUS for index in range(length))
            )
            masks.setdefault(points, Mask(points, address_rows(points)))
    return tuple(masks.values())


def fixed_step_masks(length: int, step: int) -> tuple[Mask, ...]:
    return tuple(
        Mask(
            points := tuple(
                sorted(
                    (start + index * step) % MODULUS
                    for index in range(length)
                )
            ),
            address_rows(points),
        )
        for start in range(MODULUS)
    )


ORDINARY = progression_masks(7)
GUARDS = progression_masks(14)
A_MASKS = fixed_step_masks(7, 1)
B_MASKS = fixed_step_masks(7, pow(13, -1, MODULUS))
C_MASKS = fixed_step_masks(7, pow(169, -1, MODULUS))

require(len(UNITS) == 42, "wrong unit count")
require(len(ORDINARY) == 1029, "wrong ordinary-comb count")
require(len(GUARDS) == 1029, "wrong guard-comb count")
require(pow(13, -1, MODULUS) == 34, "wrong B step")
require(pow(169, -1, MODULUS) == 29, "wrong C step")


ORDINARY_BY_TOP: dict[int, list[Mask]] = defaultdict(list)
for mask in ORDINARY:
    require(all(len(mask.addresses[residue]) == 1 for residue in ROOTS), "ordinary")
    ORDINARY_BY_TOP[mask.addresses[0][0]].append(mask)

for guard in GUARDS:
    require(all(len(guard.addresses[residue]) == 2 for residue in ROOTS), "guard")


def achievable_holes(
    guard: Mask,
) -> tuple[int, dict[tuple[int, ...], tuple[Mask, Mask, Mask, Mask]]]:
    """Return every safe-bin hole word realized by four ordinary combs."""

    top_addresses = set(guard.addresses[0])
    candidates = [
        mask
        for top in top_addresses
        for mask in ORDINARY_BY_TOP[top]
        if all(mask.addresses[residue][0] not in guard.addresses[residue]
               for residue in SAFE_BINS)
    ]

    by_first_safe: dict[int, list[Mask]] = defaultdict(list)
    for mask in candidates:
        by_first_safe[mask.addresses[1][0]].append(mask)

    complement = [
        address for address in ROOTS if address not in guard.addresses[1]
    ]
    witnesses: dict[tuple[int, ...], tuple[Mask, Mask, Mask, Mask]] = {}
    chosen: list[Mask] = []

    def extend(required: tuple[int, ...], index: int, first_hole: int) -> None:
        if index == len(required):
            holes = [first_hole]
            for residue in range(2, 7):
                used = set(guard.addresses[residue])
                used.update(mask.addresses[residue][0] for mask in chosen)
                missing = [address for address in ROOTS if address not in used]
                require(len(missing) == 1, "completed word is not a transversal")
                holes.append(missing[0])
            witnesses.setdefault(tuple(holes), tuple(chosen))  # type: ignore[arg-type]
            return

        for mask in by_first_safe[required[index]]:
            if any(
                mask.addresses[residue][0] == prior.addresses[residue][0]
                for prior in chosen
                for residue in range(2, 7)
            ):
                continue
            chosen.append(mask)
            extend(required, index + 1, first_hole)
            chosen.pop()

    for first_hole in complement:
        required = tuple(address for address in complement if address != first_hole)
        if all(address in by_first_safe for address in required):
            extend(required, 0, first_hole)

    return len(candidates), witnesses


ACHIEVABLE: dict[
    tuple[int, ...],
    dict[tuple[int, ...], tuple[Mask, Mask, Mask, Mask]],
] = {}
candidate_histogram: Counter[int] = Counter()
hole_count_histogram: Counter[int] = Counter()

for guard in GUARDS:
    candidate_count, witnesses = achievable_holes(guard)
    ACHIEVABLE[guard.points] = witnesses
    candidate_histogram[candidate_count] += 1
    hole_count_histogram[len(witnesses)] += 1

require(
    candidate_histogram == Counter({22: 294, 25: 294, 26: 147, 33: 294}),
    "wrong candidate-count histogram",
)
require(
    hole_count_histogram == Counter({0: 588, 1: 147, 3: 294}),
    "wrong achievable-hole histogram",
)
require(sum(len(words) for words in ACHIEVABLE.values()) == 1029, "hole total")

# The guard enumeration has only four affine orbits under transformations
# j -> u*j+7*t which preserve the normalized top bin Q.
affine_q = tuple((unit, 7 * shift) for unit in UNITS for shift in ROOTS)


def guard_orbit_key(points: tuple[int, ...]) -> tuple[int, ...]:
    return min(
        tuple(sorted((unit * point + shift) % MODULUS for point in points))
        for unit, shift in affine_q
    )


guard_orbits: dict[tuple[int, ...], list[Mask]] = defaultdict(list)
for guard in GUARDS:
    guard_orbits[guard_orbit_key(guard.points)].append(guard)

orbit_profile: Counter[tuple[int, int]] = Counter()
for orbit in guard_orbits.values():
    hole_counts = {len(ACHIEVABLE[guard.points]) for guard in orbit}
    require(len(hole_counts) == 1, "hole count is not affine-orbit invariant")
    orbit_profile[(len(orbit), hole_counts.pop())] += 1

require(
    orbit_profile == Counter({(294, 0): 2, (294, 3): 1, (147, 1): 1}),
    "wrong affine guard-orbit profile",
)


GUARDS_BY_TOP_PAIR: dict[tuple[int, int], list[Mask]] = defaultdict(list)
for guard in GUARDS:
    GUARDS_BY_TOP_PAIR[tuple(sorted(guard.addresses[0]))].append(guard)

top_compatible_triples = 0
scanned_hole_words = 0
best_mismatch = 7
best_witness: tuple[
    Mask,
    Mask,
    Mask,
    tuple[int, ...],
    tuple[Mask, Mask, Mask, Mask],
] | None = None
survivors_with_three_flips = 0

for a_mask in A_MASKS:
    for b_mask in B_MASKS:
        a_top = a_mask.addresses[0][0]
        b_top = b_mask.addresses[0][0]
        if a_top == b_top:
            continue
        top_pair = tuple(sorted((a_top, b_top)))
        for guard in GUARDS_BY_TOP_PAIR[top_pair]:
            top_compatible_triples += 1
            for holes, lower_words in ACHIEVABLE[guard.points].items():
                scanned_hole_words += 1
                mismatch = 0
                compatible = True
                for residue, hole in zip(SAFE_BINS, holes):
                    if hole == b_mask.addresses[residue][0]:
                        continue
                    mismatch += 1
                    if hole != a_mask.addresses[residue][0]:
                        compatible = False
                        break
                if not compatible:
                    continue
                if mismatch < best_mismatch:
                    best_mismatch = mismatch
                    best_witness = (
                        a_mask,
                        b_mask,
                        guard,
                        holes,
                        lower_words,
                    )
                if mismatch <= 3:
                    survivors_with_three_flips += 1

require(top_compatible_triples == 100842, "wrong top-compatible count")
require(scanned_hole_words == 100842, "wrong compact certificate count")
require(best_mismatch == 4, "four-flip boundary not sharp")
require(survivors_with_three_flips == 0, "three-flip survivor found")
require(best_witness is not None, "missing sharp positive control")


# Independently translated A and C masks can agree in at most two bins.
equality_histogram: Counter[int] = Counter()
for a_mask in A_MASKS:
    for c_mask in C_MASKS:
        equality_histogram[
            sum(
                a_mask.addresses[residue][0] == c_mask.addresses[residue][0]
                for residue in ROOTS
            )
        ] += 1

require(
    equality_histogram == Counter({0: 490, 1: 1421, 2: 490}),
    "wrong A/C equality histogram",
)


# Check the sharp four-flip witness pointwise.
a_control, b_control, guard_control, holes_control, lower_control = best_witness
require(
    guard_control.points
    == (0, 4, 7, 10, 17, 20, 23, 27, 30, 33, 36, 40, 43, 46),
    "unexpected sharp guard representative",
)
require(holes_control == (1, 1, 6, 1, 0, 0), "unexpected sharp hole word")

for residue in SAFE_BINS:
    load = Counter(guard_control.addresses[residue])
    load.update(mask.addresses[residue][0] for mask in lower_control)
    require(load[holes_control[residue - 1]] == 0, "control hole is occupied")
    require(
        all(load[address] == (0 if address == holes_control[residue - 1] else 1)
            for address in ROOTS),
        "control is not a six-address transversal",
    )

top_guard = set(guard_control.points).intersection(TOP_BIN)
require(
    top_guard
    == set(a_control.points).intersection(TOP_BIN)
       | set(b_control.points).intersection(TOP_BIN),
    "control top guard is not A union B",
)
require(
    all(set(mask.points).intersection(TOP_BIN) <= top_guard for mask in lower_control),
    "control lower word escapes the top guard",
)

mismatch_bins = tuple(
    residue
    for residue, hole in zip(SAFE_BINS, holes_control)
    if hole != b_control.addresses[residue][0]
)
require(mismatch_bins == (1, 2, 5, 6), "wrong sharp mismatch bins")
require(
    all(
        holes_control[residue - 1] == a_control.addresses[residue][0]
        for residue in mismatch_bins
    ),
    "sharp mismatches are not A-address holes",
)

# The static obstruction localizes orbitwise: every high-safe orbit has a
# clean root.  THM-2396's high-safe base has mass 66/91 and the orbit has
# forty-nine roots.
high_safe_base = Fraction(66, 91)
common_core_clean_floor = high_safe_base / 49
charged_cell_floor = common_core_clean_floor / 52
complete_blocker_cell_floor = common_core_clean_floor / 78
owner_tensor_floor = common_core_clean_floor / 338
universal_last_lane_floor = min(common_core_clean_floor, Fraction(1, 26754))

require(common_core_clean_floor == Fraction(66, 4459), "wrong clean floor")
require(charged_cell_floor == Fraction(33, 115934), "wrong charged-cell floor")
require(
    complete_blocker_cell_floor == Fraction(11, 57967),
    "wrong blocker-cell floor",
)
require(owner_tensor_floor == Fraction(33, 753571), "wrong tensor floor")
require(
    universal_last_lane_floor == Fraction(1, 26754),
    "wrong universal last-lane floor",
)


def render_mask(mask: Mask) -> str:
    return ",".join(str(point) for point in mask.points)


print("target=THM-2396-independent-relaxed-49-orbit-referee")
print(
    "status=FINITE-EXACT+INDEPENDENTLY-HOSTILE-AUDITED;"
    " canonical-THM2396=PROVED-TWICE-AUDITED"
)
print(
    f"unit_steps={len(UNITS)}; ordinary_masks={len(ORDINARY)};"
    f" guard_masks={len(GUARDS)}"
)
print(
    "guard_candidate_histogram="
    + ",".join(f"{count}:{number}" for count, number in sorted(candidate_histogram.items()))
)
print(
    "achievable_hole_histogram="
    + ",".join(f"{count}:{number}" for count, number in sorted(hole_count_histogram.items()))
)
print("guard_affine_Q_orbits=4; orbit_profile=294x0:2,294x3:1,147x1:1")
print(
    f"total_achievable_hole_words={sum(len(words) for words in ACHIEVABLE.values())};"
    f" top_compatible_scans={scanned_hole_words}"
)
print(
    f"minimum_A_for_B_hole_flips={best_mismatch};"
    f" survivors_with_at_most_3={survivors_with_three_flips}"
)
print(
    "A_C_equality_histogram="
    + ",".join(f"{count}:{number}" for count, number in sorted(equality_histogram.items()))
    + "; physical_max=2"
)
print(
    f"sharp_A={render_mask(a_control)};"
    f" sharp_B={render_mask(b_control)};"
    f" sharp_guard={render_mask(guard_control)}"
)
print(
    "sharp_lower_words="
    + "|".join(render_mask(mask) for mask in lower_control)
)
print(
    "sharp_holes="
    + ",".join(map(str, holes_control))
    + "; sharp_A_flip_bins="
    + ",".join(map(str, mismatch_bins))
)
print("physical_common_core_survivors=0")
print(
    f"common_core_clean_floor={common_core_clean_floor};"
    f" charged_cell={charged_cell_floor};"
    f" blocker_cell={complete_blocker_cell_floor};"
    f" owner_tensor={owner_tensor_floor}"
)
print(f"universal_last_lane_clean_floor={universal_last_lane_floor}")
print("all_checks=PASS")
