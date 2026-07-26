#!/usr/bin/env python3
"""Exact companion for THM-2435.

The theorem is primarily a Haar-disintegration argument.  This
dependency-free companion checks every finite root incidence used by
that argument, the exact mass arithmetic, deterministic label/word
assignments, sharp abstract measure controls, and explicit punctured
stalk positive and hostile examples.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import gcd


L = 91
FULL = frozenset(range(L))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def progression(start: int, step: int, length: int) -> frozenset[int]:
    return frozenset((start + index * step) % L for index in range(length))


def circle_norm(value: Fraction) -> Fraction:
    residue = value - value.numerator // value.denominator
    return min(residue, 1 - residue)


def root_mask(
    speed: int,
    parent: Fraction,
    threshold: Fraction,
) -> frozenset[int]:
    return frozenset(
        root
        for root in range(L)
        if circle_norm(Fraction(speed) * (parent + root) / L) < threshold
    )


def row_histogram(mask: frozenset[int]) -> Counter[int]:
    return Counter(point % 13 for point in mask)


def weak_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, parts - 1):
            yield (first, *tail)


# Unit ordinary words have one point in every F_13 row; the unit guard
# has two.  A top blocker has speed 13*w with w a septimal unit, hence
# is one vertical coset: one point in every F_13 row and one fixed
# residue modulo seven.
unit_residues = tuple(speed for speed in range(1, L) if gcd(speed, L) == 1)
generic_parents = (
    Fraction(1, 257),
    Fraction(19, 257),
    Fraction(123, 257),
)

ordinary_controls = 0
guard_controls = 0
for speed in unit_residues:
    for parent in generic_parents:
        ordinary = root_mask(speed, parent, Fraction(1, 14))
        guard = root_mask(speed, parent, Fraction(1, 7))
        require(len(ordinary) == 13, "ordinary root capacity changed")
        require(
            row_histogram(ordinary) == Counter({row: 1 for row in range(13)}),
            "ordinary word is not a one-point section",
        )
        require(len(guard) == 26, "guard root capacity changed")
        require(
            row_histogram(guard) == Counter({row: 2 for row in range(13)}),
            "guard word is not a two-point section",
        )
        ordinary_controls += 1
        guard_controls += 1

top_blocker_controls = 0
for reduced_speed in range(1, 7):
    for parent in generic_parents:
        blocker = root_mask(13 * reduced_speed, parent, Fraction(1, 14))
        require(len(blocker) == 13, "top-blocker root capacity changed")
        require(
            row_histogram(blocker) == Counter({row: 1 for row in range(13)}),
            "top blocker is not a one-point section",
        )
        require(
            len({point % 7 for point in blocker}) == 1,
            "top blocker is not a vertical coset",
        )
        top_blocker_controls += 1


# Exact measure transfer for the three surviving (k,b) shapes.  The
# sharpened quotient-parent floor is (4+k)/7, while the common
# six-unit exact-tiling locus has mass at most 3/7.  Every parent in
# the difference contributes at least one essential top-blocker point
# among its 91 quotient branches.
typed_shapes = ((1, 1), (2, 2), (2, 1))
unit_tiling_cap = Fraction(3, 7)
parent_floors = tuple(Fraction(4 + k, 7) for k, _ in typed_shapes)
essential_parent_floors = tuple(
    floor - unit_tiling_cap for floor in parent_floors
)
aggregate_carrier_floors = tuple(
    floor / L for floor in essential_parent_floors
)
fixed_label_floors = tuple(
    floor / b
    for floor, (_, b) in zip(aggregate_carrier_floors, typed_shapes)
)
canonical_word_floors = tuple(
    floor / (2**b - 1)
    for floor, (_, b) in zip(aggregate_carrier_floors, typed_shapes)
)

parent_floor = min(parent_floors)
essential_parent_floor = min(essential_parent_floors)
aggregate_carrier_floor = min(aggregate_carrier_floors)
fixed_label_floor = min(fixed_label_floors)
canonical_word_floor = min(canonical_word_floors)

require(parent_floor == Fraction(5, 7), "uniform parent floor changed")
require(essential_parent_floor == Fraction(2, 7), "parent gap changed")
require(
    aggregate_carrier_floor == Fraction(2, 637),
    "aggregate carrier floor changed",
)
require(
    fixed_label_floor == Fraction(3, 1274),
    "fixed-label carrier floor changed",
)
require(
    canonical_word_floor == Fraction(1, 637),
    "canonical-word carrier floor changed",
)


# Every row has unit incidence seven.  Its signed defect
# d(r)=1-m(r) is an integral zero-sum vector.  Exhausting all 1,716
# weak compositions checks the sharp nonzero cyclic Dirichlet floor.
row_multiplicity_vectors = tuple(weak_compositions(7, 7))
require(
    len(row_multiplicity_vectors) == 1716,
    "wrong seven-incidence composition count",
)
nonzero_row_energies = []
for multiplicities in row_multiplicity_vectors:
    defect = tuple(1 - value for value in multiplicities)
    require(sum(defect) == 0, "row defect stopped having zero sum")
    if any(defect):
        energy = sum(
            (defect[(residue + 1) % 7] - defect[residue]) ** 2
            for residue in range(7)
        )
        nonzero_row_energies.append(energy)
        # Over Q, vanishing of one primitive seventh-root Fourier
        # value would force all seven coefficients equal.  A zero-sum
        # equal vector is zero, contrary to this branch.
        require(
            len(set(defect)) > 1,
            "nonzero zero-sum defect became cyclotomic-constant",
        )
require(min(nonzero_row_energies) == 4, "sharp row energy floor changed")
parent_dirichlet_floor = Fraction(4) * essential_parent_floor
require(
    parent_dirichlet_floor == Fraction(8, 7),
    "integrated parent Dirichlet floor changed",
)


# With one essential blocker source, coverage off that source leaves
# exactly six directed source-to-duplicate arrows.  Adjacent arrows
# have energy six and the other four targets have energy four.
b1_arrow_energy_histogram: Counter[int] = Counter()
for source in range(7):
    admissible = []
    for multiplicities in row_multiplicity_vectors:
        if multiplicities[source] != 0:
            continue
        if any(
            multiplicities[residue] < 1
            for residue in range(7)
            if residue != source
        ):
            continue
        targets = [
            residue
            for residue, value in enumerate(multiplicities)
            if value == 2
        ]
        require(
            len(targets) == 1
            and all(value <= 2 for value in multiplicities),
            "one-source essential row is not one arrow",
        )
        defect = tuple(1 - value for value in multiplicities)
        energy = sum(
            (defect[(residue + 1) % 7] - defect[residue]) ** 2
            for residue in range(7)
        )
        b1_arrow_energy_histogram[energy] += 1
        admissible.append(targets[0])
    require(
        sorted(admissible) == [
            residue for residue in range(7) if residue != source
        ],
        "one-source arrow target bank changed",
    )
require(
    b1_arrow_energy_histogram == Counter({4: 28, 6: 14}),
    "one-source arrow energy histogram changed",
)


# Exact branch-Jacobian control for
#
#   integral sum_s 1_E((Y+s)/91) dY = 91*mu(E).
#
# A parent-length alpha_s in branch s becomes a physical interval of
# length alpha_s/91.  Nonuniform rational lengths ensure this is not
# merely a constant-function check.
parent_branch_lengths = tuple(
    Fraction((17 * branch + 11) % 101, 101)
    for branch in range(L)
)
integrated_fibre_count = sum(parent_branch_lengths, Fraction(0))
physical_mass = sum(
    (length / L for length in parent_branch_lengths),
    Fraction(0),
)
require(
    integrated_fibre_count == L * physical_mass,
    "91-branch disintegration identity failed",
)


# The typed constants are sharp from the two measure inputs alone.
# Equal abstract pieces attain each deterministic b-label and
# intrinsic (2^b-1)-word pigeonhole factor.
sharp_parent_differences = tuple(
    floor - unit_tiling_cap for floor in parent_floors
)
require(
    sharp_parent_differences == essential_parent_floors,
    "abstract parent-difference hostile changed",
)
assigned_label_masses = tuple(
    tuple(floor / b for _ in range(b))
    for floor, (_, b) in zip(aggregate_carrier_floors, typed_shapes)
)
require(
    tuple(max(masses) for masses in assigned_label_masses)
    == fixed_label_floors,
    "typed label pigeonhole hostile changed",
)
canonical_word_masses = tuple(
    tuple(floor / (2**b - 1) for _ in range(2**b - 1))
    for floor, (_, b) in zip(aggregate_carrier_floors, typed_shapes)
)
require(
    tuple(max(masses) for masses in canonical_word_masses)
    == canonical_word_floors,
    "typed word pigeonhole hostile changed",
)


# Finite punctured-stalk positive controls.  The unpunctured step-one
# tiling has guard 0..25 and ordinary intervals starting at
# 26,39,52,65,78.  Moving the last interval by one creates one hole
# and one duplicate in the same F_13 row; the vertical blocker through
# that hole repairs the cover.
guard = progression(0, 1, 26)
b1_ordinary = tuple(
    progression(start, 1, 13)
    for start in (26, 39, 52, 65, 79)
)
b1_unit_multiset = tuple(guard) + tuple(
    point for mask in b1_ordinary for point in mask
)
b1_unit_union = guard.union(*b1_ordinary)
b1_missing = FULL - b1_unit_union
b1_blocker = frozenset(point for point in range(L) if point % 7 == 1)
require(len(b1_unit_multiset) == L, "b=1 unit incidence changed")
require(b1_missing == {78}, "b=1 puncture changed")
require(b1_unit_union | b1_blocker == FULL, "b=1 blocker stopped repairing")


def row_multiplicity(
    masks: tuple[frozenset[int], ...],
    row: int,
) -> Counter[int]:
    counts: Counter[int] = Counter()
    for mask in masks:
        for point in mask:
            if point % 13 == row:
                counts[point % 7] += 1
    return counts


b1_flows = {}
for row in range(13):
    counts = row_multiplicity((guard, *b1_ordinary), row)
    missing = [residue for residue in range(7) if counts[residue] == 0]
    doubled = [residue for residue in range(7) if counts[residue] == 2]
    require(
        all(value <= 2 for value in counts.values()),
        "b=1 row has multiplicity above two",
    )
    if missing:
        require(
            missing == [1] and len(doubled) == 1,
            "b=1 row is not a source-to-double arrow",
        )
        b1_flows[row] = (1, doubled[0])
    else:
        require(
            counts == Counter({residue: 1 for residue in range(7)}),
            "b=1 flat row is not a permutation",
        )
require(b1_flows == {0: (1, 0)}, "b=1 arrow word changed")


# Two distinct blocker sources repair two holes in one row.  Repeating
# the b=1 blocker instead is also a valid cover and demonstrates that
# a coincident second source does not force a second hole.
b2_ordinary = tuple(
    progression(start, 1, 13)
    for start in (26, 40, 52, 65, 79)
)
b2_unit_union = guard.union(*b2_ordinary)
b2_missing = FULL - b2_unit_union
b2_blocker_1 = frozenset(point for point in range(L) if point % 7 == 1)
b2_blocker_4 = frozenset(point for point in range(L) if point % 7 == 4)
require(b2_missing == {39, 78}, "b=2 punctures changed")
require(
    b2_unit_union | b2_blocker_1 | b2_blocker_4 == FULL,
    "b=2 blockers stopped repairing",
)
require(
    b1_unit_union | b1_blocker | b1_blocker == FULL,
    "coincident b=2 source control changed",
)


# Scope hostile: THM-2427's strict local packet is an exact six-unit
# tiling at its displayed parent, so it has zero essential points.
# The carrier is guaranteed only after integrating over P\\Z.
hostile_parent = Fraction(11, 83)
hostile_masks = (
    root_mask(1, hostile_parent, Fraction(1, 7)),
    *(
        root_mask(speed, hostile_parent, Fraction(1, 14))
        for speed in (547, 1821, 3095, 4369, 5643)
    ),
)
require(
    sum(map(len, hostile_masks)) == L
    and set().union(*hostile_masks) == set(FULL),
    "local zero-essential hostile stopped tiling",
)


print("THM-2435 exact companion")
print(f"ordinary_row_controls={ordinary_controls}")
print(f"guard_row_controls={guard_controls}")
print(f"top_blocker_coset_controls={top_blocker_controls}")
print("typed_k_b=1,1;2,2;2,1")
print(
    "typed_essential_parent_floors="
    + ",".join(map(str, essential_parent_floors))
)
print(f"essential_parent_floor={essential_parent_floor}")
print(f"aggregate_carrier_floor={aggregate_carrier_floor}")
print(f"fixed_label_floor={fixed_label_floor}")
print(f"canonical_word_floor={canonical_word_floor}")
print(f"nonzero_row_dirichlet_floor={min(nonzero_row_energies)}")
print(f"integrated_parent_dirichlet_floor={parent_dirichlet_floor}")
print("all_six_nontrivial_F7_colors=NONZERO")
print("b1_arrow_energy_histogram=4:28,6:14")
print("disintegration_91_branch_identity=PASS")
print("b1_flow_word=h0:1->0")
print("b2_distinct_sources=holes_1,4_repaired")
print("b2_coincident_source_control=PASS")
print("local_zero_essential_hostile=PASS")
print("ALL CHECKS PASSED")
