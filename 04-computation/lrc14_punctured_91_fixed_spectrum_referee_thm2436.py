#!/usr/bin/env python3
"""Independent fixed-spectrum referee for THM-2436.

This companion does not import the production C++ atlas or any stored
solution.  It specializes the exceptional M=0 spectra to their native
geometry:

* a step-one AP13 is a cyclic interval of radius six in Z/91Z;
* every non-step-one AP13 is indexed independently by its centre; and
* one or two top blockers delete one or two residue columns modulo seven.

For a fixed non-step-one support, a dynamic missing-point recursion finds
all sets of step-one intervals completing the cover.  At each pivot it
branches on the least selected interval containing that pivot and removes
only lower-index competitors through the same pivot.  This canonical-owner
rule is compatible with a changing pivot; no increasing-index restriction
is imposed on unrelated later intervals.

The script independently reconstructs the complete one-source exceptional
table and the five critical two-label spectra with at least four step-one
labels.  In particular it replays the load-bearing maximum

    (1,1,1,1,45): 62 directed centre differences, 67 sharp cells.

The two remaining production-atlas spectra have only three and two
step-one labels and sharp banks 16 and 8; they are not needed to identify
or audit the 67-cell extremum here.
"""

from __future__ import annotations

from itertools import combinations, product


L = 91
FULL = (1 << L) - 1
GUARD = (1 << 26) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ap13(center: int, step: int) -> int:
    """The unoriented length-thirteen AP with displayed centre."""

    return sum(
        1 << ((center + offset * step) % L)
        for offset in range(-6, 7)
    )


STEP_ONE = tuple(ap13(center, 1) for center in range(L))
require(len(set(STEP_ONE)) == L, "step-one centre parametrization collided")
require(all(mask.bit_count() == 13 for mask in STEP_ONE), "wrong AP size")

THROUGH = [0] * L
for center, mask in enumerate(STEP_ONE):
    for point in range(L):
        if (mask >> point) & 1:
            THROUGH[point] |= 1 << center
require(
    all(owners.bit_count() == 13 for owners in THROUGH),
    "a point is not in exactly thirteen step-one intervals",
)


def blocker_mask(sources: tuple[int, int]) -> int:
    source_set = set(sources)
    return sum(
        1 << point
        for point in range(L)
        if point % 7 in source_set
    )


def interval_completions(
    required: int,
    base_covered: int,
    interval_count: int,
) -> tuple[set[tuple[int, ...]], int]:
    """All step-one centre sets completing one fixed partial cover."""

    solutions: set[tuple[int, ...]] = set()
    nodes = 0

    def record(selected: tuple[int, ...]) -> None:
        solutions.add(tuple(sorted(selected)))

    def visit(
        covered: int,
        allowed: int,
        depth: int,
        selected: tuple[int, ...],
    ) -> None:
        nonlocal nodes
        nodes += 1
        missing = required & ~covered
        remaining = interval_count - depth

        if missing.bit_count() > 13 * remaining:
            return
        if remaining == 0:
            if missing == 0:
                record(selected)
            return

        if missing == 0:
            available = [
                center for center in range(L)
                if (allowed >> center) & 1
            ]
            for extra in combinations(available, remaining):
                record(selected + extra)
            return

        pivot = -1
        owners = 0
        best_branching = L + 1
        pending = missing
        while pending:
            point_bit = pending & -pending
            point = point_bit.bit_length() - 1
            pending -= point_bit
            point_owners = allowed & THROUGH[point]
            branching = point_owners.bit_count()
            if branching == 0:
                return
            if branching < best_branching:
                pivot = point
                owners = point_owners
                best_branching = branching

        while owners:
            owner_bit = owners & -owners
            center = owner_bit.bit_length() - 1
            owners -= owner_bit

            # The chosen centre is the least selected owner of this pivot.
            # Lower owners through this pivot may not appear later; unrelated
            # lower centres remain available.
            lower_centers = owner_bit - 1
            next_allowed = (
                allowed
                & ~owner_bit
                & ~(THROUGH[pivot] & lower_centers)
            )
            visit(
                covered | STEP_ONE[center],
                next_allowed,
                depth + 1,
                selected + (center,),
            )

    visit(base_covered, (1 << L) - 1, 0, ())
    return solutions, nodes


def fixed_spectrum_bank(
    source_configurations: tuple[tuple[int, int], ...],
    spectrum: tuple[int, ...],
) -> tuple[int, int, int, int]:
    """Return assignments, directed differences, sharp cells, search nodes."""

    step_one_count = spectrum.count(1)
    other_steps = tuple(step for step in spectrum if step != 1)
    require(
        len(set(other_steps)) == len(other_steps),
        "referee handles only singleton non-step-one classes",
    )

    other_masks = {
        step: tuple(ap13(center, step) for center in range(L))
        for step in other_steps
    }

    assignments = 0
    directed_differences: set[int] = set()
    nodes = 0

    for sources in source_configurations:
        required = FULL & ~(GUARD | blocker_mask(sources))

        for other_centers in product(range(L), repeat=len(other_steps)):
            base_covered = 0
            for step, center in zip(other_steps, other_centers):
                base_covered |= other_masks[step][center]

            if (required & ~base_covered).bit_count() > 13 * step_one_count:
                continue

            completions, used_nodes = interval_completions(
                required,
                base_covered,
                step_one_count,
            )
            nodes += used_nodes
            assignments += len(completions)

            for centers in completions:
                for first in centers:
                    for second in centers:
                        if first != second:
                            directed_differences.add((first - second) % L)

    sharp_cells = directed_differences | {
        (difference - 1) % L
        for difference in directed_differences
    }
    return (
        assignments,
        len(directed_differences),
        len(sharp_cells),
        nodes,
    )


ONE_SOURCE_CONFIGURATIONS = tuple((source, source) for source in range(7))
TWO_LABEL_CONFIGURATIONS = tuple(
    (first, second)
    for first in range(7)
    for second in range(first, 7)
)

ONE_SOURCE_EXPECTED = {
    (1, 1, 1, 1, 1): (144, 40, 46),
    (1, 1, 1, 1, 30): (106, 36, 41),
    (1, 1, 1, 1, 45): (262, 48, 57),
}

TWO_LABEL_CRITICAL_EXPECTED = {
    (1, 1, 1, 1, 1): (2681, 60, 66),
    (1, 1, 1, 1, 15): (85, 30, 35),
    (1, 1, 1, 1, 18): (13, 26, 38),
    (1, 1, 1, 1, 30): (1611, 54, 59),
    (1, 1, 1, 1, 45): (4314, 62, 67),
}


def print_table(
    title: str,
    configurations: tuple[tuple[int, int], ...],
    expected: dict[tuple[int, ...], tuple[int, int, int]],
) -> None:
    print(title)
    print("spectrum assignments centre_differences sharp_cells nodes")
    for spectrum, target in expected.items():
        assignments, differences, sharp, nodes = fixed_spectrum_bank(
            configurations,
            spectrum,
        )
        observed = (assignments, differences, sharp)
        require(
            observed == target,
            f"{title} drift at {spectrum}: {observed} != {target}",
        )
        print(
            ",".join(map(str, spectrum)),
            assignments,
            differences,
            sharp,
            nodes,
        )


print("THM-2436 independent fixed-spectrum referee")
print("method=centred-interval-canonical-dynamic-pivot")
print(f"step_one_supports={len(STEP_ONE)}")
print(f"one_source_configurations={len(ONE_SOURCE_CONFIGURATIONS)}")
print(f"two_label_source_multisets={len(TWO_LABEL_CONFIGURATIONS)}")
print_table(
    "ONE_SOURCE_EXCEPTIONAL",
    ONE_SOURCE_CONFIGURATIONS,
    ONE_SOURCE_EXPECTED,
)
print_table(
    "TWO_LABEL_CRITICAL",
    TWO_LABEL_CONFIGURATIONS,
    TWO_LABEL_CRITICAL_EXPECTED,
)
print("one_source_max_sharp=57<65")
print("two_label_critical_max_sharp=67<78")
print("critical_extremum_spectrum=1,1,1,1,45")
print("ALL CHECKS PASSED")
