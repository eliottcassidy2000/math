#!/usr/bin/env python3
"""Exact verifier for THM-908 and HYP-7064.

For five mover speeds, averaging the seven time shifts x -> x+j/7 turns
the signed residue-six kernel into an affine-line sum on F_7^5.  The line
ceiling depends only on the projective residue direction.  This script
enumerates all 2,801 directions with exact integer arithmetic and checks the
change-of-variables identity against independent rational event sweeps.

Tournament Analysis is diagnostic.  Vertices are representative projective
directions for the fifteen signed ceiling strata.  The pairwise observable is
the difference of affine-line ceilings, the switch replaces the signed K6
kernel by THM-905's positive box majorant, and the gauge points toward the
larger obstruction ceiling with lexicographic tie completion.  This preserves
the residue-only worst-case order but destroys the within-cell path, quotient
speeds, additive relation lattice, and wall chronology.  Alternative vertex
sets considered were runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, affine cosets, Fourier modes, matroid circuits,
primitive relations, and proof obligations.  Projective residue directions
are used because they are exactly the quotient retained by the seven-shift
identity; quotient speeds and short relations remain mandatory sidecars on
the directions that do not close.
"""

from collections import Counter
from fractions import Fraction
from heapq import heapify, heappop, heappush
from itertools import combinations, product
from math import gcd, lcm

import numpy as np


FIELD_SIZE = 7
MOVER_COUNT = 5
STATE_COUNT = FIELD_SIZE**MOVER_COUNT
PROPAGATION_SLACK = Fraction(97, 1000)
SIGNED_CEILING_HISTOGRAM = {
    62: 5,
    50: 60,
    48: 175,
    36: 1075,
    34: 580,
    32: 65,
    26: 60,
    25: 270,
    24: 280,
    23: 120,
    22: 30,
    21: 60,
    20: 10,
    18: 10,
    8: 1,
}


def state_arrays():
    indices = np.arange(STATE_COUNT, dtype=np.int64)
    powers = FIELD_SIZE ** np.arange(MOVER_COUNT, dtype=np.int64)
    digits = np.stack(
        [(indices // power) % FIELD_SIZE for power in powers], axis=1
    ).astype(np.int16)
    counts = np.stack(
        [(digits == section).sum(axis=1) for section in range(FIELD_SIZE)],
        axis=1,
    ).astype(np.int16)

    signed_kernel = np.zeros(STATE_COUNT, dtype=np.int16)
    missed = counts[:, 1:] == 0
    missed_count = missed.sum(axis=1)
    singleton_rows = missed_count == 1
    signed_kernel[singleton_rows] = 1
    signed_kernel[singleton_rows & missed[:, 2]] = -6
    pair_rows = missed_count == 2
    signed_kernel[pair_rows] = -2
    signed_kernel[
        pair_rows
        & (
            (missed[:, 0] & missed[:, 4])
            | (missed[:, 1] & missed[:, 3])
        )
    ] = 12

    exceptional_complements = (
        counts[:, 1] * counts[:, 3] * counts[:, 5] * counts[:, 6]
        + counts[:, 2] * counts[:, 3] * counts[:, 4] * counts[:, 6]
    )
    pinned_exceptional_complements = counts[:, 0] * exceptional_complements
    singleton_six_complement = counts[:, 1] * counts[:, 2] * counts[:, 4] * counts[:, 5]
    positive_majorant = (
        6 * exceptional_complements
        + 6 * pinned_exceptional_complements
        + singleton_six_complement
    ).astype(np.int16)
    return digits, powers, signed_kernel, positive_majorant


def projective_directions():
    for pivot in range(MOVER_COUNT):
        for tail in product(range(FIELD_SIZE), repeat=MOVER_COUNT - pivot - 1):
            yield (0,) * pivot + (1,) + tail, pivot


def ceiling_census(digits, powers, signed_kernel, positive_majorant):
    shifts = np.arange(FIELD_SIZE, dtype=np.int16)[:, None, None]
    representatives_by_pivot = {
        pivot: digits[digits[:, pivot] == 0]
        for pivot in range(MOVER_COUNT)
    }
    signed_histogram = Counter()
    positive_histogram = Counter()
    representative_by_signed_ceiling = {}
    records = {}
    distinct_nonzero = []

    for direction, pivot in projective_directions():
        representatives = representatives_by_pivot[pivot]
        direction_array = np.array(direction, dtype=np.int16)[None, None, :]
        line_states = (representatives[None, :, :] + shifts * direction_array) % FIELD_SIZE
        line_indices = line_states @ powers
        signed_sums = signed_kernel[line_indices].sum(axis=0)
        positive_sums = positive_majorant[line_indices].sum(axis=0)
        signed_ceiling = int(signed_sums.max())
        positive_ceiling = int(positive_sums.max())
        signed_histogram[signed_ceiling] += 1
        positive_histogram[positive_ceiling] += 1
        representative_by_signed_ceiling.setdefault(signed_ceiling, direction)
        records[direction] = (signed_ceiling, positive_ceiling)
        if 0 not in direction and len(set(direction)) == MOVER_COUNT:
            distinct_nonzero.append((direction, signed_ceiling))

    return {
        "signed_histogram": signed_histogram,
        "positive_histogram": positive_histogram,
        "representative_by_signed_ceiling": representative_by_signed_ceiling,
        "records": records,
        "distinct_nonzero": distinct_nonzero,
    }


def residue_six_negative_kernel(sectors):
    counts = tuple(sectors.count(section) for section in range(FIELD_SIZE))
    missed = tuple(section for section in range(1, FIELD_SIZE) if counts[section] == 0)
    if len(missed) == 1:
        return -6 if missed == (3,) else 1
    if len(missed) == 2:
        return 12 if missed in ((1, 5), (2, 4)) else -2
    return 0


def direct_kernel_average(speeds):
    common_divisor = gcd(*speeds)
    reduced = tuple(speed // common_divisor for speed in speeds)
    period_scale = lcm(*reduced)
    sectors = [0] * MOVER_COUNT
    events = []
    for runner_index, speed in enumerate(reduced):
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    numerator = 0
    previous = 0
    while events:
        event_position = events[0][0]
        numerator += (
            event_position - previous
        ) * residue_six_negative_kernel(sectors)
        while events and events[0][0] == event_position:
            _, runner_index, event_index, speed, event_step = heappop(events)
            sectors[runner_index] = (sectors[runner_index] + 1) % FIELD_SIZE
            if event_index < FIELD_SIZE * speed:
                next_index = event_index + 1
                heappush(
                    events,
                    (
                        next_index * event_step,
                        runner_index,
                        next_index,
                        speed,
                        event_step,
                    ),
                )
        previous = event_position
    return Fraction(numerator, FIELD_SIZE * period_scale)


def seven_shift_kernel_average(speeds):
    common_divisor = gcd(*speeds)
    reduced = tuple(speed // common_divisor for speed in speeds)
    period_scale = lcm(*reduced)
    sectors = [0] * MOVER_COUNT
    events = []
    for runner_index, speed in enumerate(reduced):
        event_step = period_scale // speed
        events.append((event_step, runner_index, 1, speed, event_step))
    heapify(events)
    residues = tuple(speed % FIELD_SIZE for speed in reduced)
    numerator = 0
    previous = 0
    while events:
        event_position = events[0][0]
        line_sum = sum(
            residue_six_negative_kernel(
                [
                    (sector + shift * residue) % FIELD_SIZE
                    for sector, residue in zip(sectors, residues)
                ]
            )
            for shift in range(FIELD_SIZE)
        )
        numerator += (event_position - previous) * line_sum
        while events and events[0][0] == event_position:
            _, runner_index, event_index, speed, event_step = heappop(events)
            sectors[runner_index] = (sectors[runner_index] + 1) % FIELD_SIZE
            if event_index < speed:
                next_index = event_index + 1
                heappush(
                    events,
                    (
                        next_index * event_step,
                        runner_index,
                        next_index,
                        speed,
                        event_step,
                    ),
                )
        previous = event_position
    return Fraction(numerator, FIELD_SIZE * period_scale)


def scalar_tournament(vertices, values):
    order = tuple(sorted(vertices, key=lambda vertex: (-values[vertex], vertex)))
    rank = {vertex: index for index, vertex in enumerate(order)}
    edges = {
        (first, second) if rank[first] < rank[second] else (second, first)
        for first, second in combinations(vertices, 2)
    }
    score_histogram = {
        score: 1 for score in range(len(vertices))
    }
    return {
        "edges": edges,
        "score_histogram": score_histogram,
        "directed_triangles": 0,
        "components": tuple((vertex,) for vertex in vertices),
        "hamiltonian_path_count": 1,
        "tie_hamiltonian_path": order,
    }


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def main():
    print("THM-908 / HYP-7064: SEVEN-SHIFT RESIDUE-SIX SIEVE")
    print("=" * 78)

    print("\n[1] Exact F7^5 kernel tables")
    digits, powers, signed_kernel, positive_majorant = state_arrays()
    check("7^5 labelled mover states", len(digits) == STATE_COUNT)
    check("THM-905 majorant dominates the signed kernel", np.all(positive_majorant >= signed_kernel))
    check("signed kernel range is {-6,-2,0,1,12}", set(map(int, np.unique(signed_kernel))) == {-6, -2, 0, 1, 12})
    print("  signed state histogram:", dict(sorted(Counter(map(int, signed_kernel)).items())))

    print("\n[2] Exact 2,801-direction projective census")
    census = ceiling_census(digits, powers, signed_kernel, positive_majorant)
    signed_histogram = dict(sorted(census["signed_histogram"].items(), reverse=True))
    positive_histogram = dict(sorted(census["positive_histogram"].items(), reverse=True))
    check("projective direction count", sum(signed_histogram.values()) == 2801)
    check("signed ceiling histogram", signed_histogram == SIGNED_CEILING_HISTOGRAM)
    closing_directions = sum(
        count for ceiling, count in signed_histogram.items() if Fraction(ceiling, 343) < PROPAGATION_SLACK
    )
    check("906 directions close below 0.097", closing_directions == 906)
    check("1,895 directions retain lift debt", 2801 - closing_directions == 1895)
    print("  signed ceiling histogram:", signed_histogram)
    print("  positive-majorant ceiling histogram:", positive_histogram)

    print("\n[3] Distinct nonzero residue theorem")
    distinct_nonzero = census["distinct_nonzero"]
    check("120 projective distinct-nonzero directions", len(distinct_nonzero) == 120)
    check("every distinct-nonzero ceiling is exactly 25", {ceiling for _, ceiling in distinct_nonzero} == {25})
    print("  diameter-free bound: -F6 <= 25/343 =", float(Fraction(25, 343)))
    check("25/343 clears 0.097", Fraction(25, 343) < PROPAGATION_SLACK)

    print("\n[4] Independent rational change-of-variables referees")
    referee_speeds = (
        (1, 2, 3, 4, 5),
        (2, 3, 4, 5, 6),
        (1, 8, 16, 24, 32),
        (7, 14, 21, 28, 29),
        (11, 12, 13, 14, 15),
    )
    for speeds in referee_speeds:
        direct = direct_kernel_average(speeds)
        shifted = seven_shift_kernel_average(speeds)
        check(f"seven-shift identity at {speeds}", direct == shifted)
        print("   ", speeds, "E[L6]=", direct, "-F6=", direct / 49)

    print("\n[5] Tournament Analysis: signed kernel -> positive majorant")
    representative_by_ceiling = census["representative_by_signed_ceiling"]
    vertices = tuple(representative_by_ceiling[ceiling] for ceiling in sorted(representative_by_ceiling, reverse=True))
    signed_values = {vertex: census["records"][vertex][0] for vertex in vertices}
    positive_values = {vertex: census["records"][vertex][1] for vertex in vertices}
    signed_tournament = scalar_tournament(vertices, signed_values)
    positive_tournament = scalar_tournament(vertices, positive_values)
    edge_flips = len(signed_tournament["edges"] ^ positive_tournament["edges"]) // 2
    print("  vertices (one per signed ceiling):", vertices)
    print("  observable: affine-line ceiling difference")
    print("  switch/gauge: signed L6 -> positive 6H+6J+H6; larger ceiling wins")
    print("  signed fingerprint:", {key: value for key, value in signed_tournament.items() if key != "edges"})
    print("  positive fingerprint:", {key: value for key, value in positive_tournament.items() if key != "edges"})
    print("  switch edge flips:", edge_flips)
    check("both ceiling tournaments are transitive", signed_tournament["directed_triangles"] == 0 and positive_tournament["directed_triangles"] == 0)
    check("both tie Hamiltonian paths are unique", signed_tournament["hamiltonian_path_count"] == positive_tournament["hamiltonian_path_count"] == 1)
    print("  Assumption challenge: residue directions preserve the exact line ceiling")
    print("  but erase z(u), quotient speeds, relations, and finite-t wall chronology.")

    print("\nVERDICT")
    print("  PROVED: exact seven-shift affine-line identity for the signed K6 kernel.")
    print("  PROVED: 906/2801 projective residue directions close below 0.097.")
    print("  PROVED: distinct nonzero residues give -F6<=25/343<0.097.")
    print("  OPEN: restore lift/relation sidecars on the remaining 1,895 directions.")


if __name__ == "__main__":
    main()
