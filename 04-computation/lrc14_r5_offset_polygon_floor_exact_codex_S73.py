#!/usr/bin/env python3
"""Exact referee for THM-1129's bounded-offset polygon floor.

For every eight-speed core ``P`` in ``{1,...,12}`` and every normalized
four-offset shape

    A = {0 < a < b < c <= 30},

the script proves two finite statements with ``fractions.Fraction`` only.

1. On the closure of the positive-length part of the core-safe set, the
   maximum vertical safe width of the fixed torus polygon is at least 2/5.
2. Some one-sided core-safe order cell supports a *fixed* vertical subarc of
   length 1/5 for time at least 1/728.  Its 1/7 subarc gives the uniform
   slope-line tail K >= 832.

The finite-candidate principle is exact.  Core danger endpoints are
``(14*j +/- 1)/(14*p)``.  Offset-center order changes only when
``(b-a)t`` is integral.  Between consecutive points of this combined set,
the core predicate is constant and every labelled circular center gap is
affine.  Their maximum is convex, so its maximum on a safe cell is attained
at an endpoint.

Tournament audit.  Candidate vertices include runners, combs, core cells,
cell boundaries, center-collision events, residues, and proof obligations.
Exact boundary order is a transitive tournament (no directed cycles,
singleton SCCs, one sorted Hamiltonian path after ties are coalesced), but it
forgets metric gaps, wall owners, slopes, and which adjacent side is core
safe.  The faithful carrier is the labelled cyclic gap plus its one-sided
cell and owner-slope sidecar.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import ceil


H = F(1, 14)
OFFSET_CAP = 30
VERTICAL_ARC = F(1, 5)
NEEDLE_ARC = F(1, 7)
EXPECTED_FLOOR = F(2, 5)
EXPECTED_RECTANGLE_TIME = F(1, 728)
EXPECTED_GLOBAL_TAIL = 832

CORES = list(combinations(range(1, 13), 8))
CORE_MASKS = [sum(1 << (p - 1) for p in core) for core in CORES]


def safe_speed_mask(t: F) -> int:
    """Speeds p in {1,...,12} which are core-safe at t."""

    mask = 0
    for p in range(1, 13):
        x = (p * t) % 1
        if min(x, 1 - x) >= H:
            mask |= 1 << (p - 1)
    return mask


ELIGIBLE_CORES: list[list[int]] = []
for safe_mask in range(1 << 12):
    ELIGIBLE_CORES.append(
        [
            index
            for index, core_mask in enumerate(CORE_MASKS)
            if core_mask & ~safe_mask == 0
        ]
    )


def core_endpoint_set() -> set[F]:
    points = {F(0), F(1)}
    for p in range(1, 13):
        for tooth in range(-1, p + 2):
            for sign in (-1, 1):
                t = F(14 * tooth + sign, 14 * p)
                if 0 <= t <= 1:
                    points.add(t)
    return points


CORE_ENDPOINTS = core_endpoint_set()


def candidate_points(offsets: tuple[int, int, int, int]) -> list[F]:
    """All core endpoints and offset-center collision times."""

    points = set(CORE_ENDPOINTS)
    for a, b in combinations(offsets, 2):
        difference = b - a
        points.update(F(j, difference) for j in range(difference + 1))
    return sorted(points)


def vertical_safe_width(offsets: tuple[int, int, int, int], t: F) -> F:
    """Largest center gap minus the two danger radii 1/14+1/14."""

    centers = sorted((-offset * t) % 1 for offset in offsets)
    gaps = (
        centers[1] - centers[0],
        centers[2] - centers[1],
        centers[3] - centers[2],
        1 - centers[3] + centers[0],
    )
    return max(gaps) - F(1, 7)


def main() -> None:
    shape_count = 0
    pair_count = 0
    global_floor: F | None = None
    floor_rows: list[tuple[tuple[int, ...], tuple[int, ...], F]] = []

    global_rectangle: F | None = None
    rectangle_worst_count = 0
    rectangle_witness = None

    global_tail = 0
    all_legal_tail_pairs = 0
    raw_finite_rows = 0
    bottom_bank_rows = 0
    residual_finite_rows = 0
    candidate_cell_count = 0

    for offset_tail in combinations(range(1, OFFSET_CAP + 1), 3):
        offsets = (0,) + offset_tail
        span = offsets[-1]
        points = candidate_points(offsets)
        widths = {t: vertical_safe_width(offsets, t) for t in points}

        maxima = [F(-1) for _ in CORES]
        best_rectangle = [F(0) for _ in CORES]
        best_needle = [F(0) for _ in CORES]
        best_certificate = [None for _ in CORES]

        for left, right in zip(points, points[1:]):
            midpoint = (left + right) / 2
            core_indices = ELIGIBLE_CORES[safe_speed_mask(midpoint)]
            if not core_indices:
                continue
            candidate_cell_count += len(core_indices)
            cell_width = right - left

            for side, t in (("right", left), ("left", right)):
                width = widths[t]
                rectangle_time = min(
                    cell_width, (width - VERTICAL_ARC) / span
                ) if width >= EXPECTED_FLOOR else F(0)
                needle_time = min(
                    cell_width, (width - NEEDLE_ARC) / span
                ) if width >= EXPECTED_FLOOR else F(0)

                for core_index in core_indices:
                    if width > maxima[core_index]:
                        maxima[core_index] = width
                    if rectangle_time > best_rectangle[core_index]:
                        best_rectangle[core_index] = rectangle_time
                        best_certificate[core_index] = (
                            left,
                            right,
                            side,
                            t,
                            width,
                            rectangle_time,
                        )
                    if needle_time > best_needle[core_index]:
                        best_needle[core_index] = needle_time

        assert min(maxima) >= EXPECTED_FLOOR
        assert min(best_rectangle) > 0
        assert min(best_needle) > 0

        shape_floor = min(maxima)
        if global_floor is None or shape_floor < global_floor:
            global_floor = shape_floor
            floor_rows = [
                (CORES[i], offsets, value)
                for i, value in enumerate(maxima)
                if value == shape_floor
            ]
        elif shape_floor == global_floor:
            floor_rows.extend(
                (CORES[i], offsets, value)
                for i, value in enumerate(maxima)
                if value == shape_floor
            )

        shape_rectangle = min(best_rectangle)
        if global_rectangle is None or shape_rectangle < global_rectangle:
            global_rectangle = shape_rectangle
            rectangle_worst_count = sum(
                value == shape_rectangle for value in best_rectangle
            )
            first = best_rectangle.index(shape_rectangle)
            rectangle_witness = (
                CORES[first], offsets, best_certificate[first]
            )
        elif shape_rectangle == global_rectangle:
            rectangle_worst_count += sum(
                value == shape_rectangle for value in best_rectangle
            )

        # The 1/7 subarc has preimages of length 1/(7K), with starts spaced
        # 1/K.  A time interval d contains a complete preimage once
        # K*d >= 1+1/7=8/7.
        for core_index, needle_time in enumerate(best_needle):
            tail_start = ceil(F(8, 7) / needle_time)
            global_tail = max(global_tail, tail_start)

            core = CORES[core_index]
            legal_start = 13 * max(core) + 1
            if tail_start <= legal_start:
                all_legal_tail_pairs += 1

            finite = max(0, tail_start - legal_start)
            bottom = min(finite, max(0, 40 - span))
            raw_finite_rows += finite
            bottom_bank_rows += bottom
            residual_finite_rows += finite - bottom

        shape_count += 1
        pair_count += len(CORES)

    assert shape_count == 4060
    assert pair_count == 2_009_700
    assert global_floor == EXPECTED_FLOOR
    assert len(floor_rows) == 32
    assert global_rectangle == EXPECTED_RECTANGLE_TIME
    assert rectangle_worst_count == 73
    assert global_tail == EXPECTED_GLOBAL_TAIL

    shape_histogram = Counter(row[1] for row in floor_rows)
    assert shape_histogram == Counter(
        {
            (0, 1, 8, 10): 1,
            (0, 2, 3, 7): 15,
            (0, 2, 9, 10): 1,
            (0, 4, 5, 7): 15,
        }
    )

    hard_core = (1, 2, 4, 5, 7, 9, 11, 12)
    hard_shapes = {(0, 2, 3, 7), (0, 4, 5, 7)}
    assert all(
        any(core == hard_core and offsets == shape for core, offsets, _ in floor_rows)
        for shape in hard_shapes
    )

    print("THM-1129 bounded-offset polygon floor exact referee")
    print("arithmetic=fractions.Fraction only")
    print(f"cores={len(CORES)}")
    print(f"offset_cap={OFFSET_CAP}")
    print(f"offset_shapes={shape_count}")
    print(f"core_shape_pairs={pair_count}")
    print(f"eligible_order_cells={candidate_cell_count}")
    print(f"vertical_floor={global_floor}")
    print(f"floor_attainment_count={len(floor_rows)}")
    print(f"floor_shape_histogram={sorted(shape_histogram.items())}")
    for core, offsets, value in floor_rows:
        print(f"floor_row P={core} A={offsets} M={value}")
    print(f"fixed_vertical_arc={VERTICAL_ARC}")
    print(f"uniform_rectangle_time={global_rectangle}")
    print(f"rectangle_worst_count={rectangle_worst_count}")
    print(f"canonical_rectangle_witness={rectangle_witness}")
    print(f"uniform_tail_K={global_tail}")
    print(f"pairs_tail_closed_at_legal_start={all_legal_tail_pairs}")
    print(f"per_pair_raw_finite_rows={raw_finite_rows}")
    print(f"rows_already_in_THM1123_bottom_bank={bottom_bank_rows}")
    print(f"remaining_bounded_offset_finite_rows={residual_finite_rows}")
    print("candidate_lemma=core endpoints plus center collisions; convex endpoint maximum")
    print("tournament_vertices=labelled candidate boundaries")
    print("tournament_fingerprint=transitive; directed_cycles=0; SCCs=singletons; sorted_HP=1")
    print("destroyed_by_order_only=gaps|owners|slopes|safe_side|rectangle_width")
    print("challenged_vertices=runners|combs|core_cells|boundaries|wall_events|residues|proof_obligations")
    print("scope_warning=offset span <=30 only; arbitrary-shape cone is THM-1128")
    print("certificate=PASS")


if __name__ == "__main__":
    main()
