#!/usr/bin/env python3
"""Independent Fraction replay for THM-1133's accounting and extremal.

The full 3,539,936-row quantifier belongs to the C++ referee.  This script is
an independent guard with a deliberately different endpoint engine:

* it constructs a safe set from the complete simultaneous breakpoint
  arrangement and midpoint classification, rather than successive tooth
  subtraction;
* it recomputes the THM-1133 extremal and all of its longest intervals;
* it independently reconstructs the extremal pair's one-sided polygon tail;
* it cross-checks THM-1123's shared row and the THM-1129/1133 accounting
  identity from the two frozen result ledgers.

All arithmetic uses ``fractions.Fraction``.  This is an extremal/accounting
replay, not a second claim to enumerate all 3,539,936 rows.
"""

from __future__ import annotations

import ast
from fractions import Fraction as F
from itertools import combinations
from math import ceil, comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "05-knowledge" / "results"
H = F(1, 14)


def distance(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def arrangement_safe(speeds: tuple[int, ...]) -> list[tuple[F, F]]:
    """Complete endpoint arrangement, independent of staged subtraction."""

    points = {F(0), F(1)}
    for speed in speeds:
        for tooth in range(speed + 1):
            for sign in (-1, 1):
                point = F(14 * tooth + sign, 14 * speed)
                if 0 <= point <= 1:
                    points.add(point)
    cells: list[tuple[F, F]] = []
    ordered = sorted(points)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(distance(speed * midpoint) >= H for speed in speeds):
            if cells and cells[-1][1] == left:
                cells[-1] = (cells[-1][0], right)
            else:
                cells.append((left, right))
    return cells


def vertical_width(offsets: tuple[int, ...], t: F) -> F:
    centers = sorted((-offset * t) % 1 for offset in offsets)
    gaps = [
        centers[1] - centers[0],
        centers[2] - centers[1],
        centers[3] - centers[2],
        1 - centers[3] + centers[0],
    ]
    return max(gaps) - F(1, 7)


def direct_tail_certificate(
    core: tuple[int, ...], offsets: tuple[int, ...]
) -> tuple[F, int, tuple[F, F, str, F, F, F]]:
    """Rebuild one pair's conservative THM-1129 one-sided tail."""

    points = {F(0), F(1)}
    # Match the proved universal refinement: all possible core endpoints,
    # including harmless splits owned by speeds outside this particular core.
    for speed in range(1, 13):
        for tooth in range(-1, speed + 2):
            for sign in (-1, 1):
                point = F(14 * tooth + sign, 14 * speed)
                if 0 <= point <= 1:
                    points.add(point)
    for a, b in combinations(offsets, 2):
        difference = b - a
        points.update(F(j, difference) for j in range(difference + 1))

    best = F(0)
    witness = None
    ordered = sorted(points)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if not all(distance(p * midpoint) >= H for p in core):
            continue
        for side, t in (("right", left), ("left", right)):
            width = vertical_width(offsets, t)
            if width < F(2, 5):
                continue
            certificate = min(
                right - left, (width - F(1, 7)) / offsets[-1]
            )
            if certificate > best:
                best = certificate
                witness = (left, right, side, t, width, certificate)
    assert witness is not None
    return best, ceil(F(8, 7) / best), witness


def parse_ledger(path: Path) -> dict[str, str]:
    ledger = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            ledger[key] = value
    return ledger


def main() -> None:
    hard_core = (1, 2, 4, 5, 7, 8, 9, 11)
    hard_offsets = (0, 4, 6, 8)
    hard_k = 186
    hard_speeds = hard_core + tuple(hard_k + offset for offset in hard_offsets)
    hard_region = arrangement_safe(hard_speeds)
    hard_L = max(right - left for left, right in hard_region)
    hard_intervals = [
        interval
        for interval in hard_region
        if interval[1] - interval[0] == hard_L
    ]
    assert hard_L == F(86, 61845)
    assert 7 * (hard_k + hard_offsets[-1]) * hard_L == F(16684, 8835)
    assert hard_intervals == [
        (F(267, 2660), F(265, 2604)),
        (F(2339, 2604), F(2393, 2660)),
    ]

    tail_width, tail_start, tail_witness = direct_tail_certificate(
        hard_core, hard_offsets
    )
    assert tail_width == F(1, 308)
    assert tail_start == 352
    assert hard_k < tail_start
    legal_start = 13 * max(hard_core) + 1
    first_residual = legal_start + 40 - hard_offsets[-1]
    assert (legal_start, first_residual) == (144, 176)
    assert first_residual <= hard_k

    shared_core = (1, 2, 4, 5, 7, 9, 11, 12)
    shared_killers = (158, 160, 162, 164)
    shared_region = arrangement_safe(shared_core + shared_killers)
    shared_L = max(right - left for left, right in shared_region)
    assert shared_L == F(41, 25920)
    assert 7 * shared_killers[-1] * shared_L == F(11767, 6480)

    thm1129 = parse_ledger(
        RESULTS / "lrc14_r5_offset_polygon_floor_exact_codex_S73.out"
    )
    thm1133 = parse_ledger(
        RESULTS
        / "lrc14_r5_bounded_offset_all_scale_complement_exact_codex_S73.out"
    )
    assert int(thm1133["core_shape_pairs"]) == 495 * comb(30, 3)
    assert thm1133["maximum_individual_tail_K"] == thm1129["uniform_tail_K"]
    for key in (
        "pairs_tail_closed_at_legal_start",
        "per_pair_raw_finite_rows",
        "rows_already_in_THM1123_bottom_bank",
    ):
        assert thm1133[key] == thm1129[key]
    raw = int(thm1133["per_pair_raw_finite_rows"])
    bottom = int(thm1133["rows_already_in_THM1123_bottom_bank"])
    expected = int(thm1133["residual_rows_expected"])
    tested = int(thm1133["residual_rows_tested"])
    assert raw - bottom == expected == tested == 3_539_936
    assert int(thm1133["failures_total"]) == 0
    assert F(thm1133["hardest_longest_component"]) == hard_L
    assert ast.literal_eval(thm1133["hardest_offsets"]) == hard_offsets

    print("THM-1133 independent Fraction extremal/accounting replay")
    print("method=complete simultaneous breakpoint arrangement")
    print(f"hardest_core={hard_core}")
    print(f"hardest_offsets={hard_offsets}")
    print(f"hardest_K={hard_k}")
    print(f"hardest_L={hard_L}")
    print(f"hardest_7k4L={7*(hard_k+hard_offsets[-1])*hard_L}")
    print(f"hardest_intervals={hard_intervals}")
    print(f"hardest_pair_rectangle_time={tail_width}")
    print(f"hardest_pair_tail_start={tail_start}")
    print(f"hardest_pair_rectangle_witness={tail_witness}")
    print(f"hardest_pair_legal_start={legal_start}")
    print(f"hardest_pair_first_residual_K={first_residual}")
    print(f"THM1123_shared_L={shared_L}")
    print(f"THM1123_shared_7k4L={7*shared_killers[-1]*shared_L}")
    print(f"accounting_pairs={495*comb(30,3)}")
    print(f"accounting_raw={raw}")
    print(f"accounting_bottom={bottom}")
    print(f"accounting_residual={expected}")
    print("full_quantifier_owner=C++ referee; this replay guards extremal/accounting")
    print("tournament=endpoint order transitive; metric/owner sidecars retained")
    print("certificate=PASS")


if __name__ == "__main__":
    main()
