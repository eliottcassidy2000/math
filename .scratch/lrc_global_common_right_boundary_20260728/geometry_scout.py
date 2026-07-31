#!/usr/bin/env python3
"""Exact geometry scout for the 567-cell common/right bank.

Scratch only.  It asks whether equal-length M and R interval components are
related by the half-step lattice, before any coefficient or endpoint claim.
"""

from collections import Counter, defaultdict
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_semantic_arm_right_wing_central_digit_thm2782 as arm


P = 13
H = arm.T // (2 * P**5)
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    context = arm.build_context()
    module, _delayed, _present, source, clocks, _weight, rail = context
    cells = []
    exact_length_cover = 0
    half_lattice_cover = 0
    unique_offset_cover = 0
    offset_histogram = Counter()
    candidate_histogram = Counter()
    unmatched = []
    per_cell = []

    for clock in range(7):
        for sigma in COMMON_S:
            for target in COMMON_T:
                objects = arm.carrier_objects(
                    module, source, clocks, rail, clock, sigma, target
                )
                common = objects["M"]
                right = objects["R"]
                require(bool(common) == bool(right), "support shadows split")
                if not right:
                    continue
                cells.append((clock, sigma, target))
                common_by_length = defaultdict(list)
                for left, stop in common:
                    common_by_length[stop - left].append((left, stop))
                cell_offsets = Counter()
                cell_exact = cell_lattice = cell_unique = 0
                for rleft, rstop in right:
                    candidates = common_by_length[rstop - rleft]
                    if candidates:
                        exact_length_cover += 1
                        cell_exact += 1
                    lattice = []
                    for mleft, mstop in candidates:
                        delta = rleft - mleft
                        if delta % H == 0:
                            lattice.append(delta // H)
                    candidate_histogram[len(lattice)] += 1
                    if lattice:
                        half_lattice_cover += 1
                        cell_lattice += 1
                        cell_offsets.update(lattice)
                    else:
                        unmatched.append(
                            ((clock, sigma, target), (rleft, rstop), len(candidates))
                        )
                    if len(lattice) == 1:
                        unique_offset_cover += 1
                        cell_unique += 1
                        offset_histogram[lattice[0]] += 1
                per_cell.append((
                    (clock, sigma, target), len(common), len(right),
                    cell_exact, cell_lattice, cell_unique,
                    tuple(sorted(cell_offsets.items())),
                ))

    total_right = sum(row[2] for row in per_cell)
    print("LRC GLOBAL COMMON/RIGHT BOUNDARY GEOMETRY SCOUT")
    print(f"half_step={H};cells={len(cells)};total_right_pieces={total_right}")
    print(
        f"equal_length_cover={exact_length_cover}/{total_right};"
        f"half_lattice_cover={half_lattice_cover}/{total_right};"
        f"unique_half_lattice_offset={unique_offset_cover}/{total_right}"
    )
    print(f"lattice_candidate_count_histogram={tuple(sorted(candidate_histogram.items()))}")
    print(f"unique_offset_histogram={tuple(sorted(offset_histogram.items()))}")
    print(f"unmatched_count={len(unmatched)};first_unmatched={tuple(unmatched[:12])}")
    print(f"per_cell={tuple(per_cell)}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
