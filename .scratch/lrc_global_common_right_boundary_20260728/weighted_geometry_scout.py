#!/usr/bin/env python3
"""Full semantic weighted common/right interval matching scout."""

from collections import Counter, defaultdict
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_right_cofiber_positive_copy_stratification_thm2818 as copies


P = 13
H = copies.T // (2 * P**5)
PERIOD_STEPS = copies.T // H
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cell_objects(details, full_module, e3, clocks, clock, sigma, target):
    section = copies.physical.target.source_present_section(
        full_module, e3, clock, sigma, target, clocks
    )
    source_base, target_base = details[clock]
    source = copies.weighted_intersection(source_base, section)
    target_physical = copies.weighted_intersection(target_base, section)
    target_pullback = copies.physical.overlap.shift_weighted(
        target_physical, -copies.SHIFT
    )
    aligned = copies.physical.overlap.intersect_weighted_profiles(
        source, target_pullback
    )
    common = tuple(
        (left, right, a)
        for left, right, a, b in aligned
        if a == b
    )
    unequal = tuple(row for row in aligned if row[2] != row[3])
    require(not unequal, "common overlap acquired unequal weights")
    right = copies.physical.subtract_weighted(target_pullback, common)
    return source, target_physical, common, right


def main():
    (
        _module, _rails, _present, details, full_module, e3, clocks,
        _q_pairs, _delayed, _sw, _tw, _rail,
    ) = copies.physical_setup()

    cells = []
    total_common = total_right = 0
    exact_cover = lattice_cover = unique_cover = 0
    candidate_hist = Counter()
    unique_offsets = Counter()
    nearest_offsets = Counter()
    nearest_ties = Counter()
    nearest_collision_cells = []
    length_ratio_hist = Counter()
    unmatched = []
    cell_summary = []

    for clock in range(7):
        for sigma in COMMON_S:
            for target in COMMON_T:
                _source, _target, common, right = cell_objects(
                    details, full_module, e3, clocks, clock, sigma, target
                )
                require(bool(common) == bool(right), "support shadows split")
                if not right:
                    continue
                cells.append((clock, sigma, target))
                total_common += len(common)
                total_right += len(right)
                common_by_key = defaultdict(list)
                for left, stop, weight in common:
                    common_by_key[(stop - left, weight)].append((left, stop))
                local = Counter()
                local_match = 0
                nearest_indices = []
                for rleft, rstop, weight in right:
                    candidates = common_by_key[(rstop - rleft, weight)]
                    if candidates:
                        exact_cover += 1
                    lattice = tuple(
                        (rleft - mleft) // H
                        for mleft, _mstop in candidates
                        if (rleft - mleft) % H == 0
                    )
                    candidate_hist[len(lattice)] += 1
                    if lattice:
                        lattice_cover += 1
                        local_match += 1
                        local.update(lattice)
                    else:
                        unmatched.append(
                            ((clock, sigma, target), (rleft, rstop, weight),
                             len(candidates))
                        )
                    if len(lattice) == 1:
                        unique_cover += 1
                        unique_offsets[lattice[0]] += 1
                    circular = []
                    for index, (mleft, _mstop) in enumerate(candidates):
                        offset = ((rleft - mleft) // H) % PERIOD_STEPS
                        if offset > PERIOD_STEPS // 2:
                            offset -= PERIOD_STEPS
                        circular.append((abs(offset), offset, index))
                    minimum = min(row[0] for row in circular)
                    nearest = tuple(row for row in circular if row[0] == minimum)
                    nearest_ties[len(nearest)] += 1
                    if len(nearest) == 1:
                        nearest_offsets[nearest[0][1]] += 1
                        nearest_indices.append(nearest[0][2])
                    else:
                        nearest_indices.append(None)
                collisions = len(nearest_indices) - len(set(nearest_indices))
                if collisions:
                    nearest_collision_cells.append(
                        ((clock, sigma, target), collisions, tuple(nearest_indices))
                    )
                for _left, _stop, weight in common:
                    length_ratio_hist[weight] += 1
                cell_summary.append((
                    (clock, sigma, target), len(common), len(right),
                    local_match, tuple(sorted(local.items())),
                ))

    print("LRC FULL-SEMANTIC WEIGHTED COMMON/RIGHT GEOMETRY SCOUT")
    print(
        f"half_step={H};nonempty_cells={len(cells)};"
        f"common_pieces={total_common};right_pieces={total_right}"
    )
    print(
        f"equal_length_weight_cover={exact_cover}/{total_right};"
        f"half_lattice_cover={lattice_cover}/{total_right};"
        f"unique_half_lattice_offset={unique_cover}/{total_right}"
    )
    print(f"candidate_count_histogram={tuple(sorted(candidate_hist.items()))}")
    print(f"unique_offset_histogram={tuple(sorted(unique_offsets.items()))}")
    print(f"nearest_tie_histogram={tuple(sorted(nearest_ties.items()))}")
    print(f"nearest_signed_offset_histogram={tuple(sorted(nearest_offsets.items()))}")
    print(
        f"nearest_collision_cells={len(nearest_collision_cells)};"
        f"first_collisions={tuple(nearest_collision_cells[:12])}"
    )
    print(f"weight_piece_histogram={tuple(sorted(length_ratio_hist.items()))}")
    print(f"unmatched_count={len(unmatched)};first_unmatched={tuple(unmatched[:20])}")
    print(
        "complete_bipartite_cells="
        f"{sum(sum(count for _offset, count in row[4]) == row[1] * row[2] for row in cell_summary)}/"
        f"{len(cell_summary)}"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
