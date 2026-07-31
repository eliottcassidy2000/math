#!/usr/bin/env python3
"""Coefficient audit of the unique nearest half-step R-to-M collar."""

from collections import Counter
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
    require(all(a == b for _l, _r, a, b in aligned), "unequal common weight")
    common = tuple((left, right, a) for left, right, a, _b in aligned)
    right = copies.physical.subtract_weighted(target_pullback, common)
    return common, right


def circular_offset(rleft, mleft):
    offset = ((rleft - mleft) // H) % PERIOD_STEPS
    if offset > PERIOD_STEPS // 2:
        offset -= PERIOD_STEPS
    return offset


def semantic_value(interval, q_pair, cache):
    key = interval[:2]
    if key not in cache:
        unit = ((*key, 1),)
        target = copies.physical.overlap.shift_weighted(unit, copies.SHIFT)
        source_value = copies.physical.relative.private.delayed_carry_pair(
            unit, q_pair, {}
        )[12][1]
        target_value = copies.physical.relative.private.delayed_carry_pair(
            target, q_pair, {}
        )[6][1]
        cache[key] = (source_value, target_value)
    return cache[key]


def main():
    (
        _module, _rails, _present, details, full_module, e3, clocks,
        q_pairs, _delayed, _sw, _tw, _rail,
    ) = copies.physical_setup()

    nearest_contingency = Counter()
    right_status = Counter()
    common_status = Counter()
    same_status_offsets = Counter()
    same_status_ties = Counter()
    same_status_collision_cells = []
    dead_right = []
    relation_parity = Counter()
    relation_edges = 0
    pair_count = 0

    for clock in range(7):
        for sigma in COMMON_S:
            for target in COMMON_T:
                common, right = cell_objects(
                    details, full_module, e3, clocks, clock, sigma, target
                )
                require(bool(common) == bool(right), "support shadows split")
                if not right:
                    continue
                q_pair = q_pairs[clock]
                cache = {}
                common_by_left = {row[0]: row for row in common}
                image_lefts = []
                for rpiece in right:
                    rleft, rstop, weight = rpiece
                    nearest = common_by_left.get(rleft + H)
                    require(nearest is not None, "nearest +h common collar missing")
                    require(
                        nearest[1] - nearest[0] == rstop - rleft
                        and nearest[2] == weight,
                        "nearest collar changed length or weight",
                    )
                    rvalue = semantic_value(rpiece, q_pair, cache)
                    mvalue = semantic_value(nearest, q_pair, cache)
                    right_status[rvalue] += 1
                    if rvalue == (0, 0):
                        dead_right.append(((clock, sigma, target), rpiece[:2]))
                    common_status[mvalue] += 1
                    nearest_contingency[rvalue, mvalue] += 1
                    pair_count += 1

                    candidates = []
                    for mpiece in common:
                        if (
                            mpiece[1] - mpiece[0] != rstop - rleft
                            or mpiece[2] != weight
                        ):
                            continue
                        mcontent = semantic_value(mpiece, q_pair, cache)
                        offset = circular_offset(rleft, mpiece[0])
                        same = mcontent == rvalue
                        relation_parity[(offset % 2, same)] += 1
                        relation_edges += 1
                        if not same:
                            continue
                        candidates.append((abs(offset), offset, mpiece[0]))
                    require(candidates, "right piece has no same-content common mate")
                    minimum = min(row[0] for row in candidates)
                    closest = tuple(row for row in candidates if row[0] == minimum)
                    same_status_ties[len(closest)] += 1
                    if len(closest) == 1:
                        same_status_offsets[closest[0][1]] += 1
                        image_lefts.append(closest[0][2])
                    else:
                        image_lefts.append(None)
                collisions = len(image_lefts) - len(set(image_lefts))
                if collisions:
                    same_status_collision_cells.append(
                        ((clock, sigma, target), collisions, tuple(image_lefts))
                    )

    print("LRC NEAREST HALF-STEP SEMANTIC COLLAR SCOUT")
    print(f"half_step={H};paired_right_pieces={pair_count}")
    print(f"right_status={tuple(sorted(right_status.items()))}")
    print(f"nearest_common_status={tuple(sorted(common_status.items()))}")
    print(f"nearest_contingency={tuple(sorted(nearest_contingency.items()))}")
    require(
        relation_parity[(0, True)] + relation_parity[(1, False)]
        == relation_edges
        and relation_parity[(0, False)] == 0
        and relation_parity[(1, True)] == 0,
        "semantic equality stopped agreeing with half-step parity",
    )
    print(
        f"relation_edges={relation_edges};"
        f"relation_parity={tuple(sorted(relation_parity.items()))}"
    )
    print(f"same_status_tie_histogram={tuple(sorted(same_status_ties.items()))}")
    print(f"same_status_signed_offset_histogram={tuple(sorted(same_status_offsets.items()))}")
    print(f"dead_right={tuple(dead_right)}")
    print(
        f"same_status_collision_cells={len(same_status_collision_cells)};"
        f"first={tuple(same_status_collision_cells[:8])}"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
