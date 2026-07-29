#!/usr/bin/env python3
"""Exact all-cell companion for THM-2825.

The computation works on the complete 7 x 9 x 9 semantic common/right bank
above THM-2806.  It proves that the unscaled equal-copy relation is complete
bipartite in every nonempty cell, that the circular metric selects a unique
right-to-common neighbor at one half-step, and that delayed semantic content
is exactly the parity character of the half-step displacement.  The separate
physical companion retains native-factor, carrier, and endpoint sidecars.
No Python ``assert`` is used.
"""

from collections import Counter
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_right_cofiber_positive_copy_stratification_thm2818.py":
        "85edac9bb03f1fef198343268f4fb1226bec122d45ded79a049f8fa9a73882a8",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")


import lrc14_right_cofiber_positive_copy_stratification_thm2818 as copies


P = 13
H = copies.T // (2 * P**5)
PERIOD_STEPS = copies.T // H
SEMANTIC_PERIOD = 2 * H
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)

EXPECTED_CELL_TYPES = Counter({
    (1, 189, 3): 54,
    (1, 213, 2): 18,
    (1, 239, 2): 9,
    (2, 190, 3): 7,
    (2, 426, 4): 28,
    (2, 474, 2): 14,
    (2, 526, 2): 7,
    (3, 382, 5): 7,
    (3, 404, 4): 28,
    (3, 409, 3): 7,
    (3, 452, 2): 7,
    (3, 504, 2): 7,
})

EXPECTED_CANDIDATES = Counter({
    189: 162,
    190: 21,
    213: 36,
    239: 18,
    382: 35,
    404: 112,
    409: 21,
    426: 112,
    452: 14,
    474: 28,
    504: 14,
    526: 14,
})

EXPECTED_DEAD_RIGHT = (
    ((2, 3, 3), (145528486363260, 145528512808140)),
    ((2, 3, 4), (145528486363260, 145528512808140)),
    ((2, 3, 5), (145528486363260, 145528512808140)),
    ((2, 3, 6), (145528486363260, 145528512808140)),
    ((2, 3, 7), (145528486363260, 145528512808140)),
    ((2, 3, 8), (145528486363260, 145528512808140)),
    ((2, 3, 9), (145528486363260, 145528512808140)),
    ((3, 2, 3), (147164895537660, 147164921982540)),
    ((3, 2, 4), (147164895537660, 147164921982540)),
    ((3, 2, 5), (147164895537660, 147164921982540)),
    ((3, 2, 6), (147164895537660, 147164921982540)),
    ((3, 2, 7), (147164895537660, 147164921982540)),
    ((3, 2, 10), (147164895537660, 147164921982540)),
    ((3, 2, 11), (147164895537660, 147164921982540)),
)


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
    require(
        all(a == b for _left, _right, a, b in aligned),
        "common overlap acquired unequal weights",
    )
    common = tuple(
        (left, right, a) for left, right, a, _b in aligned
    )
    right = copies.physical.subtract_weighted(target_pullback, common)
    return source, target_physical, common, right


def circular_offset(rleft, mleft):
    difference = rleft - mleft
    require(difference % H == 0, "relation left the half-step lattice")
    offset = (difference // H) % PERIOD_STEPS
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
        q_pairs, _delayed, _source_weight, _target_weight, _rail_common,
    ) = copies.physical_setup()

    nonempty = 0
    empty = 0
    total_common = 0
    total_right = 0
    relation_edges = 0
    cell_types = Counter()
    candidate_histogram = Counter()
    relation_parity = Counter()
    nearest_offsets = Counter()
    same_content_offsets = Counter()
    right_status = Counter()
    nearest_status = Counter()
    second_status = Counter()
    nearest_contingency = Counter()
    dead_right = []
    semantic_caches = tuple({} for _clock in range(7))

    for clock in range(7):
        for sigma in COMMON_S:
            for target in COMMON_T:
                source, target_physical, common, right = cell_objects(
                    details, full_module, e3, clocks, clock, sigma, target
                )
                require(
                    bool(common) == bool(right),
                    "common/right support shadows split",
                )
                if not right:
                    empty += 1
                    continue
                nonempty += 1
                total_common += len(common)
                total_right += len(right)
                cell_types[clock, len(common), len(right)] += 1

                lengths = {
                    stop - left
                    for left, stop, _weight in (*common, *right)
                }
                weights = {
                    weight for _left, _stop, weight in (*common, *right)
                }
                expected_weight = copies.W1 if clock == 1 else copies.W2
                require(
                    lengths == {copies.LENGTH}
                    and weights == {expected_weight},
                    "a nonempty cell stopped being an equal-copy bank",
                )
                common_by_left = {piece[0]: piece for piece in common}
                require(
                    len(common_by_left) == len(common),
                    "common left endpoints collided inside a cell",
                )
                image_lefts = []
                cache = semantic_caches[clock]

                for rpiece in right:
                    rleft, rstop, weight = rpiece
                    candidates = tuple(
                        mpiece for mpiece in common
                        if mpiece[1] - mpiece[0] == rstop - rleft
                        and mpiece[2] == weight
                    )
                    require(
                        len(candidates) == len(common),
                        "equal-copy relation stopped being complete bipartite",
                    )
                    candidate_histogram[len(candidates)] += 1
                    relation_edges += len(candidates)

                    offset_rows = tuple(
                        (abs(circular_offset(rleft, mpiece[0])),
                         circular_offset(rleft, mpiece[0]), mpiece)
                        for mpiece in candidates
                    )
                    minimum = min(row[0] for row in offset_rows)
                    nearest = tuple(
                        row for row in offset_rows if row[0] == minimum
                    )
                    require(
                        len(nearest) == 1 and nearest[0][1] == -1
                        and nearest[0][2][0] == rleft + H,
                        "metric nearest neighbor stopped being the +h collar",
                    )
                    m1 = nearest[0][2]
                    nearest_offsets[-1] += 1
                    image_lefts.append(m1[0])
                    m2 = common_by_left.get(rleft + 2 * H)
                    require(m2 is not None, "+2h common collar piece missing")

                    rvalue = semantic_value(rpiece, q_pairs[clock], cache)
                    m1value = semantic_value(m1, q_pairs[clock], cache)
                    m2value = semantic_value(m2, q_pairs[clock], cache)
                    right_status[rvalue] += 1
                    nearest_status[m1value] += 1
                    second_status[m2value] += 1
                    nearest_contingency[rvalue, m1value] += 1
                    require(
                        rvalue != m1value and rvalue == m2value,
                        "one/two half-step semantic law changed",
                    )
                    if rvalue == (0, 0):
                        dead_right.append(
                            ((clock, sigma, target), rpiece[:2])
                        )

                    same_content = []
                    for _distance, offset, mpiece in offset_rows:
                        mvalue = semantic_value(
                            mpiece, q_pairs[clock], cache
                        )
                        same = mvalue == rvalue
                        relation_parity[offset % 2, same] += 1
                        if same:
                            same_content.append(
                                (abs(offset), offset, mpiece)
                            )
                    same_minimum = min(row[0] for row in same_content)
                    same_nearest = tuple(
                        row for row in same_content
                        if row[0] == same_minimum
                    )
                    require(
                        len(same_nearest) == 1
                        and same_nearest[0][1] == -2
                        and same_nearest[0][2] == m2,
                        "nearest same-content mate stopped being +2h",
                    )
                    same_content_offsets[-2] += 1

                require(
                    len(image_lefts) == len(set(image_lefts)),
                    "nearest collar collided inside a labelled cell",
                )

    expected_content = (copies.C, copies.C)
    require(
        (nonempty, empty, total_common, total_right, relation_edges)
        == (193, 374, 63308, 587, 195517),
        "global common/right census changed",
    )
    require(
        cell_types == EXPECTED_CELL_TYPES,
        "common/right cell-type census changed",
    )
    require(
        candidate_histogram == EXPECTED_CANDIDATES,
        "complete-bipartite candidate census changed",
    )
    require(
        relation_parity
        == Counter({(0, True): 97661, (1, False): 97856}),
        "semantic equality stopped being exactly half-step parity",
    )
    require(
        nearest_offsets == Counter({-1: 587})
        and same_content_offsets == Counter({-2: 587}),
        "nearest collar offsets changed",
    )
    require(
        right_status == Counter({expected_content: 573, (0, 0): 14})
        and nearest_status
        == Counter({(0, 0): 573, expected_content: 14})
        and second_status == right_status,
        "collar semantic contingency changed",
    )
    require(
        nearest_contingency == Counter({
            (expected_content, (0, 0)): 573,
            ((0, 0), expected_content): 14,
        }),
        "nearest semantic reversal table changed",
    )
    require(
        tuple(dead_right) == EXPECTED_DEAD_RIGHT,
        "dead right-copy addresses changed",
    )
    require(
        PERIOD_STEPS == 2 * P**5
        and copies.LENGTH * 91 == 6 * H
        and P * copies.SHIFT == 14 * H
        and 14 * H == 7 * SEMANTIC_PERIOD,
        "half-step scale identities changed",
    )

    print("THM-2825 NEAREST HALF-STEP COMMON/RIGHT COLLAR")
    print(
        f"cells=567;nonempty={nonempty};empty={empty};"
        f"common={total_common};right={total_right};edges={relation_edges}"
    )
    print(f"cell_types={tuple(sorted(cell_types.items()))}")
    print(
        "candidate_histogram="
        f"{tuple(sorted(candidate_histogram.items()))}"
    )
    print(
        "relation_parity="
        f"{tuple(sorted(relation_parity.items()))}"
    )
    print(
        f"nearest_offsets={tuple(sorted(nearest_offsets.items()))};"
        f"same_content_offsets={tuple(sorted(same_content_offsets.items()))}"
    )
    print(
        f"right_status={tuple(sorted(right_status.items()))};"
        f"nearest_status={tuple(sorted(nearest_status.items()))};"
        f"second_status={tuple(sorted(second_status.items()))}"
    )
    print(
        "nearest_contingency="
        f"{tuple(sorted(nearest_contingency.items()))}"
    )
    print(f"dead_right={tuple(dead_right)}")
    print(
        f"h={H};semantic_period={SEMANTIC_PERIOD};"
        f"piece_length={copies.LENGTH}=6h/91;"
        f"13SHIFT={P * copies.SHIFT}=14h=7(2h)"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
