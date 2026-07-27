#!/usr/bin/env python3
"""Exact companion for THM-2609: target-section saturation and state loss."""

from collections import Counter
from itertools import product
from math import gcd
import hashlib

import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = old.P
Q7 = old.Q7
T = old.T
R = old.R
H = 6
GLOBAL_CONTENT = 4_244_240
GRID = R * T

base = old.base
host = old.host
rail = old.rail
sat = old.sat


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rebuild_unit_qsets():
    """Rebuild THM-2600's two-rail bank and retain unit q sections by cell."""
    module = old.load_present_module()
    require(module.W == base.W, "canonical row changed")

    # Only the physically used future digit h=6 is needed here.
    word = module.build_word_Ta()
    word_prefix = []
    for ell in range(Q7):
        q_ell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        digit = sat.intersect_interval(
            q_ell, H * T // P, (H + 1) * T // P
        )
        prefix = module.make_prefix(digit)
        require(prefix[2][-1] > 0, "future digit six disappeared")
        word_prefix.append(prefix)

    e4 = base.build_set(base.PAT_E, base.ZELL)
    qb = base.build_set(host.PAT_QB, base.ZELL)
    ust, uv, vst, vv = rail.packet_profiles(e4, qb)
    _, _, k_tensor = rail.exact_tensor(e4, qb)
    owner = base.clock_cells(P, Q7, T, P * P)
    deep = rail.deep_cells()
    arrival = [(H * (T // P), (H + 1) * (T // P))]

    rail_bank = []
    for s in range(1, P):
        rvst, rvv = base.rotate_profile(vst, vv, s * (T // P), T)
        ps, pv, _, _ = base.product_cum(ust, uv, rvst, rvv, T)
        for ell in range(Q7):
            for t in (12, 0):
                if k_tensor[s][H][t][ell] == 0:
                    continue
                cell = old.intersect_sorted(
                    old.intersect_sorted(owner[ell], deep[t]), arrival
                )
                pieces = old.profile_on_intervals(cell, ps, pv)
                numerator = P * sum(
                    weight * (right - left)
                    for left, right, weight in pieces
                )
                require(
                    numerator == k_tensor[s][H][t][ell] > 0,
                    "middle-rail mass mismatch",
                )
                rail_bank.append((s, ell, t, pieces))
    require(len(rail_bank) == 162, "middle-rail census changed")

    present = {}
    starts = {}
    for ell5 in range(Q7):
        for s5 in range(P):
            intervals = module.build_F(ell5, s5)
            present[ell5, s5] = intervals
            starts[ell5, s5] = [left for left, _ in intervals]

    joint = [[[[0] * P for _ in range(Q7)] for _ in range(P)]
             for _ in rail_bank]
    content = 0
    positive_entries = 0
    for j, (_, _, _, pieces) in enumerate(rail_bank):
        for q in range(P):
            s5 = (-q) % P
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, s5], starts[ell5, s5]
                )
                for r in range(1, P):
                    probed = old.intersect_weighted_comb(
                        overlap, module.C3, 182, 14 * r - 13, 14 * r + 13
                    )
                    value = old.delayed_weighted_numerator(
                        probed, word_prefix[ell5]
                    )
                    joint[j][q][ell5][r] = value
                    if value:
                        positive_entries += 1
                        content = gcd(content, value)

    require(content == GLOBAL_CONTENT, "global primitive content changed")
    require(positive_entries == 61_248, "fine positive-entry census changed")

    unit = [[False] * P for _ in rail_bank]
    unit_slices = 0
    for j in range(len(rail_bank)):
        for q in range(P):
            y = tuple(
                sum(
                    (joint[j][q][ell][r] // content) * pow(r, -1, P)
                    for r in range(1, P)
                ) % P
                for ell in range(Q7)
            )
            reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
            determinant = sat.multiplication_determinant_7(reduced)
            if determinant:
                require(
                    any(joint[j][q][ell][r] > 0
                        for ell in range(Q7) for r in range(1, P)),
                    "unit slice has no positive physical atom",
                )
                unit[j][q] = True
                unit_slices += 1
    require(unit_slices == 1_740, "unit-slice census changed")

    by_cell = {}
    for j, (s, ell, t, _) in enumerate(rail_bank):
        by_cell.setdefault((s, ell), []).append((j, t))
    require(len(by_cell) == 84, "base-cell census changed")

    qsets = {}
    for cell, edges in sorted(by_cell.items()):
        qsets[cell] = tuple(
            q for q in range(P) if any(unit[j][q] for j, _ in edges)
        )
    return qsets, len(rail_bank), positive_entries, content


def difference_set(values):
    return {(right - left) % P for left in values for right in values}


def main():
    qsets, rail_count, positive_entries, content = rebuild_unit_qsets()

    size_hist = Counter(map(len, qsets.values()))
    pattern_hist = Counter(qsets.values())
    require(size_hist == Counter({13: 74, 12: 8, 11: 2}),
            "unit q-set size histogram changed")
    require(all(difference_set(values) == set(range(P))
                for values in qsets.values()),
            "one unit q-set lost full difference support")
    require(all(0 in values and H in values for values in qsets.values()),
            "uniform q=0 or q=6 section disappeared")

    translation_stabilizers = {
        cell: tuple(
            shift for shift in range(P)
            if {(q + shift) % P for q in values} == set(values)
        )
        for cell, values in qsets.items()
    }
    stabilizer_hist = Counter(map(len, translation_stabilizers.values()))
    require(stabilizer_hist == Counter({13: 74, 1: 10}),
            "external-section translation stabilizers changed")
    require(all(
        stabilizer == tuple(range(P)) if len(qsets[cell]) == P
        else stabilizer == (0,)
        for cell, stabilizer in translation_stabilizers.items()
    ), "a deficient q-set acquired a nontrivial C13 stabilizer")

    witnesses = []
    for cell, values in sorted(qsets.items()):
        row = tuple(
            min((left, right) for left in values for right in values
                if (right - left) % P == displacement)
            for displacement in range(P)
        )
        witnesses.append((cell, row))
    witness_digest = hashlib.sha256(repr(witnesses).encode("ascii")).hexdigest()
    require(
        witness_digest
        == "133003c6f7146394d2396a37924885882027118b967171155437d7309005c6e6",
        "endpoint witness table changed",
    )

    # A seven-edge clock loop has eight vertex events.  The two endpoints
    # revisit the same base cell and pay any requested external correction;
    # q=0 supplies the six intermediate clocks.
    itineraries = 0
    for s, ell0, marker in product(range(1, P), range(Q7), range(1, P)):
        endpoint = qsets[s, ell0]
        desired = (-Q7 * marker) % P
        q0, q7 = min(
            (left, right) for left in endpoint for right in endpoint
            if (right - left) % P == desired
        )
        labels = (q0,) + (0,) * (Q7 - 1) + (q7,)
        require(len(labels) == 8, "wrong number of itinerary vertices")
        for step, q in enumerate(labels):
            require(q in qsets[s, (ell0 + step) % Q7],
                    "itinerary left the unit-section atlas")
        require((labels[-1] - labels[0]) % P == desired,
                "external endpoint correction changed")
        itineraries += 1
    require(itineraries == 12 * 7 * 12, "itinerary census changed")

    # Every positive exact entry is a union of open atoms on the 1/(R*T)
    # grid.  An interval longer than 2/13^d contains a closed depth-d
    # base-13 cylinder.  Depth 20 is the first universal bound obtained from
    # this grid invoice, and eight blocks concatenate to depth 160.
    require(R == 13 ** 6 and T == 297_836_897_838_480,
            "clock or base grid changed")
    depth = min(d for d in range(1, 100) if 13 ** d > 2 * GRID)
    require(depth == 20, "universal cylinder depth changed")
    require(13 ** 19 < 2 * GRID < 13 ** 20,
            "depth-20 cylinder invoice changed")
    total_depth = 8 * depth
    require(total_depth == 160 and 13 ** total_depth > 0,
            "itinerary cylinder invoice changed")

    # q is an external section label; h is the physical future digit and is
    # identically six throughout this atlas.  The common-event diagonal has
    # one state only, so it supports no nonzero endpoint correction.
    state_supports = {
        cell: {(q, H) for q in values} for cell, values in qsets.items()
    }
    diagonal = {
        cell: {(q, h) for q, h in states if q == h}
        for cell, states in state_supports.items()
    }
    require(all(states == {(H, H)} for states in diagonal.values()),
            "common-event q=h state boundary changed")
    diagonal_differences = {
        (right[0] - left[0]) % P
        for states in diagonal.values() for left in states for right in states
    }
    require(diagonal_differences == {0},
            "common-event diagonal acquired a nonzero correction")

    print("THM-2609 exact external-section itinerary controls")
    print(f"inherited_rails={rail_count} fine_positive_entries={positive_entries} "
          f"global_content={content}")
    print(f"unit_qset_size_hist={sorted(size_hist.items())} "
          f"full_difference_cells={len(qsets)}/84")
    print(f"unit_qset_pattern_hist={sorted(pattern_hist.items())}")
    print(f"external_translation_stabilizer_hist="
          f"{sorted(stabilizer_hist.items())}")
    print(f"unit_endpoint_witness_digest={witness_digest}")
    print("uniform_sections=q0:84/84,q6:84/84")
    print(f"external_correction_itineraries={itineraries} "
          f"events_per_itinerary=8")
    print(f"grid_denominator={GRID} depth_test=13^19<2RT<13^20 "
          f"itinerary_measure=13^-{total_depth}")
    print("common_event_q_equals_h={(6,6)} diagonal_differences=(0,)")
    print("verdict=PASS: differences and temporal itineraries saturate; "
          "the missing datum is a common-event external-q/physical-h state map")


if __name__ == "__main__":
    main()
