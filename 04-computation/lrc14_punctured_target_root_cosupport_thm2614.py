#!/usr/bin/env python3
"""Exact controls for THM-2614's punctured target/root cosupport boundary.

The program rebuilds the complete THM-2600 constant-six two-rail bank.  It
keeps a target section q only when the globally primitive aggregate slice has
unit septimal Bockstein, then asks which deep probes r occur positively.

The answer is support-level factorization, not a charged-sheet statement:
every active fine (rail,q,future-clock) fibre contains all twelve nonzero r,
while r=0 is absent identically.  In particular no unit Bockstein is assigned
to an individual positive Boolean component.
"""

from collections import Counter
from math import gcd

import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = old.P
Q7 = old.Q7
T = old.T
H = 6
GLOBAL_CONTENT = 4_244_240


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rebuild_bank():
    """Return the exact two-rail fine bank and its aggregate unit flags."""
    module = old.load_present_module()
    require(module.W == old.base.W, "canonical THM-2600 row changed")

    word = module.build_word_Ta()
    word_prefix = []
    for ell in range(Q7):
        q_ell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        digit = old.sat.intersect_interval(
            q_ell, H * T // P, (H + 1) * T // P
        )
        prefix = module.make_prefix(digit)
        require(prefix[2][-1] > 0, "future digit six disappeared")
        word_prefix.append(prefix)

    e4 = old.base.build_set(old.base.PAT_E, old.base.ZELL)
    qb = old.base.build_set(old.host.PAT_QB, old.base.ZELL)
    ust, uv, vst, vv = old.rail.packet_profiles(e4, qb)
    _, _, k_tensor = old.rail.exact_tensor(e4, qb)
    owner = old.base.clock_cells(P, Q7, T, P * P)
    deep = old.rail.deep_cells()
    arrival = [(H * (T // P), (H + 1) * (T // P))]

    rails = []
    for s in range(1, P):
        rvst, rvv = old.base.rotate_profile(vst, vv, s * (T // P), T)
        ps, pv, _, _ = old.base.product_cum(ust, uv, rvst, rvv, T)
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
                rails.append((s, ell, t, pieces))
    require(len(rails) == 162, "middle-rail census changed")

    present = {}
    starts = {}
    for ell5 in range(Q7):
        for shift in range(P):
            intervals = module.build_F(ell5, shift)
            present[ell5, shift] = intervals
            starts[ell5, shift] = [left for left, _ in intervals]

    joint = [[[[0] * P for _ in range(Q7)] for _ in range(P)]
             for _ in rails]
    content = 0
    positive_entries = 0
    for j, (_, _, _, pieces) in enumerate(rails):
        for q in range(P):
            shift = (-q) % P
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, shift], starts[ell5, shift]
                )
                for r in range(1, P):
                    probed = old.intersect_weighted_comb(
                        overlap, module.C3, 182,
                        14 * r - 13, 14 * r + 13,
                    )
                    value = old.delayed_weighted_numerator(
                        probed, word_prefix[ell5]
                    )
                    joint[j][q][ell5][r] = value
                    if value:
                        positive_entries += 1
                        content = gcd(content, value)
    require(content == GLOBAL_CONTENT, "global primitive content changed")
    require(positive_entries == 61_248,
            "fine positive-entry census changed")

    unit = [[False] * P for _ in rails]
    unit_slices = 0
    for j in range(len(rails)):
        for q in range(P):
            y = tuple(
                sum(
                    (joint[j][q][ell][r] // content) * pow(r, -1, P)
                    for r in range(1, P)
                ) % P
                for ell in range(Q7)
            )
            reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
            determinant = old.sat.multiplication_determinant_7(reduced)
            if determinant:
                require(
                    any(joint[j][q][ell][r] > 0
                        for ell in range(Q7) for r in range(1, P)),
                    "unit aggregate slice has no positive fine entry",
                )
                unit[j][q] = True
                unit_slices += 1
    require(unit_slices == 1_740, "unit-slice census changed")

    by_cell = {}
    for j, (s, ell, t, _) in enumerate(rails):
        by_cell.setdefault((s, ell), []).append((j, t))
    require(len(by_cell) == 84, "base-cell census changed")
    return rails, joint, unit, by_cell, positive_entries, unit_slices


def translate_relation(relation, dq, dr):
    return {((q + dq) % P, (r + dr) % P) for q, r in relation}


def main():
    rails, joint, unit, by_cell, positive_entries, unit_slices = rebuild_bank()

    # Root support factors already at the finest retained resolution.
    fine_root_support_hist = Counter()
    for j in range(len(rails)):
        for q in range(P):
            for ell5 in range(Q7):
                support = tuple(
                    r for r in range(P) if joint[j][q][ell5][r] > 0
                )
                require(support in ((), tuple(range(1, P))),
                        "a fine fibre acquired partial nonzero-root support")
                fine_root_support_hist[len(support)] += 1
    require(fine_root_support_hist == Counter({0: 9_638, 12: 5_104}),
            "fine root-support histogram changed")

    relations = {}
    qsets = {}
    future_relations = {}
    for cell, edges in sorted(by_cell.items()):
        relation = set()
        future = []
        for ell5 in range(Q7):
            rel5 = {
                (q, r)
                for j, _ in edges
                for q in range(P)
                for r in range(P)
                if unit[j][q] and joint[j][q][ell5][r] > 0
            }
            active_q = {q for q, _ in rel5}
            require(rel5 == {(q, r) for q in active_q
                                      for r in range(1, P)},
                    "a future-clock cosupport failed to factor")
            future.append(rel5)
            relation |= rel5
        active_q = {q for q, _ in relation}
        require(relation == {(q, r) for q in active_q
                                      for r in range(1, P)},
                "a base-cell cosupport failed to equal Q times F13-star")
        relations[cell] = relation
        qsets[cell] = tuple(sorted(active_q))
        future_relations[cell] = tuple(future)

    qset_size_hist = Counter(map(len, qsets.values()))
    relation_size_hist = Counter(map(len, relations.values()))
    require(qset_size_hist == Counter({13: 74, 12: 8, 11: 2}),
            "unit target-section size histogram changed")
    require(relation_size_hist == Counter({156: 74, 144: 8, 132: 2}),
            "same-event relation size histogram changed")
    require(all(0 in values for values in qsets.values()),
            "uniform q=0 section disappeared")

    deficient_atlas = tuple(
        (cell, tuple(q for q in range(P) if q not in values))
        for cell, values in sorted(qsets.items()) if len(values) < P
    )

    future_active_q_hist = Counter()
    future_difference_hist = Counter()
    for future in future_relations.values():
        for relation in future:
            active_q = {q for q, _ in relation}
            future_active_q_hist[len(active_q)] += 1
            differences = {(r - q) % P for q, r in relation}
            require(not relation or differences == set(range(P)),
                    "a nonempty future-clock relation lost a difference")
            future_difference_hist[len(differences)] += 1
    require(future_active_q_hist == Counter({
        0: 204, 9: 32, 10: 30, 11: 226, 12: 76, 13: 20,
    }), "future-clock active-q histogram changed")
    require(future_difference_hist == Counter({0: 204, 13: 384}),
            "future-clock difference histogram changed")

    difference_hist = Counter()
    diagonal_count_hist = Counter()
    diagonal_stabilizer_hist = Counter()
    affine_coverage_hist = Counter()
    punctured_gauge_count_hist = Counter()
    completion_profile_hist = Counter()
    contained_total_fixed = 0
    contained_total_semilinear = 0

    for relation in relations.values():
        differences = {(r - q) % P for q, r in relation}
        difference_hist[len(differences)] += 1
        diagonal_count_hist[sum((q, q) in relation for q in range(P))] += 1
        stabilizer = tuple(
            a for a in range(P)
            if translate_relation(relation, a, a) == relation
        )
        require(stabilizer == (0,),
                "punctured relation acquired a diagonal deck symmetry")
        diagonal_stabilizer_hist[len(stabilizer)] += 1

        coverages = tuple(
            sum(((c + r) % P, r) in relation for r in range(1, P))
            for c in range(P)
        )
        affine_coverage_hist.update(coverages)
        punctured_gauge_count_hist[sum(value == P - 1
                                       for value in coverages)] += 1
        defects = tuple(P - value for value in coverages)
        completion_profile_hist[(min(defects), defects.count(min(defects)))] += 1

        for c in range(P):
            graph = {((c + r) % P, r) for r in range(P)}
            contained_total_fixed += int(graph <= relation)
        for kappa in range(1, P):
            for c in range(P):
                graph = {((kappa * r + c) % P, r) for r in range(P)}
                contained_total_semilinear += int(graph <= relation)

    require(difference_hist == Counter({13: 84}),
            "same-event relation lost full differences")
    require(diagonal_count_hist == Counter({12: 74, 11: 8, 10: 2}),
            "literal q=r diagonal census changed")
    require(diagonal_stabilizer_hist == Counter({1: 84}),
            "diagonal stabilizer census changed")
    require(affine_coverage_hist == Counter({12: 970, 11: 100, 10: 22}),
            "punctured affine-graph coverage changed")
    require(punctured_gauge_count_hist == Counter({13: 74, 1: 8, 0: 2}),
            "full punctured-gauge census changed")
    require(completion_profile_hist == Counter({
        (1, 13): 74, (1, 1): 8, (2, 2): 2,
    }), "affine-graph completion profile changed")
    require(contained_total_fixed == 0 and contained_total_semilinear == 0,
            "the missing r=0 column unexpectedly admitted a total graph")

    print("THM-2614 exact punctured target/root cosupport controls")
    print(f"inherited_rails={len(rails)} fine_positive_entries={positive_entries} "
          f"unit_aggregate_slices={unit_slices} global_content={GLOBAL_CONTENT}")
    print(f"fine_root_support_hist={sorted(fine_root_support_hist.items())}")
    print(f"unit_qset_size_hist={sorted(qset_size_hist.items())} "
          f"relation_size_hist={sorted(relation_size_hist.items())}")
    print(f"deficient_qset_atlas={deficient_atlas}")
    print(f"future_clock_active_q_hist={sorted(future_active_q_hist.items())} "
          f"difference_hist={sorted(future_difference_hist.items())}")
    print(f"base_difference_size_hist={sorted(difference_hist.items())} "
          f"literal_diagonal_count_hist={sorted(diagonal_count_hist.items())}")
    print(f"diagonal_stabilizer_size_hist="
          f"{sorted(diagonal_stabilizer_hist.items())}")
    print(f"punctured_affine_coverage_hist="
          f"{sorted(affine_coverage_hist.items())}")
    print(f"full_punctured_gauges_per_cell_hist="
          f"{sorted(punctured_gauge_count_hist.items())}")
    print(f"minimum_total_graph_completion_hist="
          f"{sorted(completion_profile_hist.items())}")
    print("contained_total_fixed_graphs=0 contained_total_semilinear_graphs=0")
    print("verdict=PASS: co-support is Q times F13-star; full differences do not "
          "supply a principal C13 connector")


if __name__ == "__main__":
    main()
