#!/usr/bin/env python3
"""Exact referee for THM-2629's fixed-deep affine graph spectrum.

The program rebuilds THM-2616's globally primitive diagonal tensor, retains
one nonzero deep label before the seven-clock unit test, exhausts all 169
affine graphs r=alpha*q+beta, and checks the puncture boundary.
"""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import gcd

import lrc14_cross_time_target_future_diagonal_thm2616 as thm


P = thm.P
Q7 = thm.Q7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fixed_r_unit(values, r):
    """Test one fixed-r seven-clock class in F13[z]/(Phi_7)."""
    require(all(values[ell][r] % thm.GLOBAL_CONTENT == 0 for ell in range(Q7)),
            "fixed-r slice is not divisible by the inherited global content")
    scalar = pow(r, -1, P)
    y = tuple(
        ((values[ell][r] // thm.GLOBAL_CONTENT) * scalar) % P
        for ell in range(Q7)
    )
    reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
    return thm.old.sat.multiplication_determinant_7(reduced) != 0


def main():
    with ProcessPoolExecutor(max_workers=len(thm.SHARDS)) as pool:
        results = list(pool.map(thm.compute_shard, thm.SHARDS))
    require(tuple(result[0] for result in results) == thm.SHARDS,
            "rail shards returned out of order")

    metadata = []
    diagonal = []
    for result in results:
        metadata.extend(result[4])
        diagonal.extend(result[5])
    require(len(metadata) == len(diagonal) == 162, "rail bank changed")

    diagonal_content = 0
    for raw in diagonal:
        for q in range(P):
            for ell in range(Q7):
                for r in range(1, P):
                    value = raw[q][ell][r]
                    if value:
                        diagonal_content = gcd(diagonal_content, value)
    require(diagonal_content == thm.GLOBAL_CONTENT == 26,
            "fixed-deep bank did not inherit the one global primitive content")

    units = [[[False] * P for _ in range(P)] for _ in diagonal]
    rail_hist = Counter()
    row_support_hist = Counter()
    for j, raw in enumerate(diagonal):
        for q in range(P):
            for ell in range(Q7):
                support = tuple(r for r in range(1, P) if raw[q][ell][r] > 0)
                if support:
                    row_support_hist[len(support)] += 1
            for r in range(1, P):
                units[j][q][r] = fixed_r_unit(raw[q], r)
        rail_hist[sum(units[j][q][r]
                      for q in range(P) for r in range(1, P))] += 1

    expected_rail_hist = Counter({
        0: 20, 76: 1, 78: 1, 102: 1, 108: 19, 110: 1,
        120: 16, 122: 7, 126: 1, 130: 2, 131: 35, 132: 58,
    })
    require(row_support_hist == Counter({12: 4_156}),
            "a nonzero labelled row lost full deep support")
    require(rail_hist == expected_rail_hist,
            "fixed-r rail unit histogram changed")

    by_cell = {}
    for j, (s, ell, _t) in enumerate(metadata):
        by_cell.setdefault((s, ell), []).append(j)
    require(len(by_cell) == 84, "base-cell bank changed")

    cell_unit = {
        cell: {
            (q, r)
            for q in range(P)
            for r in range(1, P)
            if any(units[j][q][r] for j in rails)
        }
        for cell, rails in sorted(by_cell.items())
    }
    pair_count_hist = Counter(len(relation) for relation in cell_unit.values())
    require(pair_count_hist == Counter({108: 4, 120: 2, 122: 1,
                                       130: 1, 131: 6, 132: 70}),
            "base-cell fixed-r unit histogram changed")

    common_pairs = set.intersection(*(set(v) for v in cell_unit.values()))
    common_r_by_q = {
        q: tuple(r for r in range(1, P) if (q, r) in common_pairs)
        for q in range(P)
    }
    expected_common = {
        0: (),
        1: (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
        2: (),
        3: (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
        4: (7, 8, 9, 10, 11, 12),
        5: (7, 8, 9, 10, 11, 12),
        6: (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
        7: (1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12),
        8: (1, 2, 3, 4, 5, 6),
        9: (1, 2, 3, 4, 5, 6),
        10: (1, 2),
        11: (),
        12: (),
    }
    require(len(common_pairs) == 72 and common_r_by_q == expected_common,
            "common fixed-r pair atlas changed")

    affine = []
    for alpha in range(P):
        for beta in range(P):
            qsets = []
            for relation in cell_unit.values():
                qsets.append(tuple(
                    q for q in range(P)
                    if (alpha * q + beta) % P != 0
                    and (q, (alpha * q + beta) % P) in relation
                ))
            common = tuple(sorted(set.intersection(*(set(v) for v in qsets))))
            affine.append((min(map(len, qsets)), len(common), alpha, beta,
                           common, Counter(map(len, qsets))))
    affine.sort(reverse=True)

    min_ten_maps = tuple(
        (alpha, beta, common, tuple(sorted(hist.items())))
        for minimum, _common_count, alpha, beta, common, hist in affine
        if minimum == 10
    )
    expected_min_ten = (
        (12, 12, (1, 3, 4, 5, 6, 7, 8, 9, 10), ((10, 6), (11, 78))),
        (5, 0, (1, 3, 4, 5, 6, 7, 8, 9), ((10, 7), (11, 77))),
        (2, 0, (1, 3, 4, 5, 6, 7, 8, 9), ((10, 7), (11, 77))),
        (12, 0, (3, 4, 5, 6, 7, 8, 9), ((10, 9), (11, 75))),
    )
    require(min_ten_maps == expected_min_ten,
            "sharp affine-map bank changed")
    require(affine[0][0:5] ==
            (10, 9, 12, 12, (1, 3, 4, 5, 6, 7, 8, 9, 10)),
            "r=-q-1 is no longer the unique lexicographic optimum")

    best_qsets = {
        cell: tuple(q for q in range(P)
                    if (-q - 1) % P != 0
                    and (q, (-q - 1) % P) in relation)
        for cell, relation in cell_unit.items()
    }
    best_pattern_hist = Counter(best_qsets.values())
    expected_best_patterns = Counter({
        tuple(range(1, 12)): 78,
        tuple(range(1, 11)): 4,
        (1, 3, 4, 5, 6, 7, 8, 9, 10, 11): 2,
    })
    require(best_pattern_hist == expected_best_patterns,
            "sharp graph pattern histogram changed")
    best_exceptional = tuple(
        (cell, values) for cell, values in best_qsets.items()
        if values != tuple(range(1, 12))
    )
    expected_exceptional = (
        ((2, 2), tuple(range(1, 11))),
        ((2, 3), tuple(range(1, 11))),
        ((6, 5), (1, 3, 4, 5, 6, 7, 8, 9, 10, 11)),
        ((7, 2), tuple(range(1, 11))),
        ((7, 3), tuple(range(1, 11))),
        ((11, 5), (1, 3, 4, 5, 6, 7, 8, 9, 10, 11)),
    )
    require(best_exceptional == expected_exceptional,
            "sharp graph exceptional atlas changed")

    positive_by_cell = {
        cell: tuple(
            q for q in range(P)
            if (-q - 1) % P != 0
            and any(diagonal[j][q][ell][(-q - 1) % P] > 0
                    for j in rails for ell in range(Q7))
        )
        for cell, rails in sorted(by_cell.items())
    }
    positive_pattern_hist = Counter(positive_by_cell.values())
    require(positive_pattern_hist == Counter({
        tuple(range(1, 12)): 81,
        tuple(range(1, 11)): 2,
        (1, 3, 4, 5, 6, 7, 8, 9, 10, 11): 1,
    }), "sharp graph positivity atlas changed")

    def affine_record(alpha, beta):
        return next(row for row in affine if row[2:4] == (alpha, beta))

    identity = affine_record(1, 0)
    inverse_two = affine_record(7, 0)
    constant_best = max(row for row in affine if row[2] == 0)
    require(identity == (8, 4, 1, 0, (1, 3, 6, 7),
                         Counter({11: 75, 10: 5, 8: 4})),
            "identity graph control changed")
    require(inverse_two == (9, 6, 7, 0, (1, 3, 5, 6, 7, 8),
                            Counter({11: 77, 9: 4, 10: 3})),
            "inverse-two graph control changed")
    require(constant_best == (8, 7, 0, 2, (1, 3, 6, 7, 8, 9, 10),
                              Counter({11: 77, 10: 5, 8: 2})),
            "constant-r control changed")

    # The best bijective affine graph has its unique r=0 value at q=12.
    # Hence it avoids the deep puncture exactly on the existing q=1,...,11
    # future carrier, but no bijective affine graph can avoid r=0 on all F13.
    require((-12 - 1) % P == 0 and
            {(-q - 1) % P for q in range(1, 12)} == set(range(1, 12)),
            "puncture alignment changed")
    require(all({(alpha * q + beta) % P for q in range(P)} == set(range(P))
                for alpha in range(1, P) for beta in range(P)),
            "an affine bijection avoided the deep zero sheet")

    print("THM-2629 exact fixed-deep affine-graph controls")
    print(f"inherited_diagonal_global_content={diagonal_content} divisibility=PASS")
    print(f"nonzero_labelled_row_deep_support_hist={tuple(sorted(row_support_hist.items()))}")
    print(f"rail_fixed_r_unit_count_hist={tuple(sorted(rail_hist.items()))}")
    print(f"cell_qr_unit_count_hist={tuple(sorted(pair_count_hist.items()))}")
    print(f"common_qr_pairs={len(common_pairs)} common_r_by_q={common_r_by_q}")
    print(f"affine_max_min=10 min_ten_maps={min_ten_maps}")
    print("sharp_graph=(alpha,beta)=(-1,-1) min_cell=10 common_count=9 common=(1,3,4,5,6,7,8,9,10)")
    print(f"sharp_graph_pattern_hist={tuple(sorted(best_pattern_hist.items()))}")
    print(f"sharp_graph_exceptional={best_exceptional}")
    print(f"sharp_graph_positive_pattern_hist={tuple(sorted(positive_pattern_hist.items()))}")
    print("controls=identity:(8,4) inverse_two:(9,6) best_constant_r2:(8,7)")
    print("verdict=PASS: fixed-deep unit refinement has a sharp opposite-affine graph and an unavoidable full-torsor puncture")
    print("SCOPE: coefficient/unit atlas only; rail choices vary by cell and no chronological, semantic, or ordered gluing is inferred")


if __name__ == "__main__":
    main()
