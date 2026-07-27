#!/usr/bin/env python3
"""Exact controls for THM-2616's cross-time target/future diagonal.

The full carrier retains every delayed physical digit h before imposing the
numerical diagonal h=q.  Its global content is computed over every rail,
target section, future clock, nonzero deep probe, and future digit.  The
diagonal Bockstein is then normalized by that one full-carrier content.
"""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from math import gcd

import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = old.P
Q7 = old.Q7
T = old.T
H_ARRIVAL = 6
SHARDS = ((0, 41), (41, 82), (82, 123), (123, 162))
GLOBAL_CONTENT = 26


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def build_carrier_data():
    """Rebuild the proved THM-2600 rails and all delayed digit prefixes."""
    module = old.load_present_module()
    require(module.W == old.base.W, "canonical typed row changed")

    word = module.build_word_Ta()
    prefixes = [[None] * P for _ in range(Q7)]
    whole_prefixes = [None] * Q7
    digit_masses = [[0] * P for _ in range(Q7)]
    for ell in range(Q7):
        qell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        whole_prefixes[ell] = module.make_prefix(qell)
        for h in range(P):
            digit = old.sat.intersect_interval(
                qell, h * T // P, (h + 1) * T // P
            )
            prefixes[ell][h] = module.make_prefix(digit)
            digit_masses[ell][h] = prefixes[ell][h][2][-1]
        require(
            sum(digit_masses[ell]) == whole_prefixes[ell][2][-1],
            "future digits do not partition the delayed word",
        )

    e4 = old.base.build_set(old.base.PAT_E, old.base.ZELL)
    qb = old.base.build_set(old.host.PAT_QB, old.base.ZELL)
    ust, uv, vst, vv = old.rail.packet_profiles(e4, qb)
    _, _, k_tensor = old.rail.exact_tensor(e4, qb)
    owner = old.base.clock_cells(P, Q7, T, P * P)
    deep = old.rail.deep_cells()
    arrival = [(H_ARRIVAL * T // P, (H_ARRIVAL + 1) * T // P)]

    rails = []
    for s in range(1, P):
        rvst, rvv = old.base.rotate_profile(vst, vv, s * T // P, T)
        ps, pv, _, _ = old.base.product_cum(ust, uv, rvst, rvv, T)
        for ell in range(Q7):
            for t in (12, 0):
                if k_tensor[s][H_ARRIVAL][t][ell] == 0:
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
                    numerator == k_tensor[s][H_ARRIVAL][t][ell] > 0,
                    "middle-rail mass mismatch",
                )
                rails.append((s, ell, t, pieces))
    require(len(rails) == 162, "middle-rail census changed")

    present = {}
    starts = {}
    for ell in range(Q7):
        for shift in range(P):
            intervals = module.build_F(ell, shift)
            present[ell, shift] = intervals
            starts[ell, shift] = [left for left, _ in intervals]
    return module, prefixes, whole_prefixes, digit_masses, rails, present, starts


def digit_phi_vector(x, prefixes, lengths, cache):
    """Return all thirteen digit-prefix values at one exact endpoint."""
    if x in cache:
        return cache[x]
    h = min(P - 1, P * x // T)
    values = tuple(
        lengths[j]
        if j < h
        else old.phi_at(x, *prefixes[j])
        if j == h
        else 0
        for j in range(P)
    )
    cache[x] = values
    return values


def delayed_all_digits(pieces, prefixes, cache):
    """Compute all digit-restricted delayed numerators in one exact sweep."""
    lengths = tuple(prefix[2][-1] for prefix in prefixes)
    weighted_len = 0
    acc_r = 0
    acc = [0] * P
    rred = old.R % T
    for left, right, weight in pieces:
        rleft = left * rred % T
        rright = right * rred % T
        weighted_len += weight * (right - left)
        acc_r += weight * (rright - rleft)
        vleft = digit_phi_vector(rleft, prefixes, lengths, cache)
        vright = digit_phi_vector(rright, prefixes, lengths, cache)
        for h in range(P):
            acc[h] += weight * (vright[h] - vleft[h])
    floor_numerator = old.R * weighted_len - acc_r
    require(floor_numerator % T == 0, "weighted floor count is not integral")
    quotient = floor_numerator // T
    result = tuple(lengths[h] * quotient + acc[h] for h in range(P))
    require(min(result) >= 0, "negative delayed overlap")
    return result


def compute_shard(bounds):
    """Compute one deterministic rail shard of the enlarged tensor."""
    start, stop = bounds
    (module, prefixes, whole_prefixes, digit_masses, rails,
     present, present_starts) = build_carrier_data()
    require(0 <= start < stop <= len(rails), "invalid rail shard")

    content = 0
    full_positive = 0
    diagonal_positive = 0
    scalar_controls = 0
    vector_cache = [dict() for _ in range(Q7)]
    diagonal = []
    metadata = []

    for j in range(start, stop):
        s, ell, t, pieces = rails[j]
        metadata.append((s, ell, t))
        rail_diagonal = [[[0] * P for _ in range(Q7)] for _ in range(P)]
        for q in range(P):
            shift = (-q) % P
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, shift], present_starts[ell5, shift]
                )
                for r in range(1, P):
                    probed = old.intersect_weighted_comb(
                        overlap, module.C3, 182,
                        14 * r - 13, 14 * r + 13,
                    )
                    values = delayed_all_digits(
                        probed, prefixes[ell5], vector_cache[ell5]
                    )
                    # The thirteen digit pieces are an exact partition of the
                    # unconditioned delayed word on every fine fibre.
                    require(
                        sum(values) == old.delayed_weighted_numerator(
                            probed, whole_prefixes[ell5]
                        ),
                        "digit-vector sweep lost the whole-word partition",
                    )
                    if scalar_controls < 16:
                        expected = tuple(
                            old.delayed_weighted_numerator(
                                probed, prefixes[ell5][h]
                            )
                            for h in range(P)
                        )
                        require(values == expected,
                                "vector sweep disagrees with scalar sweep")
                        scalar_controls += 1
                    for value in values:
                        if value:
                            full_positive += 1
                            content = gcd(content, value)
                    value = values[q]
                    rail_diagonal[q][ell5][r] = value
                    diagonal_positive += int(value > 0)
        diagonal.append(tuple(tuple(tuple(row) for row in qrows)
                              for qrows in rail_diagonal))

    require(content > 0 and scalar_controls == 16,
            "shard did not exercise the exact sweep")
    return (
        bounds, content, full_positive, diagonal_positive,
        tuple(metadata), tuple(diagonal), tuple(tuple(row) for row in digit_masses),
        sum(map(len, vector_cache)), scalar_controls,
    )


def translation_stabilizer(values):
    values = set(values)
    return tuple(
        a for a in range(P)
        if {(x + a) % P for x in values} == values
    )


def main():
    with ProcessPoolExecutor(max_workers=len(SHARDS)) as pool:
        results = list(pool.map(compute_shard, SHARDS))
    require(tuple(result[0] for result in results) == SHARDS,
            "rail shards returned out of order")

    full_content = 0
    full_positive = 0
    diagonal_positive = 0
    metadata = []
    diagonal = []
    digit_masses = results[0][6]
    cache_size = 0
    scalar_controls = 0
    shard_summary = []
    for (bounds, content, positives, diag_positives, meta, raw,
         masses, cache, controls) in results:
        require(masses == digit_masses, "workers disagree on digit masses")
        full_content = gcd(full_content, content)
        full_positive += positives
        diagonal_positive += diag_positives
        metadata.extend(meta)
        diagonal.extend(raw)
        cache_size += cache
        scalar_controls += controls
        shard_summary.append((bounds[0], bounds[1], positives, content))

    require(len(metadata) == len(diagonal) == 162,
            "assembled rail bank changed")
    require(full_content == GLOBAL_CONTENT,
            "enlarged full-carrier global content changed")
    require(full_positive == 649_968,
            "enlarged full-carrier positive census changed")
    require(diagonal_positive == 49_872,
            "cross-time diagonal fine census changed")

    diagonal_content = 0
    positive_pairs = 0
    units = [[False] * P for _ in diagonal]
    unit_slices = 0
    for j in range(len(diagonal)):
        for q in range(P):
            values = diagonal[j][q]
            positive_pairs += int(any(values[ell][r] > 0
                                      for ell in range(Q7)
                                      for r in range(1, P)))
            for ell in range(Q7):
                for r in range(1, P):
                    value = values[ell][r]
                    if value:
                        diagonal_content = gcd(diagonal_content, value)
            y = tuple(
                sum(
                    (values[ell][r] // full_content) * pow(r, -1, P)
                    for r in range(1, P)
                ) % P
                for ell in range(Q7)
            )
            reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
            determinant = old.sat.multiplication_determinant_7(reduced)
            if determinant:
                require(any(values[ell][r] > 0
                            for ell in range(Q7) for r in range(1, P)),
                        "unit diagonal slice has no positive fine component")
                units[j][q] = True
                unit_slices += 1
    require(diagonal_content == GLOBAL_CONTENT,
            "diagonal subset did not witness the full global content")
    require(positive_pairs == 1_703 and unit_slices == 1_483,
            "diagonal pair/unit census changed")

    by_cell = {}
    for j, (s, ell, t) in enumerate(metadata):
        by_cell.setdefault((s, ell), []).append((j, t))
    require(len(by_cell) == 84, "base-cell census changed")
    qsets = {
        cell: tuple(q for q in range(P)
                    if any(units[j][q] for j, _ in edges))
        for cell, edges in sorted(by_cell.items())
    }
    size_hist = Counter(map(len, qsets.values()))
    pattern_hist = Counter(qsets.values())
    require(size_hist == Counter({11: 76, 10: 8}),
            "diagonal unit-section size census changed")
    expected_patterns = Counter({
        tuple(range(1, 11)): 4,
        tuple(range(1, 12)): 76,
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 11): 2,
        (1, 3, 4, 5, 6, 7, 8, 9, 10, 11): 2,
    })
    require(pattern_hist == expected_patterns,
            "diagonal unit-section pattern census changed")

    differences = {
        cell: {(b - a) % P for a in values for b in values}
        for cell, values in qsets.items()
    }
    require(all(value == set(range(P)) for value in differences.values()),
            "a diagonal cell lost an additive difference")
    require(all(translation_stabilizer(values) == (0,)
                for values in qsets.values()),
            "a proper diagonal section set acquired a deck action")

    common = tuple(sorted(set.intersection(*(set(v) for v in qsets.values()))))
    require(common == (1, 3, 4, 5, 6, 7, 8, 9),
            "global constant diagonal sections changed")
    main_pattern = tuple(range(1, 12))
    exceptional = tuple(
        (cell, values) for cell, values in qsets.items()
        if values != main_pattern
    )
    expected_exceptional = (
        ((2, 2), tuple(range(1, 11))),
        ((2, 3), tuple(range(1, 11))),
        ((6, 5), (1, 3, 4, 5, 6, 7, 8, 9, 10, 11)),
        ((6, 6), (1, 2, 3, 4, 5, 6, 7, 8, 9, 11)),
        ((7, 2), tuple(range(1, 11))),
        ((7, 3), tuple(range(1, 11))),
        ((11, 5), (1, 3, 4, 5, 6, 7, 8, 9, 10, 11)),
        ((11, 6), (1, 2, 3, 4, 5, 6, 7, 8, 9, 11)),
    )
    require(exceptional == expected_exceptional,
            "exceptional diagonal atlas changed")

    digit_support_counts = tuple(
        sum(digit_masses[ell][h] > 0 for ell in range(Q7))
        for h in range(P)
    )
    future_support = tuple(
        h for h in range(P) if digit_support_counts[h] > 0
    )
    require(
        digit_support_counts == (0, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 0),
        "delayed future-digit support changed",
    )
    require(future_support == tuple(range(1, 12)),
            "future carrier is not the expected eleven-point puncture")
    require(translation_stabilizer(future_support) == (0,),
            "punctured future carrier acquired a C13 translation action")

    nonzero_root_slots = len(diagonal) * P * Q7 * (P - 1) * P
    zero_root_slots = len(diagonal) * P * Q7 * P
    require(nonzero_root_slots == 2_299_752 and zero_root_slots == 191_646,
            "enlarged tensor universe changed")

    print("THM-2616 exact cross-time target/future diagonal controls")
    print(f"full_nonzero_root_slots={nonzero_root_slots} retained_zero_root_slots={zero_root_slots}")
    print(f"full_positive={full_positive} pre_route_global_content={full_content} physical_content={P*full_content}")
    print(f"shards={tuple(shard_summary)} cache_endpoints={cache_size} scalar_route_controls={scalar_controls}")
    print(f"future_digit_phase_support={digit_support_counts} future_support={future_support}")
    print(f"diagonal_positive_fine={diagonal_positive} positive_pairs={positive_pairs} unit_slices={unit_slices}")
    print(f"diagonal_qset_size_hist={sorted(size_hist.items())}")
    print(f"diagonal_qset_pattern_hist={sorted(pattern_hist.items())}")
    print(f"common_constant_sections={common}")
    print(f"exceptional_diagonal_atlas={exceptional}")
    print("difference_full_cells=84 diagonal_stabilizer_size_hist=((1,84),)")
    print("verdict=PASS: positive cross-time sections exist; no principal C13 future action")


if __name__ == "__main__":
    main()
