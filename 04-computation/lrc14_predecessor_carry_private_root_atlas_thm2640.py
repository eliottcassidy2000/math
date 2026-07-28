#!/usr/bin/env python3
"""Exact THM-2640 companion for the missing THM-2623 predecessor carry.

The selector c=floor(R*x) mod 13 is retained before delayed integration by
descending once from R=13^6 to S=13^5.  For t={S*x}, c=floor(13*t) and
{R*x}={13*t}.  This lets one integrate all thirteen carry cells without
materialising 13^6 pullback intervals.
"""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import lrc14_cross_time_target_future_diagonal_thm2616 as core
import lrc14_successor_halfcell_carry_no_go_thm2623 as prior

old = core.old
P, Q7, T = core.P, core.Q7, core.T
R = old.R
S = R // P


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def build_pair_prefixes(module):
    """[sector][clock][h][kappa] prefixes on delayed y={R*x}."""
    result = []
    for word in prior.sector_words(module):
        by_clock = []
        for ell in range(Q7):
            qell = module.subtract_comb(
                word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
            )
            by_h = []
            for h in range(P):
                pair = []
                for kappa in range(2):
                    d = 2 * h + kappa
                    digit = old.sat.intersect_interval(
                        qell, d * T // (2 * P), (d + 1) * T // (2 * P)
                    )
                    pair.append(module.make_prefix(digit))
                by_h.append(tuple(pair))
            by_clock.append(tuple(by_h))
        result.append(tuple(by_clock))
    return tuple(result)


def phi_pair(y, prefixes, cache):
    if y in cache:
        return cache[y]
    values = tuple(old.phi_at(y, *prefix) for prefix in prefixes)
    cache[y] = values
    return values


def delayed_carry_pair(pieces, prefixes, cache):
    """Return [c][kappa] exact overlaps and check the carry partition.

    If z={S*x}, then c=floor(13*z) and y={13*z}={R*x}.  On the
    c-cell, the transformed prefix is Phi(y)/13; every earlier c-cell is
    already complete.  A carry range-difference therefore evaluates all
    thirteen cells in O(2) work per physical endpoint.
    """
    if not pieces:
        return ((0, 0),) * P
    lengths = tuple(prefix[2][-1] for prefix in prefixes)
    require(all(length % P == 0 for length in lengths),
            "delayed half-digit mass does not descend through one 13-cell")
    weighted_len = 0
    acc_s = 0
    acc_y = 0
    carry_diff = [0] * (P + 1)
    partial = [[0, 0] for _ in range(P)]
    original_phi = [0, 0]

    for left, right, weight in pieces:
        weighted_len += weight * (right - left)
        for x, sign in ((right, 1), (left, -1)):
            signed = sign * weight
            z = x * S % T
            c = P * z // T
            require(0 <= c < P, "carry cell escaped F13")
            y = P * z - c * T
            require(y == x * R % T, "one-step carry descent changed future phase")
            values = phi_pair(y, prefixes, cache)
            require(all(value % P == 0 for value in values),
                    "partial delayed prefix does not descend through carry cell")
            acc_s += signed * z
            acc_y += signed * y
            if c:
                carry_diff[0] += signed
                carry_diff[c] -= signed
            for kappa, value in enumerate(values):
                partial[c][kappa] += signed * (value // P)
                original_phi[kappa] += signed * value

    floor_s = S * weighted_len - acc_s
    floor_r = R * weighted_len - acc_y
    require(floor_s % T == 0 and floor_r % T == 0,
            "weighted carry/future floor count is not integral")
    quotient_s = floor_s // T
    quotient_r = floor_r // T

    carry_coeff = 0
    result = []
    for c in range(P):
        carry_coeff += carry_diff[c]
        # The S-scale prefix formula has denominator S*T, whereas the
        # canonical delayed carrier is normalized at R*T=P*S*T.  Multiply
        # each carry piece by P so that the thirteen pieces partition the
        # existing raw numerator, not merely its one-step descended value.
        result.append(tuple(P * (
            (lengths[kappa] // P) * (quotient_s + carry_coeff)
            + partial[c][kappa]
        ) for kappa in range(2)))
    result = tuple(result)
    require(min(value for row in result for value in row) >= 0,
            "negative carry-refined delayed overlap")
    original = tuple(
        lengths[kappa] * quotient_r + original_phi[kappa]
        for kappa in range(2)
    )
    carry_sum = tuple(sum(result[c][kappa] for c in range(P))
                      for kappa in range(2))
    require(carry_sum == original,
            f"thirteen carry cells do not partition the delayed half-digit: "
            f"carry={carry_sum} original={original} qS={quotient_s} "
            f"qR={quotient_r} accS={acc_s} accY={acc_y}")
    return result


def freeze(value):
    return tuple(freeze(item) for item in value) if isinstance(value, list) else value


def shard(bounds):
    start, stop = bounds
    module, _, _, _, rails, present, starts = core.build_carrier_data()
    prefixes = build_pair_prefixes(module)
    caches = [[[[{} for _ in range(P)] for _ in range(Q7)] for _ in range(2)]][0]
    content = 0
    positive = 0
    partition_checks = 0
    singleton_checks = 0
    metadata = []
    rows_out = []
    for j in range(start, stop):
        s, ell, theta, pieces = rails[j]
        metadata.append((s, ell, theta))
        # [sector][edge L/R][carry][kappa][h][ell5]
        raw = [[[[[[0] * Q7 for _ in range(P)] for _ in range(2)]
                  for _ in range(P)] for _ in range(2)] for _ in range(2)]
        for h in range(P):
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, (-h) % P], starts[ell5, (-h) % P]
                )
                for r in range(1, P):
                    halves = (
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182, 14 * r - 13, 14 * r
                        ),
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182, 14 * r, 14 * r + 13
                        ),
                    )
                    for sector in range(2):
                        for edge, half in enumerate(halves):
                            values = delayed_carry_pair(
                                half, prefixes[sector][ell5][h],
                                caches[sector][ell5][h]
                            )
                            partition_checks += 2
                            for c in range(P):
                                for kappa in range(2):
                                    d = 2 * h + kappa
                                    b = d // P
                                    predicted = (2 * c + b + (edge == 0)) % P
                                    value = values[c][kappa]
                                    if r == predicted:
                                        raw[sector][edge][c][kappa][h][ell5] = value
                                        if value:
                                            positive += 1
                                            content = gcd(content, value)
                                    else:
                                        require(value == 0,
                                                "carry-refined half-tooth row is not private")
                                    singleton_checks += 1
        rows_out.append(freeze(raw))
    return (bounds, content, positive, partition_checks, singleton_checks,
            tuple(metadata), tuple(rows_out))


def is_unit(values, root, content):
    if root == 0:
        return False
    scale = pow(root, -1, P)
    y = tuple((value // content) * scale % P for value in values)
    reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
    return bool(old.sat.multiplication_determinant_7(reduced))


def main():
    require(R == P ** 6 and S == P ** 5 and T % R == 0,
            "canonical carry scales changed")
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(shard, core.SHARDS))
    content = 0
    positive = 0
    partition_checks = 0
    singleton_checks = 0
    metadata = []
    rows = []
    for _, c, p, pc, sc, meta, raw in results:
        content = gcd(content, c)
        positive += p
        partition_checks += pc
        singleton_checks += sc
        metadata.extend(meta)
        rows.extend(raw)

    flags = []
    unit_count = [[[[0 for _ in range(2)] for _ in range(P)]
                   for _ in range(2)] for _ in range(2)]
    for rail_rows in rows:
        rail_flags = [[[[[False] * P for _ in range(2)] for _ in range(P)]
                       for _ in range(2)] for _ in range(2)]
        for sector in range(2):
            for edge in range(2):
                for c in range(P):
                    for kappa in range(2):
                        for h in range(P):
                            d = 2 * h + kappa
                            root = (2 * c + d // P + (edge == 0)) % P
                            values = rail_rows[sector][edge][c][kappa][h]
                            if is_unit(values, root, content):
                                rail_flags[sector][edge][c][kappa][h] = True
                                unit_count[sector][edge][c][kappa] += 1
        flags.append(freeze(rail_flags))

    # Uniformly forgetting c recovers the old THM-2623 root profile.  Its
    # zero root is literal zero, and every nonempty rational profile is
    # nonconstant; Phi_13 irreducibility then makes all twelve charged
    # Fourier modes nonzero (the companion checks the exact support census).
    aggregate_support = [[[Counter(), Counter()], [Counter(), Counter()]],
                         [[Counter(), Counter()], [Counter(), Counter()]]]
    aggregate_profiles = 0
    for rail_rows in rows:
        for sector in range(2):
            for edge in range(2):
                for kappa in range(2):
                    for h in range(P):
                        root_for_c = tuple(
                            (2 * c + (2 * h + kappa) // P + (edge == 0)) % P
                            for c in range(P)
                        )
                        require(len(set(root_for_c)) == P,
                                "fixed half-digit carry/root map is not bijective")
                        for ell in range(Q7):
                            profile = [0] * P
                            for c, root in enumerate(root_for_c):
                                profile[root] = rail_rows[sector][edge][c][kappa][h][ell]
                            support_size = sum(value != 0 for value in profile)
                            if not support_size:
                                continue
                            require(profile[0] == 0 and len(set(profile)) > 1,
                                    "uniform carry sum became flat or gained root zero")
                            aggregate_support[sector][edge][kappa][support_size] += 1
                            aggregate_profiles += 1

    bycell = {}
    for j, (s, ell, theta) in enumerate(metadata):
        bycell.setdefault((s, ell), []).append(j)

    cell_rootsets = {}
    cell_carry_rootsets = {}
    cell_carry_h_rootsets = {}
    cell_qsets = {}
    witnesses = {}
    for cell, rails in sorted(bycell.items()):
        roots = set()
        qs = set()
        per_c = {}
        per_c_h = {}
        for c in range(P):
            c_roots = set()
            h_roots = {h: set() for h in range(P)}
            for j in rails:
                for sector in range(2):
                    for edge in range(2):
                        for kappa in range(2):
                            for h in range(P):
                                if not flags[j][sector][edge][c][kappa][h]:
                                    continue
                                root = (2 * c + (2 * h + kappa) // P
                                        + (edge == 0)) % P
                                require(root != 0, "unit flag has zero root")
                                roots.add(root)
                                c_roots.add(root)
                                h_roots[h].add(root)
                                qs.add(h)
                                witnesses.setdefault((cell, c, root),
                                                     (j, sector, edge, kappa, h))
            per_c[c] = tuple(sorted(c_roots))
            per_c_h[c] = {h: tuple(sorted(values))
                          for h, values in h_roots.items()}
        cell_rootsets[cell] = tuple(sorted(roots))
        cell_carry_rootsets[cell] = per_c
        cell_carry_h_rootsets[cell] = per_c_h
        cell_qsets[cell] = tuple(sorted(qs))

    root_hist = tuple(sorted(Counter(map(len, cell_rootsets.values())).items()))
    q_hist = tuple(sorted(Counter(map(len, cell_qsets.values())).items()))
    carry_root_hist = tuple(sorted(Counter(
        len(rootset) for per_c in cell_carry_rootsets.values()
        for rootset in per_c.values()
    ).items()))
    carry_h_root_hist = tuple(sorted(Counter(
        len(rootset) for per_c in cell_carry_h_rootsets.values()
        for per_h in per_c.values() for rootset in per_h.values()
    ).items()))
    deficient = tuple(
        (cell, tuple(r for r in range(1, P) if r not in rootset))
        for cell, rootset in cell_rootsets.items()
        if len(rootset) < P - 1
    )
    carry_geometric_failures = []
    carry_h_geometric_failures = []
    for cell, per_c in cell_carry_rootsets.items():
        for c, rootset in per_c.items():
            expected = tuple(sorted({(2 * c + b + edge) % P
                                     for b in (0, 1) for edge in (0, 1)}
                                    - {0}))
            if rootset != expected:
                carry_geometric_failures.append((cell, c, rootset, expected))
            for h, h_rootset in cell_carry_h_rootsets[cell][c].items():
                bset = {(2 * h) // P, (2 * h + 1) // P}
                expected_h = tuple(sorted({(2 * c + b + edge) % P
                                           for b in bset for edge in (0, 1)}
                                          - {0}))
                if h_rootset != expected_h:
                    carry_h_geometric_failures.append(
                        (cell, c, h, h_rootset, expected_h)
                    )
    q_deficient = tuple(
        (cell, tuple(h for h in range(P) if h not in qset))
        for cell, qset in cell_qsets.items() if len(qset) < P
    )
    require(not carry_geometric_failures,
            "a carry cell lost a geometrically possible unit root")
    empty_carry_h = tuple(failure for failure in carry_h_geometric_failures
                          if not failure[3])
    partial_carry_h = tuple(failure for failure in carry_h_geometric_failures
                            if failure[3])
    transition_survival_hist = tuple(sorted(Counter(
        sum(bool(rootset) for per_h in per_c.values()
            for rootset in per_h.values())
        for per_c in cell_carry_h_rootsets.values()
    ).items()))
    empty_transition_hist = tuple(sorted(Counter(
        sum(not rootset for per_h in per_c.values()
            for rootset in per_h.values())
        for per_c in cell_carry_h_rootsets.values()
    ).items()))
    graph_metric_hist = Counter()
    for per_c in cell_carry_h_rootsets.values():
        adjacency = {
            c: {h for h, rootset in per_h.items() if rootset}
            for c, per_h in per_c.items()
        }
        indegree = [sum(h in adjacency[c] for c in range(P)) for h in range(P)]
        diameter = 0
        strongly_connected = True
        for source in range(P):
            distance = {source: 0}
            queue = [source]
            for vertex in queue:
                for target in adjacency[vertex]:
                    if target not in distance:
                        distance[target] = distance[vertex] + 1
                        queue.append(target)
            if len(distance) < P:
                strongly_connected = False
            else:
                diameter = max(diameter, max(distance.values()))
        graph_metric_hist[(sum(map(len, adjacency.values())),
                           min(map(len, adjacency.values())), min(indegree),
                           diameter if strongly_connected else -1)] += 1
    aggregate_tuple = tuple(tuple(tuple(tuple(sorted(counter.items()))
                                        for counter in edge)
                                  for edge in sector)
                            for sector in aggregate_support)
    expected_aggregate = (
        (
            (((11, 3512), (12, 356)), ((11, 3855),)),
            (((11, 3857),), ((11, 3510), (12, 356))),
        ),
        (
            (((11, 1130),), ((11, 1150), (12, 376))),
            (((11, 1130), (12, 396)), ((11, 1150),)),
        ),
    )
    expected_q_deficient = (
        ((2, 2), (11,)), ((6, 5), (2,)),
        ((7, 2), (11,)), ((11, 5), (2,)),
    )
    expected_partial = (
        ((6, 3), 8, 12, (5,), (4, 5)),
        ((6, 3), 9, 12, (7,), (6, 7)),
    )
    graph_tuple = tuple(sorted(graph_metric_hist.items()))
    require((content, positive, partition_checks, singleton_checks)
            == (26, 230_042, 1_415_232, 18_398_016),
            "carry-refined raw census changed")
    require(aggregate_profiles == 20_778 and aggregate_tuple == expected_aggregate,
            "uniform-carry charged profile census changed")
    require(root_hist == ((12, 84),) and not deficient,
            "private nonzero root atlas changed")
    require(q_hist == ((12, 4), (13, 80))
            and q_deficient == expected_q_deficient,
            "future-label atlas changed")
    require(carry_root_hist == ((2, 252), (3, 840)),
            "fixed-carry geometric root atlas changed")
    require(carry_h_root_hist == ((0, 2050), (1, 170), (2, 11136), (3, 840)),
            "fixed carry/future edge histogram changed")
    require(len(carry_h_geometric_failures) == 2052
            and len(empty_carry_h) == 2050
            and partial_carry_h == expected_partial,
            "same-labelled transition failure census changed")
    require(transition_survival_hist == ((119, 4), (142, 2), (145, 2), (146, 76))
            and empty_transition_hist == ((23, 76), (24, 2), (27, 2), (50, 4)),
            "digit-edge survival census changed")
    require(graph_tuple == (((119, 2, 0, -1), 4),
                            ((142, 2, 7, 2), 2),
                            ((145, 2, 10, 3), 2),
                            ((146, 2, 11, 2), 76)),
            "digit-edge graph metrics changed")
    require(len(witnesses) == 3024,
            "fixed carry/private root witness count changed")
    require(unit_count[0][0][6] == [0, 0]
            and unit_count[0][1][6] == [0, 0]
            and unit_count[1][0][6] == [0, 134]
            and unit_count[1][1][6] == [134, 0],
            "guard cospan middle-carry seam changed")
    # The formal slope-seven clutch has an exact physical lift on the
    # carry/root/delayed-prefix subpacket: x -> x+7*delta/R.  It is not a
    # symmetry of the full packet.  Already the speed-one guard-safe factor
    # changes on a symmetric difference of measure 2*(7*delta/R).
    slope7_guard_defects = tuple(Fraction(14 * delta, R)
                                 for delta in range(1, P))
    require(all((2 * (7 * delta) - delta) % P == 0
                for delta in range(1, P)),
            "slope-seven carry/root covariance changed")
    require(7 * (P - 1) < R // 7
            and slope7_guard_defects[0] == Fraction(14, R),
            "guard symmetric-difference boundary changed")
    print("carry-refined THM-2623 guard-cospan probe")
    print(f"R={R} S={S} global_content={content} positives={positive}")
    print(f"partition_checks={partition_checks} singleton_checks={singleton_checks}")
    print(f"unit_count={freeze(unit_count)}")
    print(f"uniform_carry_profiles={aggregate_profiles} "
          f"support_hist={aggregate_tuple}; every charged DFT mode nonzero")
    print(f"cell_root_size_hist={root_hist} deficient={deficient}")
    print(f"cell_q_size_hist={q_hist} deficient={q_deficient}")
    print(f"carry_cell_root_size_hist={carry_root_hist}")
    print(f"carry_h_root_size_hist={carry_h_root_hist}")
    print(f"carry_h_geometric_failures={len(carry_h_geometric_failures)} "
          f"first={tuple(carry_h_geometric_failures[:12])}")
    print(f"carry_h_empty={len(empty_carry_h)} partial={partial_carry_h}")
    print(f"transition_survival_per_cell_hist={transition_survival_hist} "
          f"empty_per_cell_hist={empty_transition_hist}")
    print(f"digit_edge_graph_metric_hist={graph_tuple}")
    print(f"slope7_guard_symmetric_difference_delta1="
          f"{slope7_guard_defects[0]}")
    print(f"witness_count={len(witnesses)}")


if __name__ == "__main__":
    main()
