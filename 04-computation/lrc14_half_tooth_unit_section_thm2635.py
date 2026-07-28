#!/usr/bin/env python3
"""Exact companion for THM-2635's globally primitive half-tooth section.

This refines THM-2616's complete common-x carrier by the two literal halves
of every translated deep-probe tooth.  It recomputes the global content on
that *full* refined carrier before testing the THM-2629 graph r=-h-1.
"""

from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from math import gcd

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fixed_deep_affine_graph_spectrum_thm2629 as fixed
import lrc14_old_wall_successor_sector_thm2630 as wall


P = 13
Q7 = 7
SHARDS = cross.SHARDS


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def half_shard(bounds):
    start, stop = bounds
    (module, prefixes, whole_prefixes, _digit_masses, rails,
     present, starts) = cross.build_carrier_data()
    vector_cache = [dict() for _ in range(Q7)]
    content = 0
    positives = 0
    partition_checks = 0
    graph = {}

    for s, ell4, root, pieces in rails[start:stop]:
        theta = (root - 12) % P
        require(theta in (0, 1), "rail left theta chart")
        # One row stores graph values graph[epsilon][h][ell5].
        row = [[[0] * Q7 for _h in range(P)] for _eps in range(2)]
        for q in range(P):
            shift = (-q) % P
            for ell5 in range(Q7):
                overlap = cross.old.intersect_weighted_union(
                    pieces, present[ell5, shift], starts[ell5, shift]
                )
                for probe in range(1, P):
                    split = []
                    for epsilon, left, right in (
                        (0, 14 * probe, 14 * probe + 13),
                        (1, 14 * probe - 13, 14 * probe),
                    ):
                        half = cross.old.intersect_weighted_comb(
                            overlap, module.C3, 182, left, right
                        )
                        values = cross.delayed_all_digits(
                            half, prefixes[ell5], vector_cache[ell5]
                        )
                        require(
                            sum(values) == cross.old.delayed_weighted_numerator(
                                half, whole_prefixes[ell5]
                            ),
                            "half digit vector lost whole-word partition",
                        )
                        split.append(values)
                        for value in values:
                            content = gcd(content, value)
                            positives += int(value > 0)
                        if q == (-probe - 1) % P:
                            row[epsilon][q][ell5] = values[q]

                    whole = cross.old.intersect_weighted_comb(
                        overlap, module.C3, 182,
                        14 * probe - 13, 14 * probe + 13,
                    )
                    whole_values = cross.delayed_all_digits(
                        whole, prefixes[ell5], vector_cache[ell5]
                    )
                    require(
                        all(split[0][h] + split[1][h] == whole_values[h]
                            for h in range(P)),
                        "half-tooth vectors did not add to whole tooth",
                    )
                    partition_checks += P
        key = (s, ell4, theta)
        require(key not in graph, "duplicate rail key")
        graph[key] = tuple(
            tuple(tuple(values) for values in by_h) for by_h in row
        )
    return bounds, content, positives, partition_checks, graph


def unit(values, probe, content):
    require(content > 0 and all(value % content == 0 for value in values),
            "refined global content does not divide a graph row")
    scalar = pow(probe, -1, P)
    y = tuple((value // content) * scalar % P for value in values)
    reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
    return cross.old.sat.multiplication_determinant_7(reduced) != 0


def main():
    with ProcessPoolExecutor(max_workers=len(SHARDS)) as pool:
        shards = list(pool.map(half_shard, SHARDS))
    require(tuple(x[0] for x in shards) == SHARDS,
            "shards returned out of order")

    content = 0
    positives = 0
    checks = 0
    graph = {}
    for _bounds, shard_content, shard_pos, shard_checks, rows in shards:
        content = gcd(content, shard_content)
        positives += shard_pos
        checks += shard_checks
        require(not set(graph).intersection(rows), "duplicate assembled rail")
        graph.update(rows)
    require(len(graph) == 162, "refined rail bank changed")
    require(content == 26, "full half-tooth primitive content changed")
    require(positives == 1_201_816,
            "full half-tooth positive census changed")
    require(checks == 2_299_752,
            "full half-tooth partition-check universe changed")

    digest = sha256()
    for key in sorted(graph):
        digest.update(bytes(key))
        for epsilon in range(2):
            for h in range(P):
                for value in graph[key][epsilon][h]:
                    payload = str(value).encode("ascii")
                    digest.update(len(payload).to_bytes(4, "big"))
                    digest.update(payload)
    graph_digest = digest.hexdigest()
    require(graph_digest
            == "6827925d5fd66cc96d92db80e61b7b3f0ed275cc9a45708015d697f71dba9eb9",
            "complete half-graph coefficient digest changed")

    cells = [(s, ell) for s in range(1, P) for ell in range(Q7)]
    positive = defaultdict(set)
    units = defaultdict(set)
    for (s, ell, theta), row in graph.items():
        for epsilon in range(2):
            for h in range(1, 12):
                probe = (-h - 1) % P
                values = row[epsilon][h]
                if any(values):
                    positive[s, ell, theta, epsilon].add(h)
                if unit(values, probe, content):
                    units[s, ell, theta, epsilon].add(h)

    expected_positive_common = (1, 3, 4, 5, 6, 7, 8, 9, 10)
    expected_selected_unit_common = {0: (9,), 1: (3, 8, 10)}
    expected_selected_positive_size_hist = {
        0: Counter({10: 21, 11: 63}),
        1: Counter({9: 2, 10: 82}),
    }
    expected_selected_unit_size_hist = {
        0: Counter({5: 1, 9: 18, 10: 15, 11: 50}),
        1: Counter({5: 1, 8: 2, 9: 21, 10: 60}),
    }
    expected_either_unit_common = {
        0: (3, 5, 6, 7, 9),
        1: (1, 3, 4, 6, 8, 10),
    }
    expected_either_unit_size_hist = {
        0: Counter({9: 5, 10: 6, 11: 73}),
        1: Counter({8: 2, 9: 3, 10: 79}),
    }

    selected_unit_common = {}
    combined_selected = []
    for epsilon in range(2):
        selected_positive = [
            positive[s, ell, wall.selected_theta(s, ell), epsilon]
            for s, ell in cells
        ]
        selected_units = [
            units[s, ell, wall.selected_theta(s, ell), epsilon]
            for s, ell in cells
        ]
        positive_common = tuple(sorted(
            set.intersection(*map(set, selected_positive))
        ))
        unit_common = tuple(sorted(
            set.intersection(*map(set, selected_units))
        ))
        selected_unit_common[epsilon] = unit_common
        positive_size_hist = Counter(map(len, selected_positive))
        unit_size_hist = Counter(map(len, selected_units))
        require(positive_common == expected_positive_common,
                "canonical half-edge positivity intersection changed")
        require(unit_common == expected_selected_unit_common[epsilon],
                "canonical half-edge unit intersection changed")
        require(positive_size_hist
                == expected_selected_positive_size_hist[epsilon],
                "canonical half-edge positivity size spectrum changed")
        require(unit_size_hist == expected_selected_unit_size_hist[epsilon],
                "canonical half-edge unit size spectrum changed")
        print(
            f"epsilon={epsilon} selected_positive_common="
            f"{positive_common} "
            f"hist={tuple(sorted(positive_size_hist.items()))}"
        )
        print(
            f"epsilon={epsilon} selected_unit_common="
            f"{unit_common} "
            f"hist={tuple(sorted(unit_size_hist.items()))}"
        )

        # Allow either already-proved rail in a base cell, but keep epsilon.
        cell_units = [
            units[s, ell, 0, epsilon] | units[s, ell, 1, epsilon]
            for s, ell in cells
        ]
        either_common = tuple(sorted(
            set.intersection(*map(set, cell_units))
        ))
        either_size_hist = Counter(map(len, cell_units))
        require(either_common == expected_either_unit_common[epsilon],
                "either-rail half-edge unit intersection changed")
        require(either_size_hist == expected_either_unit_size_hist[epsilon],
                "either-rail half-edge unit size spectrum changed")
        print(
            f"epsilon={epsilon} either_rail_unit_common="
            f"{either_common} "
            f"hist={tuple(sorted(either_size_hist.items()))}"
        )

    # Recombining the two half-tooth coefficient rows recovers the whole
    # fixed-r Bockstein.  This is a hostile control against accidentally
    # changing the primitive lattice when the physical edge is split.
    for s, ell in cells:
        theta = wall.selected_theta(s, ell)
        row = graph[s, ell, theta]
        combined_selected.append({
            h for h in range(1, 12)
            if unit(tuple(row[0][h][clock] + row[1][h][clock]
                          for clock in range(Q7)), (-h - 1) % P, content)
        })
    combined_common = tuple(sorted(
        set.intersection(*map(set, combined_selected))
    ))
    combined_size_hist = Counter(map(len, combined_selected))
    require(combined_common == (3, 4, 5, 8, 9, 10),
            "recombined canonical fixed-r unit intersection changed")
    require(combined_size_hist == Counter({7: 1, 10: 22, 11: 61}),
            "recombined canonical fixed-r size spectrum changed")

    # The opposite graph always has abstract outgoing coordinate -r=h+1.
    # THM-2630's literal predecessor digit is
    # j=2^{-1}(r-kappa-epsilon).  Away from the split h=6 cell, kappa is
    # fixed by h; h=6 meets both halves.  Enumerate where the abstract outgoing
    # endpoint is exactly that adjacent predecessor digit.  This is a reversed
    # chronological closure, not a successor checksum.
    chronological_closures = []
    for epsilon in range(2):
        for h in range(1, 12):
            kappas = (0, 1) if h == 6 else (int(h >= 7),)
            for kappa in kappas:
                probe = (-h - 1) % P
                predecessor = 7 * (probe - kappa - epsilon) % P
                if predecessor == (-probe) % P:
                    chronological_closures.append(
                        (epsilon, h, kappa, predecessor, probe)
                    )
    require(chronological_closures
            == [(1, 3, 0, 4, 9), (1, 7, 1, 8, 5)],
            "opposite-graph chronological closure set changed")
    unit_chronological_closures = tuple(
        row for row in chronological_closures
        if row[1] in selected_unit_common[row[0]]
    )
    require(unit_chronological_closures == ((1, 3, 0, 4, 9),),
            "uniform unit chronological closure changed")

    print(f"global_half_content={content}")
    print(f"half_positive_entries={positives}")
    print(f"half_partition_digit_checks={checks}")
    print(f"half_graph_digest={graph_digest}")
    print(f"recombined_selected_unit_common={combined_common} "
          f"hist={tuple(sorted(combined_size_hist.items()))}")
    print(f"chronological_closures={tuple(chronological_closures)}")
    print(f"unit_chronological_closures={unit_chronological_closures}")
    print("verdict=PASS: one literal half-tooth edge carries a uniform "
          "opposite-graph unit; only h=3 also closes to the adjacent "
          "predecessor digit, in reversed orientation")
    print("SCOPE: imposed coefficient section on the THM-2616 carrier; "
          "not a successor checksum, semantic two-root gluing, THM-2625 "
          "current transport, or LRC(14)")


if __name__ == "__main__":
    main()
