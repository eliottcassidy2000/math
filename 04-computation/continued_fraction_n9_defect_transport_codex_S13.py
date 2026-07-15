#!/usr/bin/env python3
"""Exact centered-CF transport of THM-828's rank-four n=9 defect.

The THM-812/813 coordinate copies are continued recursively through X10.
The resulting X9->X10 map is applied both linearly to all fifteen nonzero
defect sectors and literally to all 58 false-palindrome Q witnesses.

Assumption challenge: the proof-bearing vertices here are defect sectors and
literal complement/reflection orbits, not runners or bare tournament nodes.
The computation preserves coordinate reflection, complement, Q distinctions,
and exact tournament isomorphism within the finite witness image.  It does not
preserve or claim any metric LRC wall, owner, or loneliness predicate.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from continued_fraction_edge_descent_boundary_codex_S13 import RHO_6_7
from continued_fraction_metagraph_edge_transport_codex_S13 import (
    centered_height,
    centered_increments,
    coordinate_copy,
    rho_5_6,
)
from mobius_cech_metagraph_codec_codex_S12 import compact_partition
from tournament_tiling_metagraph_address_codex_S4 import carrier_tournament


WITNESSES = Path("05-knowledge/results/mobius_cech_n9_exact_join_witnesses_codex_S13.tsv")
OUTPUT = Path("05-knowledge/results/continued_fraction_n9_defect_transport_codex_S13.out")
JSON_OUTPUT = Path("05-knowledge/results/continued_fraction_n9_defect_transport_codex_S13.json")

DEFECT_BASIS = (0x192486, 0x8C2C0C, 0x11B4600, 0x4483414)
RHO_4_5 = (0, 1, 1, 2, 0, 2)


def tiles(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((a, b) for b in range(1, n - 1) for a in range(n, b + 1, -1))


def reflection(n: int) -> tuple[int, ...]:
    schema = tiles(n)
    position = {tile: i for i, tile in enumerate(schema)}
    return tuple(position[(n - b + 1, n - a + 1)] for a, b in schema)


def reflect(mask: int, n: int) -> int:
    return sum(((mask >> i) & 1) << j for i, j in enumerate(reflection(n)))


def phase_schedule() -> dict[int, dict[str, object]]:
    """Continue the THM-812 phase by THM-778's s'=(a-s) mod 2."""
    phase = 1
    result = {}
    for source_n in range(5, 10):
        p, q, quotient = source_n - 3, source_n - 2, 1
        increments = centered_increments(q, p, phase)
        next_phase = (quotient - phase) % 2
        result[source_n] = {
            "source_n": source_n,
            "target_n": source_n + 1,
            "p": p,
            "q": q,
            "phase": phase,
            "heights": [centered_height(q, p, phase, i) for i in range(p + 1)],
            "increments": list(increments),
            "next_phase": next_phase,
        }
        phase = next_phase
    return result


def recursive_rho(source_n: int, core_rho: tuple[int, ...], increments: tuple[int, ...]) -> tuple[int, ...]:
    """Build X_n->X_(n+1): copied legs plus the X_(n-2)->X_(n-1) core map."""
    source, target = tiles(source_n), tiles(source_n + 1)
    source_position = {tile: i for i, tile in enumerate(source)}
    core_source, core_target = tiles(source_n - 2), tiles(source_n - 1)
    assert len(core_rho) == len(core_target)
    mapping: dict[tuple[int, int], tuple[int, int]] = {(source_n + 1, 1): (source_n, 1)}

    target_high = [(a, 1) for a in range(source_n, 2, -1)]
    source_high = [(a, 1) for a in range(source_n - 1, 2, -1)]
    target_low = [(source_n + 1, b) for b in range(2, source_n)]
    source_low = [(source_n, b) for b in range(2, source_n - 1)]
    for target_leg, source_leg in ((target_high, source_high), (target_low, source_low)):
        cursor = 0
        for source_tile, copies in zip(source_leg, increments, strict=True):
            for _ in range(copies):
                mapping[target_leg[cursor]] = source_tile
                cursor += 1
        assert cursor == len(target_leg)

    core_target_position = {tile: i for i, tile in enumerate(core_target)}
    for a, b in target:
        if a < source_n + 1 and b >= 2:
            local_target = (a - 1, b - 1)
            local_source = core_source[core_rho[core_target_position[local_target]]]
            mapping[(a, b)] = (local_source[0] + 1, local_source[1] + 1)
    assert len(mapping) == len(target)
    return tuple(source_position[mapping[tile]] for tile in target)


def rho_audit(source_n: int, rho: tuple[int, ...]) -> dict[str, object]:
    source_sigma, target_sigma = reflection(source_n), reflection(source_n + 1)
    return {
        "source_n": source_n,
        "target_n": source_n + 1,
        "rho": list(rho),
        "source_coordinates": len(tiles(source_n)),
        "target_coordinates": len(tiles(source_n + 1)),
        "surjective": len(set(rho)) == len(tiles(source_n)),
        "reflection_intertwining_failures": sum(
            rho[target_sigma[j]] != source_sigma[rho[j]] for j in range(len(rho))
        ),
    }


def gf2_rank(words: list[int] | tuple[int, ...]) -> int:
    pivots: list[int] = []
    for word in words:
        value = word
        for pivot in pivots:
            value = min(value, value ^ pivot)
        if value:
            pivots.append(value)
            pivots.sort(reverse=True)
    return len(pivots)


def copy_mask(mask: int, source_n: int, rho: tuple[int, ...]) -> int:
    return coordinate_copy(mask, source_n, source_n + 1, rho)


def syndrome_vector(mask: int, n: int) -> tuple[int, ...]:
    schema, sigma = tiles(n), reflection(n)
    result = [mask & 1]
    for tau in range(3, n):
        low = [i for i, tile in enumerate(schema) if sum(tile) - 1 == tau]
        result.append(sum((mask >> i) & 1 for i in low) & 1)
        result.append(sum((mask >> sigma[i]) & 1 for i in low) & 1)
    fixed = [i for i, tile in enumerate(schema) if sum(tile) - 1 == n]
    result.append(sum((mask >> i) & 1 for i in fixed) & 1)
    return tuple(result)


def syndrome_failures(mask: int, n: int) -> list[str]:
    vector = syndrome_vector(mask, n)
    labels = ["apex"]
    for tau in range(3, n):
        labels.extend((f"tau{tau}_low", f"tau{tau}_high"))
    labels.append(f"tau{n}_fixed")
    return [label for label, value in zip(labels, vector, strict=True) if value]


def s2_word(mask: int, n: int) -> tuple[object, ...]:
    schema, sigma = tiles(n), reflection(n)
    result: list[object] = []
    for tau in range(3, n):
        counts = [0, 0, 0, 0]
        for i, tile in enumerate(schema):
            if sum(tile) - 1 == tau:
                counts[2 * ((mask >> i) & 1) + ((mask >> sigma[i]) & 1)] += 1
        result.append(tuple(counts))
    result.append(sum((mask >> i) & 1 for i, tile in enumerate(schema) if sum(tile) - 1 == n))
    return tuple(result)


def s2_first_moments(mask: int, n: int) -> tuple[object, ...]:
    schema, sigma = tiles(n), reflection(n)
    result: list[object] = []
    for tau in range(3, n):
        moments = [0, 0, 0, 0]
        positions = [i for i, tile in enumerate(schema) if sum(tile) - 1 == tau]
        for position, i in enumerate(positions):
            state = 2 * ((mask >> i) & 1) + ((mask >> sigma[i]) & 1)
            moments[state] += position
        result.append(tuple(moments))
    fixed = [i for i, tile in enumerate(schema) if sum(tile) - 1 == n]
    result.append(sum(position * ((mask >> i) & 1) for position, i in enumerate(fixed)))
    return tuple(result)


def skew(mask: int, n: int) -> tuple[int, ...]:
    schema, sigma = tiles(n), reflection(n)
    result = []
    for tau in range(3, n):
        positions = [i for i, tile in enumerate(schema) if sum(tile) - 1 == tau]
        result.append(
            sum(
                position * (((mask >> i) & 1) - ((mask >> sigma[i]) & 1))
                for position, i in enumerate(positions)
            )
        )
    return tuple(result)


def eta(moment: tuple[int, ...]) -> int:
    return next((1 if value > 0 else -1 for value in moment if value), 0)


def tournament(mask: int, n: int) -> tuple[int, ...]:
    """Tournament in reversed vertex labels and the complement bit convention.

    Relative to the direct 1..n chart used by the exact subset-DP metagraph
    address scripts, this is the tournament of ``mask ^ FULL`` after reversing
    the vertex labels.  Merged-node equality and the unordered projected cell
    are unchanged, but the two ordered endpoint projections are interchanged.
    """
    out = [0] * n
    for i in range(n - 1):
        out[i] |= 1 << (i + 1)
    for bit, (a, b) in enumerate(tiles(n)):
        i, j = n - a, n - b
        if (mask >> bit) & 1:
            out[j] |= 1 << i
        else:
            out[i] |= 1 << j
    return tuple(out)


def tournament_converse(value: tuple[int, ...]) -> tuple[int, ...]:
    full = (1 << len(value)) - 1
    return tuple((~row) & (full ^ (1 << i)) for i, row in enumerate(value))


def isomorphic(left: tuple[int, ...], right: tuple[int, ...], counter: dict[str, int]) -> bool:
    """Exact degree-refined backtracking; fingerprints only prune."""
    counter["calls"] += 1
    n = len(left)
    left_degree = [row.bit_count() for row in left]
    right_degree = [row.bit_count() for row in right]
    if sorted(left_degree) != sorted(right_degree):
        return False

    def signature(graph: tuple[int, ...], degree: list[int], u: int) -> tuple[object, ...]:
        out_degrees = tuple(sorted(degree[v] for v in range(n) if (graph[u] >> v) & 1))
        in_degrees = tuple(
            sorted(degree[v] for v in range(n) if v != u and not ((graph[u] >> v) & 1))
        )
        return degree[u], out_degrees, in_degrees

    left_signature = [signature(left, left_degree, u) for u in range(n)]
    right_signature = [signature(right, right_degree, v) for v in range(n)]
    if Counter(left_signature) != Counter(right_signature):
        return False
    candidates = {
        u: [v for v in range(n) if left_signature[u] == right_signature[v]] for u in range(n)
    }
    order = sorted(range(n), key=lambda u: (len(candidates[u]), left_signature[u], u))
    image: dict[int, int] = {}
    used: set[int] = set()

    def extend(position: int) -> bool:
        counter["states"] += 1
        if position == n:
            return True
        u = order[position]
        for v in candidates[u]:
            if v in used:
                continue
            if all(((left[u] >> w) & 1) == ((right[v] >> image[w]) & 1) for w in image):
                image[u] = v
                used.add(v)
                if extend(position + 1):
                    return True
                used.remove(v)
                del image[u]
        return False

    return extend(0)


def merged_node_partition(masks: set[int], n: int) -> tuple[dict[int, int], dict[str, int]]:
    representatives: list[tuple[int, ...]] = []
    rank: dict[int, int] = {}
    counter = {"calls": 0, "states": 0}
    for mask in sorted(masks):
        value = tournament(mask, n)
        for node, representative in enumerate(representatives):
            if isomorphic(value, representative, counter) or isomorphic(
                value, tournament_converse(representative), counter
            ):
                rank[mask] = node
                break
        else:
            rank[mask] = len(representatives)
            representatives.append(value)
    return rank, {"merged_nodes": len(representatives), **counter}


def classifier_self_test() -> dict[str, object]:
    rows = {}
    for n, expected in ((5, 10), (6, 34)):
        masks = set(range(1 << len(tiles(n))))
        _, stats = merged_node_partition(masks, n)
        assert stats["merged_nodes"] == expected
        rows[str(n)] = stats
    return rows


def read_witnesses(path: Path) -> list[dict[str, int]]:
    rows = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({key: int(row[key], 16) for key in ("D", "u", "v")})
    assert len(rows) == 58
    return rows


def projected_cells(rows: list[dict[str, int]], rho: tuple[int, ...]) -> dict[str, object]:
    source_full = (1 << len(tiles(9))) - 1
    target_full = (1 << len(tiles(10))) - 1
    source_masks: set[int] = set()
    target_masks: set[int] = set()
    for row in rows:
        for mask in (row["u"], row["v"]):
            source_masks.update((mask, mask ^ source_full))
            target = copy_mask(mask, 9, rho)
            target_masks.update((target, target ^ target_full))
    source_node, source_stats = merged_node_partition(source_masks, 9)
    target_node, target_stats = merged_node_partition(target_masks, 10)

    source_p, target_p = [], []
    reflection_failures = 0
    for row in rows:
        u, v = row["u"], row["v"]
        pu = tuple(sorted((source_node[u], source_node[u ^ source_full])))
        pv = tuple(sorted((source_node[v], source_node[v ^ source_full])))
        tu, tv = copy_mask(u, 9, rho), copy_mask(v, 9, rho)
        ptu = tuple(sorted((target_node[tu], target_node[tu ^ target_full])))
        ptv = tuple(sorted((target_node[tv], target_node[tv ^ target_full])))
        reflection_failures += pu != pv or ptu != ptv
        source_p.append(pu)
        target_p.append(ptu)

    source_fibres: dict[tuple[int, int], list[int]] = defaultdict(list)
    target_fibres: dict[tuple[int, int], list[int]] = defaultdict(list)
    images: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
    for index, (source, target) in enumerate(zip(source_p, target_p, strict=True)):
        source_fibres[source].append(index)
        target_fibres[target].append(index)
        images[source].add(target)

    def fibre_histogram(values: list[int]) -> dict[int, int]:
        return dict(sorted(Counter(Counter(values).values()).items()))

    source_direct = [source_node[row["u"]] for row in rows]
    source_partner = [source_node[row["u"] ^ source_full] for row in rows]
    target_direct = [target_node[copy_mask(row["u"], 9, rho)] for row in rows]
    target_partner = [
        target_node[copy_mask(row["u"], 9, rho) ^ target_full] for row in rows
    ]
    # Independently replayed with the subset-DP lex canonical codes used by the
    # metagraph-address exporter: direct-chart counts are 53/54, exchanged here.
    assert (len(set(source_direct)), fibre_histogram(source_direct)) == (
        54, {1: 51, 2: 2, 3: 1}
    )
    assert (len(set(source_partner)), fibre_histogram(source_partner)) == (
        53, {1: 48, 2: 5}
    )
    assert len(set(zip(source_direct, source_partner, strict=True))) == 58
    assert len(source_fibres) == 58 and all(len(fibre) == 1 for fibre in source_fibres.values())
    return {
        "source_masks_classified": len(source_masks),
        "target_masks_classified": len(target_masks),
        "source_classifier": source_stats,
        "target_classifier": target_stats,
        "reflection_P_failures": reflection_failures,
        "source_P_cells": len(source_fibres),
        "source_direct_projection_nodes": len(set(source_direct)),
        "source_direct_projection_fibre_histogram": fibre_histogram(source_direct),
        "source_complement_partner_nodes": len(set(source_partner)),
        "source_complement_partner_fibre_histogram": fibre_histogram(source_partner),
        "source_P_loops": sum(left == right for left, right in source_p),
        "source_P_fibre_histogram": dict(sorted(Counter(map(len, source_fibres.values())).items())),
        "target_P_cells": len(target_fibres),
        "target_direct_projection_nodes": len(set(target_direct)),
        "target_direct_projection_fibre_histogram": fibre_histogram(target_direct),
        "target_complement_partner_nodes": len(set(target_partner)),
        "target_complement_partner_fibre_histogram": fibre_histogram(target_partner),
        "target_P_loops": sum(left == right for left, right in target_p),
        "target_P_fibre_histogram": dict(sorted(Counter(map(len, target_fibres.values())).items())),
        "P_descent_failures": sum(len(values) != 1 for values in images.values()),
        "source_P": source_p,
        "target_P": target_p,
        "convention_reconciliation": (
            "tournament(mask) here is the direct-chart tournament(mask^FULL) up to label "
            "reversal; hence the 54-node direct count here is the exact subset-DP chart's "
            "complement-partner count, while the 53-node partner count here is its direct count"
        ),
    }


def audit(witness_path: Path) -> dict[str, object]:
    schedules = phase_schedule()
    assert [schedules[n]["phase"] for n in range(5, 10)] == [1, 0, 1, 0, 1]
    assert schedules[5]["increments"] == [1, 2]
    assert schedules[6]["increments"] == [1, 2, 1]

    rebuilt_6_7 = recursive_rho(6, RHO_4_5, tuple(schedules[6]["increments"]))
    assert rebuilt_6_7 == RHO_6_7
    rhos: dict[int, tuple[int, ...]] = {5: rho_5_6(), 6: RHO_6_7}
    rhos[7] = recursive_rho(7, rhos[5], tuple(schedules[7]["increments"]))
    rhos[8] = recursive_rho(8, rhos[6], tuple(schedules[8]["increments"]))
    rhos[9] = recursive_rho(9, rhos[7], tuple(schedules[9]["increments"]))
    rho_rows = [rho_audit(n, rhos[n]) for n in range(5, 10)]
    assert all(row["surjective"] and row["reflection_intertwining_failures"] == 0 for row in rho_rows)

    rho = rhos[9]
    assert rho == (
        0, 1, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12, 12, 13, 14, 15,
        16, 17, 17, 18, 19, 20, 15, 21, 22, 23, 24, 21, 25, 26, 24, 27, 26, 27,
    )
    image_basis = tuple(copy_mask(word, 9, rho) for word in DEFECT_BASIS)
    assert gf2_rank(DEFECT_BASIS) == gf2_rank(image_basis) == 4
    sectors = []
    source_by_coordinate: dict[int, int] = {}
    for coordinate in range(16):
        source = 0
        for bit, basis_word in enumerate(DEFECT_BASIS):
            if (coordinate >> bit) & 1:
                source ^= basis_word
        source_by_coordinate[coordinate] = source
    source_coordinate = {word: coordinate for coordinate, word in source_by_coordinate.items()}
    assert len(source_coordinate) == 16

    rows = read_witnesses(witness_path)
    multiplicity = Counter(row["D"] for row in rows)
    for coordinate in range(1, 16):
        source = source_by_coordinate[coordinate]
        target = copy_mask(source, 9, rho)
        sectors.append(
            {
                "coordinate": f"{coordinate:04b}",
                "source_D": f"0x{source:07x}",
                "target_D": f"0x{target:09x}",
                "target_weight": target.bit_count(),
                "target_syndrome": not any(syndrome_vector(target, 10)),
                "target_syndrome_failures": syndrome_failures(target, 10),
                "witness_pairs": multiplicity[source],
            }
        )
    good_coordinates = [
        coordinate
        for coordinate, source in source_by_coordinate.items()
        if not any(syndrome_vector(copy_mask(source, 9, rho), 10))
    ]
    assert good_coordinates == [0, 1, 2, 3]
    assert gf2_rank([source_by_coordinate[c] for c in good_coordinates]) == 2
    assert syndrome_failures(image_basis[2], 10) == ["tau7_low", "tau7_high"]
    assert syndrome_failures(image_basis[3], 10) == ["tau5_low", "tau5_high"]
    assert gf2_rank([sum(value << bit for bit, value in enumerate(syndrome_vector(word, 10))) for word in image_basis]) == 2

    projected = projected_cells(rows, rho)
    assert projected["reflection_P_failures"] == 0
    assert projected["source_P_cells"] == projected["target_P_cells"] == 58
    assert projected["source_P_loops"] == projected["target_P_loops"] == 0
    assert projected["source_classifier"]["merged_nodes"] == 97
    assert projected["target_classifier"]["merged_nodes"] == 116
    assert projected["source_P_fibre_histogram"] == projected["target_P_fibre_histogram"] == {1: 58}
    assert projected["P_descent_failures"] == 0

    prefix = [0] * 8
    first_s2_failure: Counter[int] = Counter()
    first_m1_separator: Counter[int] = Counter()
    eta_transport: Counter[str] = Counter()
    surviving_sector_pairs: Counter[int] = Counter()
    target_q_representatives = set()
    carrier_values: dict[str, dict[int, object]] = defaultdict(dict)
    for index, row in enumerate(rows):
        source_u, source_v = row["u"], row["v"]
        target_u, target_v = copy_mask(source_u, 9, rho), copy_mask(source_v, 9, rho)
        assert target_v == reflect(target_u, 10)
        target_q = min(target_u, target_v)
        target_q_representatives.add(target_q)

        left, right = s2_word(target_u, 10), s2_word(target_v, 10)
        equal_so_far = True
        for layer in range(8):
            if equal_so_far and left[layer] == right[layer]:
                prefix[layer] += 1
            elif equal_so_far:
                first_s2_failure[layer + 3] += 1
                equal_so_far = False
        if equal_so_far:
            surviving_sector_pairs[row["D"]] += 1
            left_m1, right_m1 = s2_first_moments(target_u, 10), s2_first_moments(target_v, 10)
            first = next(layer for layer in range(8) if left_m1[layer] != right_m1[layer])
            first_m1_separator[first + 3] += 1

        source_eta, target_eta = eta(skew(source_u, 9)), eta(skew(target_u, 10))
        eta_transport["preserved" if source_eta == target_eta else "reversed"] += 1
        target_s2_orbit = tuple(sorted((left, right)))
        carrier_values["source_sector"][index] = row["D"]
        carrier_values["target_difference"][index] = copy_mask(row["D"], 9, rho)
        carrier_values["target_syndrome"][index] = syndrome_vector(copy_mask(row["D"], 9, rho), 10)
        carrier_values["target_S2_orbit"][index] = target_s2_orbit
        carrier_values["source_projected_P"][index] = projected["source_P"][index]
        carrier_values["target_projected_P"][index] = projected["target_P"][index]
        carrier_values["target_literal_Q"][index] = target_q

    assert len(target_q_representatives) == 58
    assert prefix == [58, 58, 20, 20, 10, 10, 10, 10]
    assert first_s2_failure == {5: 38, 7: 10}
    assert surviving_sector_pairs == {0x192486: 6, 0x8C2C0C: 2, 0x95088A: 2}
    assert first_m1_separator == {7: 4, 8: 6}
    assert eta_transport == {"preserved": 50, "reversed": 8}
    reversed_eta_sectors = Counter()
    for row in rows:
        if eta(skew(row["u"], 9)) != eta(skew(copy_mask(row["u"], 9, rho), 10)):
            reversed_eta_sectors[row["D"]] += 1
    assert reversed_eta_sectors == {0x18E4E8A: 4, 0x1976A0C: 4}

    carrier_order = (
        "source_sector",
        "target_difference",
        "target_syndrome",
        "target_S2_orbit",
        "source_projected_P",
        "target_projected_P",
        "target_literal_Q",
    )
    carrier_stats = {
        name: compact_partition(carrier_values[name]) for name in carrier_order
    }
    assert {name: carrier_stats[name]["cells"] for name in carrier_order} == {
        "source_sector": 11,
        "target_difference": 11,
        "target_syndrome": 4,
        "target_S2_orbit": 58,
        "source_projected_P": 58,
        "target_projected_P": 58,
        "target_literal_Q": 58,
    }
    retention = carrier_tournament(carrier_stats, "retention")
    economy = carrier_tournament(carrier_stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carrier_order))
        for j in range(i + 1, len(carrier_order))
    )

    classifier_regression = classifier_self_test()
    return {
        "schema_version": 1,
        "source_theorems": ["THM-778", "THM-812", "THM-813", "THM-828"],
        "phase_schedules": schedules,
        "rho_audits": rho_rows,
        "rho_9_10": list(rho),
        "linear_defect_action": {
            "source_basis": [f"0x{word:07x}" for word in DEFECT_BASIS],
            "target_basis": [f"0x{word:09x}" for word in image_basis],
            "source_rank": gf2_rank(DEFECT_BASIS),
            "target_rank": gf2_rank(image_basis),
            "target_syndrome_intersection_rank": 2,
            "target_syndrome_coordinates": [f"{c:04b}" for c in good_coordinates],
            "parity_quotient": {
                "basis_1": [],
                "basis_2": [],
                "basis_3": ["tau7_low", "tau7_high"],
                "basis_4": ["tau5_low", "tau5_high"],
            },
            "sectors": sectors,
        },
        "literal_witness_action": {
            "source_Q_orbits": len(rows),
            "distinct_target_Q_orbits": len(target_q_representatives),
            "target_S2_prefix_tau3_to_fixed10": prefix,
            "first_target_S2_failure": dict(sorted(first_s2_failure.items())),
            "target_raw_S2_equal_pairs": sum(surviving_sector_pairs.values()),
            "target_raw_S2_equal_sectors": {
                f"0x{key:07x}": value for key, value in sorted(surviving_sector_pairs.items())
            },
            "first_target_M1_separator": dict(sorted(first_m1_separator.items())),
            "eta_transport": dict(sorted(eta_transport.items())),
            "eta_reversed_sectors": {
                f"0x{key:07x}": value for key, value in sorted(reversed_eta_sectors.items())
            },
        },
        "projected_cell_action": {key: value for key, value in projected.items() if key not in ("source_P", "target_P")},
        "classifier_regression": classifier_regression,
        "tournament_analysis": {
            "vertices": list(carrier_order),
            "pairwise_observable": "number of unordered THM-828 Q-orbit pairs separated by the carrier",
            "switches": ["retention", "retention_per_log2_cells"],
            "tie_hamiltonian_path": list(carrier_order),
            "carrier_stats": carrier_stats,
            "retention": retention,
            "economy": economy,
            "edge_flips_between_gauges": flips,
        },
        "preservation_boundary": {
            "preserves": [
                "literal masks and complement lines",
                "staircase reflection and Q orbits",
                "rank-four linear defect distinctions",
                "exact projected P equality within the 58-orbit image",
            ],
            "destroys_or_does_not_claim": [
                "48 of 58 source raw-S2 reflection equalities",
                "the eta orientation on 8 witnesses",
                "a global n=10 metagraph atlas or node ordering",
                "metric LRC gaps, owners, walls, or loneliness",
            ],
        },
    }


def render(result: dict[str, object]) -> str:
    linear = result["linear_defect_action"]
    literal = result["literal_witness_action"]
    projected = result["projected_cell_action"]
    lines = [
        "CENTERED-CF ACTION ON THM-828'S N=9 RANK-FOUR DEFECT",
        "=" * 78,
        "phase schedule n=5..9: "
        + ", ".join(
            f"{n}->{int(n)+1}:s={row['phase']} d={tuple(row['increments'])}"
            for n, row in result["phase_schedules"].items()
        ),
        f"rho_9_10={tuple(result['rho_9_10'])}",
        "rho audits: "
        + ", ".join(
            f"{row['source_n']}->{row['target_n']} surj={row['surjective']} "
            f"sigma_fail={row['reflection_intertwining_failures']}"
            for row in result["rho_audits"]
        ),
        "",
        "LINEAR DEFECT ACTION",
        f"  source basis={tuple(linear['source_basis'])}",
        f"  target basis={tuple(linear['target_basis'])}",
        f"  rank source/target={linear['source_rank']}/{linear['target_rank']}",
        f"  target S2-syndrome intersection rank={linear['target_syndrome_intersection_rank']} "
        f"coordinates={tuple(linear['target_syndrome_coordinates'])}",
        f"  parity quotient={linear['parity_quotient']}",
        "  verdict: rank four is retained literally but not syndrome-invariant; "
        "the last two basis directions become visible at target tau7 and tau5.",
        "",
        "LITERAL 58-Q-ORBIT ACTION",
        f"  source/target distinct Q={literal['source_Q_orbits']}/{literal['distinct_target_Q_orbits']}",
        f"  target S2 prefix tau3..fixed10={tuple(literal['target_S2_prefix_tau3_to_fixed10'])}",
        f"  first S2 failure={literal['first_target_S2_failure']}",
        f"  target raw-S2 equal pairs/sectors={literal['target_raw_S2_equal_pairs']}/"
        f"{literal['target_raw_S2_equal_sectors']}",
        f"  first target M1 separator={literal['first_target_M1_separator']}; "
        f"eta transport={literal['eta_transport']} reversed sectors={literal['eta_reversed_sectors']}",
        "",
        "EXACT PROJECTED-CELL AUDIT INSIDE THE WITNESS IMAGE",
        f"  classified masks source/target={projected['source_masks_classified']}/"
        f"{projected['target_masks_classified']}; local merged nodes="
        f"{projected['source_classifier']['merged_nodes']}/"
        f"{projected['target_classifier']['merged_nodes']}",
        f"  P cells source/target={projected['source_P_cells']}/{projected['target_P_cells']}; "
        f"descent/reflection failures={projected['P_descent_failures']}/"
        f"{projected['reflection_P_failures']}",
        "  ordered endpoint projections source direct/partner="
        f"{projected['source_direct_projection_nodes']}/"
        f"{projected['source_complement_partner_nodes']} with fibre histograms "
        f"{projected['source_direct_projection_fibre_histogram']}/"
        f"{projected['source_complement_partner_fibre_histogram']}",
        "  convention audit: local tournament(mask) is direct-chart tournament(mask^FULL) "
        "up to vertex-label reversal; this exchanges the 54/53 ordered projections",
        f"  P fibre histograms={projected['source_P_fibre_histogram']}/"
        f"{projected['target_P_fibre_histogram']}",
        f"  source/target P loops={projected['source_P_loops']}/{projected['target_P_loops']}",
    ]
    ta = result["tournament_analysis"]
    lines.extend(
        [
            "",
            "TOURNAMENT ANALYSIS (information carriers as vertices)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; flips={ta['edge_flips_between_gauges']}",
            f"  retention scores={ta['retention']['score_hist']} C3="
            f"{ta['retention']['directed_3cycles']} SCC={ta['retention']['scc_sizes']} "
            f"Hpaths={ta['retention']['hamiltonian_paths']}",
            "",
            "Boundary: this is an exact coordinate/Q/finite-image P result, not metric LRC transport "
            "and not a global n=10 metagraph classification.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--witnesses", type=Path, default=WITNESSES)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--json", type=Path, default=JSON_OUTPUT)
    args = parser.parse_args()
    result = audit(args.witnesses)
    text = render(result)
    print(text, end="")
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
