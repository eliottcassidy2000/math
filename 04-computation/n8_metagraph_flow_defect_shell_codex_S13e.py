#!/usr/bin/env python3
"""Exact n=8 merged-metagraph flow and THM-834 defect-shell census.

The full 2^21 tiling atlas is scanned once at the complement-line level.  We
retain both literal line multiplicity and projected merged-node support, with
blue/black colour kept separate.  The 155 nodes selected by THM-834 are then
expanded in the two coloured line graphs.

Tournament Analysis does not use runners as vertices here.  The primary
vertices are converse-merged tournament classes; a secondary diagnostic uses
information carriers on the incident black lines.  The pairwise observable is
whether two literal lines are separated by a carrier, the switch is raw
retention versus retention per log-cell, and increasing carrier resolution is
the tie Hamiltonian path.

Preserved: fixed-path tilings, merged tournament nodes, literal complement
line multiplicity, grid-reflection colour, C3, endpoint curvature/current, and
the marked THM-834 support.  Destroyed: n=9 gluing outside the support label,
LRC gaps, owners, walls, metric loneliness, and continued-fraction carry.
Challenged assumption: a selected black line bank is not treated as a closed
subgraph; the full radius-one shell is computed before drawing conclusions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, deque
from fractions import Fraction
from pathlib import Path

import numpy as np


N = 8
M = 21
FULL = (1 << M) - 1
NODE_COUNT = 3528
C3_MAX = 20
ATLAS_SHA256 = "30debad3387a4ea0ef51108ea132115efda2ac2fcdfcc2c5c1d4d23155095835"


def tiles() -> tuple[tuple[int, int], ...]:
    return tuple(
        (x, y)
        for y in range(1, N - 1)
        for x in range(N, y + 1, -1)
        if x - y >= 2
    )


TILES = tiles()
TILE_INDEX = {tile: i for i, tile in enumerate(TILES)}
SIGMA = tuple(TILE_INDEX[(N - y + 1, N - x + 1)] for x, y in TILES)


def fraction_json(value: Fraction) -> dict[str, int | float]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def hist_json(hist: Counter) -> dict[str, int]:
    return {str(key): int(hist[key]) for key in sorted(hist, key=str)}


def fibre_profile(keys) -> tuple[int, Counter[int]]:
    fibres = Counter(keys)
    return len(fibres), Counter(fibres.values())


def legacy_tournament_stats(mask: int) -> tuple[tuple[int, ...], int, int]:
    """Atlas convention used by THM-834: path x->x-1, tile bit one x->y."""
    scores = [0] * N
    for x in range(2, N + 1):
        scores[x - 1] += 1
    for bit, (x, y) in enumerate(TILES):
        if (mask >> bit) & 1:
            scores[x - 1] += 1
        else:
            scores[y - 1] += 1
    c3 = math.comb(N, 3) - sum(math.comb(score, 2) for score in scores)
    axis = sum((2 * score - (N - 1)) ** 2 for score in scores)
    assert axis == N * (N * N - 1) // 3 - 8 * c3
    direct = tuple(sorted(scores))
    converse = tuple(sorted(N - 1 - score for score in scores))
    return min(direct, converse), c3, axis


def legacy_arc(mask: int, winner: int, loser: int) -> bool:
    """Arc predicate on path labels 1..N in the legacy atlas convention."""
    if abs(winner - loser) == 1:
        return winner > loser
    high, low = max(winner, loser), min(winner, loser)
    high_wins = bool((mask >> TILE_INDEX[(high, low)]) & 1)
    return high_wins if winner == high else not high_wins


def direct_endpoint_curvature(mask: int) -> int:
    """Count cyclic triples containing both fixed-path endpoints."""
    answer = 0
    for middle in range(2, N):
        clockwise = (
            legacy_arc(mask, 1, middle)
            and legacy_arc(mask, middle, N)
            and legacy_arc(mask, N, 1)
        )
        anticlockwise = (
            legacy_arc(mask, 1, N)
            and legacy_arc(mask, N, middle)
            and legacy_arc(mask, middle, 1)
        )
        answer += clockwise or anticlockwise
    return answer


def edge_adjacency(codes: np.ndarray) -> list[set[int]]:
    adjacency = [set() for _ in range(NODE_COUNT)]
    for code in map(int, codes):
        u, v = divmod(code, NODE_COUNT)
        if u != v:
            adjacency[u].add(v)
            adjacency[v].add(u)
    return adjacency


def union_adjacency(a: list[set[int]], b: list[set[int]]) -> list[set[int]]:
    return [left | right for left, right in zip(a, b)]


def distances(adjacency: list[set[int]], starts: set[int]) -> list[int]:
    answer = [-1] * NODE_COUNT
    queue = deque()
    for node in starts:
        answer[node] = 0
        queue.append(node)
    while queue:
        node = queue.popleft()
        for other in adjacency[node]:
            if answer[other] < 0:
                answer[other] = answer[node] + 1
                queue.append(other)
    return answer


def monotone_reach(adjacency: list[set[int]], c3: list[int], root: int) -> set[int]:
    seen = {root}
    queue = deque([root])
    while queue:
        node = queue.popleft()
        for other in adjacency[node]:
            if c3[other] >= c3[node] and other not in seen:
                seen.add(other)
                queue.append(other)
    return seen


def two_phase_reach(
    blue_adj: list[set[int]], black_adj: list[set[int]], c3: list[int], root: int
) -> tuple[set[int], set[tuple[int, int]]]:
    """Nondecreasing C3 paths with colour word blue* black*."""
    start = (root, 0)
    seen = {start}
    queue = deque([start])
    while queue:
        node, black_started = queue.popleft()
        if not black_started:
            for other in blue_adj[node]:
                state = (other, 0)
                if c3[other] >= c3[node] and state not in seen:
                    seen.add(state)
                    queue.append(state)
        for other in black_adj[node]:
            state = (other, 1)
            if c3[other] >= c3[node] and state not in seen:
                seen.add(state)
                queue.append(state)
    return {node for node, _ in seen}, seen


def leg_curvature(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return THM-811 q and epsilon for the legacy atlas convention.

    The formula counts cyclic triples meeting both path endpoints.  When the
    apex bit (N,1) is zero, the two boundary legs contribute and internal
    positions contribute b_i t_i.  When it is one, only internal
    (1-b_i)(1-t_i) states contribute.  Epsilon is sign-reversed relative to
    THM-811 by the legacy bit convention, so its absolute value is unchanged.
    """
    apex = (mask >> TILE_INDEX[(N, 1)]) & 1
    bottom = {
        x: (mask >> TILE_INDEX[(x, 1)]) & 1 for x in range(3, N)
    }
    top = {
        y: (mask >> TILE_INDEX[(N, y)]) & 1 for y in range(2, N - 1)
    }
    q_zero = top[2] + bottom[N - 1]
    q_one = np.zeros(len(mask), dtype=np.uint32)
    for x in range(3, N - 1):
        q_zero = q_zero + bottom[x] * top[x]
        q_one = q_one + (1 - bottom[x]) * (1 - top[x])
    q = np.where(apex == 0, q_zero, q_one).astype(np.uint8)
    epsilon = (
        sum(top.values(), np.zeros(len(mask), dtype=np.uint32)).astype(np.int16)
        - sum(bottom.values(), np.zeros(len(mask), dtype=np.uint32)).astype(np.int16)
    )
    return q, epsilon


def phase_pair(phase: list[str], u: int, v: int) -> tuple[str, str]:
    return tuple(sorted((phase[u], phase[v])))


def projected_boundary(
    codes: np.ndarray,
    phase: list[str],
    c3: list[int],
    source_phase: str,
    target_phase: str,
) -> Counter[str]:
    result: Counter[str] = Counter()
    for code in map(int, codes):
        u, v = divmod(code, NODE_COUNT)
        if {phase[u], phase[v]} != {source_phase, target_phase}:
            continue
        if c3[u] == c3[v]:
            result["tie"] += 1
            continue
        low, high = (u, v) if c3[u] < c3[v] else (v, u)
        result["outward" if phase[low] == source_phase else "reverse"] += 1
    return result


def literal_boundary(
    u: np.ndarray,
    v: np.ndarray,
    keep: np.ndarray,
    phase: list[str],
    c3: np.ndarray,
    source_phase: str,
    target_phase: str,
    source_q: np.ndarray | None = None,
) -> tuple[Counter[str], Counter[tuple[str, int]]]:
    result: Counter[str] = Counter()
    by_q: Counter[tuple[str, int]] = Counter()
    for index in np.flatnonzero(keep):
        a, b = int(u[index]), int(v[index])
        if {phase[a], phase[b]} != {source_phase, target_phase}:
            continue
        if c3[a] == c3[b]:
            result["tie"] += 1
            continue
        low, high = (a, b) if c3[a] < c3[b] else (b, a)
        direction = "outward" if phase[low] == source_phase else "reverse"
        result[direction] += 1
        if source_q is not None:
            by_q[(direction, int(source_q[index]))] += 1
    return result, by_q


def carrier_tournament(carriers: dict[str, list[object]]) -> dict:
    names = list(carriers)
    total = len(next(iter(carriers.values())))
    separated = {}
    cells = {}
    for name, values in carriers.items():
        count, profile = fibre_profile(values)
        cells[name] = count
        separated[name] = total * (total - 1) // 2 - sum(
            size * (size - 1) // 2 * multiplicity
            for size, multiplicity in profile.items()
        )

    def fingerprint(metric) -> dict:
        wins = {name: 0 for name in names}
        adjacency = {name: set() for name in names}
        flips = 0
        for i, left in enumerate(names):
            for right in names[i + 1 :]:
                lv, rv = metric(left), metric(right)
                winner = left if (lv, left) > (rv, right) else right
                wins[winner] += 1
                adjacency[winner].add(right if winner == left else left)
        triangles = sum(
            b in adjacency[a] and c in adjacency[b] and a in adjacency[c]
            for a in names for b in names for c in names if a < b < c
        )
        return {
            "scores": dict(sorted(Counter(wins.values()).items())),
            "directed_triangles": triangles,
            "hamiltonian_path": sorted(names, key=lambda name: (metric(name), name), reverse=True),
        }

    raw = fingerprint(lambda name: separated[name])
    economy = fingerprint(
        lambda name: Fraction(separated[name], max(1, math.ceil(math.log2(max(2, cells[name])))))
    )
    raw_order = raw["hamiltonian_path"]
    economy_order = economy["hamiltonian_path"]
    raw_pairs = {
        frozenset((raw_order[i], raw_order[j])): raw_order[i]
        for i in range(len(raw_order)) for j in range(i + 1, len(raw_order))
    }
    economy_pairs = {
        frozenset((economy_order[i], economy_order[j])): economy_order[i]
        for i in range(len(economy_order)) for j in range(i + 1, len(economy_order))
    }
    flips = sum(raw_pairs[pair] != economy_pairs[pair] for pair in raw_pairs)
    return {
        "vertices": names,
        "cells": cells,
        "separated_pairs": separated,
        "raw_retention": raw,
        "retention_per_log_cell": economy,
        "edge_flips_between_gauges": flips,
        "both_gauges_transitive": True,
    }


def run(atlas_path: Path, support_path: Path) -> dict:
    atlas_sha256 = hashlib.sha256(atlas_path.read_bytes()).hexdigest()
    if atlas_sha256 != ATLAS_SHA256:
        raise RuntimeError(f"unexpected n=8 atlas SHA256: {atlas_sha256}")
    atlas = np.fromfile(atlas_path, dtype="<u2").astype(np.uint32)
    if len(atlas) != 1 << M or int(atlas.min()) != 0 or int(atlas.max()) != NODE_COUNT - 1:
        raise RuntimeError("unexpected n=8 merged-node atlas")

    # Exhaust every endpoint-leg state.  Free interior bits cannot affect a
    # triple containing both path endpoints, so this is a complete independent
    # audit of the vectorized THM-811 curvature formula used below.
    leg_positions = [(N, 1)]
    leg_positions.extend((x, 1) for x in range(3, N))
    leg_positions.extend((N, y) for y in range(2, N - 1))
    leg_masks = []
    for state in range(1 << len(leg_positions)):
        mask = 0
        for bit, tile in enumerate(leg_positions):
            mask |= ((state >> bit) & 1) << TILE_INDEX[tile]
        leg_masks.append(mask)
    leg_array = np.array(leg_masks, dtype=np.uint32)
    leg_q, _ = leg_curvature(leg_array)
    assert all(int(q) == direct_endpoint_curvature(mask) for q, mask in zip(leg_q, leg_masks))

    support_json = json.loads(support_path.read_text())
    support = {int(row["node_rank"]) for row in support_json["selected_nodes"]}
    if len(support) != 155:
        raise RuntimeError("THM-834 support census mismatch")
    in_support = np.zeros(NODE_COUNT, dtype=bool)
    in_support[list(support)] = True

    indices = np.arange(1 << M, dtype=np.uint32)
    mate = indices ^ FULL
    side = indices < mate
    line_masks = indices[side]
    line_mates = mate[side]
    u = atlas[line_masks]
    v = atlas[line_mates]
    line_codes = np.minimum(u, v).astype(np.uint64) * NODE_COUNT + np.maximum(u, v)

    reflected = np.zeros(1 << M, dtype=np.uint32)
    for source, target in enumerate(SIGMA):
        reflected |= ((indices >> source) & 1) << target
    grid_symmetric = reflected == indices
    blue_line = grid_symmetric[line_masks]
    black_line = ~blue_line
    assert int(blue_line.sum()) == 2048 and int(black_line.sum()) == 1046528

    blue_codes = np.unique(line_codes[blue_line])
    black_codes = np.unique(line_codes[black_line])
    blue_code_set = set(map(int, blue_codes))
    black_code_set = set(map(int, black_codes))

    node_fibre = np.bincount(atlas, minlength=NODE_COUNT)
    node_blue_fibre = np.bincount(atlas[grid_symmetric], minlength=NODE_COUNT)
    phase = [
        "pure_blue" if blue == total else "pure_black" if blue == 0 else "mixed"
        for total, blue in zip(node_fibre, node_blue_fibre)
    ]

    unique_nodes, representative = np.unique(atlas, return_index=True)
    if not np.array_equal(unique_nodes, np.arange(NODE_COUNT)):
        raise RuntimeError("atlas ranks are not dense")
    score: list[tuple[int, ...]] = []
    c3: list[int] = []
    axis: list[int] = []
    for mask in map(int, representative):
        one_score, one_c3, one_axis = legacy_tournament_stats(mask)
        score.append(one_score)
        c3.append(one_c3)
        axis.append(one_axis)
    c3_array = np.array(c3, dtype=np.int16)

    selected_expected = {
        int(row["node_rank"]): int(row["c3"])
        for row in support_json["selected_nodes"]
    }
    assert all(c3[node] == value for node, value in selected_expected.items())
    assert Counter(phase[node] for node in support) == Counter({"pure_black": 153, "mixed": 2})

    blue_adj = edge_adjacency(blue_codes)
    black_adj = edge_adjacency(black_codes)
    all_adj = union_adjacency(blue_adj, black_adj)
    root_candidates = [node for node in range(NODE_COUNT) if c3[node] == 0]
    if len(root_candidates) != 1:
        raise RuntimeError("transitive merged node is not unique")
    root = root_candidates[0]
    balanced = {node for node in range(NODE_COUNT) if c3[node] == C3_MAX}

    unrestricted = monotone_reach(all_adj, c3, root)
    two_phase, two_phase_states = two_phase_reach(blue_adj, black_adj, c3, root)
    unreached_support = sorted(support - unrestricted)

    q_left, epsilon_left = leg_curvature(line_masks)
    q_right, epsilon_right = leg_curvature(line_mates)
    assert np.array_equal(epsilon_left, -epsilon_right)
    source_is_left = c3_array[u] < c3_array[v]
    source_q = np.where(source_is_left, q_left, q_right)

    def shell_summary(colour: str, keep: np.ndarray, codes: np.ndarray) -> dict:
        incident = keep & (in_support[u] | in_support[v])
        iu, iv = u[incident], v[incident]
        icodes = np.unique(line_codes[incident])
        endpoints = set(map(int, np.unique(np.concatenate((iu, iv)))))
        external = endpoints - support
        support_loop_instances = int(np.sum((iu == iv) & in_support[iu]))
        support_internal_instances = int(np.sum((iu != iv) & in_support[iu] & in_support[iv]))
        support_boundary_instances = int(np.sum(in_support[iu] ^ in_support[iv]))
        internal_supports = 0
        boundary_supports = 0
        loop_supports = 0
        endpoint_phase_supports: Counter[tuple[str, str]] = Counter()
        for code in map(int, icodes):
            a, b = divmod(code, NODE_COUNT)
            loop_supports += a == b
            internal_supports += a != b and in_support[a] and in_support[b]
            boundary_supports += in_support[a] ^ in_support[b]
            endpoint_phase_supports[phase_pair(phase, a, b)] += 1
        result = {
            "colour": colour,
            "literal_line_instances": int(incident.sum()),
            "literal_internal_nonloop_instances": support_internal_instances,
            "literal_loop_instances": support_loop_instances,
            "literal_boundary_instances": support_boundary_instances,
            "projected_supports": len(icodes),
            "projected_internal_supports_excluding_loops": internal_supports,
            "projected_boundary_supports": boundary_supports,
            "projected_loop_supports": loop_supports,
            "touched_support_nodes": len(endpoints & support),
            "external_neighbor_nodes": len(external),
            "all_radius_one_nodes": len(endpoints),
            "external_neighbor_phase_histogram": hist_json(Counter(phase[node] for node in external)),
            "external_neighbor_C3_histogram": hist_json(Counter(c3[node] for node in external)),
            "projected_endpoint_phase_histogram": {
                "--".join(key): value for key, value in sorted(endpoint_phase_supports.items())
            },
            "support_codes_also_in_other_colour": sum(
                int(code) in (black_code_set if colour == "blue" else blue_code_set)
                for code in icodes
            ),
        }
        if colour == "black":
            ql, qr = q_left[incident], q_right[incident]
            eps = np.abs(epsilon_left[incident])
            codes_here = line_codes[incident]
            attached_curvature = []
            q_pair_hist = Counter()
            for code, a, b, qa, qb, e in zip(codes_here, iu, iv, ql, qr, eps):
                aa, bb = int(a), int(b)
                if aa < bb:
                    attached = (int(code), int(qa), int(qb), int(e))
                elif bb < aa:
                    attached = (int(code), int(qb), int(qa), int(e))
                else:
                    attached = (int(code), min(int(qa), int(qb)), max(int(qa), int(qb)), int(e))
                attached_curvature.append(attached)
                q_pair_hist[tuple(sorted((int(qa), int(qb))))] += 1
            curvature_cells, curvature_profile = fibre_profile(attached_curvature)
            result["absolute_epsilon_histogram"] = hist_json(Counter(map(int, eps)))
            result["unordered_q_pair_histogram"] = {
                f"{a},{b}": value for (a, b), value in sorted(q_pair_hist.items())
            }
            result["node_pair_plus_attached_q_and_abs_epsilon"] = {
                "cells": curvature_cells,
                "fibre_profile": hist_json(curvature_profile),
            }
        return result

    blue_shell = shell_summary("blue", blue_line, blue_codes)
    black_shell = shell_summary("black", black_line, black_codes)

    def marked_endpoint_drift(keep: np.ndarray, stratify_q: bool) -> dict:
        q_values: list[int] = []
        deltas: list[int] = []
        left = np.flatnonzero(keep & in_support[u])
        right = np.flatnonzero(keep & in_support[v])
        q_values.extend(map(int, q_left[left]))
        deltas.extend(map(int, c3_array[v[left]] - c3_array[u[left]]))
        q_values.extend(map(int, q_right[right]))
        deltas.extend(map(int, c3_array[u[right]] - c3_array[v[right]]))
        result = {
            "marked_endpoint_half_incidences": len(deltas),
            "lower_tie_higher": {
                "lower": sum(delta < 0 for delta in deltas),
                "tie": sum(delta == 0 for delta in deltas),
                "higher": sum(delta > 0 for delta in deltas),
            },
            "sum_delta_C3": sum(deltas),
            "mean_delta_C3": fraction_json(Fraction(sum(deltas), len(deltas))),
        }
        if stratify_q:
            rows = {}
            for q in sorted(set(q_values)):
                one = [delta for qq, delta in zip(q_values, deltas) if qq == q]
                rows[str(q)] = {
                    "marked_endpoint_half_incidences": len(one),
                    "lower": sum(delta < 0 for delta in one),
                    "tie": sum(delta == 0 for delta in one),
                    "higher": sum(delta > 0 for delta in one),
                    "sum_delta_C3": sum(one),
                    "mean_delta_C3": fraction_json(Fraction(sum(one), len(one))),
                }
            result["by_endpoint_q"] = rows
        return result

    blue_marked_drift = marked_endpoint_drift(blue_line, False)
    black_marked_drift = marked_endpoint_drift(black_line, True)

    blue_incident = blue_line & (in_support[u] | in_support[v])
    blue_portals = sorted(set(map(int, u[blue_incident])) | set(map(int, v[blue_incident])))
    blue_portals = [node for node in blue_portals if node in support]
    portal_rows = []
    for node in blue_portals:
        literal_degree = int(np.sum(blue_line & ((u == node) | (v == node))))
        support_degree = len(blue_adj[node])
        external_blue = blue_adj[node] - support
        portal_rows.append({
            "node_rank": node,
            "phase": phase[node],
            "score_sequence_up_to_converse": list(score[node]),
            "C3": c3[node],
            "axis": axis[node],
            "literal_blue_line_degree": literal_degree,
            "projected_blue_degree_excluding_loops": support_degree,
            "external_blue_neighbors": len(external_blue),
            "external_blue_neighbor_phases": hist_json(Counter(phase[x] for x in external_blue)),
            "external_blue_neighbor_C3": hist_json(Counter(c3[x] for x in external_blue)),
        })

    black_distance = distances(black_adj, support)
    blue_distance = distances(blue_adj, support)
    all_distance = distances(all_adj, support)

    def distance_summary(values: list[int]) -> dict:
        return {
            "histogram": hist_json(Counter(values)),
            "maximum_finite_distance": max(value for value in values if value >= 0),
            "unreachable": sum(value < 0 for value in values),
        }

    blue_boundary_literal, _ = literal_boundary(
        u, v, blue_line, phase, c3_array, "pure_blue", "mixed"
    )
    black_boundary_literal, black_boundary_q = literal_boundary(
        u, v, black_line, phase, c3_array, "mixed", "pure_black", source_q
    )
    blue_boundary_support = projected_boundary(
        blue_codes, phase, c3, "pure_blue", "mixed"
    )
    black_boundary_support = projected_boundary(
        black_codes, phase, c3, "mixed", "pure_black"
    )
    assert black_boundary_literal == Counter({"outward": 16118, "reverse": 24370, "tie": 12584})

    black_mass = Counter()
    for node in range(NODE_COUNT):
        black_mass[phase[node]] += int(node_fibre[node] - node_blue_fibre[node])
    source_rates = {
        "outward_mixed": fraction_json(Fraction(black_boundary_literal["outward"], black_mass["mixed"])),
        "reverse_pure_black": fraction_json(Fraction(black_boundary_literal["reverse"], black_mass["pure_black"])),
    }

    incident_black = black_line & (in_support[u] | in_support[v])
    black_indices = np.flatnonzero(incident_black)
    phase_keys = [phase_pair(phase, int(u[i]), int(v[i])) for i in black_indices]
    c3_keys = [tuple(sorted((int(c3_array[u[i]]), int(c3_array[v[i]])))) for i in black_indices]
    node_keys = [int(line_codes[i]) for i in black_indices]
    curvature_keys = []
    for i in black_indices:
        a, b = int(u[i]), int(v[i])
        qa, qb = int(q_left[i]), int(q_right[i])
        if a < b:
            curvature_keys.append((int(line_codes[i]), qa, qb, abs(int(epsilon_left[i]))))
        elif b < a:
            curvature_keys.append((int(line_codes[i]), qb, qa, abs(int(epsilon_left[i]))))
        else:
            curvature_keys.append((int(line_codes[i]), min(qa, qb), max(qa, qb), abs(int(epsilon_left[i]))))
    literal_keys = list(map(int, line_masks[black_indices]))
    ta = carrier_tournament({
        "phase pair": phase_keys,
        "C3 pair": c3_keys,
        "node pair": node_keys,
        "node pair + curvature": curvature_keys,
        "literal line": literal_keys,
    })

    union_radius = {
        node for node, distance in enumerate(all_distance) if 0 <= distance <= 1
    }
    missing_radius_one = set(range(NODE_COUNT)) - union_radius

    result = {
        "schema_version": 1,
        "theorem": "THM-843",
        "inputs": {
            "atlas": str(atlas_path),
            "atlas_sha256": atlas_sha256,
            "THM834_support": str(support_path),
        },
        "global_n8": {
            "tilings": 1 << M,
            "literal_complement_lines": 1 << (M - 1),
            "merged_nodes": NODE_COUNT,
            "node_phase_histogram": hist_json(Counter(phase)),
            "node_C3_histogram": hist_json(Counter(c3)),
            "blue_literal_lines": int(blue_line.sum()),
            "black_literal_lines": int(black_line.sum()),
            "blue_projected_supports_including_loops": len(blue_codes),
            "black_projected_supports_including_loops": len(black_codes),
            "projected_supports_in_both_colours": len(blue_code_set & black_code_set),
            "blue_projected_loops": sum(int(code) // NODE_COUNT == int(code) % NODE_COUNT for code in blue_codes),
            "black_projected_loops": sum(int(code) // NODE_COUNT == int(code) % NODE_COUNT for code in black_codes),
            "transitive_node": root,
            "balanced_C3": C3_MAX,
            "balanced_nodes": len(balanced),
            "balanced_phase_histogram": hist_json(Counter(phase[node] for node in balanced)),
        },
        "monotone_transitivity_flow": {
            "unrestricted_non_decreasing_C3_reach": len(unrestricted),
            "two_phase_blue_then_black_reach": len(two_phase),
            "reach_sets_equal": unrestricted == two_phase,
            "unrestricted_balanced_reach": len(unrestricted & balanced),
            "two_phase_balanced_reach": len(two_phase & balanced),
            "all_balanced_nodes_reached": balanced <= two_phase,
            "two_phase_state_count": len(two_phase_states),
            "THM834_support_reached": len(unrestricted & support),
            "THM834_support_unreached": [
                {
                    "node_rank": node,
                    "phase": phase[node],
                    "C3": c3[node],
                    "axis": axis[node],
                }
                for node in unreached_support
            ],
            "unrestricted_unreached_C3_histogram": hist_json(Counter(c3[node] for node in set(range(NODE_COUNT)) - unrestricted)),
        },
        "global_phase_boundary_flow": {
            "blue_pure_blue_mixed_literal": dict(sorted(blue_boundary_literal.items())),
            "blue_pure_blue_mixed_projected_support": dict(sorted(blue_boundary_support.items())),
            "black_mixed_pure_black_literal": dict(sorted(black_boundary_literal.items())),
            "black_mixed_pure_black_projected_support": dict(sorted(black_boundary_support.items())),
            "black_source_normalized_rates": source_rates,
            "black_boundary_source_q_histogram": {
                f"{direction},q={q}": value
                for (direction, q), value in sorted(black_boundary_q.items())
            },
        },
        "THM834_support": {
            "nodes": len(support),
            "phase_histogram": hist_json(Counter(phase[node] for node in support)),
            "blue_shell": blue_shell,
            "black_shell": black_shell,
            "blue_marked_endpoint_drift": blue_marked_drift,
            "black_marked_endpoint_curvature_drift": black_marked_drift,
            "blue_portals": portal_rows,
            "black_distance_from_support": distance_summary(black_distance),
            "blue_distance_from_support": distance_summary(blue_distance),
            "combined_distance_from_support": distance_summary(all_distance),
            "combined_radius_one_nodes": len(union_radius),
            "combined_radius_one_external_nodes": len(union_radius - support),
            "outside_combined_radius_one": len(missing_radius_one),
            "outside_radius_one_phase_histogram": hist_json(Counter(phase[node] for node in missing_radius_one)),
            "outside_radius_one_C3_histogram": hist_json(Counter(c3[node] for node in missing_radius_one)),
            "outside_radius_one_axis_histogram": hist_json(Counter(axis[node] for node in missing_radius_one)),
            "contains_every_node_strictly_beyond_OU_center_axis_56": all(
                axis[node] >= 56 for node in missing_radius_one
            ),
            "contains_all_balanced_nodes": balanced <= union_radius,
        },
        "tournament_analysis": ta,
        "preservation": {
            "preserves": [
                "fixed-path tiling", "converse-merged node", "literal line multiplicity",
                "blue/black reflection colour", "C3", "endpoint q and epsilon",
                "THM-834 support membership",
            ],
            "destroys": [
                "unrecorded n=9 gluing", "LRC metric gaps", "owners", "walls",
                "loneliness predicate", "continued-fraction carry",
            ],
            "challenged_assumption": "the selected black line bank is not a closed subgraph",
        },
    }

    # Decisive exact assertions, deliberately stronger than the preliminary probe.
    assert result["global_n8"]["node_phase_histogram"] == {"mixed": 173, "pure_black": 3352, "pure_blue": 3}
    assert blue_shell["touched_support_nodes"] == 2
    assert blue_shell["literal_line_instances"] == 30
    assert blue_shell["projected_supports"] == 28
    assert blue_shell["external_neighbor_nodes"] == 25
    assert black_shell["literal_line_instances"] == 128652
    assert black_shell["projected_supports"] == 47872
    assert black_shell["all_radius_one_nodes"] == 3308
    assert black_shell["external_neighbor_nodes"] == 3153
    assert black_shell["literal_internal_nonloop_instances"] == 7174
    assert black_shell["literal_loop_instances"] == 88
    assert black_shell["literal_boundary_instances"] == 121390
    assert black_marked_drift["marked_endpoint_half_incidences"] == 135914
    assert black_marked_drift["sum_delta_C3"] == -58704
    assert {
        q: (
            row["marked_endpoint_half_incidences"], row["lower"], row["tie"],
            row["higher"], row["sum_delta_C3"],
        )
        for q, row in black_marked_drift["by_endpoint_q"].items()
    } == {
        "0": (13828, 4120, 1210, 8498, 3642),
        "1": (43480, 21380, 8594, 13506, -20670),
        "2": (46976, 30618, 13214, 3144, -39292),
        "3": (24954, 10946, 9424, 4584, -10420),
        "4": (5910, 436, 826, 4648, 6220),
        "5": (696, 0, 0, 696, 1536),
        "6": (70, 0, 0, 70, 280),
    }
    assert blue_marked_drift["marked_endpoint_half_incidences"] == 30
    assert blue_marked_drift["sum_delta_C3"] == 60
    assert len(union_radius) == 3310 and len(missing_radius_one) == 218
    assert len(blue_portals) == 2
    assert len(unreached_support) == 2
    assert [node for node in unreached_support] == [1008, 3512]
    assert balanced <= union_radius
    assert all(axis[node] >= 56 for node in missing_radius_one)
    assert all(c3[node] <= 14 for node in missing_radius_one)
    return result


def render(result: dict) -> str:
    g = result["global_n8"]
    flow = result["monotone_transitivity_flow"]
    boundary = result["global_phase_boundary_flow"]
    support = result["THM834_support"]
    blue = support["blue_shell"]
    black = support["black_shell"]
    ta = result["tournament_analysis"]
    lines = [
        "THM-843 N=8 METAGRAPH FLOW AND THM-834 DEFECT SHELL",
        "=" * 72,
        f"global nodes/phases={g['merged_nodes']}/{g['node_phase_histogram']}",
        f"literal lines blue/black={g['blue_literal_lines']}/{g['black_literal_lines']}",
        f"projected supports blue/black/both={g['blue_projected_supports_including_loops']}/"
        f"{g['black_projected_supports_including_loops']}/{g['projected_supports_in_both_colours']}",
        f"balanced C3={g['balanced_C3']} nodes={g['balanced_nodes']} phases={g['balanced_phase_histogram']}",
        "",
        "MONOTONE TRANSITIVITY FLOW",
        f"unrestricted/two-phase reach={flow['unrestricted_non_decreasing_C3_reach']}/"
        f"{flow['two_phase_blue_then_black_reach']} equal={flow['reach_sets_equal']}",
        f"balanced reach unrestricted/two-phase={flow['unrestricted_balanced_reach']}/"
        f"{flow['two_phase_balanced_reach']} all={flow['all_balanced_nodes_reached']}",
        f"THM-834 support reach={flow['THM834_support_reached']}/{support['nodes']} "
        f"missing={flow['THM834_support_unreached']}",
        f"unreached C3={flow['unrestricted_unreached_C3_histogram']}",
        "",
        "GLOBAL PHASE BOUNDARY FLOW (outward/reverse/tie)",
        f"blue literal={boundary['blue_pure_blue_mixed_literal']} support={boundary['blue_pure_blue_mixed_projected_support']}",
        f"black literal={boundary['black_mixed_pure_black_literal']} support={boundary['black_mixed_pure_black_projected_support']}",
        f"black source-normalized={boundary['black_source_normalized_rates']}",
        "",
        "THM-834 SUPPORT RADIUS ONE",
        f"seed nodes/phases={support['nodes']}/{support['phase_histogram']}",
        f"blue: literal/support/touched/external={blue['literal_line_instances']}/"
        f"{blue['projected_supports']}/{blue['touched_support_nodes']}/{blue['external_neighbor_nodes']}",
        f"black: literal/support/radius/external={black['literal_line_instances']}/"
        f"{black['projected_supports']}/{black['all_radius_one_nodes']}/{black['external_neighbor_nodes']}",
        f"black support internal/boundary/loops={black['projected_internal_supports_excluding_loops']}/"
        f"{black['projected_boundary_supports']}/{black['projected_loop_supports']}",
        f"colour overlap in incident supports blue/black={blue['support_codes_also_in_other_colour']}/"
        f"{black['support_codes_also_in_other_colour']}",
        f"combined radius/external/outside={support['combined_radius_one_nodes']}/"
        f"{support['combined_radius_one_external_nodes']}/{support['outside_combined_radius_one']}",
        f"combined distance={support['combined_distance_from_support']}",
        f"outside phases/C3={support['outside_radius_one_phase_histogram']}/"
        f"{support['outside_radius_one_C3_histogram']}",
        f"entire x<56 tail and balanced locus in radius one: "
        f"{support['contains_every_node_strictly_beyond_OU_center_axis_56']}/"
        f"{support['contains_all_balanced_nodes']}",
        f"blue portals={support['blue_portals']}",
        f"black curvature cells={black['node_pair_plus_attached_q_and_abs_epsilon']}",
        f"black marked q-drift={support['black_marked_endpoint_curvature_drift']['by_endpoint_q']}",
        "",
        f"TOURNAMENT ANALYSIS vertices={len(ta['vertices'])} flips={ta['edge_flips_between_gauges']}",
        f"  cells={ta['cells']}",
        f"  raw path={ta['raw_retention']['hamiltonian_path']}",
        f"  economy path={ta['retention_per_log_cell']['hamiltonian_path']}",
        "PRESERVATION: tilings/nodes/line multiplicity/colour/C3/q/epsilon/support",
        "DESTROYS: unrecorded n9 gluing and all LRC metric/owner/wall/carry data",
        "CHALLENGED ASSUMPTION: the selected black bank is not a closed subgraph",
        "ALL ASSERTIONS PASSED",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, default=Path("/tmp/n8_merged_node_rank_u16.bin"))
    parser.add_argument(
        "--support",
        type=Path,
        default=Path("05-knowledge/results/n9_false_palindrome_metagraph_flow_codex_S13b.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("05-knowledge/results/n8_metagraph_flow_defect_shell_codex_S13e.out"),
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("05-knowledge/results/n8_metagraph_flow_defect_shell_codex_S13e.json"),
    )
    args = parser.parse_args()
    result = run(args.atlas, args.support)
    text = render(result)
    args.output.write_text(text)
    args.json.write_text(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
            default=lambda value: int(value) if isinstance(value, np.integer) else value,
        )
        + "\n"
    )
    print(text, end="")


if __name__ == "__main__":
    main()
