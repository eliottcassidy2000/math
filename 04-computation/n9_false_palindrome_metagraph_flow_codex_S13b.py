#!/usr/bin/env python3
"""Tournament/metagraph coordinates of THM-828's 58 n=9 doubletons.

The input rows carry three marked n=8 faces.  Each marked face ``x`` selects
the complement line joining the merged nodes of ``x`` and ``x xor FULL``.
This audit keeps three relations rigorously separate:

* ``D``-cube adjacency: XOR by one of THM-828's four defect basis words;
* local/wiggly adjacency: one free staircase tile is flipped;
* complement-line adjacency: all free tiles are flipped, coloured blue/black
  by anti-diagonal grid symmetry.

For every face endpoint tournament we compute the score sequence, score axis,
cyclic triples, Hamiltonian paths, automorphism order, THM-810 stationary
numerator H/|Aut|, full-arc flip margins, and marked-wiggly flip margins.  The
2^21-entry merged-node atlas also gives exact local depth from the transitive
node and pure-blue/mixed/pure-black node category.

No tournament vertices are runners here.  The vertices are the eight entries
of the marked face tournament; the quotient preserves its isomorphism class,
converse orbit, score/triangle/Hamiltonian data, and both metagraph edge
channels.  It destroys the upper n=9 gluing and the defect word D; those are
retained explicitly as a separate carrier.  Challenged assumption: a
syndrome-cube edge is not presumed to be a metagraph edge, and is tested
against both actual edge layers below.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from collections import Counter, defaultdict, deque
from fractions import Fraction
from pathlib import Path

import numpy as np


N = 8
M = 21
FULL = (1 << M) - 1
NODE_COUNT = 3528
DEFECT_BASIS = (0x192486, 0x8C2C0C, 0x11B4600, 0x4483414)
EXPECTED_SECTORS = {
    0x192486: 6,
    0x8C2C0C: 2,
    0x95088A: 2,
    0x11B4600: 2,
    0x18E4E8A: 4,
    0x1976A0C: 4,
    0x4483414: 4,
    0x4511092: 2,
    0x4C41818: 26,
    0x5C67A9E: 2,
    0x5DF5E18: 4,
}


def tiles(n: int = N) -> tuple[tuple[int, int], ...]:
    return tuple((x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2)


TILES = tiles()
TILE_INDEX = {tile: i for i, tile in enumerate(TILES)}
SIGMA = tuple(TILE_INDEX[(N - y + 1, N - x + 1)] for x, y in TILES)


def reflect_mask(mask: int) -> int:
    ans = 0
    for i, j in enumerate(SIGMA):
        ans |= ((mask >> i) & 1) << j
    return ans


def adjacency(mask: int) -> tuple[int, ...]:
    """Legacy classifier convention: vertices 1..8 and path v -> v-1."""
    adj = [0] * N
    for v in range(2, N + 1):
        adj[v - 1] |= 1 << (v - 2)
    for i, (x, y) in enumerate(TILES):
        if (mask >> i) & 1:
            adj[x - 1] |= 1 << (y - 1)
        else:
            adj[y - 1] |= 1 << (x - 1)
    return tuple(adj)


def hamiltonian_paths(adj: tuple[int, ...]) -> int:
    dp = {(1 << v, v): 1 for v in range(N)}
    for size in range(1, N):
        nxt: dict[tuple[int, int], int] = defaultdict(int)
        for (used, v), count in dp.items():
            for w in range(N):
                bit = 1 << w
                if not used & bit and adj[v] & bit:
                    nxt[(used | bit, w)] += count
        dp = nxt
    return sum(dp.values())


def automorphism_order(adj: tuple[int, ...], scores: tuple[int, ...]) -> int:
    """Exact S_8 stabilizer, reduced to permutations inside score cells."""
    cells: list[tuple[int, ...]] = []
    for score in sorted(set(scores)):
        cells.append(tuple(v for v, s in enumerate(scores) if s == score))
    ans = 0
    for blocks in itertools.product(*(itertools.permutations(cell) for cell in cells)):
        p = list(range(N))
        for cell, image in zip(cells, blocks):
            for v, w in zip(cell, image):
                p[v] = w
        if all(bool(adj[i] & (1 << j)) == bool(adj[p[i]] & (1 << p[j]))
               for i in range(N) for j in range(i + 1, N)):
            ans += 1
    return ans


def delta_histogram(adj: tuple[int, ...], allowed: set[tuple[int, int]] | None = None) -> Counter[int]:
    scores = tuple(x.bit_count() for x in adj)
    hist: Counter[int] = Counter()
    for i in range(N):
        for j in range(i + 1, N):
            if allowed is not None and (i, j) not in allowed:
                continue
            if adj[i] & (1 << j):
                winner, loser = i, j
            else:
                winner, loser = j, i
            # THM-810 F3: 4(d_l-d_w)+8 = 8-8(s_w-s_l).
            hist[8 - 8 * (scores[winner] - scores[loser])] += 1
    return hist


FREE_ARCS = {(min(x - 1, y - 1), max(x - 1, y - 1)) for x, y in TILES}


def normalized_score(scores: tuple[int, ...]) -> tuple[int, ...]:
    direct = tuple(sorted(scores))
    converse = tuple(sorted(N - 1 - s for s in scores))
    return min(direct, converse)


def tournament_invariant(mask: int) -> dict:
    adj = adjacency(mask)
    raw_scores = tuple(x.bit_count() for x in adj)
    score = normalized_score(raw_scores)
    axis = sum((2 * s - (N - 1)) ** 2 for s in raw_scores)
    c3 = (N * (N - 1) * (N - 2) // 6) - sum(s * (s - 1) // 2 for s in raw_scores)
    h = hamiltonian_paths(adj)
    aut = automorphism_order(adj, raw_scores)
    assert h % aut == 0
    full_delta = delta_histogram(adj)
    marked_delta = delta_histogram(adj, FREE_ARCS)
    assert sum(full_delta.values()) == 28 and sum(marked_delta.values()) == 21
    assert axis == 168 - 8 * c3
    return {
        "score": score,
        "axis": axis,
        "c3": c3,
        "H": h,
        "aut": aut,
        "class_stationary_numerator": h // aut,
        "arc_delta": full_delta,
        "wiggly_delta": marked_delta,
    }


def quotient_edges(atlas: np.ndarray) -> set[int]:
    idx = np.arange(1 << M, dtype=np.uint32)
    chunks = []
    for bit in range(M):
        lo = idx[(idx & (1 << bit)) == 0]
        a, b = atlas[lo], atlas[lo | (1 << bit)]
        chunks.append(np.unique(np.minimum(a, b) * NODE_COUNT + np.maximum(a, b)))
    return set(map(int, np.unique(np.concatenate(chunks))))


def edge_adjacency(codes: set[int], include_loops: bool = False) -> list[set[int]]:
    adj = [set() for _ in range(NODE_COUNT)]
    for code in codes:
        u, v = divmod(code, NODE_COUNT)
        if u == v and not include_loops:
            continue
        adj[u].add(v)
        adj[v].add(u)
    return adj


def bfs(adj: list[set[int]], starts: set[int]) -> list[int]:
    dist = [-1] * len(adj)
    q = deque()
    for v in starts:
        dist[v] = 0
        q.append(v)
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] < 0:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def fraction_json(x: Fraction) -> dict:
    return {"numerator": x.numerator, "denominator": x.denominator, "decimal": float(x)}


def hist_json(hist: Counter) -> dict[str, int]:
    return {str(k): hist[k] for k in sorted(hist)}


def summarize_signs(hist: Counter[int]) -> dict[str, int]:
    return {
        "climb": sum(v for k, v in hist.items() if k > 0),
        "level": hist[0],
        "descend": sum(v for k, v in hist.items() if k < 0),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, default=Path("/tmp/n8_merged_node_rank_u16.bin"))
    parser.add_argument("--witnesses", type=Path,
                        default=Path("05-knowledge/results/mobius_cech_n9_exact_join_witnesses_codex_S13.tsv"))
    parser.add_argument("--json", type=Path,
                        default=Path("05-knowledge/results/n9_false_palindrome_metagraph_flow_codex_S13b.json"))
    parser.add_argument("--output", type=Path,
                        default=Path("05-knowledge/results/n9_false_palindrome_metagraph_flow_codex_S13b.out"))
    args = parser.parse_args()

    atlas = np.fromfile(args.atlas, dtype="<u2").astype(np.uint32)
    if len(atlas) != 1 << M or int(atlas.max()) != NODE_COUNT - 1:
        raise RuntimeError("unexpected n=8 merged-node atlas")

    rows = []
    with args.witnesses.open() as f:
        for raw in csv.DictReader(f, delimiter="\t"):
            row = {key: int(raw[key], 16) for key in ("D", "u", "v", "xA", "xB", "xC", "dA", "dB", "dC")}
            rows.append(row)
    sector_rows: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        sector_rows[row["D"]].append(row)
    if len(rows) != 58 or {d: len(rs) for d, rs in sector_rows.items()} != EXPECTED_SECTORS:
        raise RuntimeError("THM-828 witness census mismatch")

    # Exact global n=8 quotient structure.
    local_codes = quotient_edges(atlas)
    local_adj = edge_adjacency(local_codes)
    root = int(atlas[0])
    local_depth = bfs(local_adj, {root})
    if min(local_depth) < 0:
        raise RuntimeError("local quotient is disconnected")

    idx = np.arange(1 << M, dtype=np.uint32)
    mate = idx ^ FULL
    side = idx < mate
    sigma_image = np.zeros(1 << M, dtype=np.uint32)
    for source, target in enumerate(SIGMA):
        sigma_image |= ((idx >> source) & 1) << target
    grid_symmetric = sigma_image == idx
    line_code = np.minimum(atlas[side], atlas[mate[side]]) * NODE_COUNT + np.maximum(atlas[side], atlas[mate[side]])
    blue_codes = set(map(int, np.unique(line_code[grid_symmetric[side]])))
    black_codes = set(map(int, np.unique(line_code[~grid_symmetric[side]])))
    blue_adj, black_adj = edge_adjacency(blue_codes), edge_adjacency(black_codes)

    node_fibre = np.bincount(atlas, minlength=NODE_COUNT)
    node_blue_fibre = np.bincount(atlas[grid_symmetric], minlength=NODE_COUNT)
    node_category = []
    for total, blue in zip(node_fibre, node_blue_fibre):
        node_category.append("pure_blue" if blue == total else "pure_black" if blue == 0 else "mixed")

    # Node invariants are evaluated only on the 155 nodes selected by the 58 rows.
    masks = []
    for row in rows:
        for role in ("A", "B", "C"):
            x = row["x" + role]
            masks.extend((x, x ^ FULL))
    invariant_by_node: dict[int, dict] = {}
    marked_by_mask: dict[int, Counter[int]] = {}
    for mask in masks:
        node = int(atlas[mask])
        inv = tournament_invariant(mask)
        marked_by_mask[mask] = inv.pop("wiggly_delta")
        if node in invariant_by_node:
            old = invariant_by_node[node]
            for key in ("score", "axis", "c3", "H", "aut", "class_stationary_numerator", "arc_delta"):
                if old[key] != inv[key]:
                    raise RuntimeError(f"merged-node invariant mismatch at node {node}: {key}")
        else:
            invariant_by_node[node] = inv
    if len(set(masks)) != 348 or len(invariant_by_node) != 155:
        raise RuntimeError("face endpoint support mismatch")

    # F2 identity after converse merging.
    for node, inv in invariant_by_node.items():
        self_converse = bool(node_fibre[node] & 1)
        expected = inv["class_stationary_numerator"] * (1 if self_converse else 2)
        if int(node_fibre[node]) != expected:
            raise RuntimeError(f"stationary fibre identity failed at node {node}")
        inv["self_converse"] = self_converse
        inv["merged_stationary_numerator"] = int(node_fibre[node])
        inv["category"] = node_category[node]
        inv["local_depth"] = local_depth[node]

    sector_nodes: dict[int, set[int]] = defaultdict(set)
    sector_lines: dict[int, set[int]] = defaultdict(set)
    sector_line_occurrences: dict[int, list[int]] = defaultdict(list)
    sector_endpoint_occurrences: dict[int, list[tuple[int, int]]] = defaultdict(list)
    all_line_occurrences: list[int] = []
    for d, rs in sector_rows.items():
        for row in rs:
            for role in ("A", "B", "C"):
                x = row["x" + role]
                u, v = int(atlas[x]), int(atlas[x ^ FULL])
                code = min(u, v) * NODE_COUNT + max(u, v)
                if code not in black_codes or code in blue_codes:
                    raise RuntimeError("selected face line is not exclusively black")
                sector_lines[d].add(code)
                sector_line_occurrences[d].append(code)
                all_line_occurrences.append(code)
                for mask in (x, x ^ FULL):
                    node = int(atlas[mask])
                    sector_nodes[d].add(node)
                    sector_endpoint_occurrences[d].append((node, mask))

    line_mult = Counter(all_line_occurrences)
    if len(line_mult) != 87 or Counter(line_mult.values()) != Counter({1: 58, 4: 29}):
        raise RuntimeError("black-line recurrence mismatch")

    # Defect basis coordinates.
    coordinate_of = {}
    for z in range(16):
        d = 0
        for bit, basis in enumerate(DEFECT_BASIS):
            if (z >> bit) & 1:
                d ^= basis
        coordinate_of[d] = z
    if any(d not in coordinate_of for d in sector_rows):
        raise RuntimeError("sector outside rank-four defect span")

    sector_json = []
    sector_coordinate: dict[int, tuple[Fraction, Fraction, Fraction, int]] = {}
    for d in sorted(sector_rows):
        endpoint_occ = sector_endpoint_occurrences[d]
        line_occ = sector_line_occurrences[d]
        axis_values = [invariant_by_node[node]["axis"] for node, _ in endpoint_occ]
        depth_values = [invariant_by_node[node]["local_depth"] for node, _ in endpoint_occ]
        stat_values = [invariant_by_node[node]["class_stationary_numerator"] for node, _ in endpoint_occ]
        line_delta = Counter()
        line_types = Counter()
        for code in line_occ:
            u, v = divmod(code, NODE_COUNT)
            iu, iv = invariant_by_node[u], invariant_by_node[v]
            line_delta[abs(iu["axis"] - iv["axis"])] += 1
            line_types[tuple(sorted((iu["category"], iv["category"]))) ] += 1
        arc_delta = Counter()
        marked_delta = Counter()
        for node, mask in endpoint_occ:
            arc_delta.update(invariant_by_node[node]["arc_delta"])
            marked_delta.update(marked_by_mask[mask])
        axis_mean = Fraction(sum(axis_values), len(axis_values))
        depth_mean = Fraction(sum(depth_values), len(depth_values))
        stat_mean = Fraction(sum(stat_values), len(stat_values))
        sector_coordinate[d] = (axis_mean, depth_mean, stat_mean, coordinate_of[d])
        sector_json.append({
            "D": f"0x{d:x}",
            "cube_coordinate_lsb_first": format(coordinate_of[d], "04b")[::-1],
            "witnesses": len(sector_rows[d]),
            "support_nodes": len(sector_nodes[d]),
            "selected_black_lines": len(sector_lines[d]),
            "endpoint_axis_mean": fraction_json(axis_mean),
            "endpoint_axis_range": [min(axis_values), max(axis_values)],
            "endpoint_C3_mean": fraction_json(Fraction(sum(168 - x for x in axis_values), 8 * len(axis_values))),
            "endpoint_depth_mean": fraction_json(depth_mean),
            "endpoint_depth_range": [min(depth_values), max(depth_values)],
            "class_stationary_numerator_mean": fraction_json(stat_mean),
            "node_category_occurrences": dict(sorted(Counter(invariant_by_node[node]["category"] for node, _ in endpoint_occ).items())),
            "black_line_axis_transport_histogram": hist_json(line_delta),
            "black_line_endpoint_category_histogram": {"--".join(k): v for k, v in sorted(line_types.items())},
            "full_arc_flip_signs": summarize_signs(arc_delta),
            "marked_wiggly_flip_signs": summarize_signs(marked_delta),
        })

    # Exact transitivity-to-distribution order of sector barycentres.
    ordered_sectors = sorted(sector_rows, key=lambda d: (
        -sector_coordinate[d][0], sector_coordinate[d][1], sector_coordinate[d][2], coordinate_of[d]
    ))
    order_rank = {d: i for i, d in enumerate(ordered_sectors)}
    for item in sector_json:
        item["transitivity_order_rank"] = order_rank[int(item["D"], 16)]

    # Pairwise audit: D-cube edges versus actual metagraph relations.
    local_pair_codes = {code for code in local_codes if code // NODE_COUNT != code % NODE_COUNT}
    blue_pair_codes = {code for code in blue_codes if code // NODE_COUNT != code % NODE_COUNT}
    black_pair_codes = {code for code in black_codes if code // NODE_COUNT != code % NODE_COUNT}

    def bridge_count(a: set[int], b: set[int], codes: set[int]) -> int:
        a_only, b_only = a - b, b - a
        return sum(1 for code in codes
                   if ((code // NODE_COUNT in a_only and code % NODE_COUNT in b_only)
                       or (code % NODE_COUNT in a_only and code // NODE_COUNT in b_only)))

    pair_rows = []
    cube_summary = Counter()
    noncube_summary = Counter()
    cube_axis_gaps: list[Fraction] = []
    noncube_axis_gaps: list[Fraction] = []
    for da, db in itertools.combinations(sorted(sector_rows), 2):
        cube = (coordinate_of[da] ^ coordinate_of[db]).bit_count() == 1
        common_nodes = len(sector_nodes[da] & sector_nodes[db])
        common_lines = len(sector_lines[da] & sector_lines[db])
        local_bridges = bridge_count(sector_nodes[da], sector_nodes[db], local_pair_codes)
        blue_bridges = bridge_count(sector_nodes[da], sector_nodes[db], blue_pair_codes)
        black_bridges = bridge_count(sector_nodes[da], sector_nodes[db], black_pair_codes)
        distance = min(bfs(local_adj, sector_nodes[da])[v] for v in sector_nodes[db])
        gap = abs(sector_coordinate[da][0] - sector_coordinate[db][0])
        rec = {
            "D1": f"0x{da:x}", "D2": f"0x{db:x}", "cube_adjacent": cube,
            "common_nodes": common_nodes, "common_selected_black_lines": common_lines,
            "local_bridges_excluding_overlap": local_bridges,
            "blue_bridges_excluding_overlap": blue_bridges,
            "black_bridges_excluding_overlap": black_bridges,
            "minimum_local_node_distance": distance,
            "axis_barycentre_gap": fraction_json(gap),
            "transitivity_order_gap": abs(order_rank[da] - order_rank[db]),
        }
        pair_rows.append(rec)
        summary = cube_summary if cube else noncube_summary
        summary["pairs"] += 1
        summary["node_overlap_pairs"] += common_nodes > 0
        summary["selected_line_overlap_pairs"] += common_lines > 0
        summary["local_bridge_pairs"] += local_bridges > 0
        summary["blue_bridge_pairs"] += blue_bridges > 0
        summary["black_bridge_pairs"] += black_bridges > 0
        summary[f"distance_{distance}"] += 1
        (cube_axis_gaps if cube else noncube_axis_gaps).append(gap)

    # Cube adjacency has 14 edges by THM-828, but test its relation to the
    # transitivity order without promoting the order to an edge structure.
    cube_order_gaps = [r["transitivity_order_gap"] for r in pair_rows if r["cube_adjacent"]]
    if cube_summary["pairs"] != 14 or noncube_summary["pairs"] != 41:
        raise RuntimeError("occupied-cube pair census mismatch")

    global_endpoint_nodes = set().union(*sector_nodes.values())
    global_node_category = Counter(node_category[v] for v in global_endpoint_nodes)
    global_line_delta = Counter()
    global_line_type = Counter()
    for code in all_line_occurrences:
        u, v = divmod(code, NODE_COUNT)
        global_line_delta[abs(invariant_by_node[u]["axis"] - invariant_by_node[v]["axis"])] += 1
        global_line_type[tuple(sorted((node_category[u], node_category[v])))] += 1

    global_arc_delta = Counter()
    global_marked_delta = Counter()
    for d in sector_rows:
        for node, mask in sector_endpoint_occurrences[d]:
            global_arc_delta.update(invariant_by_node[node]["arc_delta"])
            global_marked_delta.update(marked_by_mask[mask])

    # Orient occupied-cube edges by decreasing exact axis barycentre and test
    # whether this scalar flow radiates from its most-transitive sector.
    cube_adj: dict[int, set[int]] = {d: set() for d in sector_rows}
    for da, db in itertools.combinations(sector_rows, 2):
        if (coordinate_of[da] ^ coordinate_of[db]).bit_count() == 1:
            cube_adj[da].add(db)
            cube_adj[db].add(da)
    transitive_root = max(sector_rows, key=lambda d: (sector_coordinate[d][0], -coordinate_of[d]))
    cube_distance = {transitive_root: 0}
    queue = deque([transitive_root])
    while queue:
        u = queue.popleft()
        for v in cube_adj[u]:
            if v not in cube_distance:
                cube_distance[v] = cube_distance[u] + 1
                queue.append(v)
    directed_cube: dict[int, set[int]] = {d: set() for d in sector_rows}
    radial_change = Counter()
    for da, db in itertools.combinations(sector_rows, 2):
        if db not in cube_adj[da]:
            continue
        high, low = (da, db) if sector_coordinate[da][0] > sector_coordinate[db][0] else (db, da)
        directed_cube[high].add(low)
        radial_change[cube_distance[low] - cube_distance[high]] += 1
    directed_reach = {transitive_root}
    queue = deque([transitive_root])
    while queue:
        u = queue.popleft()
        for v in directed_cube[u]:
            if v not in directed_reach:
                directed_reach.add(v)
                queue.append(v)
    cube_local_maxima = sorted(d for d in sector_rows
                               if all(sector_coordinate[d][0] > sector_coordinate[e][0] for e in cube_adj[d]))
    cube_local_minima = sorted(d for d in sector_rows
                               if all(sector_coordinate[d][0] < sector_coordinate[e][0] for e in cube_adj[d]))

    node_rows = []
    for node in sorted(global_endpoint_nodes):
        inv = invariant_by_node[node]
        node_rows.append({
            "node_rank": node,
            "score_sequence_up_to_converse": list(inv["score"]),
            "score_histogram": {str(k): v for k, v in sorted(Counter(inv["score"]).items())},
            "axis": inv["axis"],
            "c3": inv["c3"],
            "H": inv["H"],
            "aut": inv["aut"],
            "ordinary_class_stationary_numerator_H_over_aut": inv["class_stationary_numerator"],
            "ordinary_class_stationary_probability": {
                "numerator": inv["class_stationary_numerator"], "denominator": 1 << M,
            },
            "merged_stationary_numerator": inv["merged_stationary_numerator"],
            "merged_stationary_probability": {
                "numerator": inv["merged_stationary_numerator"], "denominator": 1 << M,
            },
            "self_converse": inv["self_converse"],
            "blueblack_category": inv["category"],
            "local_depth_from_transitive": inv["local_depth"],
            "full_arc_delta_x_histogram": hist_json(inv["arc_delta"]),
            "full_arc_flip_signs": summarize_signs(inv["arc_delta"]),
        })

    face_line_rows = []
    for witness_index, row in enumerate(rows):
        for role in ("A", "B", "C"):
            x, y = row["x" + role], row["x" + role] ^ FULL
            u, v = int(atlas[x]), int(atlas[y])
            code = min(u, v) * NODE_COUNT + max(u, v)
            face_line_rows.append({
                "witness_index": witness_index,
                "D": f"0x{row['D']:x}",
                "face_role": role,
                "endpoint_masks": [f"0x{x:x}", f"0x{y:x}"],
                "endpoint_node_ranks": [u, v],
                "endpoint_axis": [invariant_by_node[u]["axis"], invariant_by_node[v]["axis"]],
                "endpoint_local_depth": [local_depth[u], local_depth[v]],
                "endpoint_categories": [node_category[u], node_category[v]],
                "line_color": "black",
                "line_is_exclusively_black_at_node_pair_level": code not in blue_codes,
                "absolute_axis_transport": abs(invariant_by_node[u]["axis"] - invariant_by_node[v]["axis"]),
                "marked_wiggly_delta_x_histograms": [hist_json(marked_by_mask[x]), hist_json(marked_by_mask[y])],
                "marked_wiggly_flip_signs": [summarize_signs(marked_by_mask[x]), summarize_signs(marked_by_mask[y])],
                "global_selected_line_multiplicity": line_mult[code],
            })

    result = {
        "schema_version": 1,
        "theorems": ["THM-810", "THM-828"],
        "inputs": {"atlas": str(args.atlas), "witnesses": str(args.witnesses)},
        "global_n8": {
            "merged_nodes": NODE_COUNT,
            "local_edge_pairs_including_loops": len(local_codes),
            "blue_edge_pairs_including_loops": len(blue_codes),
            "black_edge_pairs_including_loops": len(black_codes),
            "local_depth_maximum": max(local_depth),
        },
        "face_support": {
            "endpoint_tiling_occurrences": 348,
            "distinct_endpoint_tilings": len(set(masks)),
            "distinct_merged_nodes": len(global_endpoint_nodes),
            "node_categories": dict(sorted(global_node_category.items())),
            "selected_black_line_occurrences": len(all_line_occurrences),
            "distinct_selected_black_lines": len(line_mult),
            "selected_black_line_multiplicity_histogram": hist_json(Counter(line_mult.values())),
            "selected_black_line_axis_transport_histogram": hist_json(global_line_delta),
            "selected_black_line_endpoint_category_histogram": {"--".join(k): v for k, v in sorted(global_line_type.items())},
            "all_selected_node_pairs_exclusively_black": True,
            "full_arc_flip_signs_over_348_endpoint_occurrences": summarize_signs(global_arc_delta),
            "marked_wiggly_flip_signs_over_348_endpoint_occurrences": summarize_signs(global_marked_delta),
        },
        "sector_order_definition": [
            "decreasing endpoint-axis barycentre",
            "increasing endpoint local-depth barycentre",
            "increasing class stationary-numerator barycentre H/|Aut|",
            "defect-basis coordinate tie-break",
        ],
        "sector_order": [f"0x{d:x}" for d in ordered_sectors],
        "sectors": sorted(sector_json, key=lambda x: x["transitivity_order_rank"]),
        "selected_nodes": node_rows,
        "witness_face_lines": face_line_rows,
        "cube_vs_metagraph": {
            "cube_adjacent": dict(sorted(cube_summary.items())),
            "not_cube_adjacent": dict(sorted(noncube_summary.items())),
            "mean_axis_barycentre_gap_cube": fraction_json(sum(cube_axis_gaps, Fraction()) / len(cube_axis_gaps)),
            "mean_axis_barycentre_gap_noncube": fraction_json(sum(noncube_axis_gaps, Fraction()) / len(noncube_axis_gaps)),
            "cube_transitivity_order_gap_histogram": hist_json(Counter(cube_order_gaps)),
            "axis_oriented_cube_flow": {
                "most_transitive_sector": f"0x{transitive_root:x}",
                "distance_from_root": {f"0x{d:x}": cube_distance[d] for d in sorted(cube_distance)},
                "radial_distance_change_histogram": hist_json(radial_change),
                "directed_reach_from_root": len(directed_reach),
                "unreached_sectors": [f"0x{d:x}" for d in sorted(set(sector_rows) - directed_reach)],
                "local_maxima": [f"0x{d:x}" for d in cube_local_maxima],
                "local_minima": [f"0x{d:x}" for d in cube_local_minima],
            },
            "pair_table": pair_rows,
        },
    }

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n")

    lines = [
        "THM-810 x THM-828: N=9 FALSE-PALINDROME METAGRAPH FLOW",
        "exact on 58 witnesses, 348 distinct n=8 endpoint tilings, and the full 2^21 atlas",
        "",
        f"global n=8 quotient: nodes={NODE_COUNT} local-pairs={len(local_codes)} "
        f"blue-pairs={len(blue_codes)} black-pairs={len(black_codes)} local-depth-max={max(local_depth)}",
        f"face support: merged-nodes={len(global_endpoint_nodes)} categories={dict(sorted(global_node_category.items()))}",
        f"selected lines: occurrences={len(all_line_occurrences)} distinct={len(line_mult)} "
        f"multiplicity={dict(sorted(Counter(line_mult.values()).items()))}",
        f"selected black-line |Delta x|={dict(sorted(global_line_delta.items()))}",
        f"selected line endpoint types={dict(sorted(global_line_type.items()))}",
        f"full-arc flip signs={summarize_signs(global_arc_delta)}; marked-wiggly signs={summarize_signs(global_marked_delta)}",
        "",
        "OBJECTIVE SECTOR ORDER",
        "rank D cube witnesses nodes lines xbar depthbar H/Aut-bar categories K-|Dx|",
    ]
    sector_lookup = {int(x["D"], 16): x for x in sector_json}
    for rank, d in enumerate(ordered_sectors):
        x = sector_lookup[d]
        fm = lambda q: f"{q['numerator']}/{q['denominator']}" if q["denominator"] != 1 else str(q["numerator"])
        lines.append(
            f"{rank:2d} {d:07x} {x['cube_coordinate_lsb_first']} {x['witnesses']:2d} "
            f"{x['support_nodes']:3d} {x['selected_black_lines']:2d} "
            f"{fm(x['endpoint_axis_mean']):>7s} {fm(x['endpoint_depth_mean']):>7s} "
            f"{fm(x['class_stationary_numerator_mean']):>9s} "
            f"{x['node_category_occurrences']} {x['black_line_axis_transport_histogram']}"
        )
    lines += [
        "",
        "DEFECT-CUBE ADJACENCY VERSUS ACTUAL METAGRAPH RELATIONS",
        f"cube-adjacent: {dict(sorted(cube_summary.items()))}",
        f"not-cube-adjacent: {dict(sorted(noncube_summary.items()))}",
        f"mean |axis barycentre gap| cube={fraction_json(sum(cube_axis_gaps, Fraction()) / len(cube_axis_gaps))} "
        f"noncube={fraction_json(sum(noncube_axis_gaps, Fraction()) / len(noncube_axis_gaps))}",
        f"cube edge gaps in sector order={dict(sorted(Counter(cube_order_gaps).items()))}",
        f"axis-oriented cube root={transitive_root:07x} radial changes={dict(sorted(radial_change.items()))} "
        f"directed reach={len(directed_reach)}/11",
        f"cube local maxima={[f'{d:07x}' for d in cube_local_maxima]} "
        f"local minima={[f'{d:07x}' for d in cube_local_minima]}",
        "",
        "INTERPRETATION",
        "Every witness selects three exclusively black complement lines; its reflected mate selects the same lines.",
        "The 11 D sectors are therefore carriers of finite sets of black lines, not metagraph nodes.",
        "Cube adjacency is compared with, but never identified with, local, blue, or black metagraph adjacency.",
    ]
    args.output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
