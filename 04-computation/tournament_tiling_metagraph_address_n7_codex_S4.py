#!/usr/bin/env python3
"""Exact n=7 continuation of the HYP-6825 metagraph address audit.

The explorer stops at n=6 because naive canonicalization is expensive.  This
script uses klein-S161's exact refinement strategy: canonically colour vertices
by iterated in/out-neighbour profiles, then minimize over every permutation
inside the resulting colour cells.  It then applies the same two-channel
address audit as the n<=6 script.

The committed JSON is compact but bidirectional: arrays indexed by tiling mask
give the canonical class, merged node, global index, and fibre-local index;
each node lists its complete inverse mask fibre.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from merged_metagraph_lines_n7_klein_S161 import canon_tournament
from tournament_tiling_metagraph_address_codex_S4 import (
    carrier_tournament,
    deletion_profile,
    edge_type_local,
    incidence_profiles,
    lex_shortest_paths,
    partition_stats,
    stable_refinement,
    tile_schema,
    tiling_tournament,
)


def adjacency_from_tiling(mask: int, n: int, tiles: tuple[tuple[int, int], ...]) -> list[list[int]]:
    adjacency = [[0] * n for _ in range(n)]
    for label in range(n, 1, -1):
        adjacency[label - 1][label - 2] = 1
    for bit, (x, y) in enumerate(tiles):
        if (mask >> bit) & 1:
            adjacency[y - 1][x - 1] = 1
        else:
            adjacency[x - 1][y - 1] = 1
    return adjacency


def exact_classes() -> tuple[list[int], list[int], tuple[int, ...], tuple[tuple[int, int], ...]]:
    n = 7
    tiles, sigma = tile_schema(n)
    m = len(tiles)
    nt = 1 << m
    codes = []
    trans_masks = []
    for mask in range(nt):
        key = canon_tournament(adjacency_from_tiling(mask, n, tiles), n)
        codes.append(int.from_bytes(key, "big"))
        trans = 0
        for bit, target in enumerate(sigma):
            trans |= ((mask >> bit) & 1) << target
        trans_masks.append(trans)
    return codes, trans_masks, sigma, tiles


def analyze(parent_json: Path) -> dict:
    n = 7
    canon_arr, trans_masks, sigma, tiles = exact_classes()
    m = len(tiles)
    nt = 1 << m
    full = nt - 1
    converse_arr = [canon_arr[trans_masks[t]] for t in range(nt)]

    classes = sorted(set(canon_arr))
    assert len(classes) == 456
    class_representative = {}
    for t, cls in enumerate(canon_arr):
        class_representative.setdefault(cls, t)
    converse_by_class: dict[int, int] = {}
    for cls in classes:
        vals = {converse_arr[t] for t in range(nt) if canon_arr[t] == cls}
        assert len(vals) == 1
        converse_by_class[cls] = next(iter(vals))
    assert all(converse_by_class[converse_by_class[c]] == c for c in converse_by_class)

    merged_arr = [min(canon_arr[t], converse_arr[t]) for t in range(nt)]
    node_codes = sorted(set(merged_arr))
    assert len(node_codes) == 272
    root = int(merged_arr[0])
    node_sc = {u: converse_by_class[u] == u for u in node_codes}
    node_classes = {
        u: tuple(sorted(c for c in converse_by_class if min(c, converse_by_class[c]) == u))
        for u in node_codes
    }
    node_fibres: dict[int, list[int]] = defaultdict(list)
    for t, u in enumerate(merged_arr):
        node_fibres[u].append(t)

    parent_data = json.loads(parent_json.read_text())
    n6 = next(result for result in parent_data["sizes"] if result["n"] == 6)
    parent_rank = {
        int(code, 16): node["rank"]
        for node in n6["nodes"]
        for code in node["canonical_orbit_codes"]
    }
    class_parent = {}
    for cls, representative in class_representative.items():
        pairmask = tiling_tournament(representative, n, tiles)
        profile = deletion_profile(pairmask, n)
        assert all(code in parent_rank for code, _ in profile)
        class_parent[cls] = tuple((parent_rank[code], multiplicity) for code, multiplicity in profile)
    node_parent = {u: tuple(sorted(class_parent[c] for c in node_classes[u])) for u in node_codes}

    edges: dict[tuple[int, int], Counter[str]] = defaultdict(Counter)
    local_edges = 0
    for t in range(nt):
        u = int(merged_arr[t])
        for bit in range(m):
            if (t >> bit) & 1:
                continue
            v = int(merged_arr[t ^ (1 << bit)])
            edges[(min(u, v), max(u, v))][edge_type_local(node_sc[u], node_sc[v])] += 1
            local_edges += 1
    assert local_edges == m * (1 << (m - 1))

    grid_symmetric = [
        all(((t >> bit) & 1) == ((t >> target) & 1) for bit, target in enumerate(sigma))
        for t in range(nt)
    ]
    blue = black = 0
    for t in range(nt):
        other = t ^ full
        if t > other:
            continue
        u, v = int(merged_arr[t]), int(merged_arr[other])
        channel = "L_B" if grid_symmetric[t] else "L_K"
        edges[(min(u, v), max(u, v))][channel] += 1
        blue += channel == "L_B"
        black += channel == "L_K"
    assert blue == 256 and black == 16128

    local_channels = {"F_S", "F_R", "F_G"}
    line_channels = {"L_B", "L_K"}
    all_channels = local_channels | line_channels
    line_incidence = incidence_profiles(node_codes, edges, line_channels)
    local_incidence = incidence_profiles(node_codes, edges, local_channels)
    line_initial = {u: (int(u == root), int(node_sc[u]), line_incidence[u]) for u in node_codes}
    line_wl, line_rounds = stable_refinement(node_codes, line_initial, edges, line_channels)
    combined_initial = {
        u: (
            int(u == root),
            int(node_sc[u]),
            len(node_fibres[u]),
            node_parent[u],
            line_incidence[u],
            local_incidence[u],
        )
        for u in node_codes
    }
    combined_wl, combined_rounds = stable_refinement(node_codes, combined_initial, edges, all_channels)

    local_dist, local_word, local_components = lex_shortest_paths(
        node_codes,
        root,
        edges,
        local_channels,
        {"F_S": "S", "F_R": "R", "F_G": "G"},
        {"S": 0, "R": 1, "G": 2},
    )
    line_dist, line_word, line_components = lex_shortest_paths(
        node_codes, root, edges, line_channels, {"L_B": "B", "L_K": "K"}
    )
    fibre_depth = {u: min(t.bit_count() for t in node_fibres[u]) for u in node_codes}

    structural = {
        u: (
            local_dist[u],
            tuple({"S": 0, "R": 1, "G": 2}[letter] for letter in local_word[u]),
            (line_dist[u], line_word[u]) if u in line_dist else (10**9, "~"),
            line_wl[u],
            combined_wl[u],
            node_parent[u],
        )
        for u in node_codes
    }
    exact = {u: structural[u] + (u,) for u in node_codes}
    ordered_nodes = sorted(node_codes, key=lambda u: exact[u])
    node_rank = {u: rank for rank, u in enumerate(ordered_nodes)}

    global_tilings = sorted(
        range(nt), key=lambda t: (t.bit_count(), node_rank[int(merged_arr[t])], int(canon_arr[t]), t)
    )
    global_index = [0] * nt
    for i, t in enumerate(global_tilings):
        global_index[t] = i
    fibre_index = [0] * nt
    for u in node_codes:
        for i, t in enumerate(sorted(node_fibres[u], key=lambda t: (t.bit_count(), int(canon_arr[t]), t))):
            fibre_index[t] = i

    carriers = {
        "raw_sc_type": {u: int(node_sc[u]) for u in node_codes},
        "radial_depth": {u: local_dist[u] for u in node_codes},
        "blueblack_incidence": {u: line_incidence[u] for u in node_codes},
        "rooted_blueblack_wl": {u: line_wl[u] for u in node_codes},
        "local_parent_atlas": {u: (local_dist[u], local_word[u], node_parent[u]) for u in node_codes},
        "combined_colored_wl": {u: combined_wl[u] for u in node_codes},
        "structural_address": structural,
        "exact_address": exact,
    }
    carrier_stats = {name: partition_stats(values) for name, values in carriers.items()}
    retention = carrier_tournament(carrier_stats, "retention")
    economy = carrier_tournament(carrier_stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carriers))
        for j in range(i + 1, len(carriers))
    )

    assert local_components == 1
    assert all(local_dist[u] <= fibre_depth[u] for u in node_codes)
    assert all(node_sc[u] == (len(node_fibres[u]) % 2 == 1) for u in node_codes)
    assert sorted(global_index) == list(range(nt))

    node_records = []
    for u in ordered_nodes:
        fibre = sorted(node_fibres[u], key=lambda t: fibre_index[t])
        node_records.append(
            {
                "rank": node_rank[u],
                "id": f"n7-a{node_rank[u]:03d}",
                "canonical_orbit_codes": [hex(c) for c in node_classes[u]],
                "self_converse": node_sc[u],
                "local_depth": local_dist[u],
                "minimum_labelled_tiling_weight": fibre_depth[u],
                "local_path_word": local_word[u],
                "blueblack_root_distance": line_dist.get(u),
                "blueblack_root_word": line_word.get(u),
                "blueblack_wl_color": line_wl[u],
                "combined_wl_color": combined_wl[u],
                "recursive_parent_address": repr(node_parent[u]),
                "line_incidence": [list(item) for item in line_incidence[u]],
                "tiling_count": len(fibre),
                "tiling_masks": fibre,
            }
        )

    return {
        "schema_version": 1,
        "n": 7,
        "tile_count": m,
        "tile_order": [list(tile) for tile in tiles],
        "tilings": nt,
        "classes": len(classes),
        "merged_nodes": len(node_codes),
        "blue_lines": blue,
        "black_lines": black,
        "local_components": local_components,
        "blueblack_components": line_components,
        "blueblack_root_reach": len(line_dist),
        "local_depth_equals_min_tiling_weight": sum(local_dist[u] == fibre_depth[u] for u in node_codes),
        "line_wl_rounds": line_rounds,
        "combined_wl_rounds": combined_rounds,
        "carrier_stats": carrier_stats,
        "tournament_analysis": {
            "pair_observable": "number of unordered node pairs separated by the carrier partition",
            "retention": retention,
            "economy": economy,
            "edge_flips_between_gauges": flips,
        },
        "nodes": node_records,
        "class_code_by_mask": canon_arr,
        "node_rank_by_mask": [node_rank[int(x)] for x in merged_arr],
        "global_index_by_mask": global_index,
        "fibre_index_by_mask": fibre_index,
    }


def render(result: dict) -> str:
    lines = [
        "TOURNAMENT-TILING METAGRAPH ADDRESS: EXACT n=7 CONTINUATION",
        "exact canonical colour refinement plus exhaustive within-cell permutations",
        "",
        f"tilings={result['tilings']} classes={result['classes']} merged_nodes={result['merged_nodes']}",
        f"complement lines: blue={result['blue_lines']} black={result['black_lines']}",
        f"blue/black components={result['blueblack_components']}; root reaches "
        f"{result['blueblack_root_reach']}/{result['merged_nodes']}",
        f"local components={result['local_components']}; quotient depth equals minimum tiling weight at "
        f"{result['local_depth_equals_min_tiling_weight']}/{result['merged_nodes']} nodes",
        f"WL rounds: line={result['line_wl_rounds']} combined={result['combined_wl_rounds']}",
        "",
        "controlled-forgetting audit (cells / max collision / separated pairs):",
    ]
    for name, stats in result["carrier_stats"].items():
        lines.append(f"  {name:24s} {stats['cells']:3d} / {stats['max_cell']:3d} / {stats['separated_pairs']:5d}")
    for gauge in ("retention", "economy"):
        data = result["tournament_analysis"][gauge]
        lines.append(
            f"Tournament Analysis {gauge}: scores={data['score_hist']} c3={data['directed_3cycles']} "
            f"SCC={data['scc_sizes']} Hpaths={data['hamiltonian_paths']}"
        )
    lines.append(f"gauge edge flips={result['tournament_analysis']['edge_flips_between_gauges']}")
    lines.append("")
    lines.append("first 20 objective node addresses:")
    lines.append("rank code(s) SC depth/path line-root line-WL combined-WL fibre")
    for node in result["nodes"][:20]:
        root_path = f"{node['blueblack_root_distance']}:{node['blueblack_root_word'] or '-'}"
        lines.append(
            f"{node['rank']:3d} {','.join(node['canonical_orbit_codes']):>18s} "
            f"{'Y' if node['self_converse'] else 'N'} {node['local_depth']}:{node['local_path_word'] or '-':<6s} "
            f"{root_path:<14s} {node['blueblack_wl_color']:3d} {node['combined_wl_color']:3d} "
            f"x{node['tiling_count']}"
        )
    lines.extend(
        [
            "",
            "ROUND TRIP: mask-indexed arrays give exact class/node/global/fibre indices; every node lists its inverse mask fibre.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--parent-json",
        type=Path,
        default=Path("05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json"),
    )
    parser.add_argument("--json", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.parent_json)
    text = render(result)
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")
    if args.json:
        # The mask-indexed arrays contain 32768 entries each.  Compact JSON
        # keeps the exact atlas reviewable as one generated payload instead of
        # creating roughly 170,000 mechanically wrapped lines.
        args.json.write_text(json.dumps(result, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
