#!/usr/bin/env python3
"""Exact addresses and fibres for the tournament-tiling merged metagraph.

This is deliberately an audit of two different graphs which older notes often
draw with overlapping colour names.

* LOCAL FLIP GRAPH: one staircase tile is flipped.  Its quotient supplies the
  radial coordinate from the transitive class.  We call its edge types
  S/R/G = SC--SC / SC--NS / NS--NS.
* COMPLEMENT-LINE GRAPH: all staircase tiles are flipped.  Each antipodal line
  is B (blue) when its tiling is fixed by the anti-diagonal grid reflection and
  K (black) otherwise.  This is the blue/black structure drawn by
  tournament-tiling-explorer.html.

The two graphs are not interchangeable.  We combine them only after retaining
their edge-channel names.  The resulting address is objective relative to the
explorer's declared fixed Hamiltonian path and tile order:

  (local depth, local path word, rooted blue/black path, blue/black WL colour,
   combined WL colour, deletion-parent profile, canonical converse-orbit code)

The last coordinate is an explicit canonical tie-breaker.  Before it is used,
the program reports all collisions, so no graph profile is silently promoted
to a complete invariant.

The JSON map is two-sorted.  A tiling maps to one canonical tournament class
and one converse-merged node; a node maps back to its complete fibre of
tilings.  It never chooses a fictitious unique tiling for a class.

Tournament Analysis treats *candidate information carriers*, not runners, as
vertices.  The pair observable is the number of node pairs separated by each
carrier.  The switch is either retention or retention-per-description-bit;
the listed carrier order is the tie Hamiltonian path.

Usage:
  python3 tournament_tiling_metagraph_address_codex_S4.py
  python3 tournament_tiling_metagraph_address_codex_S4.py \
      --json path.json --output path.out
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path


EXPECTED_CLASSES = {3: 2, 4: 4, 5: 12, 6: 56}
EXPECTED_MERGED = {3: 2, 4: 3, 5: 10, 6: 34}


@lru_cache(maxsize=None)
def pairs(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


@lru_cache(maxsize=None)
def perms(n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(permutations(range(n)))


@lru_cache(maxsize=None)
def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {pair: k for k, pair in enumerate(pairs(n))}


@lru_cache(maxsize=None)
def permutation_arc_maps(n: int) -> tuple[tuple[tuple[int, bool], ...], ...]:
    """For each permutation/new pair: old pair-bit index and reversal flag."""
    index = pair_index(n)
    maps = []
    for p in perms(n):
        one = []
        for i, j in pairs(n):
            a, b = p[i], p[j]
            one.append((index[(min(a, b), max(a, b))], a > b))
        maps.append(tuple(one))
    return tuple(maps)


def arc(mask: int, n: int, i: int, j: int) -> int:
    """Return 1 iff i -> j in the labelled tournament mask."""
    if i == j:
        return 0
    if i < j:
        k = pair_index(n)[(i, j)]
        return (mask >> k) & 1
    k = pair_index(n)[(j, i)]
    return 1 - ((mask >> k) & 1)


@lru_cache(maxsize=None)
def canonical(mask: int, n: int) -> int:
    """Lexicographically least upper-triangle adjacency mask under S_n."""
    best = None
    for pmap in permutation_arc_maps(n):
        value = 0
        for k, (old_k, reversed_arc) in enumerate(pmap):
            bit = (mask >> old_k) & 1
            if bit ^ reversed_arc:
                value |= 1 << k
        if best is None or value < best:
            best = value
    assert best is not None
    return best


def converse(mask: int, n: int) -> int:
    return mask ^ ((1 << len(pairs(n))) - 1)


def score_sequence(mask: int, n: int) -> tuple[int, ...]:
    return tuple(sorted((sum(arc(mask, n, i, j) for j in range(n)) for i in range(n)), reverse=True))


def c3_count(mask: int, n: int) -> int:
    ans = 0
    for i, j, k in combinations(range(n), 3):
        ans += int(
            (arc(mask, n, i, j) and arc(mask, n, j, k) and arc(mask, n, k, i))
            or (arc(mask, n, i, k) and arc(mask, n, k, j) and arc(mask, n, j, i))
        )
    return ans


def hamiltonian_paths(mask: int, n: int) -> int:
    dp: dict[tuple[int, int], int] = {(1 << v, v): 1 for v in range(n)}
    for subset in range(1, 1 << n):
        for v in range(n):
            count = dp.get((subset, v), 0)
            if not count:
                continue
            for w in range(n):
                if not (subset >> w) & 1 and arc(mask, n, v, w):
                    key = (subset | (1 << w), w)
                    dp[key] = dp.get(key, 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def automorphisms(mask: int, n: int) -> int:
    ans = 0
    for pmap in permutation_arc_maps(n):
        image = 0
        for k, (old_k, reversed_arc) in enumerate(pmap):
            if ((mask >> old_k) & 1) ^ reversed_arc:
                image |= 1 << k
        ans += image == mask
    return ans


def delete_vertex(mask: int, n: int, dead: int) -> int:
    live = [v for v in range(n) if v != dead]
    value = 0
    for k, (i, j) in enumerate(pairs(n - 1)):
        if arc(mask, n, live[i], live[j]):
            value |= 1 << k
    return value


def deletion_profile(mask: int, n: int) -> tuple[tuple[int, int], ...]:
    counts = Counter(canonical(delete_vertex(mask, n, v), n - 1) for v in range(n))
    return tuple(sorted(counts.items()))


def tile_schema(n: int) -> tuple[tuple[tuple[int, int], ...], tuple[int, ...]]:
    """Return explorer tile order and anti-diagonal reflection permutation."""
    tiles = tuple((x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1))
    index = {tile: k for k, tile in enumerate(tiles)}
    sigma = tuple(index[(n - y + 1, n - x + 1)] for x, y in tiles)
    return tiles, sigma


def tiling_tournament(tile_mask: int, n: int, tiles: tuple[tuple[int, int], ...]) -> int:
    """Explorer convention: vertices n,...,1; all tiles off is transitive."""
    vertex_index = {n - i: i for i in range(n)}
    value = 0
    pindex = pair_index(n)
    for i in range(n - 1):
        value |= 1 << pindex[(i, i + 1)]
    for k, (x, y) in enumerate(tiles):
        i, j = vertex_index[x], vertex_index[y]
        assert i < j
        if not ((tile_mask >> k) & 1):
            value |= 1 << pindex[(i, j)]
    return value


def is_grid_symmetric(tile_mask: int, sigma: tuple[int, ...]) -> bool:
    return all(((tile_mask >> i) & 1) == ((tile_mask >> sigma[i]) & 1) for i in range(len(sigma)))


def incidence_profiles(
    nodes: list[int],
    edges: dict[tuple[int, int], Counter[str]],
    channels: set[str],
) -> dict[int, tuple[tuple[str, int, int], ...]]:
    out: dict[int, Counter[tuple[str, int]]] = {u: Counter() for u in nodes}
    for (u, v), counts in edges.items():
        for channel, count in counts.items():
            if channel not in channels:
                continue
            if u == v:
                out[u][(channel, 1)] += 2 * count
            else:
                out[u][(channel, 0)] += count
                out[v][(channel, 0)] += count
    return {u: tuple(sorted((ch, loop, count) for (ch, loop), count in out[u].items())) for u in nodes}


def stable_refinement(
    nodes: list[int],
    initial: dict[int, tuple],
    edges: dict[tuple[int, int], Counter[str]],
    channels: set[str],
) -> tuple[dict[int, int], int]:
    """Canonical 1-WL on an edge-coloured weighted multigraph."""
    initial_values = sorted(set(initial.values()), key=repr)
    rank = {value: i for i, value in enumerate(initial_values)}
    colors = {u: rank[initial[u]] for u in nodes}
    rounds = 0
    while True:
        rounds += 1
        neighborhoods: dict[int, Counter[tuple[str, int]]] = {u: Counter() for u in nodes}
        for (u, v), counts in edges.items():
            for channel, count in counts.items():
                if channel not in channels:
                    continue
                if u == v:
                    neighborhoods[u][(channel, colors[u])] += 2 * count
                else:
                    neighborhoods[u][(channel, colors[v])] += count
                    neighborhoods[v][(channel, colors[u])] += count
        signatures = {
            u: (colors[u], tuple(sorted((ch, col, count) for (ch, col), count in neighborhoods[u].items())))
            for u in nodes
        }
        values = sorted(set(signatures.values()), key=repr)
        sig_rank = {value: i for i, value in enumerate(values)}
        new_colors = {u: sig_rank[signatures[u]] for u in nodes}
        if len(set(new_colors.values())) == len(set(colors.values())):
            return new_colors, rounds
        colors = new_colors


def edge_type_local(sc_u: bool, sc_v: bool) -> str:
    if sc_u and sc_v:
        return "F_S"  # local SC--SC (spine)
    if sc_u != sc_v:
        return "F_R"  # local SC--NS (rib)
    return "F_G"      # local NS--NS (sea)


def lex_shortest_paths(
    nodes: list[int],
    root: int,
    edges: dict[tuple[int, int], Counter[str]],
    allowed: set[str],
    letters: dict[str, str],
    letter_order: dict[str, int] | None = None,
) -> tuple[dict[int, int], dict[int, str], int]:
    adjacency: dict[int, list[tuple[int, str]]] = {u: [] for u in nodes}
    for (u, v), counts in edges.items():
        if u == v:
            continue
        for channel, count in counts.items():
            if channel in allowed and count:
                adjacency[u].append((v, letters[channel]))
                adjacency[v].append((u, letters[channel]))
    dist = {root: 0}
    queue = deque([root])
    while queue:
        u = queue.popleft()
        for v, _ in adjacency[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)
    words = {root: ""}
    letter_order = letter_order or {letter: i for i, letter in enumerate(sorted(set(letters.values())))}
    word_key = lambda word: tuple(letter_order[letter] for letter in word)
    for d in range(1, max(dist.values(), default=0) + 1):
        for v in sorted(u for u in nodes if dist.get(u) == d):
            candidates = [words[u] + letter for u, letter in adjacency[v] if dist.get(u) == d - 1]
            words[v] = min(candidates, key=word_key)
    components = 0
    seen: set[int] = set()
    for start in nodes:
        if start in seen:
            continue
        components += 1
        todo = [start]
        seen.add(start)
        while todo:
            u = todo.pop()
            for v, _ in adjacency[u]:
                if v not in seen:
                    seen.add(v)
                    todo.append(v)
    return dist, words, components


def partition_stats(values: dict[int, object]) -> dict[str, int]:
    cells = Counter(values.values())
    return {
        "cells": len(cells),
        "colliding_cells": sum(size > 1 for size in cells.values()),
        "max_cell": max(cells.values()),
        "separated_pairs": sum(a * b for i, a in enumerate(cells.values()) for b in list(cells.values())[i + 1 :]),
    }


def count_directed_triangles(adj: list[list[int]]) -> int:
    n = len(adj)
    return sum(adj[i][j] and adj[j][k] and adj[k][i] for i, j, k in combinations(range(n), 3)) + sum(
        adj[i][k] and adj[k][j] and adj[j][i] for i, j, k in combinations(range(n), 3)
    )


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)
    order: list[int] = []
    seen: set[int] = set()

    def dfs(u: int) -> None:
        seen.add(u)
        for v in range(n):
            if adj[u][v] and v not in seen:
                dfs(v)
        order.append(u)

    for u in range(n):
        if u not in seen:
            dfs(u)
    rev_seen: set[int] = set()
    sizes: list[int] = []

    def rdfs(u: int) -> int:
        rev_seen.add(u)
        return 1 + sum(rdfs(v) for v in range(n) if adj[v][u] and v not in rev_seen)

    for u in reversed(order):
        if u not in rev_seen:
            sizes.append(rdfs(u))
    return sorted(sizes, reverse=True)


def count_tournament_hamilton_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = {(1 << v, v): 1 for v in range(n)}
    for subset in range(1, 1 << n):
        for v in range(n):
            count = dp.get((subset, v), 0)
            if not count:
                continue
            for w in range(n):
                if not (subset >> w) & 1 and adj[v][w]:
                    dp[(subset | (1 << w), w)] = dp.get((subset | (1 << w), w), 0) + count
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))


def carrier_tournament(stats: dict[str, dict[str, int]], gauge: str) -> dict:
    names = list(stats)
    n = len(names)

    def value(name: str) -> float:
        separated = stats[name]["separated_pairs"]
        if gauge == "retention":
            return float(separated)
        bits = math.log2(max(2, stats[name]["cells"]))
        return separated / bits

    adjacency = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            vi, vj = value(names[i]), value(names[j])
            if vi > vj or (vi == vj and i < j):
                adjacency[i][j] = 1
            else:
                adjacency[j][i] = 1
    scores = [sum(row) for row in adjacency]
    return {
        "gauge": gauge,
        "vertices": names,
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3cycles": count_directed_triangles(adjacency),
        "scc_sizes": scc_sizes(adjacency),
        "hamiltonian_paths": count_tournament_hamilton_paths(adjacency),
        "adjacency": adjacency,
    }


def analyze(n: int, parent_rank_by_class: dict[int, int] | None = None) -> dict:
    tiles, sigma = tile_schema(n)
    m = len(tiles)
    full_tile = (1 << m) - 1
    tiling_pairmask = {t: tiling_tournament(t, n, tiles) for t in range(1 << m)}
    tiling_class = {t: canonical(mask, n) for t, mask in tiling_pairmask.items()}
    class_fibres: dict[int, list[int]] = defaultdict(list)
    for t, cls in tiling_class.items():
        class_fibres[cls].append(t)
    classes = sorted(class_fibres)

    class_info = {}
    for cls in classes:
        rep = tiling_pairmask[class_fibres[cls][0]]
        opp = canonical(converse(rep, n), n)
        class_info[cls] = {
            "converse": opp,
            "sc": cls == opp,
            "score": score_sequence(rep, n),
            "c3": c3_count(rep, n),
            "H": hamiltonian_paths(rep, n),
            "aut": automorphisms(rep, n),
            "deletion": deletion_profile(rep, n),
            "fibre": len(class_fibres[cls]),
        }

    merged_key_of_class = {cls: min(cls, class_info[cls]["converse"]) for cls in classes}
    tiling_node = {t: merged_key_of_class[cls] for t, cls in tiling_class.items()}
    node_classes: dict[int, list[int]] = defaultdict(list)
    for cls in classes:
        node_classes[merged_key_of_class[cls]].append(cls)
    nodes = sorted(node_classes)
    root = tiling_node[0]
    node_fibres: dict[int, list[int]] = defaultdict(list)
    for t, node in tiling_node.items():
        node_fibres[node].append(t)

    node_sc = {u: len(node_classes[u]) == 1 for u in nodes}
    node_parent = {
        u: tuple(sorted(tuple((code, multiplicity) for code, multiplicity in class_info[cls]["deletion"]) for cls in node_classes[u]))
        for u in nodes
    }
    node_parent_address = {
        u: tuple(
            sorted(
                tuple(
                    (
                        parent_rank_by_class.get(code, code) if parent_rank_by_class is not None else code,
                        multiplicity,
                    )
                    for code, multiplicity in class_info[cls]["deletion"]
                )
                for cls in node_classes[u]
            )
        )
        for u in nodes
    }
    node_base = {
        u: (
            int(u == root),
            int(node_sc[u]),
            len(node_fibres[u]),
            tuple(sorted(class_info[c]["fibre"] for c in node_classes[u])),
            tuple(sorted(class_info[c]["H"] for c in node_classes[u])),
            tuple(sorted(class_info[c]["score"] for c in node_classes[u])),
            tuple(sorted(class_info[c]["c3"] for c in node_classes[u])),
            tuple(sorted(class_info[c]["aut"] for c in node_classes[u])),
            node_parent_address[u],
        )
        for u in nodes
    }

    edges: dict[tuple[int, int], Counter[str]] = defaultdict(Counter)
    local_cube_edges = 0
    for t in range(1 << m):
        for bit in range(m):
            other = t ^ (1 << bit)
            if t > other:
                continue
            u, v = tiling_node[t], tiling_node[other]
            key = (min(u, v), max(u, v))
            edges[key][edge_type_local(node_sc[u], node_sc[v])] += 1
            local_cube_edges += 1

    blue_lines = black_lines = 0
    for t in range(1 << m):
        other = t ^ full_tile
        if t > other:
            continue
        u, v = tiling_node[t], tiling_node[other]
        key = (min(u, v), max(u, v))
        channel = "L_B" if is_grid_symmetric(t, sigma) else "L_K"
        edges[key][channel] += 1
        blue_lines += channel == "L_B"
        black_lines += channel == "L_K"

    line_channels = {"L_B", "L_K"}
    local_channels = {"F_S", "F_R", "F_G"}
    all_channels = line_channels | local_channels
    line_incidence = incidence_profiles(nodes, edges, line_channels)
    local_incidence = incidence_profiles(nodes, edges, local_channels)

    line_initial = {u: (int(u == root), int(node_sc[u]), line_incidence[u]) for u in nodes}
    line_wl, line_rounds = stable_refinement(nodes, line_initial, edges, line_channels)
    combined_initial = {u: (node_base[u], line_incidence[u], local_incidence[u]) for u in nodes}
    combined_wl, combined_rounds = stable_refinement(nodes, combined_initial, edges, all_channels)

    local_dist, local_word, local_components = lex_shortest_paths(
        nodes,
        root,
        edges,
        local_channels,
        {"F_S": "S", "F_R": "R", "F_G": "G"},
        {"S": 0, "R": 1, "G": 2},
    )
    line_dist, line_word, line_components = lex_shortest_paths(
        nodes, root, edges, line_channels, {"L_B": "B", "L_K": "K"}
    )
    fibre_depth = {u: min(t.bit_count() for t in node_fibres[u]) for u in nodes}

    def structure_address(u: int) -> tuple:
        # Line distance is a coordinate only where the blue/black graph reaches
        # the node.  The sentinel makes disconnection explicit, not invisible.
        line_coordinate = (0, line_dist[u], line_word[u]) if u in line_dist else (1, math.inf, "~")
        return (
            local_dist[u],
            tuple({"S": 0, "R": 1, "G": 2}[letter] for letter in local_word[u]),
            line_coordinate,
            line_wl[u],
            combined_wl[u],
            node_parent_address[u],
        )

    structural = {u: structure_address(u) for u in nodes}
    full_address = {u: structural[u] + (u,) for u in nodes}
    ordered_nodes = sorted(nodes, key=lambda u: full_address[u])
    node_rank = {u: rank for rank, u in enumerate(ordered_nodes)}

    global_tilings = sorted(range(1 << m), key=lambda t: (t.bit_count(), node_rank[tiling_node[t]], tiling_class[t], t))
    global_index = {t: i for i, t in enumerate(global_tilings)}
    fibre_index = {}
    for u in nodes:
        for i, t in enumerate(sorted(node_fibres[u], key=lambda t: (t.bit_count(), tiling_class[t], t))):
            fibre_index[t] = i

    carriers = {
        "raw_sc_type": {u: int(node_sc[u]) for u in nodes},
        "radial_depth": {u: local_dist[u] for u in nodes},
        "blueblack_incidence": {u: line_incidence[u] for u in nodes},
        "rooted_blueblack_wl": {u: line_wl[u] for u in nodes},
        "local_parent_atlas": {u: (local_dist[u], local_word[u], node_parent_address[u]) for u in nodes},
        "combined_colored_wl": {u: combined_wl[u] for u in nodes},
        "structural_address": structural,
        "exact_address": full_address,
    }
    carrier_stats = {name: partition_stats(values) for name, values in carriers.items()}
    retention_tournament = carrier_tournament(carrier_stats, "retention")
    economy_tournament = carrier_tournament(carrier_stats, "economy")
    flips = 0
    A, B = retention_tournament["adjacency"], economy_tournament["adjacency"]
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            flips += A[i][j] != B[i][j]

    # Exact round-trip and structural laws.
    assert len(classes) == EXPECTED_CLASSES[n]
    assert len(nodes) == EXPECTED_MERGED[n]
    assert sum(len(v) for v in class_fibres.values()) == 1 << m
    assert sum(len(v) for v in node_fibres.values()) == 1 << m
    assert local_cube_edges == m * (1 << (m - 1))
    assert blue_lines + black_lines == 1 << (m - 1)
    assert len(set(global_index.values())) == 1 << m
    assert global_tilings[0] == 0
    assert all(t in node_fibres[tiling_node[t]] for t in range(1 << m))
    assert all(node_sc[u] == (len(node_fibres[u]) % 2 == 1) for u in nodes)
    assert all(local_dist[u] <= fibre_depth[u] for u in nodes)

    node_json = []
    for u in ordered_nodes:
        classes_u = sorted(node_classes[u])
        masks = sorted(node_fibres[u], key=lambda t: fibre_index[t])
        node_json.append(
            {
                "rank": node_rank[u],
                "id": f"n{n}-a{node_rank[u]:02d}",
                "canonical_orbit_codes": [hex(c) for c in classes_u],
                "self_converse": node_sc[u],
                "local_depth": local_dist[u],
                "minimum_labelled_tiling_weight": fibre_depth[u],
                "local_path_word": local_word[u],
                "blueblack_root_distance": line_dist.get(u),
                "blueblack_root_word": line_word.get(u),
                "blueblack_wl_color": line_wl[u],
                "combined_wl_color": combined_wl[u],
                "line_incidence": [list(item) for item in line_incidence[u]],
                "deletion_parent_profile": repr(node_parent[u]),
                "recursive_parent_address": repr(node_parent_address[u]),
                "tiling_count": len(masks),
                "tiling_masks": masks,
                "global_tiling_indices": [global_index[t] for t in masks],
            }
        )

    tiling_json = []
    for t in global_tilings:
        u = tiling_node[t]
        tiling_json.append(
            {
                "global_index": global_index[t],
                "mask": t,
                "bit_word": "".join(str((t >> k) & 1) for k in range(m)),
                "popcount": t.bit_count(),
                "class_code": hex(tiling_class[t]),
                "node_rank": node_rank[u],
                "node_id": f"n{n}-a{node_rank[u]:02d}",
                "fibre_index": fibre_index[t],
                "grid_symmetric": is_grid_symmetric(t, sigma),
                "complement_mask": t ^ full_tile,
                "line_color": "blue" if is_grid_symmetric(t, sigma) else "black",
            }
        )

    return {
        "n": n,
        "tile_count": m,
        "tile_order": [list(tile) for tile in tiles],
        "tilings": 1 << m,
        "classes": len(classes),
        "merged_nodes": len(nodes),
        "root_node_rank": node_rank[root],
        "local_cube_edges": local_cube_edges,
        "blue_lines": blue_lines,
        "black_lines": black_lines,
        "local_components": local_components,
        "blueblack_components": line_components,
        "blueblack_root_reach": len(line_dist),
        "local_depth_equals_min_tiling_weight": sum(local_dist[u] == fibre_depth[u] for u in nodes),
        "line_wl_rounds": line_rounds,
        "combined_wl_rounds": combined_rounds,
        "carrier_stats": carrier_stats,
        "tournament_analysis": {
            "pair_observable": "number of unordered node pairs separated by the carrier partition",
            "switches": ["retention", "retention per log2(number of cells)"],
            "tie_hamiltonian_path": list(carriers),
            "retention": retention_tournament,
            "economy": economy_tournament,
            "edge_flips_between_gauges": flips,
        },
        "nodes": node_json,
        "tiling_map": tiling_json,
    }


def render(results: list[dict]) -> str:
    lines = []
    lines.append("TOURNAMENT-TILING METAGRAPH ADDRESS AND FIBRE AUDIT")
    lines.append("exact enumeration; explorer convention; n=3..6")
    lines.append("")
    lines.append("CHANNELS: local F_S/F_R/F_G are one-tile spine/rib/sea moves;")
    lines.append("          line L_B/L_K are all-tile complement blue/black lines.")
    lines.append("These channels are kept distinct throughout.")
    lines.append("")
    for r in results:
        lines.append(
            f"n={r['n']} m={r['tile_count']}: tilings={r['tilings']} classes={r['classes']} "
            f"merged={r['merged_nodes']} local_edges={r['local_cube_edges']}"
        )
        lines.append(
            f"  complement lines: blue={r['blue_lines']} black={r['black_lines']}; "
            f"blue/black components={r['blueblack_components']}, root reaches "
            f"{r['blueblack_root_reach']}/{r['merged_nodes']} nodes"
        )
        lines.append(
            f"  local graph: components={r['local_components']}; quotient depth equals minimum labelled "
            f"tiling weight at {r['local_depth_equals_min_tiling_weight']}/{r['merged_nodes']} nodes"
        )
        lines.append("  controlled-forgetting audit (cells / max collision / separated pairs):")
        for name, s in r["carrier_stats"].items():
            lines.append(f"    {name:24s} {s['cells']:3d} / {s['max_cell']:2d} / {s['separated_pairs']:4d}")
        ta = r["tournament_analysis"]
        for gauge in ("retention", "economy"):
            g = ta[gauge]
            lines.append(
                f"  Tournament Analysis {gauge}: scores={g['score_hist']} c3={g['directed_3cycles']} "
                f"SCC={g['scc_sizes']} Hpaths={g['hamiltonian_paths']}"
            )
        lines.append(f"  gauge edge flips={ta['edge_flips_between_gauges']}")
        lines.append("")
    r = results[-1]
    lines.append("n=6 OBJECTIVE NODE ORDER")
    lines.append("rank code(s) SC depth/path line-root line-WL combined-WL fibre")
    for node in r["nodes"]:
        root_path = (
            f"{node['blueblack_root_distance']}:{node['blueblack_root_word']}"
            if node["blueblack_root_distance"] is not None
            else "disconnected"
        )
        lines.append(
            f"{node['rank']:2d} {','.join(node['canonical_orbit_codes']):>13s} "
            f"{'Y' if node['self_converse'] else 'N'} {node['local_depth']}:{node['local_path_word'] or '-':<5s} "
            f"{root_path:<13s} {node['blueblack_wl_color']:2d} {node['combined_wl_color']:2d} "
            f"x{node['tiling_count']}"
        )
    lines.append("")
    lines.append("ROUND TRIP: every tiling mask occurs once in tiling_map; every node lists exactly its inverse fibre.")
    lines.append("ORDER: depth and coloured structure precede the canonical converse-orbit code; the code only breaks residual ties.")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--sizes", default="3,4,5,6")
    args = parser.parse_args()
    sizes = [int(x) for x in args.sizes.split(",")]
    results = []
    prior_by_size: dict[int, dict[int, int]] = {}
    for n in sizes:
        result = analyze(n, prior_by_size.get(n - 1))
        results.append(result)
        prior_by_size[n] = {
            int(code, 16): node["rank"]
            for node in result["nodes"]
            for code in node["canonical_orbit_codes"]
        }
    text = render(results)
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")
    if args.json:
        args.json.write_text(json.dumps({"schema_version": 1, "sizes": results}, indent=2) + "\n")


if __name__ == "__main__":
    main()
