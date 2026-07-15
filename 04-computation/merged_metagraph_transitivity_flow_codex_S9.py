#!/usr/bin/env python3
"""Exact transitivity flow on the tournament-tiling merged metagraph.

Nodes are placed on the objective cyclic-triangle spectrum

    transitive C3=0  --->  regular/near-regular C3=C3_max.

The script projects every explorer complement line to merged nodes, retains
its blue/grid-symmetric or black/grid-asymmetric colour, and orients it toward
larger C3.  It separates line instances (the actual tiling complement pairs)
from projected node-pair support, since their multiplicities carry different
information.

The theorem-facing identities are:

  4 * score_variance(T) = n(n^2-1)/3 - 8 C3(T),
  C3(complement(t))-C3(t) = d_0(t)-d_(n-1)(t)-1.

On a blue line, anti-diagonal grid symmetry forces
d_i+d_(n-1-i)=n-1, so the line flux is 2d_0-n.  Black lines carry the
additional left/right defect epsilon=d_0+d_(n-1)-(n-1).
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict, deque
from pathlib import Path

from tournament_tiling_metagraph_address_codex_S4 import (
    arc,
    c3_count,
    carrier_tournament,
    is_grid_symmetric,
    partition_stats,
    score_sequence,
    tile_schema,
    tiling_tournament,
)


SMALL_ATLAS = Path("05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json")
N7_ATLAS = Path("05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json")
CATEGORY_ORDER = {"pure_blue": 0, "mixed": 1, "pure_black": 2}


def c3_maximum(n: int) -> int:
    return n * (n * n - 1) // 24 if n % 2 else n * (n * n - 4) // 24


def balanced_score(n: int) -> tuple[int, ...]:
    if n % 2:
        return ((n - 1) // 2,) * n
    return (n // 2,) * (n // 2) + (n // 2 - 1,) * (n // 2)


def blue_step_prediction(n: int) -> Counter[int]:
    """Closed-form line counts by |Delta C3| on the blue locus."""
    m = math.comb(n - 1, 2)
    fixed_tiles = (n - 1) // 2
    reflection_orbits = (m + fixed_tiles) // 2
    top_variables = n - 2
    free_factor = 1 << (reflection_orbits - top_variables)
    result: Counter[int] = Counter()
    for wins in range((top_variables // 2) + 1):
        opposite = top_variables - wins
        step = abs(2 * wins + 2 - n)
        count = math.comb(top_variables, wins) * free_factor
        if wins == opposite:
            count //= 2
        result[step] += count
    return result


def all_line_step_prediction(n: int) -> Counter[int]:
    """Closed-form complement-line counts by |Delta C3| before colour split."""
    m = math.comb(n - 1, 2)
    middle = n - 3
    free_factor = 1 << (m - (2 * n - 5))
    signed: Counter[int] = Counter()
    for shared in (0, 1):
        for left_wins in range(middle + 1):
            for right_wins in range(middle + 1):
                step = 2 * shared - 1 + left_wins - right_wins
                signed[step] += (
                    math.comb(middle, left_wins)
                    * math.comb(middle, right_wins)
                    * free_factor
                )
    result = Counter({step: signed[step] for step in signed if step > 0})
    result[0] = signed[0] // 2
    return +result


def black_defect_prediction(n: int) -> Counter[int]:
    """Closed-form black-line counts by |epsilon|, epsilon=d0+dlast-(n-1)."""
    m = math.comb(n - 1, 2)
    middle = n - 3
    free_factor = 1 << (m - 2 * middle)
    blue_lines = sum(blue_step_prediction(n).values())
    result: Counter[int] = Counter()
    for defect in range(middle + 1):
        count = math.comb(2 * middle, middle + defect) * free_factor
        if defect == 0:
            count //= 2
            count -= blue_lines
        result[defect] = count
    return +result


def audit_closed_forms(max_n: int = 30) -> dict:
    """Cheap algebraic checks beyond the exhaustively enumerated sizes."""
    failures = 0
    for n in range(3, max_n + 1):
        m = math.comb(n - 1, 2)
        fixed_tiles = (n - 1) // 2
        failures += (m + fixed_tiles) % 2 != 0
        reflection_orbits = (m + fixed_tiles) // 2
        blue = blue_step_prediction(n)
        all_lines = all_line_step_prediction(n)
        black = all_lines - blue
        defects = black_defect_prediction(n)
        failures += sum(all_lines.values()) != 1 << (m - 1)
        failures += sum(blue.values()) != 1 << (reflection_orbits - 1)
        failures += sum(black.values()) != (1 << (m - 1)) - (1 << (reflection_orbits - 1))
        failures += sum(defects.values()) != sum(black.values())
        failures += max(blue) != n - 2
        failures += any(step % 2 != n % 2 for step in blue)
    assert failures == 0
    return {"n_min": 3, "n_max": max_n, "failures": failures}


def landau_slack(degrees: tuple[int, ...]) -> tuple[int, ...]:
    ascending = tuple(sorted(degrees))
    return tuple(sum(ascending[:k]) - k * (k - 1) // 2 for k in range(1, len(degrees)))


def pearson(xs: list[int], ys: list[int]) -> float:
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    numerator = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denominator = math.sqrt(
        sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)
    )
    return numerator / denominator if denominator else 0.0


def order_comparison(values: dict[int, int], depths: dict[int, int]) -> dict:
    nodes = sorted(values)
    counts = Counter()
    for i, u in enumerate(nodes):
        for v in nodes[i + 1 :]:
            dx = values[u] - values[v]
            dy = depths[u] - depths[v]
            if dx == 0 and dy == 0:
                counts["both_tie"] += 1
            elif dx == 0:
                counts["coordinate_tie"] += 1
            elif dy == 0:
                counts["triangle_tie"] += 1
            elif dx * dy > 0:
                counts["concordant"] += 1
            else:
                counts["discordant"] += 1
    counts["pairs"] = len(nodes) * (len(nodes) - 1) // 2
    return dict(counts)


def counter_json(counter: Counter) -> dict[str, int]:
    def key_text(key: object) -> str:
        if isinstance(key, tuple):
            return "|".join(map(str, key))
        return str(key)

    return {key_text(key): value for key, value in sorted(counter.items(), key=lambda kv: repr(kv[0]))}


def classify_node(masks: list[int], grid_symmetric: list[bool]) -> str:
    values = {grid_symmetric[mask] for mask in masks}
    if values == {True}:
        return "pure_blue"
    if values == {False}:
        return "pure_black"
    return "mixed"


def direction_label(source: str, target: str) -> str:
    return f"{source}->{target}"


def phase_direction(source: str, target: str) -> str:
    if CATEGORY_ORDER[source] < CATEGORY_ORDER[target]:
        return "categorical_outward"
    if CATEGORY_ORDER[source] > CATEGORY_ORDER[target]:
        return "categorical_inward"
    return "same_category"


def bfs_reach(
    root: int,
    adjacency: dict[int, set[tuple[int, str]]],
    blue_then_black: bool,
) -> tuple[set[int], dict[tuple[int, int], tuple[tuple[int, int], str] | None]]:
    start = (root, 0)
    seen = {start}
    parent: dict[tuple[int, int], tuple[tuple[int, int], str] | None] = {start: None}
    queue = deque([start])
    while queue:
        u, phase = queue.popleft()
        for v, color in sorted(adjacency[u], key=lambda item: (item[1], item[0])):
            if blue_then_black and phase == 1 and color == "blue":
                continue
            next_phase = int(phase == 1 or color == "black") if blue_then_black else 0
            state = (v, next_phase)
            if state not in seen:
                seen.add(state)
                parent[state] = ((u, phase), color)
                queue.append(state)
    return {u for u, _ in seen}, parent


def recover_path(
    target: int,
    parent: dict[tuple[int, int], tuple[tuple[int, int], str] | None],
) -> tuple[str, list[int]] | None:
    candidates = [state for state in parent if state[0] == target]
    if not candidates:
        return None
    state = min(candidates, key=lambda item: item[1])
    colors = []
    nodes = [state[0]]
    while parent[state] is not None:
        prior, color = parent[state]
        colors.append("B" if color == "blue" else "K")
        nodes.append(prior[0])
        state = prior
    return "".join(reversed(colors)), list(reversed(nodes))


def analyze(n: int, size: dict) -> dict:
    tiles, sigma = tile_schema(n)
    m = len(tiles)
    tiling_count = 1 << m
    full = tiling_count - 1
    node_records = {int(node["rank"]): node for node in size["nodes"]}
    if n <= 6:
        by_mask = {int(record["mask"]): record for record in size["tiling_map"]}
        node_by_mask = [int(by_mask[mask]["node_rank"]) for mask in range(tiling_count)]
    else:
        node_by_mask = list(map(int, size["node_rank_by_mask"]))

    grid_symmetric = [False] * tiling_count
    triangle_depth = [0] * tiling_count
    labelled_degrees: list[tuple[int, ...]] = [()] * tiling_count
    sorted_scores: list[tuple[int, ...]] = [()] * tiling_count
    energy4 = [0] * tiling_count
    energy4_transitive = n * (n * n - 1) // 3
    energy_failures = 0
    for mask in range(tiling_count):
        tournament = tiling_tournament(mask, n, tiles)
        degrees = tuple(sum(arc(tournament, n, i, j) for j in range(n)) for i in range(n))
        c3 = c3_count(tournament, n)
        value4 = sum((2 * degree - (n - 1)) ** 2 for degree in degrees)
        energy_failures += value4 != energy4_transitive - 8 * c3
        grid_symmetric[mask] = is_grid_symmetric(mask, sigma)
        triangle_depth[mask] = c3
        labelled_degrees[mask] = degrees
        sorted_scores[mask] = score_sequence(tournament, n)
        energy4[mask] = value4
    assert energy_failures == 0

    masks_by_node: dict[int, list[int]] = defaultdict(list)
    for mask, node in enumerate(node_by_mask):
        masks_by_node[node].append(mask)
    node_category = {
        node: classify_node(masks, grid_symmetric) for node, masks in masks_by_node.items()
    }
    node_depth = {node: triangle_depth[masks[0]] for node, masks in masks_by_node.items()}
    # Converse merging can join two different score sequences.  C3 and the
    # quadratic energy survive because d -> n-1-d, but the full score datum is
    # an unordered orbit of one or two sequences.
    node_score = {
        node: tuple(sorted({sorted_scores[mask] for mask in masks}))
        for node, masks in masks_by_node.items()
    }
    node_energy4 = {node: energy4[masks[0]] for node, masks in masks_by_node.items()}
    node_slack = {
        node: tuple(
            sorted({landau_slack(labelled_degrees[mask]) for mask in masks})
        )
        for node, masks in masks_by_node.items()
    }
    assert all(
        len({triangle_depth[mask] for mask in masks}) == 1
        for masks in masks_by_node.values()
    )

    observed_max = max(node_depth.values())
    expected_max = c3_maximum(n)
    target_score = balanced_score(n)
    assert observed_max == expected_max
    assert all(
        node_score[node] == (target_score,)
        for node in node_depth
        if node_depth[node] == observed_max
    )
    assert all(
        node_depth[node] < observed_max
        for node in node_depth
        if node_score[node] != (target_score,)
    )
    root = next(node for node, depth in node_depth.items() if depth == 0)
    assert list(node_depth.values()).count(0) == 1

    # A line support key retains colour and the unordered projected node pair.
    line_support: Counter[tuple[str, int, int]] = Counter()
    line_category_counts: Counter[tuple] = Counter()
    line_step_histogram: Counter[tuple] = Counter()
    line_direction_counts: Counter[tuple] = Counter()
    line_phase_counts: Counter[tuple] = Counter()
    self_loop_counts: Counter[tuple] = Counter()
    black_defect_histogram: Counter[int] = Counter()
    black_boundary_source_defect: Counter[tuple] = Counter()
    black_depth_pair_counts: Counter[tuple[int, int]] = Counter()
    black_cut_current: Counter[int] = Counter()
    flux_failures = blue_symmetry_failures = blue_parity_failures = 0
    color_endpoint_failures = 0
    for mask in range(tiling_count):
        other = mask ^ full
        if mask > other:
            continue
        color = "blue" if grid_symmetric[mask] else "black"
        u, v = node_by_mask[mask], node_by_mask[other]
        a, b = min(u, v), max(u, v)
        line_support[(color, a, b)] += 1
        category_pair = tuple(sorted((node_category[u], node_category[v])))
        line_category_counts[(color,) + category_pair] += 1
        delta = triangle_depth[other] - triangle_depth[mask]
        degrees = labelled_degrees[mask]
        flux = degrees[0] - degrees[-1] - 1
        epsilon = degrees[0] + degrees[-1] - (n - 1)
        flux_failures += delta != flux
        line_step_histogram[(color, abs(delta))] += 1
        if color == "blue":
            blue_symmetry_failures += any(
                degrees[i] + degrees[n - 1 - i] != n - 1 for i in range(n)
            )
            blue_symmetry_failures += epsilon != 0
            blue_symmetry_failures += delta != 2 * degrees[0] - n
            blue_parity_failures += abs(delta) % 2 != n % 2
            color_endpoint_failures += "pure_black" in category_pair
        else:
            black_defect_histogram[abs(epsilon)] += 1
            color_endpoint_failures += "pure_blue" in category_pair
            low_depth, high_depth = sorted((triangle_depth[mask], triangle_depth[other]))
            if low_depth != high_depth:
                black_depth_pair_counts[(low_depth, high_depth)] += 1
                for cut in range(low_depth, high_depth):
                    black_cut_current[cut] += 1
        if u == v:
            self_loop_counts[(color, node_category[u])] += 1

        if triangle_depth[mask] < triangle_depth[other]:
            source_mask, target_mask = mask, other
        elif triangle_depth[other] < triangle_depth[mask]:
            source_mask, target_mask = other, mask
        else:
            line_direction_counts[(color, "tie") + category_pair] += 1
            line_phase_counts[(color, "triangle_tie")] += 1
            if color == "black" and set(category_pair) == {"mixed", "pure_black"}:
                black_boundary_source_defect[("tie", abs(epsilon))] += 1
            continue
        source = node_by_mask[source_mask]
        target = node_by_mask[target_mask]
        source_category = node_category[source]
        target_category = node_category[target]
        line_direction_counts[(color, direction_label(source_category, target_category))] += 1
        line_phase_counts[(color, phase_direction(source_category, target_category))] += 1
        if color == "black" and {source_category, target_category} == {"mixed", "pure_black"}:
            source_degrees = labelled_degrees[source_mask]
            source_epsilon = source_degrees[0] + source_degrees[-1] - (n - 1)
            black_boundary_source_defect[
                (direction_label(source_category, target_category), source_epsilon)
            ] += 1

    assert flux_failures == blue_symmetry_failures == blue_parity_failures == 0
    assert color_endpoint_failures == 0

    observed_blue_steps = Counter(
        {step: count for (color, step), count in line_step_histogram.items() if color == "blue"}
    )
    observed_black_steps = Counter(
        {step: count for (color, step), count in line_step_histogram.items() if color == "black"}
    )
    predicted_blue_steps = blue_step_prediction(n)
    predicted_all_steps = all_line_step_prediction(n)
    predicted_black_steps = predicted_all_steps - predicted_blue_steps
    predicted_black_defects = black_defect_prediction(n)
    distribution_failures = sum(
        (
            observed_blue_steps != predicted_blue_steps,
            observed_black_steps != predicted_black_steps,
            black_defect_histogram != predicted_black_defects,
        )
    )
    assert distribution_failures == 0

    black_boundary_absolute_defect: Counter[tuple[str, int]] = Counter()
    for (direction, defect), count in black_boundary_source_defect.items():
        black_boundary_absolute_defect[(direction, abs(defect))] += count
    boundary_sign_symmetry_failures = 0
    for direction in {key[0] for key in black_boundary_source_defect if key[0] != "tie"}:
        maximum = max(abs(key[1]) for key in black_boundary_source_defect if key[0] == direction)
        for defect in range(1, maximum + 1):
            boundary_sign_symmetry_failures += (
                black_boundary_source_defect[(direction, defect)]
                != black_boundary_source_defect[(direction, -defect)]
            )
    assert boundary_sign_symmetry_failures == 0

    # Projected support has one edge regardless of tiling-line multiplicity.
    support_category_counts: Counter[tuple] = Counter()
    support_direction_counts: Counter[tuple] = Counter()
    support_phase_counts: Counter[tuple] = Counter()
    directed_adjacency: dict[int, set[tuple[int, str]]] = defaultdict(set)
    for color, u, v in line_support:
        pair = tuple(sorted((node_category[u], node_category[v])))
        support_category_counts[(color,) + pair] += 1
        if node_depth[u] < node_depth[v]:
            source, target = u, v
        elif node_depth[v] < node_depth[u]:
            source, target = v, u
        else:
            support_direction_counts[(color, "tie") + pair] += 1
            support_phase_counts[(color, "triangle_tie")] += 1
            directed_adjacency[u].add((v, color))
            directed_adjacency[v].add((u, color))
            continue
        source_category = node_category[source]
        target_category = node_category[target]
        support_direction_counts[(color, direction_label(source_category, target_category))] += 1
        support_phase_counts[(color, phase_direction(source_category, target_category))] += 1
        directed_adjacency[source].add((target, color))

    ordinary_reach, _ = bfs_reach(root, directed_adjacency, blue_then_black=False)
    phased_reach, phased_parent = bfs_reach(root, directed_adjacency, blue_then_black=True)
    balanced_nodes = sorted(node for node, depth in node_depth.items() if depth == observed_max)
    balanced_paths = {}
    for node in balanced_nodes:
        path = recover_path(node, phased_parent)
        assert path is not None
        balanced_paths[node_records[node]["id"]] = {
            "word": path[0],
            "node_ranks": path[1],
        }
    assert set(balanced_nodes) <= ordinary_reach
    assert set(balanced_nodes) <= phased_reach
    assert ordinary_reach == phased_reach

    # Incidence and parity laws at the category interface.
    node_incidence: dict[int, Counter] = {node: Counter() for node in masks_by_node}
    for (color, u, v), multiplicity in line_support.items():
        if u == v:
            node_incidence[u][f"{color}_self"] += multiplicity
        else:
            node_incidence[u][f"{color}_cross"] += multiplicity
            node_incidence[v][f"{color}_cross"] += multiplicity
    incidence_failures = 0
    for node, masks in masks_by_node.items():
        inc = node_incidence[node]
        incidence_failures += len(masks) != (
            inc["blue_cross"]
            + inc["black_cross"]
            + 2 * inc["blue_self"]
            + 2 * inc["black_self"]
        )
        if node_category[node] == "pure_blue":
            incidence_failures += bool(inc["black_cross"] or inc["black_self"])
        if node_category[node] == "pure_black":
            incidence_failures += bool(inc["blue_cross"] or inc["blue_self"])
    assert incidence_failures == 0

    category_by_depth: dict[int, Counter[str]] = defaultdict(Counter)
    for node, depth in node_depth.items():
        category_by_depth[depth][node_category[node]] += 1

    local_depth = {node: int(node_records[node]["local_depth"]) for node in masks_by_node}
    line_distance = {
        node: int(node_records[node]["blueblack_root_distance"]) for node in masks_by_node
    }
    flow_address = {
        node: (
            node_depth[node],
            CATEGORY_ORDER[node_category[node]],
            line_distance[node],
            node_records[node]["blueblack_root_word"],
            node_slack[node],
            node_score[node],
            len(masks_by_node[node]),
            int(node_records[node]["blueblack_wl_color"]),
            node,
        )
        for node in masks_by_node
    }
    flow_order = sorted(
        masks_by_node,
        key=flow_address.__getitem__,
    )
    flow_rank = {node: rank for rank, node in enumerate(flow_order)}

    carriers = {
        "pure_mixed_black_category": node_category,
        "cyclic_triangle_depth": node_depth,
        "score_sequence": node_score,
        "landau_slack_vector": node_slack,
        "local_flip_depth": local_depth,
        "blueblack_root_distance": line_distance,
        "flow_packet": {
            node: (node_depth[node], node_score[node], node_category[node]) for node in masks_by_node
        },
        "exact_HYP6825_address": {node: node for node in masks_by_node},
    }
    carrier_stats = {name: partition_stats(values) for name, values in carriers.items()}
    retention = carrier_tournament(carrier_stats, "retention")
    economy = carrier_tournament(carrier_stats, "economy")
    edge_flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carriers))
        for j in range(i + 1, len(carriers))
    )

    nodes_json = []
    for node in flow_order:
        record = node_records[node]
        inc = node_incidence[node]
        nodes_json.append(
            {
                "flow_rank": flow_rank[node],
                "flow_address": {
                    "cyclic_triangle_depth": node_depth[node],
                    "phase_rank": CATEGORY_ORDER[node_category[node]],
                    "blueblack_root_distance": line_distance[node],
                    "blueblack_root_word": record["blueblack_root_word"],
                    "landau_slack_orbit": [list(slack) for slack in node_slack[node]],
                    "score_sequence_orbit": [list(score) for score in node_score[node]],
                    "tiling_fibre_size": len(masks_by_node[node]),
                    "blueblack_wl_color": int(record["blueblack_wl_color"]),
                    "HYP6825_rank": node,
                },
                "node_id": record["id"],
                "HYP6825_rank": node,
                "category": node_category[node],
                "cyclic_triangles": node_depth[node],
                "cyclic_triangle_maximum": expected_max,
                "score_sequence_orbit": [list(score) for score in node_score[node]],
                "score_energy_times_4": node_energy4[node],
                "landau_slack_orbit": [list(slack) for slack in node_slack[node]],
                "local_flip_depth": local_depth[node],
                "blueblack_root_distance": line_distance[node],
                "blueblack_root_word": record["blueblack_root_word"],
                "tiling_count": len(masks_by_node[node]),
                "line_incidence": dict(inc),
            }
        )

    comparison_local = order_comparison(local_depth, node_depth)
    comparison_line = order_comparison(line_distance, node_depth)
    comparison_local["pearson"] = round(
        pearson([local_depth[node] for node in masks_by_node], [node_depth[node] for node in masks_by_node]),
        9,
    )
    comparison_line["pearson"] = round(
        pearson([line_distance[node] for node in masks_by_node], [node_depth[node] for node in masks_by_node]),
        9,
    )

    return {
        "n": n,
        "tile_count": m,
        "tilings": tiling_count,
        "classes": int(size["classes"]),
        "merged_nodes": len(masks_by_node),
        "root_node_id": node_records[root]["id"],
        "energy4_transitive": energy4_transitive,
        "cyclic_triangle_maximum": expected_max,
        "balanced_score_sequence": list(target_score),
        "balanced_nodes": [node_records[node]["id"] for node in balanced_nodes],
        "balanced_categories": counter_json(Counter(node_category[node] for node in balanced_nodes)),
        "balanced_monotone_blue_then_black_paths": balanced_paths,
        "node_category_counts": counter_json(Counter(node_category.values())),
        "node_categories_by_triangle_depth": {
            str(depth): counter_json(counts) for depth, counts in sorted(category_by_depth.items())
        },
        "line_instances": sum(line_support.values()),
        "line_support_edges": len(line_support),
        "line_instance_category_counts": counter_json(line_category_counts),
        "line_support_category_counts": counter_json(support_category_counts),
        "line_instance_direction_counts": counter_json(line_direction_counts),
        "line_support_direction_counts": counter_json(support_direction_counts),
        "line_instance_phase_counts": counter_json(line_phase_counts),
        "line_support_phase_counts": counter_json(support_phase_counts),
        "line_step_histogram": counter_json(line_step_histogram),
        "closed_form_blue_step_histogram": counter_json(predicted_blue_steps),
        "closed_form_black_step_histogram": counter_json(predicted_black_steps),
        "closed_form_all_line_step_histogram": counter_json(predicted_all_steps),
        "self_loop_counts": counter_json(self_loop_counts),
        "black_absolute_reflection_defect_histogram": counter_json(black_defect_histogram),
        "closed_form_black_absolute_reflection_defect_histogram": counter_json(
            predicted_black_defects
        ),
        "black_boundary_source_defect": counter_json(black_boundary_source_defect),
        "black_boundary_absolute_source_defect": counter_json(
            black_boundary_absolute_defect
        ),
        "black_boundary_sign_symmetry_failures": boundary_sign_symmetry_failures,
        "black_cyclic_triangle_level_pair_counts": counter_json(black_depth_pair_counts),
        "black_cyclic_triangle_cut_current": counter_json(black_cut_current),
        "strict_or_tie_monotone_reachable_nodes": len(ordinary_reach),
        "blue_then_black_monotone_reachable_nodes": len(phased_reach),
        "reachable_category_counts": counter_json(Counter(node_category[node] for node in ordinary_reach)),
        "energy_identity_failures": energy_failures,
        "complement_flux_failures": flux_failures,
        "blue_reflection_failures": blue_symmetry_failures,
        "blue_flux_parity_failures": blue_parity_failures,
        "line_color_endpoint_failures": color_endpoint_failures,
        "incidence_failures": incidence_failures,
        "closed_form_distribution_failures": distribution_failures,
        "local_depth_vs_triangle_depth": comparison_local,
        "line_distance_vs_triangle_depth": comparison_line,
        "carrier_stats": carrier_stats,
        "tournament_analysis": {
            "pair_observable": "number of merged-node pairs separated by the position carrier",
            "tie_hamiltonian_path": list(carriers),
            "retention": retention,
            "economy": economy,
            "edge_flips_between_gauges": edge_flips,
        },
        "nodes_in_flow_order": nodes_json,
    }


def boundary_count(size: dict, color: str, direction: str, support: bool = False) -> int:
    field = "line_support_direction_counts" if support else "line_instance_direction_counts"
    return int(size[field].get(f"{color}|{direction}", 0))


def render(result: dict) -> str:
    lines = [
        "THM-785 CYCLIC-TRIANGLE FLOW ON THE BLUE/BLACK MERGED METAGRAPH",
        "=" * 82,
        "score spectrum: sum_i(d_i-(n-1)/2)^2 = n(n^2-1)/12 - 2 C3",
        "line flux: C3(complement(t))-C3(t) = d0-dlast-1",
        "blue specialization: flux = 2 d0-n; black defect: epsilon=d0+dlast-(n-1)",
        f"closed-form aggregate audit: n={result['closed_form_general_audit']['n_min']}.."
        f"{result['closed_form_general_audit']['n_max']}, failures="
        f"{result['closed_form_general_audit']['failures']}",
        "",
    ]
    for size in result["sizes"]:
        n = size["n"]
        lines.extend(
            [
                f"n={n}: nodes={size['merged_nodes']} tilings={size['tilings']} "
                f"C3=0..{size['cyclic_triangle_maximum']} balanced={tuple(size['balanced_nodes'])}",
                f"  categories={size['node_category_counts']} balanced_categories={size['balanced_categories']}",
                f"  lines={size['line_instances']} support_edges={size['line_support_edges']} "
                f"self={size['self_loop_counts']}",
                f"  line steps={size['line_step_histogram']}",
                f"  blue phase={size['line_instance_phase_counts'].get('blue|categorical_outward',0)} outward, "
                f"{size['line_instance_phase_counts'].get('blue|categorical_inward',0)} inward, "
                f"{size['line_instance_phase_counts'].get('blue|same_category',0)} same, "
                f"{size['line_instance_phase_counts'].get('blue|triangle_tie',0)} tie",
                f"  black phase={size['line_instance_phase_counts'].get('black|categorical_outward',0)} outward, "
                f"{size['line_instance_phase_counts'].get('black|categorical_inward',0)} inward, "
                f"{size['line_instance_phase_counts'].get('black|same_category',0)} same, "
                f"{size['line_instance_phase_counts'].get('black|triangle_tie',0)} tie",
                f"  black reflection |epsilon|={size['black_absolute_reflection_defect_histogram']}",
                f"  monotone reach={size['strict_or_tie_monotone_reachable_nodes']}/{size['merged_nodes']}; "
                f"blue-then-black reach={size['blue_then_black_monotone_reachable_nodes']}",
                f"  balanced monotone paths="
                + ", ".join(
                    f"{node}:{path['word']}" for node, path in size["balanced_monotone_blue_then_black_paths"].items()
                ),
                f"  local-depth/C3 Pearson={size['local_depth_vs_triangle_depth']['pearson']}; "
                f"line-distance/C3 Pearson={size['line_distance_vs_triangle_depth']['pearson']}",
                f"  identity failures energy/flux/blue/parity/endpoints/incidence="
                f"{size['energy_identity_failures']}/{size['complement_flux_failures']}/"
                f"{size['blue_reflection_failures']}/{size['blue_flux_parity_failures']}/"
                f"{size['line_color_endpoint_failures']}/{size['incidence_failures']}",
                f"  closed-form line-distribution failures="
                f"{size['closed_form_distribution_failures']}; "
                f"black boundary sign-symmetry failures="
                f"{size['black_boundary_sign_symmetry_failures']}",
                "",
            ]
        )

    n7 = result["sizes"][-1]
    lines.extend(
        [
            "N=7 CATEGORY-BOUNDARY FLOW (line multiplicity / projected support)",
            f"  blue pure_blue->mixed: "
            f"{boundary_count(n7,'blue','pure_blue->mixed')}/"
            f"{boundary_count(n7,'blue','pure_blue->mixed',True)}",
            f"  blue mixed->pure_blue: "
            f"{boundary_count(n7,'blue','mixed->pure_blue')}/"
            f"{boundary_count(n7,'blue','mixed->pure_blue',True)}",
            f"  black mixed->pure_black: "
            f"{boundary_count(n7,'black','mixed->pure_black')}/"
            f"{boundary_count(n7,'black','mixed->pure_black',True)}",
            f"  black pure_black->mixed: "
            f"{boundary_count(n7,'black','pure_black->mixed')}/"
            f"{boundary_count(n7,'black','pure_black->mixed',True)}",
            f"  black mixed--pure_black ties: "
            f"{n7['line_instance_direction_counts'].get('black|tie|mixed|pure_black',0)}/"
            f"{n7['line_support_direction_counts'].get('black|tie|mixed|pure_black',0)}",
            f"  black boundary direction by |epsilon|: "
            f"{n7['black_boundary_absolute_source_defect']}",
            f"  black C3 cut current (k means k|k+1): "
            f"{n7['black_cyclic_triangle_cut_current']}",
            "",
            "INTERPRETATION",
            "  topology: pure-blue --blue-- mixed --black-- pure-black (no shortcut)",
            "  C3 orientation: the blue boundary is outward at n=7, but the black",
            "  boundary has a net reverse drift from pure-black back into mixed",
            "  blue left/right symmetry is exact; black |epsilon| is the transverse defect coordinate",
            "",
            "TOURNAMENT ANALYSIS",
            "  vertices: node-position carriers (category, C3, score, Landau slack, graph depths)",
            "  observable: number of unordered node pairs separated by each carrier",
            f"  n=7 carrier edge flips retention/economy={n7['tournament_analysis']['edge_flips_between_gauges']}",
            f"  tie Hamiltonian path={tuple(n7['tournament_analysis']['tie_hamiltonian_path'])}",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=SMALL_ATLAS)
    parser.add_argument("--n7-atlas", type=Path, default=N7_ATLAS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    small = json.loads(args.small_atlas.read_text())
    sizes = {int(size["n"]): size for size in small["sizes"]}
    sizes[7] = json.loads(args.n7_atlas.read_text())
    result = {
        "schema_version": 1,
        "theorem": "THM-785",
        "closed_form_general_audit": audit_closed_forms(),
        "sizes": [analyze(n, sizes[n]) for n in range(3, 8)],
    }
    output = render(result)
    print(output, end="")
    if args.output:
        args.output.write_text(output)
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
