#!/usr/bin/env python3
"""Exact audit for THM-848: the Walsh needle stalk of Hamilton paths.

The full tournament arc-flip walk is used, even though representatives come
from fixed-Hamiltonian-path tiling atlases.  Every converse-merged class has a
representative in those atlases, so one representative per node is enough.

The audit independently computes the homogeneous Walsh layers from their
even-path-forest support, the actual one-arc drift, and the complete multiset
of directional H derivatives.  It checks all merged nodes at n=4,5,6,7 and
optionally checks the gradient classifier at n=8 against a hash-guarded raw
rank atlas.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from array import array
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import ceil, comb, factorial, log2
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SMALL_ATLAS = ROOT / "05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json"
DEFAULT_N7_ATLAS = ROOT / "05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json"
DEFAULT_FLOW = ROOT / "05-knowledge/results/merged_metagraph_transitivity_flow_codex_S9.json"
N8_ATLAS_SHA256 = "30debad3387a4ea0ef51108ea132115efda2ac2fcdfcc2c5c1d4d23155095835"
N8_NODE_COUNT = 3528


def tile_schema(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (x, y)
        for y in range(1, n - 1)
        for x in range(n, y + 1, -1)
        if x - y >= 2
    )


def adjacency_from_tiling(mask: int, n: int) -> list[list[int]]:
    adjacency = [[0] * n for _ in range(n)]
    for label in range(n, 1, -1):
        adjacency[label - 1][label - 2] = 1
    for bit, (x, y) in enumerate(tile_schema(n)):
        if (mask >> bit) & 1:
            adjacency[y - 1][x - 1] = 1
        else:
            adjacency[x - 1][y - 1] = 1
    return adjacency


def adjacency_from_legacy_tiling(mask: int, n: int) -> list[list[int]]:
    """Decode the certified n=8 atlas convention: a set bit means high wins."""
    adjacency = [[0] * n for _ in range(n)]
    for label in range(n, 1, -1):
        adjacency[label - 1][label - 2] = 1
    for bit, (high, low) in enumerate(tile_schema(n)):
        if (mask >> bit) & 1:
            adjacency[high - 1][low - 1] = 1
        else:
            adjacency[low - 1][high - 1] = 1
    return adjacency


def certified_n8_representatives(path: Path) -> tuple[list[int], str]:
    """Load one mask per merged rank from THM-828/843's hash-guarded atlas."""
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != N8_ATLAS_SHA256:
        raise RuntimeError(
            f"n=8 atlas SHA256 mismatch: expected {N8_ATLAS_SHA256}, got {digest}"
        )
    ranks = array("H")
    with path.open("rb") as atlas_file:
        ranks.fromfile(atlas_file, 1 << 21)
    if ranks.itemsize != 2 or len(ranks) != 1 << 21:
        raise RuntimeError("n=8 atlas must contain exactly 2^21 uint16 ranks")
    if sys.byteorder != "little":
        ranks.byteswap()
    representatives = [-1] * N8_NODE_COUNT
    for mask, rank in enumerate(ranks):
        if rank >= N8_NODE_COUNT:
            raise RuntimeError(f"invalid n=8 merged rank {rank}")
        if representatives[rank] < 0:
            representatives[rank] = mask
    if any(mask < 0 for mask in representatives):
        raise RuntimeError("n=8 atlas does not contain all 3,528 merged ranks")
    return representatives, digest


def hamiltonian_paths(adjacency: list[list[int]]) -> int:
    n = len(adjacency)
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for support in range(1 << n):
        for last, count in enumerate(dp[support]):
            if not count:
                continue
            for nxt in range(n):
                if not (support >> nxt) & 1 and adjacency[last][nxt]:
                    dp[support | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def cyclic_triangles(adjacency: list[list[int]]) -> int:
    n = len(adjacency)
    degrees = [sum(row) for row in adjacency]
    return comb(n, 3) - sum(comb(degree, 2) for degree in degrees)


def score_sequence(adjacency: list[list[int]]) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adjacency))


def flip_gradient(adjacency: list[list[int]]) -> tuple[int, tuple[int, ...]]:
    """Return H and sorted values H(T^e)-H(T) over all unordered arcs e."""
    n = len(adjacency)
    h_value = hamiltonian_paths(adjacency)
    derivatives = []
    for left in range(n):
        for right in range(left + 1, n):
            adjacency[left][right], adjacency[right][left] = (
                adjacency[right][left],
                adjacency[left][right],
            )
            derivatives.append(hamiltonian_paths(adjacency) - h_value)
            adjacency[left][right], adjacency[right][left] = (
                adjacency[right][left],
                adjacency[left][right],
            )
    return h_value, tuple(sorted(derivatives))


def even_path_forests(n: int, level: int) -> tuple[tuple[int, tuple[tuple[int, ...], ...]], ...]:
    """Supports for Walsh degree 2*level, with coefficient numerators.

    Each retained edge set is a vertex-disjoint union of paths of positive
    even length.  The stored multiplier is 2^(number of components)*(n-2r)!;
    division by 2^(n-1) is performed during evaluation.
    """
    edge_count = 2 * level
    edges = tuple(combinations(range(n), 2))
    supports = []
    for chosen in combinations(range(len(edges)), edge_count):
        neighbours: dict[int, list[int]] = defaultdict(list)
        for edge_index in chosen:
            left, right = edges[edge_index]
            neighbours[left].append(right)
            neighbours[right].append(left)
        if any(len(adjacent) > 2 for adjacent in neighbours.values()):
            continue

        unseen = set(neighbours)
        paths = []
        valid = True
        while unseen:
            seed = min(unseen)
            component = {seed}
            frontier = [seed]
            while frontier:
                vertex = frontier.pop()
                for nxt in neighbours[vertex]:
                    if nxt not in component:
                        component.add(nxt)
                        frontier.append(nxt)
            unseen.difference_update(component)
            component_edges = sum(len(neighbours[v]) for v in component) // 2
            endpoints = sorted(v for v in component if len(neighbours[v]) == 1)
            if component_edges % 2 or len(endpoints) != 2 or component_edges != len(component) - 1:
                valid = False
                break
            path = [endpoints[0]]
            previous = -1
            current = endpoints[0]
            while current != endpoints[1]:
                nxt = next(v for v in neighbours[current] if v != previous)
                path.append(nxt)
                previous, current = current, nxt
            paths.append(tuple(path))
        if valid:
            multiplier = (1 << len(paths)) * factorial(n - edge_count)
            supports.append((multiplier, tuple(paths)))
    return tuple(supports)


def needle_layer(
    adjacency: list[list[int]],
    supports: tuple[tuple[int, tuple[tuple[int, ...], ...]], ...],
) -> Fraction:
    total = 0
    for multiplier, paths in supports:
        orientation = 1
        for path in paths:
            for left, right in zip(path, path[1:]):
                orientation *= 1 if adjacency[left][right] else -1
        total += multiplier * orientation
    return Fraction(total, 1 << (len(adjacency) - 1))


def one_defect_needles(adjacency: list[list[int]], h_value: int, drift_sum: int) -> int:
    """U_1 from the exact boundary identity U_1=(n-1)H+sum_e Delta_e H."""
    return (len(adjacency) - 1) * h_value + drift_sum


def q_coefficient(m: int, level: int) -> int:
    """Return [w^level] ((1+w)/(1-w))^m by finite binomial extraction."""
    return sum(
        comb(m, chosen) * comb(m + level - chosen - 1, level - chosen)
        for chosen in range(min(m, level) + 1)
    )


def falling(n: int, length: int) -> int:
    return factorial(n) // factorial(n - length)


def strongly_connected_components(vertices: int, arcs: set[tuple[int, int]]) -> list[int]:
    graph = [[] for _ in range(vertices)]
    reverse = [[] for _ in range(vertices)]
    for left, right in arcs:
        graph[left].append(right)
        reverse[right].append(left)
    seen: set[int] = set()
    order = []

    def visit(vertex: int) -> None:
        seen.add(vertex)
        for nxt in graph[vertex]:
            if nxt not in seen:
                visit(nxt)
        order.append(vertex)

    for vertex in range(vertices):
        if vertex not in seen:
            visit(vertex)
    seen.clear()
    sizes = []

    def collect(vertex: int) -> int:
        seen.add(vertex)
        size = 1
        for nxt in reverse[vertex]:
            if nxt not in seen:
                size += collect(nxt)
        return size

    for vertex in reversed(order):
        if vertex not in seen:
            sizes.append(collect(vertex))
    return sorted(sizes, reverse=True)


def tournament_fingerprint(values: list[Fraction], tie_order: list[int]) -> dict[str, object]:
    vertices = len(values)
    priority = {vertex: rank for rank, vertex in enumerate(tie_order)}
    arcs = set()
    for left in range(vertices):
        for right in range(left + 1, vertices):
            if values[left] > values[right] or (
                values[left] == values[right] and priority[left] < priority[right]
            ):
                arcs.add((left, right))
            else:
                arcs.add((right, left))
    scores = [sum((vertex, other) in arcs for other in range(vertices) if other != vertex) for vertex in range(vertices)]
    triangles = 0
    for a, b, c in combinations(range(vertices), 3):
        triangles += int(
            ((a, b) in arcs and (b, c) in arcs and (c, a) in arcs)
            or ((a, c) in arcs and (c, b) in arcs and (b, a) in arcs)
        )
    paths = [[0] * vertices for _ in range(1 << vertices)]
    for vertex in range(vertices):
        paths[1 << vertex][vertex] = 1
    for support in range(1 << vertices):
        for last, count in enumerate(paths[support]):
            if not count:
                continue
            for nxt in range(vertices):
                if not (support >> nxt) & 1 and (last, nxt) in arcs:
                    paths[support | (1 << nxt)][nxt] += count
    return {
        "arcs": arcs,
        "score_histogram": dict(sorted(Counter(scores).items())),
        "directed_triangles": triangles,
        "scc_sizes": strongly_connected_components(vertices, arcs),
        "hamiltonian_paths": sum(paths[-1]),
    }


def partition_stats(rows: list[dict[str, object]], key) -> tuple[int, int, int]:
    fibres: dict[object, int] = defaultdict(int)
    for row in rows:
        fibres[key(row)] += 1
    total_pairs = comb(len(rows), 2)
    separated = total_pairs - sum(comb(size, 2) for size in fibres.values())
    return len(fibres), max(fibres.values()), separated


def main() -> None:
    if not __debug__:
        raise RuntimeError("THM-848 verification requires assertions; do not run with python -O")
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=DEFAULT_SMALL_ATLAS)
    parser.add_argument("--n7-atlas", type=Path, default=DEFAULT_N7_ATLAS)
    parser.add_argument("--flow", type=Path, default=DEFAULT_FLOW)
    parser.add_argument(
        "--n8-atlas",
        type=Path,
        help="optional THM-828/843 raw uint16 rank atlas; SHA256 is mandatory",
    )
    args = parser.parse_args()

    small = json.loads(args.small_atlas.read_text())
    n7 = json.loads(args.n7_atlas.read_text())
    flow = json.loads(args.flow.read_text())
    atlas_by_n = {entry["n"]: entry for entry in small["sizes"]}
    atlas_by_n[7] = n7
    flow_by_n = {
        entry["n"]: {node["node_id"]: node for node in entry["nodes_in_flow_order"]}
        for entry in flow["sizes"]
        if 4 <= entry["n"] <= 7
    }

    print("THM-848: H-DRIFT WALSH NEEDLE STALK - EXACT AUDIT")
    print("=" * 78)
    all_rows: dict[int, list[dict[str, object]]] = {}
    support_bank = {
        n: {level: even_path_forests(n, level) for level in range(1, (n - 1) // 2 + 1)}
        for n in range(4, 8)
    }
    for n, levels in support_bank.items():
        for level, supports in levels.items():
            isolated = n - 2 * level
            weighted_supports = sum(
                (multiplier // factorial(isolated)) ** 2
                for multiplier, _ in supports
            )
            expected = (
                factorial(n)
                // factorial(isolated)
                * q_coefficient(isolated, level)
            )
            assert weighted_supports == expected

    for n in range(4, 8):
        arc_count = comb(n, 2)
        mean_h = Fraction(factorial(n), 1 << (n - 1))
        rows = []
        for node in atlas_by_n[n]["nodes"]:
            adjacency = adjacency_from_tiling(node["tiling_masks"][0], n)
            h_value, gradient = flip_gradient(adjacency)
            c3 = cyclic_triangles(adjacency)
            layers = tuple(
                needle_layer(adjacency, support_bank[n][level])
                for level in range(1, (n - 1) // 2 + 1)
            )
            assert Fraction(h_value) == mean_h + sum(layers)
            drift_sum = sum(gradient)
            predicted_sum = -4 * sum((level + 1) * value for level, value in enumerate(layers))
            assert drift_sum == predicted_sum
            drift = Fraction(drift_sum, arc_count)
            diffusion = Fraction(sum(value * value for value in gradient), arc_count)
            u1 = one_defect_needles(adjacency, h_value, drift_sum)
            assert u1 >= 0
            rows.append(
                {
                    "id": node["id"],
                    "mask": node["tiling_masks"][0],
                    "phase": flow_by_n[n][node["id"]]["category"],
                    "depth": node["local_depth"],
                    "H": h_value,
                    "c3": c3,
                    "score": score_sequence(adjacency),
                    "layers": layers,
                    "gradient": gradient,
                    "drift": drift,
                    "diffusion": diffusion,
                    "U1": u1,
                }
            )
        all_rows[n] = rows
        gradient_cells = len({row["gradient"] for row in rows})
        away = sum((Fraction(row["H"]) - mean_h) * row["drift"] > 0 for row in rows)
        print(
            f"n={n}: nodes={len(rows)}, layers={len(rows[0]['layers'])}, "
            f"forest supports={[len(support_bank[n][j]) for j in support_bank[n]]}, "
            f"gradient cells={gradient_cells}, away-from-mean nodes={away}"
        )
        assert gradient_cells == len(rows)

    print("\nLOW-SIZE CLOSED FORMS")
    print("-" * 78)
    for row in all_rows[4]:
        assert row["drift"] == Fraction(2, 3) * (3 - row["H"])
        assert row["U1"] == 12 - row["H"]
    for row in all_rows[5]:
        assert row["drift"] == Fraction(6 * row["c3"] - 4 * row["H"] + 15, 5)
        assert row["U1"] == 12 * row["c3"] - 4 * row["H"] + 30
    for row in all_rows[6]:
        assert row["drift"] == Fraction(24 * row["c3"] - 8 * row["H"] + 60, 15)
        assert row["U1"] == 24 * row["c3"] - 3 * row["H"] + 60
    print("n=4: b_H=(2/3)(3-H), U1=12-H: PASS")
    print("n=5: b_H=(6c3-4H+15)/5, U1=12c3-4H+30: PASS")
    print("n=6: b_H=(24c3-8H+60)/15, U1=24c3-3H+60: PASS")
    for n in (4, 5, 6):
        profile: dict[int, set[Fraction]] = defaultdict(set)
        for row in all_rows[n]:
            profile[row["H"]].add(row["drift"])
        rendered = ", ".join(
            f"{h}:{'/'.join(map(str, sorted(values)))}" for h, values in sorted(profile.items())
        )
        print(f"n={n} H-fibre drift profile: {rendered}")

    print("\nEXACT WALSH-LEVEL ENERGY EGF AND N=8 CORRECTION")
    print("-" * 78)
    n8_ratios = tuple(
        Fraction(q_coefficient(8 - 2 * level, level), falling(8, 2 * level))
        for level in range(1, 4)
    )
    assert n8_ratios == (Fraction(3, 14), Fraction(2, 105), Fraction(1, 1680))
    assert sum(n8_ratios) == Fraction(131, 560)
    assert Fraction(49752, factorial(8)) - 1 == Fraction(131, 560)
    assert Fraction(59, 252) != Fraction(131, 560)

    # At level six on eight vertices the forests are P7+isolate or P5 disjoint P3.
    weighted_p7 = 4 * 8 * factorial(7) // 2
    weighted_p5_p3 = 16 * comb(8, 5) * factorial(5) // 2 * factorial(3) // 2
    weighted_level6 = weighted_p7 + weighted_p5_p3
    assert weighted_level6 == falling(8, 6) * q_coefficient(2, 3) == 241920
    print("weighted forest EGF: exp(x+2x^3 z/(1-x^2 z)); Q(w)=(1+w)/(1-w)")
    print("E_(2r)/mu_n^2=[w^r]Q(w)^(n-2r)/(n)_(2r): PASS for all supports n=4..7")
    print("n=8 coefficients [w^r]Q(w)^(8-2r): (12,32,12)")
    print("n=8 level ratios (E2,E4,E6)/mu^2: (3/14,2/105,1/1680)")
    print("n=8 variance ratio: 131/560 = 49752/8!-1; old 59/252: REFUTED")

    print("\nN=7 FIRST DEGREE-SIX SPLIT")
    print("-" * 78)
    rows = all_rows[7]
    for row in rows:
        h6 = row["layers"][2]
        assert row["drift"] == Fraction(105 - 8 * row["H"] + 60 * row["c3"] - 4 * h6, 21)

    refinements = [
        ("H", lambda row: row["H"]),
        ("(H,c3)", lambda row: (row["H"], row["c3"])),
        ("(H,c3,phase)", lambda row: (row["H"], row["c3"], row["phase"])),
        ("(H,score)", lambda row: (row["H"], row["score"])),
        (
            "(H,c3,score,phase,depth)",
            lambda row: (row["H"], row["c3"], row["score"], row["phase"], row["depth"]),
        ),
        ("Walsh needle stalk", lambda row: row["layers"]),
        ("directional gradient multiset", lambda row: row["gradient"]),
    ]
    refinement_rows = []
    for name, key in refinements:
        cells, max_fibre, separated = partition_stats(rows, key)
        fibres: dict[object, list[dict[str, object]]] = defaultdict(list)
        for row in rows:
            fibres[key(row)].append(row)
        drift_splits = sum(len({item["drift"] for item in fibre}) > 1 for fibre in fibres.values())
        diffusion_splits = sum(len({item["diffusion"] for item in fibre}) > 1 for fibre in fibres.values())
        print(
            f"{name:31s}: cells={cells:3d}, max fibre={max_fibre}, "
            f"drift-split={drift_splits:2d}, diffusion-split={diffusion_splits:2d}, "
            f"separated pairs={separated}"
        )
        refinement_rows.append((name, cells, separated))

    witness_ids = {"n7-a106", "n7-a118"}
    witnesses = sorted((row for row in rows if row["id"] in witness_ids), key=lambda row: row["id"])
    assert len(witnesses) == 2
    common_fields = ("H", "c3", "score", "phase", "depth")
    assert all(witnesses[0][field] == witnesses[1][field] for field in common_fields)
    assert witnesses[0]["drift"] != witnesses[1]["drift"]
    print("\nstrong same-address witness:")
    for row in witnesses:
        h2, h4, h6 = row["layers"]
        print(
            f"  {row['id']}: mask={row['mask']}, H={row['H']}, c3={row['c3']}, "
            f"score={row['score']}, phase={row['phase']}, depth={row['depth']}, "
            f"(H2,H4,H6)=({h2},{h4},{h6}), b_H={row['drift']}"
        )

    layer_fibres: dict[object, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        layer_fibres[row["layers"]].append(row)
    first_diffusion = next(
        fibre for fibre in layer_fibres.values() if len({row["diffusion"] for row in fibre}) > 1
    )
    print("\nfirst-moment boundary of the Walsh stalk:")
    print(f"  common layers={first_diffusion[0]['layers']}")
    print(
        "  diffusion values="
        + str([(row["id"], str(row["diffusion"])) for row in first_diffusion])
    )

    print("\nTOURNAMENT ANALYSIS OF CANDIDATE NODE COLOURS")
    print("-" * 78)
    tie_order = list(range(len(refinement_rows)))
    raw_values = [Fraction(separated) for _, _, separated in refinement_rows]
    economy_values = [
        Fraction(separated, max(1, ceil(log2(cells))))
        for _, cells, separated in refinement_rows
    ]
    raw = tournament_fingerprint(raw_values, tie_order)
    economy = tournament_fingerprint(economy_values, tie_order)
    flips = len(raw["arcs"].symmetric_difference(economy["arcs"])) // 2
    for index, (name, cells, separated) in enumerate(refinement_rows):
        print(
            f"  v{index} {name:31s}: cells={cells:3d}, separated={separated:5d}, "
            f"economy={economy_values[index]}"
        )
    print(
        "raw fingerprint: score=%s, C3=%s, SCC=%s, HP=%s"
        % (
            raw["score_histogram"],
            raw["directed_triangles"],
            raw["scc_sizes"],
            raw["hamiltonian_paths"],
        )
    )
    print(
        "economy fingerprint: score=%s, C3=%s, SCC=%s, HP=%s, edge flips=%s"
        % (
            economy["score_histogram"],
            economy["directed_triangles"],
            economy["scc_sizes"],
            economy["hamiltonian_paths"],
            flips,
        )
    )

    print("\nOPTIONAL HASH-GUARDED N=8 GRADIENT CLASSIFIER")
    print("-" * 78)
    if args.n8_atlas is None:
        print("SKIPPED: pass --n8-atlas /tmp/n8_merged_node_rank_u16.bin")
    else:
        representatives, digest = certified_n8_representatives(args.n8_atlas)
        gradients = {
            flip_gradient(adjacency_from_legacy_tiling(mask, 8))[1]
            for mask in representatives
        }
        assert len(representatives) == N8_NODE_COUNT
        assert len(gradients) == N8_NODE_COUNT
        print(f"atlas SHA256={digest}")
        print("merged ranks=3528, directional-gradient cells=3528, max fibre=1")
        print("one H-seeded multiplicity-aware refinement round: COMPLETE at n=8")

    print("\nALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
