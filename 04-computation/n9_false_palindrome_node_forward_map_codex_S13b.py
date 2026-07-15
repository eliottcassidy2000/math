#!/usr/bin/env python3
"""Exact tiling -> ordinary class -> converse-merged node map for THM-828.

For each of the 116 size-nine collision endpoints this script constructs the
labelled nine-vertex tournament in the repository's fixed-Hamilton-path
chart, computes the lexicographically least upper-triangle adjacency mask
under all vertex permutations, computes the same code for the converse, and
uses their sorted pair as the merged-node key.

The canonicalizer is an exact subset dynamic program.  If vertex ``v`` is
placed first, the induced code on the remaining vertices occupies all higher
bits and the arcs from ``v`` occupy the first ``k-1`` bits.  Retaining every
order attaining the minimum handles automorphism ties exactly.  This is the
same upper-triangle integer convention as
``tournament_tiling_metagraph_address_codex_S4.py`` but avoids materializing
all ``9!`` permutations.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path


N = 9
TILES = tuple((x, y) for y in range(1, N - 1) for x in range(N, y + 1, -1) if x - y >= 2)
assert len(TILES) == 28
FULL_ARCS = (1 << (N * (N - 1) // 2)) - 1
FULL_TILES = (1 << len(TILES)) - 1


def adjacency(tile_mask: int) -> tuple[int, ...]:
    """Vertices 1..9, with the fixed Hamilton path v -> v-1."""
    adj = [0] * N
    for v in range(2, N + 1):
        adj[v - 1] |= 1 << (v - 2)
    for bit, (x, y) in enumerate(TILES):
        if (tile_mask >> bit) & 1:
            adj[x - 1] |= 1 << (y - 1)
        else:
            adj[y - 1] |= 1 << (x - 1)
    return tuple(adj)


def converse(adj: tuple[int, ...]) -> tuple[int, ...]:
    all_vertices = (1 << N) - 1
    return tuple(all_vertices ^ (1 << v) ^ adj[v] for v in range(N))


def code_for_order(adj: tuple[int, ...], order: tuple[int, ...]) -> int:
    value = 0
    bit = 0
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            value |= int(bool(adj[order[i]] & (1 << order[j]))) << bit
            bit += 1
    return value


def canonical(adj: tuple[int, ...]) -> tuple[int, tuple[tuple[int, ...], ...]]:
    """Return exact least code and all orders attaining it."""

    @lru_cache(maxsize=None)
    def solve(subset: int) -> tuple[int, tuple[tuple[int, ...], ...]]:
        size = subset.bit_count()
        if size == 0:
            return 0, ((),)
        best = None
        best_orders: list[tuple[int, ...]] = []
        for v in range(N):
            if not (subset >> v) & 1:
                continue
            child_code, tails = solve(subset ^ (1 << v))
            for tail in tails:
                low = 0
                for j, w in enumerate(tail):
                    low |= int(bool(adj[v] & (1 << w))) << j
                value = low | (child_code << (size - 1))
                order = (v,) + tail
                if best is None or value < best:
                    best = value
                    best_orders = [order]
                elif value == best:
                    best_orders.append(order)
        assert best is not None
        return best, tuple(best_orders)

    code, orders = solve((1 << N) - 1)
    assert all(code_for_order(adj, order) == code for order in orders)
    return code, orders


def hamiltonian_paths(adj: tuple[int, ...]) -> int:
    dp = {(1 << v, v): 1 for v in range(N)}
    for _ in range(1, N):
        nxt: dict[tuple[int, int], int] = defaultdict(int)
        for (used, v), count in dp.items():
            for w in range(N):
                bit = 1 << w
                if not used & bit and adj[v] & bit:
                    nxt[(used | bit, w)] += count
        dp = nxt
    return sum(dp.values())


def tournament_data(tile_mask: int) -> dict:
    adj = adjacency(tile_mask)
    code, orders = canonical(adj)
    converse_code, converse_orders = canonical(converse(adj))
    scores = tuple(sorted((row.bit_count() for row in adj), reverse=True))
    axis = sum((2 * s - (N - 1)) ** 2 for s in scores)
    c3 = N * (N - 1) * (N - 2) // 6 - sum(s * (s - 1) // 2 for s in scores)
    h = hamiltonian_paths(adj)
    aut = len(orders)
    assert len(converse_orders) == aut and h % aut == 0
    assert axis == 240 - 8 * c3
    return {
        "tile_mask": tile_mask,
        "canonical_code": code,
        "converse_code": converse_code,
        "merged_key": (min(code, converse_code), max(code, converse_code)),
        "self_converse": code == converse_code,
        "score_sequence": scores,
        "axis": axis,
        "c3": c3,
        "H": h,
        "aut": aut,
        "class_fibre_size_H_over_aut": h // aut,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--witnesses", type=Path,
                        default=Path("05-knowledge/results/mobius_cech_n9_exact_join_witnesses_codex_S13.tsv"))
    parser.add_argument("--json", type=Path,
                        default=Path("05-knowledge/results/n9_false_palindrome_node_forward_map_codex_S13b.json"))
    parser.add_argument("--output", type=Path,
                        default=Path("05-knowledge/results/n9_false_palindrome_node_forward_map_codex_S13b.out"))
    args = parser.parse_args()

    witness_rows = []
    with args.witnesses.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            witness_rows.append({key: int(row[key], 16) for key in ("D", "u", "v")})
    if len(witness_rows) != 58:
        raise RuntimeError("expected 58 THM-828 witnesses")

    endpoint_data: dict[int, dict] = {}
    pair_rows = []
    groups: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for pair_index, row in enumerate(witness_rows):
        u, v = row["u"], row["v"]
        if u in endpoint_data or v in endpoint_data:
            raise RuntimeError("THM-828 endpoints are not disjoint")
        du, dv = tournament_data(u), tournament_data(v)
        endpoint_data[u], endpoint_data[v] = du, dv
        # Grid reflection of a tiling realizes tournament converse.
        if du["canonical_code"] != dv["converse_code"] or dv["canonical_code"] != du["converse_code"]:
            raise RuntimeError("false-palindrome endpoints are not converse classes")
        if du["merged_key"] != dv["merged_key"]:
            raise RuntimeError("pair does not map to one merged node")
        pair = {
            "pair_index": pair_index,
            "D": f"0x{row['D']:x}",
            "u": f"0x{u:x}",
            "v": f"0x{v:x}",
            "u_canonical_code": f"0x{du['canonical_code']:09x}",
            "v_canonical_code": f"0x{dv['canonical_code']:09x}",
            "merged_node_key": [f"0x{z:09x}" for z in du["merged_key"]],
            "self_converse_node": du["self_converse"],
        }
        pair_rows.append(pair)
        groups[du["merged_key"]].append(pair)

    if len(endpoint_data) != 116:
        raise RuntimeError("endpoint map is incomplete")

    # The oriented P projection keeps both the node of the selected tiling and
    # the node of its all-free-tile complement.  This reconciles the two
    # convention-dependent one-coordinate counts: the source and partner
    # projections are not equinumerous, while their pairs are collision-free.
    complement_data: dict[int, dict] = {}
    projection_rows = []
    source_nodes = Counter()
    partner_nodes = Counter()
    ordered_p = Counter()
    unordered_p = Counter()
    for pair_index, row in enumerate(witness_rows):
        source_mask = row["u"]
        partner_mask = source_mask ^ FULL_TILES
        source = endpoint_data[source_mask]
        partner = complement_data.setdefault(partner_mask, tournament_data(partner_mask))
        source_key, partner_key = source["merged_key"], partner["merged_key"]
        ordered_key = (source_key, partner_key)
        unordered_key = tuple(sorted(ordered_key))
        source_nodes[source_key] += 1
        partner_nodes[partner_key] += 1
        ordered_p[ordered_key] += 1
        unordered_p[unordered_key] += 1
        projection_rows.append({
            "pair_index": pair_index,
            "D": f"0x{row['D']:x}",
            "source_mask": f"0x{source_mask:x}",
            "complement_partner_mask": f"0x{partner_mask:x}",
            "source_merged_node": [f"0x{z:09x}" for z in source_key],
            "partner_merged_node": [f"0x{z:09x}" for z in partner_key],
            "ordered_P": [[f"0x{z:09x}" for z in key] for key in ordered_key],
            "unordered_projected_P": [[f"0x{z:09x}" for z in key] for key in unordered_key],
        })
    source_multiplicity = Counter(source_nodes.values())
    partner_multiplicity = Counter(partner_nodes.values())
    if (len(source_nodes), source_multiplicity, len(partner_nodes), partner_multiplicity,
            len(ordered_p), len(unordered_p)) != (
            53, Counter({1: 48, 2: 5}), 54, Counter({1: 51, 2: 2, 3: 1}), 58, 58):
        raise RuntimeError("source/partner P projection census mismatch")
    if any(source == partner for source, partner in ordered_p):
        raise RuntimeError("a selected P cell is a node loop")
    if set(ordered_p.values()) != {1} or set(unordered_p.values()) != {1}:
        raise RuntimeError("P projection is not injective on the 58 source presentations")

    group_rows = []
    for key, pairs in sorted(groups.items()):
        masks = [int(p[e], 16) for p in pairs for e in ("u", "v")]
        class_counts = Counter(endpoint_data[mask]["canonical_code"] for mask in masks)
        representative = endpoint_data[masks[0]]
        expected_codes = {key[0]} if key[0] == key[1] else set(key)
        if set(class_counts) != expected_codes:
            raise RuntimeError("merged-node ordinary-class support mismatch")
        group_rows.append({
            "merged_node_key": [f"0x{z:09x}" for z in key],
            "pair_count": len(pairs),
            "endpoint_presentations": len(masks),
            "distinct_endpoint_tilings": len(set(masks)),
            "D_histogram": {f"0x{d:x}": count for d, count in sorted(Counter(int(p["D"], 16) for p in pairs).items())},
            "ordinary_class_presentation_counts": {f"0x{code:09x}": count for code, count in sorted(class_counts.items())},
            "self_converse": key[0] == key[1],
            "score_sequence": list(representative["score_sequence"]),
            "axis": representative["axis"],
            "c3": representative["c3"],
            "H": representative["H"],
            "aut": representative["aut"],
            "ordinary_class_fibre_size_H_over_aut": representative["class_fibre_size_H_over_aut"],
            "observed_fraction_of_ordinary_presentations": {
                f"0x{code:09x}": {
                    "numerator": count,
                    "denominator": representative["class_fibre_size_H_over_aut"],
                }
                for code, count in sorted(class_counts.items())
            },
            "pair_indices": [p["pair_index"] for p in pairs],
        })

    multiplicity = Counter(len(pairs) for pairs in groups.values())
    repeated_groups = sum(count for size, count in multiplicity.items() if size > 1)
    self_converse_nodes = sum(key[0] == key[1] for key in groups)
    self_converse_pairs = sum(len(pairs) for key, pairs in groups.items() if key[0] == key[1])
    self_converse_endpoints = 2 * self_converse_pairs
    repeated = [(key, pairs) for key, pairs in groups.items() if len(pairs) > 1]
    if not all(key[0] == key[1] and len(pairs) == 2 and
               {int(pair["D"], 16) for pair in pairs} == {0x4C41818}
               for key, pairs in repeated):
        raise RuntimeError("repeated-node classification mismatch")
    endpoints_json = []
    for mask, data in sorted(endpoint_data.items()):
        endpoints_json.append({
            "tile_mask": f"0x{mask:x}",
            "canonical_code": f"0x{data['canonical_code']:09x}",
            "converse_code": f"0x{data['converse_code']:09x}",
            "merged_node_key": [f"0x{z:09x}" for z in data["merged_key"]],
            "self_converse": data["self_converse"],
            "score_sequence": list(data["score_sequence"]),
            "axis": data["axis"],
            "c3": data["c3"],
            "H": data["H"],
            "aut": data["aut"],
            "class_fibre_size_H_over_aut": data["class_fibre_size_H_over_aut"],
        })

    result = {
        "schema_version": 1,
        "theorem": "THM-828",
        "canonical_code_convention": "least upper-triangle adjacency integer under S9",
        "endpoint_tilings": len(endpoint_data),
        "collision_pairs": len(pair_rows),
        "converse_merged_nodes": len(groups),
        "pair_multiplicity_by_merged_node": {str(size): count for size, count in sorted(multiplicity.items())},
        "merged_nodes_with_multiple_pairs": repeated_groups,
        "self_converse_merged_nodes": self_converse_nodes,
        "self_converse_pairs": self_converse_pairs,
        "self_converse_endpoints": self_converse_endpoints,
        "non_self_converse_merged_nodes": len(groups) - self_converse_nodes,
        "repeated_nodes_are_two_pairs_in_D_0x4c41818_and_self_converse": True,
        "all_pairs_are_converse_class_presentations": True,
        "source_complement_projection": {
            "source_merged_nodes": len(source_nodes),
            "source_node_multiplicity_histogram": {str(k): v for k, v in sorted(source_multiplicity.items())},
            "complement_partner_merged_nodes": len(partner_nodes),
            "partner_node_multiplicity_histogram": {str(k): v for k, v in sorted(partner_multiplicity.items())},
            "ordered_P_cells": len(ordered_p),
            "unordered_projected_P_cells": len(unordered_p),
            "P_cell_multiplicity_histogram": {"1": 58},
            "node_loops": 0,
            "rows": projection_rows,
        },
        "pairs": pair_rows,
        "merged_node_groups": sorted(group_rows, key=lambda x: (-x["pair_count"], x["merged_node_key"])),
        "endpoints": endpoints_json,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n")

    lines = [
        "THM-828 SIZE-NINE FORWARD MAP: TILING -> CLASS -> MERGED NODE",
        "canonical code = least upper-triangle adjacency integer under all S9 relabellings",
        "",
        f"endpoints={len(endpoint_data)} pairs={len(pair_rows)} merged-nodes={len(groups)}",
        f"pair multiplicity per merged node={dict(sorted(multiplicity.items()))}",
        f"merged nodes carrying multiple false-palindrome pairs={repeated_groups}",
        f"self-converse nodes/pairs/endpoints={self_converse_nodes}/{self_converse_pairs}/{self_converse_endpoints}; "
        f"non-self-converse nodes={len(groups) - self_converse_nodes}",
        "all 5 repeated nodes are self-converse, carry exactly 2 pairs, and lie in D=0x4c41818: YES",
        "all 58 pairs join converse ordinary classes and map to one merged node: YES",
        f"source-node projection={len(source_nodes)} multiplicity={dict(sorted(source_multiplicity.items()))}",
        f"complement-partner projection={len(partner_nodes)} multiplicity={dict(sorted(partner_multiplicity.items()))}",
        f"ordered/unordered projected P cells={len(ordered_p)}/{len(unordered_p)}; all singleton; loops=0",
        "",
        "count node-key pair-count endpoints H Aut H/Aut D-hist ordinary-class-presentations",
    ]
    for rank, group in enumerate(sorted(group_rows, key=lambda x: (-x["pair_count"], x["merged_node_key"]))):
        lines.append(
            f"{rank:2d} {'/'.join(group['merged_node_key'])} {group['pair_count']:2d} "
            f"{group['endpoint_presentations']:2d} {group['H']:5d} {group['aut']:2d} "
            f"{group['ordinary_class_fibre_size_H_over_aut']:5d} {group['D_histogram']} "
            f"{group['ordinary_class_presentation_counts']}"
        )
    lines += [
        "",
        "INTERPRETATION",
        "A repeated group means several of the 58 false-palindrome pairs are distinct fixed-Hamilton-path",
        "presentations inside the same converse-merged tournament node; it is not a new isomorphism class.",
        "For 54 pairs the two reflected tilings are presentations of the SAME self-converse ordinary class;",
        "only 4 pairs have endpoints in two distinct ordinary classes exchanged by converse.",
    ]
    args.output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
