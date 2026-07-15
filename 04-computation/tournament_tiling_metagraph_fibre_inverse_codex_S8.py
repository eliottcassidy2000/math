#!/usr/bin/env python3
"""Bidirectional API for fixed-path tilings and merged metagraph nodes.

The forward map is the literal quotient

    tiling -> labelled tournament -> canonical class -> converse orbit -> node.

The inverse does not scan the tiling cube.  For every unmerged class carried by
a node, it decodes one canonical tournament representative, enumerates its
directed Hamiltonian paths, relabels each path to the explorer's fixed path,
and reads the remaining arcs as a staircase tiling.  Automorphism-related paths
give the same tiling, so a class contributes exactly H(T)/|Aut(T)| tilings.

The committed HYP-6825 atlases remain the source of the objective node order.
This module turns their stored relation into explicit functions and audits the
direct Hamiltonian-path inverse against every fibre at n=3,...,7.

Examples:

  python3 04-computation/tournament_tiling_metagraph_fibre_inverse_codex_S8.py \
      --n 7 --tiling 615
  python3 04-computation/tournament_tiling_metagraph_fibre_inverse_codex_S8.py \
      --n 7 --node n7-a267
  python3 04-computation/tournament_tiling_metagraph_fibre_inverse_codex_S8.py \
      --audit-all --output result.out --json result.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import permutations
from pathlib import Path
from typing import Iterable

from merged_metagraph_lines_n7_klein_S161 import canon_tournament
from tournament_tiling_metagraph_address_codex_S4 import (
    arc,
    automorphisms,
    canonical,
    converse,
    hamiltonian_paths,
    pair_index,
    pairs,
    tile_schema,
    tiling_tournament,
)


SMALL_ATLAS = Path("05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json")
N7_ATLAS = Path("05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json")


def pairmask_to_adjacency(mask: int, n: int) -> list[list[int]]:
    adjacency = [[0] * n for _ in range(n)]
    for i, j in pairs(n):
        value = (mask >> pair_index(n)[(i, j)]) & 1
        adjacency[i][j] = value
        adjacency[j][i] = 1 - value
    return adjacency


def adjacency_to_pairmask(adjacency: list[list[int]]) -> int:
    n = len(adjacency)
    mask = 0
    for i, j in pairs(n):
        assert adjacency[i][j] + adjacency[j][i] == 1
        mask |= int(adjacency[i][j]) << pair_index(n)[(i, j)]
    return mask


def canonical_class_code(tournament: int, n: int) -> int:
    """Return the class code in the committed atlas convention for this size."""
    if n <= 6:
        return canonical(tournament, n)
    if n == 7:
        key = canon_tournament(pairmask_to_adjacency(tournament, n), n)
        return int.from_bytes(key, "big")
    raise ValueError("the committed merged-node atlas currently covers 3 <= n <= 7")


def decode_class_code(code: int, n: int) -> int:
    """Decode one stored canonical code to a representative tournament mask."""
    if n <= 6:
        return code
    if n == 7:
        row_bytes = (n + 7) // 8
        raw = code.to_bytes(n * row_bytes, "big")
        adjacency = []
        for i in range(n):
            value = int.from_bytes(raw[i * row_bytes : (i + 1) * row_bytes], "big")
            adjacency.append([(value >> (n - 1 - j)) & 1 for j in range(n)])
        assert all(adjacency[i][i] == 0 for i in range(n))
        return adjacency_to_pairmask(adjacency)
    raise ValueError("the committed merged-node atlas currently covers 3 <= n <= 7")


def hamiltonian_path_words(tournament: int, n: int) -> Iterable[tuple[int, ...]]:
    """Yield every directed Hamiltonian vertex word of a labelled tournament."""
    for path in permutations(range(n)):
        if all(arc(tournament, n, path[i], path[i + 1]) for i in range(n - 1)):
            yield path


def path_to_tiling(tournament: int, n: int, path: tuple[int, ...]) -> int:
    """Relabel ``path`` to 0->...->n-1 and read the off-path tile bits."""
    if len(path) != n or set(path) != set(range(n)):
        raise ValueError("path must be a permutation of the tournament vertices")
    if not all(arc(tournament, n, path[i], path[i + 1]) for i in range(n - 1)):
        raise ValueError("path is not a directed Hamiltonian path")
    tiles, _ = tile_schema(n)
    tile_mask = 0
    for bit, (x, y) in enumerate(tiles):
        i, j = n - x, n - y
        assert j - i >= 2
        # Tile off means the earlier fixed-path vertex points to the later one.
        tile_mask |= int(not arc(tournament, n, path[i], path[j])) << bit
    return tile_mask


@dataclass(frozen=True)
class ClassFibre:
    class_code: int
    hamiltonian_paths: int
    automorphisms: int
    tiling_masks: tuple[int, ...]

    @property
    def predicted_size(self) -> int:
        return self.hamiltonian_paths // self.automorphisms


def class_to_tilings(class_code: int, n: int) -> ClassFibre:
    """Generate a complete unmerged class fibre from Hamiltonian paths."""
    representative = decode_class_code(class_code, n)
    multiplicity: Counter[int] = Counter()
    path_count = 0
    for path in hamiltonian_path_words(representative, n):
        path_count += 1
        multiplicity[path_to_tiling(representative, n, path)] += 1
    aut = automorphisms(representative, n)
    assert path_count == hamiltonian_paths(representative, n)
    assert path_count % aut == 0
    assert multiplicity
    assert set(multiplicity.values()) == {aut}
    tilings = tuple(sorted(multiplicity))
    assert len(tilings) == path_count // aut
    return ClassFibre(class_code, path_count, aut, tilings)


class MetagraphFibreAtlas:
    """Read-only exact correspondence for the committed n=3,...,7 atlases."""

    def __init__(self, small_path: Path = SMALL_ATLAS, n7_path: Path = N7_ATLAS):
        small = json.loads(small_path.read_text())
        n7 = json.loads(n7_path.read_text())
        self._sizes = {int(record["n"]): record for record in small["sizes"]}
        self._sizes[7] = n7
        self._node_by_rank: dict[int, dict[int, dict]] = {}
        self._node_by_id: dict[int, dict[str, dict]] = {}
        self._node_by_merged_code: dict[int, dict[int, dict]] = {}
        self._tiling_record: dict[int, dict[int, dict]] = {}
        self._class_fibre: dict[int, dict[int, tuple[int, ...]]] = {}
        for n, size in self._sizes.items():
            nodes = {int(node["rank"]): node for node in size["nodes"]}
            self._node_by_rank[n] = nodes
            self._node_by_id[n] = {node["id"]: node for node in nodes.values()}
            by_code = {}
            for node in nodes.values():
                codes = tuple(int(code, 16) for code in node["canonical_orbit_codes"])
                for code in codes:
                    by_code[code] = node
                by_code[min(codes)] = node
            self._node_by_merged_code[n] = by_code

            if n <= 6:
                records = {int(record["mask"]): record for record in size["tiling_map"]}
            else:
                records = {
                    mask: {
                        "mask": mask,
                        "class_code": hex(int(size["class_code_by_mask"][mask])),
                        "node_rank": int(size["node_rank_by_mask"][mask]),
                        "node_id": nodes[int(size["node_rank_by_mask"][mask])]["id"],
                        "global_index": int(size["global_index_by_mask"][mask]),
                        "fibre_index": int(size["fibre_index_by_mask"][mask]),
                        "popcount": mask.bit_count(),
                    }
                    for mask in range(int(size["tilings"]))
                }
            self._tiling_record[n] = records
            class_fibres: dict[int, list[int]] = defaultdict(list)
            for mask, record in records.items():
                class_fibres[int(record["class_code"], 16)].append(mask)
            self._class_fibre[n] = {
                code: tuple(sorted(masks)) for code, masks in class_fibres.items()
            }

    def size(self, n: int) -> dict:
        if n not in self._sizes:
            raise ValueError("the committed merged-node atlas currently covers 3 <= n <= 7")
        return self._sizes[n]

    def node(self, n: int, node: str | int) -> dict:
        if isinstance(node, str):
            try:
                return self._node_by_id[n][node]
            except KeyError as error:
                raise KeyError(f"unknown n={n} node id {node!r}") from error
        try:
            return self._node_by_rank[n][int(node)]
        except KeyError as error:
            raise KeyError(f"unknown n={n} node rank {node}") from error

    def tiling_to_node(self, n: int, tile_mask: int) -> dict:
        """Compute the quotient map, independently checking the stored lookup."""
        size = self.size(n)
        if not 0 <= tile_mask < int(size["tilings"]):
            raise ValueError(f"n={n} tiling mask must be in [0,{int(size['tilings']) - 1}]")
        tiles, _ = tile_schema(n)
        tournament = tiling_tournament(tile_mask, n, tiles)
        cls = canonical_class_code(tournament, n)
        opp = canonical_class_code(converse(tournament, n), n)
        merged = min(cls, opp)
        node = self._node_by_merged_code[n][merged]
        stored = self._tiling_record[n][tile_mask]
        assert int(stored["class_code"], 16) == cls
        assert int(stored["node_rank"]) == int(node["rank"])
        return {
            "n": n,
            "tile_mask": tile_mask,
            "bit_word_lsb_first": "".join(
                str((tile_mask >> bit) & 1) for bit in range(int(size["tile_count"]))
            ),
            "class_code": hex(cls),
            "converse_class_code": hex(opp),
            "merged_orbit_code": hex(merged),
            "node_id": node["id"],
            "node_rank": int(node["rank"]),
            "global_index": int(stored["global_index"]),
            "fibre_index": int(stored["fibre_index"]),
        }

    def node_to_tilings(self, n: int, node: str | int) -> list[dict]:
        """Generate the inverse fibre directly from class Hamiltonian paths."""
        node_record = self.node(n, node)
        masks = set()
        class_by_mask = {}
        for text in node_record["canonical_orbit_codes"]:
            code = int(text, 16)
            fibre = class_to_tilings(code, n)
            for mask in fibre.tiling_masks:
                assert mask not in masks
                masks.add(mask)
                class_by_mask[mask] = code
        ordered = sorted(masks, key=lambda mask: int(self._tiling_record[n][mask]["fibre_index"]))
        assert ordered == list(node_record["tiling_masks"])
        return [
            {
                "n": n,
                "node_id": node_record["id"],
                "node_rank": int(node_record["rank"]),
                "fibre_index": int(self._tiling_record[n][mask]["fibre_index"]),
                "global_index": int(self._tiling_record[n][mask]["global_index"]),
                "tile_mask": mask,
                "bit_word_lsb_first": "".join(
                    str((mask >> bit) & 1) for bit in range(int(self._sizes[n]["tile_count"]))
                ),
                "class_code": hex(class_by_mask[mask]),
            }
            for mask in ordered
        ]

    def stored_class_fibre(self, n: int, class_code: int) -> tuple[int, ...]:
        return self._class_fibre[n][class_code]


def digest_masks(masks: Iterable[int]) -> str:
    text = ",".join(str(mask) for mask in masks).encode()
    return hashlib.sha256(text).hexdigest()


def audit(atlas: MetagraphFibreAtlas, exact_forward_n7: bool = True) -> dict:
    result = {"schema_version": 1, "theorem": "THM-781", "sizes": []}
    for n in range(3, 8):
        size = atlas.size(n)
        class_records = []
        node_records = []
        duplicate_multiplicity_failures = 0
        class_fibre_failures = 0
        size_formula_failures = 0
        all_codes = sorted(atlas._class_fibre[n])
        generated_by_class = {}
        for code in all_codes:
            direct = class_to_tilings(code, n)
            generated_by_class[code] = direct.tiling_masks
            stored = atlas.stored_class_fibre(n, code)
            class_fibre_failures += direct.tiling_masks != stored
            size_formula_failures += len(stored) != direct.predicted_size
            duplicate_multiplicity_failures += direct.hamiltonian_paths != (
                direct.automorphisms * len(direct.tiling_masks)
            )
            class_records.append(
                {
                    "class_code": hex(code),
                    "hamiltonian_paths": direct.hamiltonian_paths,
                    "automorphisms": direct.automorphisms,
                    "predicted_tilings": direct.predicted_size,
                    "generated_digest": digest_masks(direct.tiling_masks),
                }
            )

        merged_formula_failures = 0
        inverse_roundtrip_failures = 0
        for node in size["nodes"]:
            codes = tuple(int(code, 16) for code in node["canonical_orbit_codes"])
            generated = tuple(
                sorted(mask for code in codes for mask in generated_by_class[code])
            )
            stored = tuple(sorted(int(mask) for mask in node["tiling_masks"]))
            merged_formula_failures += len(stored) != sum(
                len(generated_by_class[code]) for code in codes
            )
            inverse_roundtrip_failures += generated != stored
            node_records.append(
                {
                    "node_id": node["id"],
                    "self_converse": bool(node["self_converse"]),
                    "class_codes": [hex(code) for code in codes],
                    "predicted_tilings": sum(len(generated_by_class[code]) for code in codes),
                    "stored_tilings": len(stored),
                    "fibre_digest": digest_masks(stored),
                }
            )

        forward_failures = 0
        forward_checked = int(size["tilings"]) if n < 7 or exact_forward_n7 else 0
        if forward_checked:
            for mask in range(forward_checked):
                try:
                    atlas.tiling_to_node(n, mask)
                except (AssertionError, KeyError):
                    forward_failures += 1

        assert not any(
            (
                duplicate_multiplicity_failures,
                class_fibre_failures,
                size_formula_failures,
                merged_formula_failures,
                inverse_roundtrip_failures,
                forward_failures,
            )
        )
        class_sizes = [record["predicted_tilings"] for record in class_records]
        node_sizes = [record["stored_tilings"] for record in node_records]
        result["sizes"].append(
            {
                "n": n,
                "tilings": int(size["tilings"]),
                "classes": int(size["classes"]),
                "merged_nodes": int(size["merged_nodes"]),
                "forward_checked": forward_checked,
                "forward_failures": forward_failures,
                "class_fibre_failures": class_fibre_failures,
                "automorphism_multiplicity_failures": duplicate_multiplicity_failures,
                "class_size_formula_failures": size_formula_failures,
                "merged_size_formula_failures": merged_formula_failures,
                "inverse_roundtrip_failures": inverse_roundtrip_failures,
                "class_fibre_size_range": [min(class_sizes), max(class_sizes)],
                "node_fibre_size_range": [min(node_sizes), max(node_sizes)],
                "class_records": class_records,
                "node_records": node_records,
            }
        )
    return result


def render_audit(result: dict) -> str:
    lines = [
        "THM-781 HAMILTONIAN-PATH INVERSE OF MERGED METAGRAPH FIBRES",
        "=" * 78,
        "forward: tiling -> tournament -> canonical class -> converse orbit -> node",
        "inverse: node classes -> Hamiltonian paths / automorphisms -> tilings",
        "",
    ]
    for size in result["sizes"]:
        lines.extend(
            [
                f"n={size['n']}: tilings={size['tilings']} classes={size['classes']} "
                f"merged_nodes={size['merged_nodes']}",
                f"  exact forward checks={size['forward_checked']} failures={size['forward_failures']}",
                f"  class inverse failures={size['class_fibre_failures']} "
                f"automorphism multiplicity failures={size['automorphism_multiplicity_failures']}",
                f"  H/Aut formula failures={size['class_size_formula_failures']} "
                f"merged-sum failures={size['merged_size_formula_failures']} "
                f"inverse round-trip failures={size['inverse_roundtrip_failures']}",
                f"  class fibre size range={tuple(size['class_fibre_size_range'])} "
                f"node fibre size range={tuple(size['node_fibre_size_range'])}",
            ]
        )
    n7 = result["sizes"][-1]
    a267 = next(node for node in n7["node_records"] if node["node_id"] == "n7-a267")
    lines.extend(
        [
            "",
            "DISTINGUISHED PRIME-SEVEN PULLBACK",
            f"  n7-a267 classes={tuple(a267['class_codes'])} self_converse={a267['self_converse']}",
            f"  inverse fibre size={a267['stored_tilings']} digest={a267['fibre_digest']}",
            "",
            "TOURNAMENT ANALYSIS",
            "  vertices: fixed Hamiltonian paths (inverse generator), not runners",
            "  observable: whether two paths normalize to the same tiling",
            "  switch: quotient by Aut(T); each equivalence block has size |Aut(T)|",
            "  tie Hamiltonian path: lexicographic path order",
            "  isomorphism-node loss: path representative and fixed-cut tiling, restored as the inverse fibre",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=SMALL_ATLAS)
    parser.add_argument("--n7-atlas", type=Path, default=N7_ATLAS)
    parser.add_argument("--n", type=int)
    parser.add_argument("--tiling", type=int)
    parser.add_argument("--node")
    parser.add_argument("--audit-all", action="store_true")
    parser.add_argument("--skip-forward-n7", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    atlas = MetagraphFibreAtlas(args.small_atlas, args.n7_atlas)

    if args.tiling is not None:
        if args.n is None:
            parser.error("--tiling requires --n")
        print(json.dumps(atlas.tiling_to_node(args.n, args.tiling), indent=2))
    if args.node is not None:
        if args.n is None:
            parser.error("--node requires --n")
        node: str | int = int(args.node) if args.node.isdigit() else args.node
        print(json.dumps(atlas.node_to_tilings(args.n, node), indent=2))
    if args.audit_all or (args.tiling is None and args.node is None):
        result = audit(atlas, exact_forward_n7=not args.skip_forward_n7)
        output = render_audit(result)
        print(output, end="")
        if args.output:
            args.output.write_text(output)
        if args.json:
            args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
