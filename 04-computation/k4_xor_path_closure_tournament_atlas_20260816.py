#!/usr/bin/env python3
"""Exact K4 XOR atlas for the pointed-P4/Boolean-closure sidecar.

A fixed transitive tournament supplies one orientation bit on each of the six
unordered pairs.  XOR with an undirected edge mask reverses precisely the
selected pairs, so the 64 masks biject with the 64 labelled tournaments.
The canonical pointed P4 and its C4 closure are audited together with the
coefficient-change obstruction from Boolean parity to oriented F_13 H^1.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json


EXPECTED_SEMANTIC_SHA256 = (
    "1d0dfe67dc8537728c3a71bed89c3109e9604130c4dcd24d56412f37d9fadf7d"
)
VERTICES = tuple(range(4))
EDGES = tuple(combinations(VERTICES, 2))
TRIANGLES = tuple(combinations(VERTICES, 3))
PATH_EDGES = ((0, 1), (1, 2), (2, 3))
CLOSURE_EDGE = (0, 3)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def edge_mask(edges: tuple[tuple[int, int], ...]) -> int:
    normalized = {tuple(sorted(edge)) for edge in edges}
    return sum(1 << index for index, edge in enumerate(EDGES) if edge in normalized)


def tournament(mask: int) -> frozenset[tuple[int, int]]:
    arcs = set()
    for index, (left, right) in enumerate(EDGES):
        arcs.add((right, left) if (mask >> index) & 1 else (left, right))
    return frozenset(arcs)


def scores(arcs: frozenset[tuple[int, int]]) -> tuple[int, ...]:
    return tuple(sum((vertex, other) in arcs for other in VERTICES if other != vertex)
                 for vertex in VERTICES)


def cyclic_triangle_count(arcs: frozenset[tuple[int, int]]) -> int:
    total = 0
    for a, b, c in TRIANGLES:
        forward = (a, b) in arcs and (b, c) in arcs and (c, a) in arcs
        backward = (a, c) in arcs and (c, b) in arcs and (b, a) in arcs
        total += int(forward or backward)
    return total


def strongly_connected(arcs: frozenset[tuple[int, int]]) -> bool:
    for source in VERTICES:
        seen = {source}
        while True:
            enlarged = seen | {
                target
                for start in seen
                for target in VERTICES
                if (start, target) in arcs
            }
            if enlarged == seen:
                break
            seen = enlarged
        if len(seen) != len(VERTICES):
            return False
    return True


def triangle_parity(mask: int) -> tuple[int, ...]:
    return tuple(
        sum((mask >> EDGES.index(tuple(sorted(edge)))) & 1
            for edge in combinations(triangle, 2)) % 2
        for triangle in TRIANGLES
    )


def cut_mask(subset: frozenset[int]) -> int:
    return edge_mask(tuple(
        edge for edge in EDGES if (edge[0] in subset) != (edge[1] in subset)
    ))


def canonical_record(mask: int) -> tuple[object, ...]:
    arcs = tournament(mask)
    return (
        mask,
        tuple(sorted(arcs)),
        scores(arcs),
        cyclic_triangle_count(arcs),
        strongly_connected(arcs),
        triangle_parity(mask),
    )


def main() -> None:
    # For a fixed orientation gauge, every six-bit reversal mask gives exactly
    # one labelled tournament, and XORing back recovers the mask.
    tournaments = tuple(tournament(mask) for mask in range(1 << len(EDGES)))
    require(len(set(tournaments)) == 64, "mask-to-tournament bijection")
    require(all(len(arcs) == 6 for arcs in tournaments), "one arc per pair")

    type_census = Counter()
    for arcs in tournaments:
        key = (tuple(sorted(scores(arcs))), cyclic_triangle_count(arcs),
               strongly_connected(arcs))
        type_census[key] += 1

    path_masks = {
        edge_mask(tuple(sorted((order[index], order[index + 1])))
                  for index in range(3))
        for order in permutations(VERTICES)
    }
    require(len(path_masks) == 12, ("Hamiltonian path masks", len(path_masks)))
    require({sum(triangle_parity(mask)) for mask in path_masks} == {2},
            "every P4 mask has two odd K4 triangles")

    cuts = {cut_mask(frozenset(subset))
            for size in range(5)
            for subset in combinations(VERTICES, size)}
    parity_kernel = {mask for mask in range(64)
                     if triangle_parity(mask) == (0, 0, 0, 0)}
    require(len(cuts) == 8 and cuts == parity_kernel,
            ("K4 cut kernel", len(cuts), len(parity_kernel)))

    path_mask = edge_mask(PATH_EDGES)
    closure_mask = path_mask ^ edge_mask((CLOSURE_EDGE,))
    path_record = canonical_record(path_mask)
    closure_record = canonical_record(closure_mask)

    require(path_record[2:] == ((2, 2, 1, 1), 2, True, (0, 1, 1, 0)),
            ("canonical path", path_record))
    require(closure_record[2:] == ((1, 2, 1, 2), 2, True, (0, 0, 0, 0)),
            ("canonical closure", closure_record))
    require(closure_mask == cut_mask(frozenset((0, 2))),
            "C4 closure mask is the alternating K4 cut")

    # Generalized-tournament typing before XOR: bidirected on mask edges and
    # missing off the mask.  The path and closure have different pair-status
    # censuses even though XOR sends both to honest tournaments.
    generalized_census = (
        ("P4", len(PATH_EDGES), 0, len(EDGES) - len(PATH_EDGES)),
        ("C4", len(PATH_EDGES) + 1, 0, len(EDGES) - len(PATH_EDGES) - 1),
        ("transitive_T4", 0, len(EDGES), 0),
    )

    # On the oriented C4 0->1->2->3->0, the all-one cochain has circulation
    # 4 mod 13, hence is nonzero in H^1(C4;F13).  Its 0/1 parity reduction has
    # circulation 0 mod 2 and is the alternating cut; XOR cannot retain the
    # characteristic-13 class.
    circulation = (4 % 13, 4 % 2)
    require(circulation == (4, 0), circulation)

    record = (
        tuple(sorted(type_census.items())),
        tuple(sorted(path_masks)),
        tuple(sorted(cuts)),
        path_record,
        closure_record,
        generalized_census,
        circulation,
        ("F13_H1_to_XOR_F2", "NO CANONICAL COEFFICIENT MAP; PARITY KILLS CIRCULATION"),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== K4 XOR/path-closure tournament atlas ==")
    print(f"labelled_masks=64;labelled_tournaments={len(set(tournaments))};edges={EDGES}")
    print(f"tournament_type_census={tuple(sorted(type_census.items()))}")
    print(f"hamiltonian_P4_masks={len(path_masks)};triangle_parity_weight=2")
    print(f"K4_cut_space={len(cuts)};triangle_parity_kernel={len(parity_kernel)}: EQUAL")
    print(f"canonical_P4={path_record}")
    print(f"canonical_C4_closure={closure_record};alternating_cut={{0,2}}|{{1,3}}")
    print(f"generalized_pair_census=(name,both,one,missing)={generalized_census}")
    print("xor_interpretation=bidirected mask reverses selected arcs of a fixed transitive T4 gauge")
    print(f"coefficient_boundary=(oriented_C4_sum_mod13,parity_sum_mod2)={circulation}")
    print("D5_boundary=no canonical F13-to-F2 coefficient transport; XOR forgets the nonzero oriented circulation")
    print(f"semantic_sha256={semantic}")
    print("scope=finite K4 Boolean/tournament atlas;no LRC current,D5 bridge,Keller flux,or row exclusion")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
