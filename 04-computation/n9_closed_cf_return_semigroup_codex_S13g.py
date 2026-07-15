#!/usr/bin/env python3
"""Exact closed-return transformation semigroup for THM-853.

The letters A,B,C are THM-846's endomorphisms F_R=R_10 Phi_9,10 of X_9.
A word is applied left-to-right, so AB means F_B after F_A.  This is a
repeated fixed-size 9->10->9 return algebra, not the increasing-size
continued-fraction cocycle; the latter also needs ambient size and phase.

The exact coordinate semigroup and its right Cayley graph are enumerated.
On THM-828's 58 rows, literal Q equality, old-pair descent, defect rank, and
THM-840 operation kernels are audited for every semigroup element.  A
partition-output Moore minimization records future Q equality under arbitrary
further closed returns.

Tournament Analysis uses information carriers rather than runners or arcs as
vertices.  Its pairwise observable is separated semigroup elements, its
switch is raw retention versus retention per logarithmic cost, and the input
carrier order is the tie Hamiltonian path.

Preserved: exact coordinate-copy maps, complement, reflection conjugation,
literal Q equality on the finite bank, and the rank-four defect images.
Destroyed/not computed: merged-node identity for all 8,940 output masks, LRC
metrics/owners/walls/loneliness, and genuine increasing-size CF evolution.
Challenged assumption: a single seam bit cannot be iterated as a universal
state; the closed returns form a 46,126-element contraction semigroup.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict, deque
from pathlib import Path

from n9_closed_cf_lift_face_kernels_codex_S13f import (
    DEFECT_BASIS,
    build_rho,
    read_witnesses,
    role_coordinate_maps,
)
from n9_defect_b3_continuation_purity_codex_S13d import (
    gf2_rank,
    q_cell,
    reflect,
    ta_fingerprint,
    tiles,
)


LETTERS = "ABC"
NCOORD = 28


def then(mapping: tuple[int, ...], generator: tuple[int, ...]) -> tuple[int, ...]:
    """Apply ``generator`` after ``mapping``."""
    return tuple(mapping[generator[target]] for target in range(NCOORD))


def enumerate_semigroup(generators: dict[str, tuple[int, ...]]):
    identity = tuple(range(NCOORD))
    word = {identity: ""}
    order = [identity]
    queue = deque([identity])
    while queue:
        mapping = queue.popleft()
        for letter in LETTERS:
            successor = then(mapping, generators[letter])
            if successor not in word:
                word[successor] = word[mapping] + letter
                order.append(successor)
                queue.append(successor)
    index = {mapping: i for i, mapping in enumerate(order)}
    adjacency = [
        tuple(index[then(mapping, generators[letter])] for letter in LETTERS)
        for mapping in order
    ]
    return order, word, adjacency


def mask_luts(mapping: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    patterns = [0] * NCOORD
    for target, source in enumerate(mapping):
        patterns[source] |= 1 << target
    tables = []
    for block in range(4):
        base = 7 * block
        table = [0] * 128
        for value in range(1, 128):
            low = value & -value
            table[value] = table[value ^ low] | patterns[base + low.bit_length() - 1]
        tables.append(tuple(table))
    return tuple(tables)


def transform_mask(mask: int, tables: tuple[tuple[int, ...], ...]) -> int:
    return (
        tables[0][mask & 127]
        | tables[1][(mask >> 7) & 127]
        | tables[2][(mask >> 14) & 127]
        | tables[3][(mask >> 21) & 127]
    )


def normalized_partition(values: list[object]) -> tuple[int, ...]:
    labels: dict[object, int] = {}
    result = []
    for value in values:
        if value not in labels:
            labels[value] = len(labels)
        result.append(labels[value])
    return tuple(result)


def fibre_profile(partition: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(Counter(partition).values()).items()))


def split_source_cells(source: tuple[int, ...], target: tuple[int, ...]) -> int:
    images: dict[int, set[int]] = defaultdict(set)
    for left, right in zip(source, target, strict=True):
        images[left].add(right)
    return sum(len(values) > 1 for values in images.values())


def strongly_connected_components(adjacency: list[tuple[int, ...]]):
    reverse = [[] for _ in adjacency]
    for source, row in enumerate(adjacency):
        for target in row:
            reverse[target].append(source)

    visited: set[int] = set()
    finish: list[int] = []
    for root in range(len(adjacency)):
        if root in visited:
            continue
        stack = [(root, 0)]
        visited.add(root)
        while stack:
            vertex, position = stack[-1]
            if position < len(LETTERS):
                target = adjacency[vertex][position]
                stack[-1] = (vertex, position + 1)
                if target not in visited:
                    visited.add(target)
                    stack.append((target, 0))
            else:
                finish.append(vertex)
                stack.pop()

    component = [-1] * len(adjacency)
    members: list[list[int]] = []
    for root in reversed(finish):
        if component[root] >= 0:
            continue
        label = len(members)
        stack = [root]
        component[root] = label
        group = []
        while stack:
            vertex = stack.pop()
            group.append(vertex)
            for source in reverse[vertex]:
                if component[source] < 0:
                    component[source] = label
                    stack.append(source)
        members.append(group)
    return component, members


def condensation(adjacency: list[tuple[int, ...]], component: list[int], count: int):
    outgoing = [set() for _ in range(count)]
    incoming = [set() for _ in range(count)]
    for source, row in enumerate(adjacency):
        for target in row:
            left, right = component[source], component[target]
            if left != right:
                outgoing[left].add(right)
                incoming[right].add(left)
    indegree = [len(row) for row in incoming]
    queue = deque(i for i, degree in enumerate(indegree) if degree == 0)
    topological = []
    while queue:
        left = queue.popleft()
        topological.append(left)
        for right in outgoing[left]:
            indegree[right] -= 1
            if indegree[right] == 0:
                queue.append(right)
    assert len(topological) == count
    return outgoing, topological


def counter_json(counter: Counter) -> dict[str, int]:
    return {str(key): counter[key] for key in sorted(counter)}


def profile_json(counter: Counter) -> dict[str, int]:
    return {
        ",".join(f"{size}x{number}" for size, number in profile): count
        for profile, count in sorted(counter.items())
    }


def run(witness_path: Path) -> dict[str, object]:
    rows = read_witnesses(witness_path)
    generators = role_coordinate_maps(build_rho())
    order, words, adjacency = enumerate_semigroup(generators)
    assert len(order) == 46126
    assert max(map(len, words.values())) == 22

    ranks = [len(set(mapping)) for mapping in order]
    rank_histogram = Counter(ranks)
    assert rank_histogram == {
        1: 6, 2: 19186, 3: 25146, 4: 1038, 5: 318, 6: 180, 7: 78,
        8: 36, 9: 32, 10: 54, 11: 12, 14: 14, 15: 16, 20: 4,
        21: 4, 27: 1, 28: 1,
    }
    length_histogram = Counter(map(len, words.values()))
    idempotent = [then(mapping, mapping) == mapping for mapping in order]
    idempotent_rank_histogram = Counter(
        rank for rank, truth in zip(ranks, idempotent, strict=True) if truth
    )
    assert sum(idempotent) == 3149
    assert idempotent_rank_histogram == {1: 6, 2: 3009, 3: 130, 4: 2, 27: 1, 28: 1}

    supports_by_rank: dict[int, set[tuple[int, ...]]] = defaultdict(set)
    for mapping, rank in zip(order, ranks, strict=True):
        supports_by_rank[rank].add(tuple(sorted(set(mapping))))
    support_histogram = {rank: len(values) for rank, values in supports_by_rank.items()}
    assert support_histogram == {
        1: 6, 2: 13, 3: 22, 4: 32, 5: 19, 6: 7, 7: 6, 8: 6, 9: 2,
        10: 4, 11: 2, 14: 3, 15: 3, 20: 2, 21: 2, 27: 1, 28: 1,
    }

    sigma = tuple(reflect(1 << bit, 9).bit_length() - 1 for bit in range(NCOORD))
    conjugate = lambda mapping: tuple(sigma[mapping[sigma[target]]] for target in range(NCOORD))
    reflection_fixed = [mapping for mapping in order if conjugate(mapping) == mapping]
    assert all(conjugate(mapping) in words for mapping in order)
    assert len(reflection_fixed) == 12
    reflection_orbits = (len(order) + len(reflection_fixed)) // 2
    assert reflection_orbits == 23069

    # Pure word laws and the complete two-letter rank table.
    def map_for_word(word: str) -> tuple[int, ...]:
        mapping = order[0]
        for letter in word:
            mapping = then(mapping, generators[letter])
        return mapping

    assert [len(set(map_for_word("A" * power))) for power in range(1, 8)] == [21, 15, 11, 8, 5, 4, 4]
    assert [len(set(map_for_word("C" * power))) for power in range(1, 8)] == [21, 15, 11, 8, 5, 4, 4]
    assert map_for_word("A" * 6) == map_for_word("A" * 7)
    assert map_for_word("C" * 6) == map_for_word("C" * 7)
    assert map_for_word("BB") == generators["B"]
    assert [len(set(map_for_word("AC" * power))) for power in range(1, 4)] == [15, 6, 3]
    assert [len(set(map_for_word("CA" * power))) for power in range(1, 4)] == [15, 6, 3]
    assert map_for_word("AC" * 3) == map_for_word("AC" * 4)
    assert map_for_word("CA" * 3) == map_for_word("CA" * 4)
    two_step_ranks = {
        left + right: len(set(map_for_word(left + right)))
        for left in LETTERS for right in LETTERS
    }
    assert two_step_ranks == {
        "AA": 15, "AB": 21, "AC": 15,
        "BA": 20, "BB": 27, "BC": 20,
        "CA": 15, "CB": 21, "CC": 15,
    }

    # Right Cayley flow, SCCs, and the six rank-one information horizons.
    component, members = strongly_connected_components(adjacency)
    outgoing, topological = condensation(adjacency, component, len(members))
    scc_size_histogram = Counter(map(len, members))
    assert len(members) == 31075
    assert scc_size_histogram == {1: 31058, 216: 6, 1252: 11}
    terminal = [label for label, targets in enumerate(outgoing) if not targets]
    assert len(terminal) == 6 and all(len(members[label]) == 1 for label in terminal)
    sink_mask = [0] * len(members)
    for bit, label in enumerate(terminal):
        sink_mask[label] = 1 << bit
    for label in reversed(topological):
        for target in outgoing[label]:
            sink_mask[label] |= sink_mask[target]
    reachable_sink_count = [sink_mask[component[i]].bit_count() for i in range(len(order))]
    reachable_sink_histogram = Counter(reachable_sink_count)
    assert reachable_sink_histogram == {1: 13032, 2: 31342, 3: 1746, 4: 4, 5: 1, 6: 1}
    assert reachable_sink_count[0] == 6

    sink_rows = []
    for bit, label in enumerate(terminal):
        vertex = members[label][0]
        mapping = order[vertex]
        assert ranks[vertex] == 1 and idempotent[vertex]
        source = mapping[0]
        sink_rows.append({
            "sink_bit": bit,
            "source_coordinate_index": source,
            "source_tile": list(tiles(9)[source]),
            "shortest_word": words[mapping],
        })
    sink_rows.sort(key=lambda row: row["source_tile"])
    assert [tuple(row["source_tile"]) for row in sink_rows] == [
        (5, 3), (6, 3), (6, 4), (7, 3), (7, 4), (7, 5),
    ]

    nontrivial_scc_profiles = Counter()
    for label, group in enumerate(members):
        if len(group) > 1:
            rank_set = tuple(sorted({ranks[vertex] for vertex in group}))
            nontrivial_scc_profiles[(len(group), rank_set, sink_mask[label].bit_count())] += 1
    assert nontrivial_scc_profiles == {(216, (3,), 3): 6, (1252, (2,), 2): 11}

    # Full finite-bank Q audit, accelerated by four seven-bit lookup tables.
    source_u = [row["u"] for row in rows]
    source_v = [row["v"] for row in rows]
    sectors = sorted({row["D"] for row in rows})
    q_cache: dict[int, int] = {}
    q_partitions: list[tuple[int, ...]] = []
    q_cell_counts = []
    q_profiles: Counter = Counter()
    old_pair_descent = []
    defect_ranks = []
    sector_image_counts = []
    output_masks: set[int] = set()
    full = (1 << NCOORD) - 1

    def q_of(mask: int) -> int:
        if mask not in q_cache:
            q_cache[mask] = q_cell(mask, 9)
        return q_cache[mask]

    for mapping in order:
        tables = mask_luts(mapping)
        images_u = [transform_mask(mask, tables) for mask in source_u]
        images_v = [transform_mask(mask, tables) for mask in source_v]
        output_masks.update(images_u)
        q_u = [q_of(mask) for mask in images_u]
        partition = normalized_partition(q_u)
        q_partitions.append(partition)
        q_cell_counts.append(len(set(partition)))
        q_profiles[fibre_profile(partition)] += 1
        old_pair_descent.append(sum(q_of(right) == left for left, right in zip(q_u, images_v, strict=True)))
        defect_images = [transform_mask(word, tables) for word in DEFECT_BASIS]
        defect_ranks.append(gf2_rank(defect_images))
        sector_image_counts.append(len({transform_mask(word, tables) for word in sectors}))

    assert len(output_masks) == 6138
    output_with_complements = output_masks | {mask ^ full for mask in output_masks}
    assert len(output_with_complements) == 8940
    full_q_words = sorted(
        (words[order[i]] or "I" for i, cells in enumerate(q_cell_counts) if cells == 58),
        key=lambda word: (0 if word == "I" else len(word), word),
    )
    full_q_and_descent_words = sorted(
        (
            words[order[i]] or "I"
            for i, (cells, descent) in enumerate(zip(q_cell_counts, old_pair_descent, strict=True))
            if cells == 58 and descent == 58
        ),
        key=lambda word: (0 if word == "I" else len(word), word),
    )
    assert len(full_q_words) == 10

    # Exact two-step cross-check against the independent scout.
    expected_two_step = {
        "AA": (37, ((1, 24), (2, 9), (4, 4)), 0),
        "AB": (58, ((1, 58),), 0),
        "AC": (48, ((1, 38), (2, 10)), 4),
        "BA": (58, ((1, 58),), 0),
        "BB": (58, ((1, 58),), 58),
        "BC": (58, ((1, 58),), 0),
        "CA": (48, ((1, 38), (2, 10)), 4),
        "CB": (58, ((1, 58),), 0),
        "CC": (37, ((1, 24), (2, 9), (4, 4)), 0),
    }
    two_step_q = {}
    index = {mapping: i for i, mapping in enumerate(order)}
    for word, expected in expected_two_step.items():
        vertex = index[map_for_word(word)]
        actual = (q_cell_counts[vertex], fibre_profile(q_partitions[vertex]), old_pair_descent[vertex])
        assert actual == expected
        two_step_q[word] = {
            "coordinate_rank": ranks[vertex],
            "Q_cells": actual[0],
            "Q_fibre_profile": {str(size): count for size, count in actual[1]},
            "old_pair_descent": actual[2],
        }

    # THM-840 operation kernels between Q partitions on every Cayley edge.
    operation_kernel_histograms = {}
    for position, letter in enumerate(LETTERS):
        histogram = Counter(
            split_source_cells(q_partitions[source], q_partitions[target])
            for source, row in enumerate(adjacency)
            for target in (row[position],)
        )
        operation_kernel_histograms[letter] = counter_json(histogram)

    # Moore minimization: states are equal only when their present Q equality
    # partition and all future partition outputs under A/B/C agree.
    initial_ids: dict[tuple[int, ...], int] = {}
    colours = []
    for partition in q_partitions:
        if partition not in initial_ids:
            initial_ids[partition] = len(initial_ids)
        colours.append(initial_ids[partition])
    refinement_counts = [len(initial_ids)]
    while True:
        identifiers: dict[tuple[int, int, int, int], int] = {}
        refined = []
        for vertex, row in enumerate(adjacency):
            signature = (colours[vertex], *(colours[target] for target in row))
            if signature not in identifiers:
                identifiers[signature] = len(identifiers)
            refined.append(identifiers[signature])
        refinement_counts.append(len(identifiers))
        if refined == colours:
            break
        colours = refined
    future_class_histogram = Counter(Counter(colours).values())

    loops_by_letter = {
        letter: sum(row[position] == source for source, row in enumerate(adjacency))
        for position, letter in enumerate(LETTERS)
    }
    rank_drop_histogram = Counter(
        ranks[source] - ranks[target]
        for source, row in enumerate(adjacency) for target in row
    )

    mapping_hash = hashlib.sha256(
        b"".join(bytes(mapping) for mapping in sorted(order))
    ).hexdigest()
    word_hash = hashlib.sha256(
        b"".join(bytes(mapping) + b"\0" + words[mapping].encode() + b"\n"
                 for mapping in sorted(order))
    ).hexdigest()

    tournament_analysis = ta_fingerprint([
        ("coordinate rank", ranks),
        ("visible support", [tuple(sorted(set(mapping))) for mapping in order]),
        ("reachable rank-one horizon", [sink_mask[component[i]] for i in range(len(order))]),
        ("defect image rank", defect_ranks),
        ("Q cell count", q_cell_counts),
        ("Q equality partition", q_partitions),
        ("old-pair Q descent", old_pair_descent),
    ])

    return {
        "schema_version": 1,
        "theorem": "THM-853",
        "word_convention": "letters applied left-to-right; AB = F_B after F_A",
        "semigroup": {
            "maps": len(order),
            "shortest_word_maximum": max(map(len, words.values())),
            "shortest_word_length_histogram": counter_json(length_histogram),
            "coordinate_rank_histogram": counter_json(rank_histogram),
            "distinct_supports_by_rank": {str(rank): support_histogram[rank] for rank in sorted(support_histogram)},
            "idempotents": sum(idempotent),
            "idempotent_rank_histogram": counter_json(idempotent_rank_histogram),
            "reflection_fixed_maps": len(reflection_fixed),
            "reflection_orbits": reflection_orbits,
            "mapping_sha256": mapping_hash,
            "shortest_word_sha256": word_hash,
            "two_step_coordinate_ranks": two_step_ranks,
        },
        "right_cayley_graph": {
            "vertices": len(order),
            "labelled_edges": len(order) * len(LETTERS),
            "loops_by_letter": loops_by_letter,
            "rank_drop_histogram": counter_json(rank_drop_histogram),
            "sccs": len(members),
            "scc_size_histogram": counter_json(scc_size_histogram),
            "nontrivial_scc_profiles": [
                {
                    "size": size,
                    "coordinate_ranks": list(rank_set),
                    "reachable_sinks": sinks,
                    "components": count,
                }
                for (size, rank_set, sinks), count in sorted(nontrivial_scc_profiles.items())
            ],
            "terminal_sinks": sink_rows,
            "reachable_sink_count_histogram": counter_json(reachable_sink_histogram),
        },
        "finite_bank": {
            "rows": len(rows),
            "distinct_output_masks": len(output_masks),
            "distinct_output_masks_with_complements": len(output_with_complements),
            "Q_cell_count_histogram": counter_json(Counter(q_cell_counts)),
            "Q_fibre_profile_histogram": profile_json(q_profiles),
            "old_pair_Q_descent_histogram": counter_json(Counter(old_pair_descent)),
            "Q_injective_shortest_words": full_q_words,
            "Q_injective_and_full_old_pair_descent_shortest_words": full_q_and_descent_words,
            "defect_image_rank_histogram": counter_json(Counter(defect_ranks)),
            "occupied_sector_image_count_histogram": counter_json(Counter(sector_image_counts)),
            "two_step_Q_census": two_step_q,
            "Q_operation_kernel_split_histograms": operation_kernel_histograms,
            "Q_partition_future_minimization": {
                "round_class_counts": refinement_counts,
                "final_classes": len(set(colours)),
                "final_class_size_histogram": counter_json(future_class_histogram),
            },
        },
        "tournament_analysis": tournament_analysis,
        "preservation": {
            "preserves": [
                "literal coordinate-copy composition", "complement",
                "reflection conjugation", "finite-bank Q equality",
                "finite-bank defect images",
            ],
            "destroys_or_defers": [
                "all-output merged-node identity", "LRC gaps", "owners",
                "walls", "metric loneliness", "increasing-size CF phase state",
            ],
            "challenged_assumption": (
                "iterating one seam bit does not close the return dynamics; the exact state is a contraction semigroup"
            ),
        },
    }


def render(result: dict[str, object]) -> str:
    semigroup = result["semigroup"]
    cayley = result["right_cayley_graph"]
    bank = result["finite_bank"]
    future = bank["Q_partition_future_minimization"]
    ta = result["tournament_analysis"]
    lines = [
        "THM-853 N=9 CLOSED CF RETURN SEMIGROUP",
        "=" * 70,
        f"maps={semigroup['maps']} max-shortest-word={semigroup['shortest_word_maximum']} idempotents={semigroup['idempotents']}",
        f"coordinate rank histogram={semigroup['coordinate_rank_histogram']}",
        f"idempotent rank histogram={semigroup['idempotent_rank_histogram']}",
        f"reflection fixed/orbits={semigroup['reflection_fixed_maps']}/{semigroup['reflection_orbits']}",
        f"two-step ranks={semigroup['two_step_coordinate_ranks']}",
        f"mapping/word SHA256={semigroup['mapping_sha256']}/{semigroup['shortest_word_sha256']}",
        "",
        "RIGHT CAYLEY FLOW",
        f"vertices/edges/SCCs={cayley['vertices']}/{cayley['labelled_edges']}/{cayley['sccs']}",
        f"SCC sizes={cayley['scc_size_histogram']} nontrivial={cayley['nontrivial_scc_profiles']}",
        f"terminal rank-one sinks={cayley['terminal_sinks']}",
        f"reachable-sink count histogram={cayley['reachable_sink_count_histogram']}",
        f"loops={cayley['loops_by_letter']} rank drops={cayley['rank_drop_histogram']}",
        "",
        "FINITE 58-ROW Q BANK",
        f"output masks/with complements={bank['distinct_output_masks']}/{bank['distinct_output_masks_with_complements']}",
        f"Q cell count histogram={bank['Q_cell_count_histogram']}",
        f"old-pair Q descent histogram={bank['old_pair_Q_descent_histogram']}",
        f"Q-injective words={bank['Q_injective_shortest_words']}",
        f"Q-injective + full-descent words={bank['Q_injective_and_full_old_pair_descent_shortest_words']}",
        f"defect rank/sector images={bank['defect_image_rank_histogram']}/{bank['occupied_sector_image_count_histogram']}",
        f"two-step Q={bank['two_step_Q_census']}",
        f"Q operation-kernel splits={bank['Q_operation_kernel_split_histograms']}",
        f"Q future refinement rounds={future['round_class_counts']}",
        f"Q future classes/profile={future['final_classes']}/{future['final_class_size_histogram']}",
        "",
        f"TOURNAMENT ANALYSIS vertices={len(ta['vertices'])} edge-flips={ta['edge_flips']}",
        f"  retention={ta['retention']}",
        "  both gauges transitive: C3=0, singleton SCCs, one Hamiltonian path",
        "PRESERVATION: literal return maps/reflection/complement/Q equality/defect images",
        "DEFERS: 8,940-mask node classification, every LRC metric field, increasing-size CF phase",
        "CHALLENGED ASSUMPTION: one iterated seam bit is not a closed recursive state",
        "ALL ASSERTIONS PASSED",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--witnesses", type=Path,
        default=Path("05-knowledge/results/mobius_cech_n9_exact_join_witnesses_codex_S13.tsv"),
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("05-knowledge/results/n9_closed_cf_return_semigroup_codex_S13g.out"),
    )
    parser.add_argument(
        "--json", type=Path,
        default=Path("05-knowledge/results/n9_closed_cf_return_semigroup_codex_S13g.json"),
    )
    args = parser.parse_args()
    result = run(args.witnesses)
    text = render(result)
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(text, end="")


if __name__ == "__main__":
    main()
