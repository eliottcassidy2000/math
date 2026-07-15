#!/usr/bin/env python3
"""Minimize bounded recursive observations on merged-metagraph line instances.

THM-796 gives deterministic maps on literal complement lines: delete either
Hamiltonian-path endpoint or reflect the staircase.  This script computes the
coarsest stable refinement of several current-state observations under those
operations over all line instances at 3 <= n <= 7.

Two gauges are kept separate.  The labelled-face gauge distinguishes top from
bottom deletion.  The unordered-face gauge retains only the multiset of face
successors, matching the symmetry that exchanges the two faces.  Both retain
reflection as a unary operation and retain the carrier's current observation;
the second gauge is therefore not a quotient of current states by reflection.
These are different continuation predicates, not alternative encodings of one
answer.

The residual audit then lifts line equality back to the marked-tournament
isomorphism groupoid.  It records chord-dual witnesses, endpoint stabilizers,
face restrictions, core ancestry, and the THM-801 gap face that resolves the
remaining two-face phase square.

Tournament Analysis uses information carriers as vertices.  The pairwise
observable is the sign of their difference in the number of n=7 sample-line
pairs separated after stable labelled-face refinement.  The switches are total
separation and separation per cell; the declared carrier order breaks ties.
Runners, arcs, gaps, wall events, and proof obligations were considered as
alternate vertex sets, but carrier partitions are required by the predicate
tested here.  The finite line state does not encode metric phase, owners,
reduced token holonomy, or LRC component translation.
"""

from __future__ import annotations

import argparse
import itertools
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from merged_metagraph_recursive_three_sort_audit_codex_S2 import (
    DEFAULT_N7,
    DEFAULT_SMALL,
    endpoint_delete,
    is_grid_symmetric,
    line_sheet_status,
    load_atlases,
    permute_mask,
    reflection_defect,
    reflection_permutation,
    tile_schema,
)
from mobius_cech_metagraph_codec_codex_S12 import (
    b3_face_mask,
    b3_line_face,
    b3_signature,
    boundary_c3_curvature,
    endpoint_epsilon,
)
from tournament_tiling_metagraph_address_codex_S4 import (
    permutation_arc_maps,
    perms,
    tile_schema as address_tile_schema,
    tiling_tournament,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge/results/three_sorted_metagraph_continuation_minimization_codex_S11.out"
)

State = tuple[int, int]


def full_mask(n: int) -> int:
    return (1 << len(tile_schema(n))) - 1


def line_id(mask: int, n: int) -> int:
    return min(mask, mask ^ full_mask(n))


def apex_zero(mask: int, n: int) -> int:
    other = mask ^ full_mask(n)
    apex = tile_schema(n).index((n, 1))
    return mask if ((mask >> apex) & 1) == 0 else other


def compact_labels(values: list[object]) -> tuple[int, ...]:
    labels: dict[object, int] = {}
    answer = []
    for value in values:
        answer.append(labels.setdefault(value, len(labels)))
    return tuple(answer)


def partition_stats(states: list[State], labels: tuple[int, ...], n: int) -> tuple[int, int, int, int]:
    cells: Counter[int] = Counter(
        label for state, label in zip(states, labels, strict=True) if state[0] == n
    )
    collisions = sum(size - 1 for size in cells.values())
    collision_pairs = sum(size * (size - 1) // 2 for size in cells.values())
    return len(cells), max(cells.values()), collisions, collision_pairs


def scc_sizes(adjacency: list[list[bool]]) -> list[int]:
    n = len(adjacency)
    unseen = set(range(n))
    sizes = []
    while unseen:
        seed = min(unseen)
        forward = {seed}
        backward = {seed}
        changed = True
        while changed:
            changed = False
            for u in range(n):
                if u not in forward and any(adjacency[v][u] for v in forward):
                    forward.add(u)
                    changed = True
            for u in range(n):
                if u not in backward and any(adjacency[u][v] for v in backward):
                    backward.add(u)
                    changed = True
        component = forward & backward
        unseen -= component
        sizes.append(len(component))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adjacency: list[list[bool]]) -> int:
    n = len(adjacency)
    dp: dict[tuple[int, int], int] = {}
    for vertex in range(n):
        dp[(1 << vertex, vertex)] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1 or not adjacency[last][nxt]:
                    continue
                key = (mask | (1 << nxt), nxt)
                dp[key] = dp.get(key, 0) + count
    return sum(dp.get(((1 << n) - 1, last), 0) for last in range(n))


def transitive_hamiltonian_path(adjacency: list[list[bool]]) -> tuple[int, ...] | None:
    scores = [sum(row) for row in adjacency]
    if sorted(scores) != list(range(len(adjacency))):
        return None
    path = tuple(sorted(range(len(adjacency)), key=lambda vertex: scores[vertex], reverse=True))
    assert all(adjacency[left][right] for left, right in itertools.pairwise(path))
    return path


def collision_groups(
    states: list[State], labels: tuple[int, ...], n: int
) -> list[list[State]]:
    groups: dict[int, list[State]] = defaultdict(list)
    for state, label in zip(states, labels, strict=True):
        if state[0] == n:
            groups[label].append(state)
    return sorted(
        (group for group in groups.values() if len(group) > 1),
        key=lambda group: (len(group), group),
    )


def same_partition(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    left_to_right: dict[int, int] = {}
    right_to_left: dict[int, int] = {}
    for left_label, right_label in zip(left, right, strict=True):
        if left_label in left_to_right and left_to_right[left_label] != right_label:
            return False
        if right_label in right_to_left and right_to_left[right_label] != left_label:
            return False
        left_to_right[left_label] = right_label
        right_to_left[right_label] = left_label
    return True


def tournament(values: list[Fraction], tie_order: list[int]) -> list[list[bool]]:
    n = len(values)
    tie_rank = {vertex: rank for rank, vertex in enumerate(tie_order)}
    adjacency = [[False] * n for _ in range(n)]
    for u in range(n):
        for v in range(u + 1, n):
            u_wins = values[u] > values[v] or (
                values[u] == values[v] and tie_rank[u] < tie_rank[v]
            )
            adjacency[u][v] = u_wins
            adjacency[v][u] = not u_wins
    return adjacency


def fingerprint(adjacency: list[list[bool]]) -> str:
    n = len(adjacency)
    scores = Counter(sum(row) for row in adjacency)
    triangles = 0
    for a, b, c in itertools.combinations(range(n), 3):
        if adjacency[a][b] and adjacency[b][c] and adjacency[c][a]:
            triangles += 1
        if adjacency[a][c] and adjacency[c][b] and adjacency[b][a]:
            triangles += 1
    return (
        f"score_hist={dict(sorted(scores.items()))} directed_3cycles={triangles} "
        f"SCCs={scc_sizes(adjacency)} HP={hamiltonian_path_count(adjacency)}"
    )


def tournament_image(mask: int, permutation_map: tuple[tuple[int, bool], ...]) -> int:
    image = 0
    for bit, (old_bit, reversed_arc) in enumerate(permutation_map):
        if ((mask >> old_bit) & 1) ^ reversed_arc:
            image |= 1 << bit
    return image


def analyze(atlases: dict[int, dict]) -> str:
    states = [
        (n, mask)
        for n in range(3, 8)
        for mask in range(1 << len(tile_schema(n)))
        if mask < (mask ^ full_mask(n))
    ]
    state_index = {state: index for index, state in enumerate(states)}
    assert len(states) == 1 + 4 + 32 + 512 + 16384

    transitions: dict[str, list[int | None]] = {name: [] for name in ("top", "bottom", "reflect")}
    records: dict[State, dict[str, object]] = {}
    for n, mask in states:
        other = mask ^ full_mask(n)
        q = list(map(int, atlases[n]["q"]))
        classes = list(map(int, atlases[n]["class_codes"]))
        u, v = sorted((q[mask], q[other]))
        colour = "B" if is_grid_symmetric(mask, n) else "K"
        sheet = line_sheet_status(q, classes, mask, other)
        defect = reflection_defect(mask, n)
        reflected = line_id(permute_mask(mask, reflection_permutation(n)), n)
        transitions["reflect"].append(state_index[(n, reflected)])

        if n == 3:
            transitions["top"].append(None)
            transitions["bottom"].append(None)
            xi = (n, u, v, colour, sheet)
            phase = (n, mask)
        else:
            top_mask = endpoint_delete(mask, n, "top")
            bottom_mask = endpoint_delete(mask, n, "bottom")
            top_line = line_id(top_mask, n - 1)
            bottom_line = line_id(bottom_mask, n - 1)
            transitions["top"].append(state_index[(n - 1, top_line)])
            transitions["bottom"].append(state_index[(n - 1, bottom_line)])

            oriented = apex_zero(mask, n)
            oriented_other = oriented ^ full_mask(n)
            top = endpoint_delete(oriented, n, "top")
            top_other = endpoint_delete(oriented_other, n, "top")
            bottom = endpoint_delete(oriented, n, "bottom")
            bottom_other = endpoint_delete(oriented_other, n, "bottom")
            low_q = list(map(int, atlases[n - 1]["q"]))
            top_colour = "B" if is_grid_symmetric(top, n - 1) else "K"
            bottom_colour = "B" if is_grid_symmetric(bottom, n - 1) else "K"
            xi = (
                n,
                q[oriented],
                q[oriented_other],
                low_q[top],
                low_q[top_other],
                low_q[bottom],
                low_q[bottom_other],
                colour + top_colour + bottom_colour,
            )
            if n == 4:
                phase = (n, mask)
            else:
                core = endpoint_delete(top, n - 1, "bottom")
                phase = (n, line_id(top, n - 1), line_id(bottom, n - 1), core & 1)

        records[(n, mask)] = {
            "plain": (n, u, v),
            "colour": (n, u, v, colour),
            "sheet": (n, u, v, colour, sheet),
            "defect": (n, u, v, colour, sheet, defect),
            "xi": xi,
            "xi_defect": (xi, sheet, defect),
            "phase": phase,
            "exact": (n, mask),
        }

    assert all(value is not None for value in transitions["reflect"])
    assert len({records[state]["phase"] for state in states}) == len(states)

    carriers = ["plain", "colour", "sheet", "defect", "xi", "xi_defect", "phase", "exact"]
    results: dict[tuple[str, str], tuple[tuple[int, ...], int]] = {}
    for carrier in carriers:
        initial = compact_labels([records[state][carrier] for state in states])
        for gauge in ("labelled", "unordered"):
            labels = initial
            rounds = 0
            while True:
                keys = []
                for index, _state in enumerate(states):
                    top = transitions["top"][index]
                    bottom = transitions["bottom"][index]
                    reflected = transitions["reflect"][index]
                    top_label = -1 if top is None else labels[top]
                    bottom_label = -1 if bottom is None else labels[bottom]
                    if gauge == "labelled":
                        future = (top_label, bottom_label, labels[reflected])
                    else:
                        future = (tuple(sorted((top_label, bottom_label))), labels[reflected])
                    keys.append((labels[index], future))
                refined = compact_labels(keys)
                if refined == labels:
                    break
                labels = refined
                rounds += 1
                assert rounds <= 10
            results[(carrier, gauge)] = (labels, rounds)

    for carrier in ("colour", "sheet", "defect"):
        assert same_partition(
            results[("plain", "labelled")][0], results[(carrier, "labelled")][0]
        )
    assert same_partition(
        results[("plain", "unordered")][0], results[("colour", "unordered")][0]
    )
    assert same_partition(
        results[("plain", "unordered")][0], results[("sheet", "unordered")][0]
    )
    for gauge in ("labelled", "unordered"):
        assert same_partition(results[("xi", gauge)][0], results[("xi_defect", gauge)][0])
        assert same_partition(results[("phase", gauge)][0], results[("exact", gauge)][0])

    lines = [
        "THM-796 BOUNDED CONTINUATION MINIMIZATION",
        "states are literal complement lines at n=3..7; operations are top/bottom deletion and reflection",
        "labelled gauge distinguishes faces; unordered-successor gauge retains their cell multiset",
        "carriers: plain=node boundary; colour=plain+blue/black; sheet=colour+loop sheet;",
        "defect=sheet+exact reflection defect; xi=joint upper/two-face node-colour cell;",
        "xi_defect=xi+sheet+defect; phase=two lower lines+coherent bit for n>=5 (exact fallback below);",
        "exact=literal line; unordered mode unorders successor cells, not ordered fields already in O",
        "",
        "N=7 CARRIER REFINEMENT",
        "carrier initial[cells,max,excess-lines,pairs] labelled[rounds,cells,max,excess-lines,pairs] "
        "unordered[rounds,cells,max,excess-lines,pairs]",
    ]
    for carrier in carriers:
        initial = compact_labels([records[state][carrier] for state in states])
        initial_stats = partition_stats(states, initial, 7)
        labelled, labelled_rounds = results[(carrier, "labelled")]
        unordered, unordered_rounds = results[(carrier, "unordered")]
        labelled_stats = partition_stats(states, labelled, 7)
        unordered_stats = partition_stats(states, unordered, 7)
        lines.append(
            f"{carrier} {list(initial_stats)} "
            f"[{labelled_rounds},{','.join(map(str, labelled_stats))}] "
            f"[{unordered_rounds},{','.join(map(str, unordered_stats))}]"
        )

    lines.extend(
        [
            "",
            "BASE OBSERVATION BY SIZE (initial -> stable labelled / stable unordered-successor cells)",
            "n plain colour sheet defect xi xi_defect phase exact",
        ]
    )
    for n in range(3, 8):
        entries = []
        for carrier in carriers:
            initial = compact_labels([records[state][carrier] for state in states])
            initial_cells = partition_stats(states, initial, n)[0]
            labelled_cells = partition_stats(states, results[(carrier, "labelled")][0], n)[0]
            unordered_cells = partition_stats(states, results[(carrier, "unordered")][0], n)[0]
            entries.append(f"{initial_cells}->{labelled_cells}/{unordered_cells}")
        lines.append(f"{n} " + " ".join(entries))

    plain_groups = collision_groups(states, results[("plain", "labelled")][0], 7)
    plain_histogram = Counter(map(len, plain_groups))
    lines.extend(
        [
            "",
            "N=7 COLLISION RESIDUA",
            f"stable labelled plain collision_cell_histogram={dict(sorted(plain_histogram.items()))}",
            "partition equalities: labelled plain=colour=sheet=defect; unordered plain=colour=sheet;",
            "both gauges xi=xi_defect and phase=exact",
        ]
    )
    for gauge in ("labelled", "unordered"):
        groups = collision_groups(states, results[("xi", gauge)][0], 7)
        lines.append(f"xi stable {gauge} collision_cells={len(groups)}")
        for group_index, group in enumerate(groups, 1):
            entries = []
            for state in group:
                record = records[state]
                _n, mask = state
                entries.append(
                    f"mask=0x{mask:04x} phase={record['phase']} "
                    f"defect=0x{int(record['defect'][-1]):04x}"
                )
            lines.append(f"  cell {group_index}: " + " ; ".join(entries))

    labelled_xi_groups = collision_groups(states, results[("xi", "labelled")][0], 7)
    residual_masks = {mask for group in labelled_xi_groups for _n, mask in group}
    assert residual_masks == {0x12CA, 0x12CB, 0x146C, 0x146D}
    rho = {
        mask: line_id(apex_zero(mask, 7) ^ (1 << tile_schema(7).index((7, 1))), 7)
        for mask in residual_masks
    }
    reflect = {
        mask: line_id(permute_mask(mask, reflection_permutation(7)), 7)
        for mask in residual_masks
    }
    assert rho == {0x12CA: 0x12CB, 0x12CB: 0x12CA, 0x146C: 0x146D, 0x146D: 0x146C}
    assert reflect == {
        0x12CA: 0x146C,
        0x12CB: 0x146D,
        0x146C: 0x12CA,
        0x146D: 0x12CB,
    }
    assert all(rho[reflect[mask]] == reflect[rho[mask]] for mask in residual_masks)
    square_records = [records[(7, mask)] for mask in sorted(residual_masks)]
    assert {record["plain"] for record in square_records} == {(7, 264, 270)}
    assert {record["colour"][-1] for record in square_records} == {"K"}
    assert {record["sheet"][-1] for record in square_records} == {"E"}
    assert {record["defect"][-1] for record in square_records} == {0x06A6}
    assert {apex_zero(mask, 7) for mask in residual_masks} == {
        0x12CA,
        0x6D34,
        0x146C,
        0x6B92,
    }
    curvature_data = {}
    for mask in residual_masks:
        oriented = apex_zero(mask, 7)
        other = oriented ^ full_mask(7)
        curvature_data[mask] = (
            boundary_c3_curvature(oriented, 7),
            boundary_c3_curvature(other, 7),
            endpoint_epsilon(oriented, 7),
            b3_signature(oriented, 7),
        )
    phase_one_b3 = (0, 0, 0, 2, 0, 2, 2)
    phase_zero_b3 = (1, 0, 1, 1, 3, 1, 1)
    assert curvature_data == {
        0x12CA: (1, 2, 0, phase_one_b3),
        0x12CB: (1, 2, 0, phase_zero_b3),
        0x146C: (1, 2, 0, phase_one_b3),
        0x146D: (1, 2, 0, phase_zero_b3),
    }
    assert all(
        curvature_data[mask][3]
        == (phase_one_b3 if records[(7, mask)]["phase"][-1] else phase_zero_b3)
        for mask in residual_masks
    )
    q7 = list(map(int, atlases[7]["q"]))
    classes7 = list(map(int, atlases[7]["class_codes"]))
    class_fibres = Counter(
        classes7[mask] for mask, node in enumerate(q7) if node in {264, 270}
    )
    assert class_fibres == {
        0x6544A611A2974: 151,
        0x65846512A2974: 151,
        0x7436470291C1A: 57,
    }
    class_pairs = {
        mask: (
            classes7[apex_zero(mask, 7)],
            classes7[apex_zero(mask, 7) ^ full_mask(7)],
        )
        for mask in residual_masks
    }
    assert class_pairs[0x12CA] == class_pairs[0x12CB] == (
        0x6544A611A2974,
        0x7436470291C1A,
    )
    assert class_pairs[0x146C] == class_pairs[0x146D] == (
        0x65846512A2974,
        0x7436470291C1A,
    )
    address_tiles, _sigma = address_tile_schema(7)
    phase_presentations = (
        (0x12CA, apex_zero(0x12CB, 7)),
        (0x146C, apex_zero(0x146D, 7)),
    )
    phase_isomorphisms = []
    opposite_isomorphisms = []
    for left, right in phase_presentations:
        left_tournament = tiling_tournament(left, 7, address_tiles)
        right_tournament = tiling_tournament(right, 7, address_tiles)
        witnesses = tuple(
            permutation
            for permutation, permutation_map in zip(perms(7), permutation_arc_maps(7), strict=True)
            if tournament_image(left_tournament, permutation_map) == right_tournament
        )
        phase_isomorphisms.append(witnesses)
        left_opposite = tiling_tournament(left ^ full_mask(7), 7, address_tiles)
        right_opposite = tiling_tournament(right ^ full_mask(7), 7, address_tiles)
        opposite_isomorphisms.append(
            tuple(
                permutation
                for permutation, permutation_map in zip(
                    perms(7), permutation_arc_maps(7), strict=True
                )
                if tournament_image(left_opposite, permutation_map) == right_opposite
            )
        )
    assert phase_isomorphisms == [
        ((0, 2, 3, 5, 1, 4, 6),),
        ((0, 2, 5, 1, 3, 4, 6),),
    ]
    assert opposite_isomorphisms == [
        (
            (0, 4, 1, 2, 5, 3, 6),
            (2, 0, 3, 4, 1, 5, 6),
            (4, 2, 5, 0, 3, 1, 6),
        ),
        (
            (0, 1, 5, 2, 3, 6, 4),
            (0, 3, 1, 4, 5, 2, 6),
            (0, 5, 3, 6, 1, 4, 2),
        ),
    ]
    for rigid, opposite in zip(phase_isomorphisms, opposite_isomorphisms, strict=True):
        endpoint_fixed = tuple(p for p in opposite if p[0] == 0 and p[-1] == 6)
        inverse_rigid = tuple(rigid[0].index(vertex) for vertex in range(7))
        assert endpoint_fixed == (inverse_rigid,)
        assert set(rigid).isdisjoint(opposite)
    path_reflection = tuple(reversed(range(7)))
    first_isomorphism = phase_isomorphisms[0][0]
    reflected_isomorphism = tuple(
        path_reflection[first_isomorphism[path_reflection[index]]] for index in range(7)
    )
    assert reflected_isomorphism == phase_isomorphisms[1][0]
    node_records = {int(node["rank"]): node for node in atlases[7]["record"]["nodes"]}
    assert not node_records[264]["self_converse"] and node_records[270]["self_converse"]
    q6 = list(map(int, atlases[6]["q"]))
    for record in square_records:
        assert all(q6[line] == q6[line ^ full_mask(6)] for line in record["phase"][1:3])
    core_lines = {
        mask: line_id(
            endpoint_delete(endpoint_delete(apex_zero(mask, 7), 7, "top"), 6, "bottom"),
            5,
        )
        for mask in residual_masks
    }
    assert core_lines == {
        0x12CA: 0x3,
        0x12CB: 0x3,
        0x146C: 0x9,
        0x146D: 0x9,
    }
    q5 = list(map(int, atlases[5]["q"]))
    classes5 = list(map(int, atlases[5]["class_codes"]))
    for core_line in {0x3, 0x9}:
        core_other = core_line ^ full_mask(5)
        assert q5[core_line] == q5[core_other] == 6
        assert classes5[core_line] == classes5[core_other] == 0xC
        assert not is_grid_symmetric(core_line, 5)
    gap_data = {}
    for mask in residual_masks:
        gap_endpoint = b3_face_mask(apex_zero(mask, 7), 7, "B")
        gap_other = gap_endpoint ^ full_mask(6)
        gap_data[mask] = (
            b3_line_face(mask, 7, "B"),
            (q6[gap_endpoint], q6[gap_other]),
            "B" if is_grid_symmetric(gap_endpoint, 6) else "K",
        )
    assert gap_data == {
        0x12CA: (0x115, (33, 33), "K"),
        0x12CB: (0x114, (21, 33), "K"),
        0x146C: (0x0C3, (33, 33), "K"),
        0x146D: (0x0C2, (21, 33), "K"),
    }
    assert all(
        int(gap_data[mask][1][0] == gap_data[mask][1][1])
        == int(records[(7, mask)]["phase"][-1])
        for mask in residual_masks
    )
    lines.extend(
        [
            "residual labelled Xi square:",
            "  rho=(12ca 12cb)(146c 146d), reflection=(12ca 146c)(12cb 146d); involutions commute",
            "  all four are black cross-lines on nodes (264 non-SC,270 SC), defect=0x06a6",
            "  THM-811: all four have oriented (q0,q1,epsilon_Smith)=(1,2,0),",
            "  B3 signatures split coherent phase exactly: phase 1=(0,0,0,2,0,2,2), phase 0=(1,0,1,1,3,1,1)",
            "  hence node pair+(q0,q1,|epsilon_Smith|) collapses the square; B3 gap/ownership data do not",
            "  endpoint class-sheet fibre sizes at those nodes are 151,151 and 57",
            "  rho changes the marked-path presentation inside fixed ordinary endpoint classes;",
            "  reflection exchanges the two 151-element converse sheets at node 264",
            "  rigid-endpoint rho isomorphisms are unique, endpoint-fixed, reflection-conjugate 5-cycles",
            f"  interior permutations={tuple(witnesses[0] for witnesses in phase_isomorphisms)}",
            "  SC endpoint has three isomorphisms; its endpoint-fixed choice is the inverse 5-cycle,",
            "  so no one relabelling maps both line endpoints (relative holonomy has order 5)",
            "  ancestry cores are reflected n=5 black same-class loops 0x3 and 0x9 at node 6",
            "  both n=6 endpoint faces are loops; the n=7 residuals are cross-lines,",
            "  so this branch is born at n=7 and cannot directly face-lift to an n=8 Xi collision",
            "  THM-801 gap lines=(0x115,0x114,0xc3,0xc2): phase 1 is a black loop at node 33,",
            "  phase 0 is a black cross-line 21--33; Omega resolves the whole square",
        ]
    )

    n7_count = sum(state[0] == 7 for state in states)
    separated: list[Fraction] = []
    economy: list[Fraction] = []
    for carrier in carriers:
        labels = results[(carrier, "labelled")][0]
        cells, _, _, collision_pairs = partition_stats(states, labels, 7)
        separated_pairs = n7_count * (n7_count - 1) // 2 - collision_pairs
        separated.append(Fraction(separated_pairs))
        economy.append(Fraction(separated_pairs, cells))
    tie = list(range(len(carriers)))
    retention_tournament = tournament(separated, tie)
    economy_tournament = tournament(economy, tie)
    retention_path = transitive_hamiltonian_path(retention_tournament)
    economy_path = transitive_hamiltonian_path(economy_tournament)
    flips = sum(
        retention_tournament[u][v] != economy_tournament[u][v]
        for u in range(len(carriers))
        for v in range(u + 1, len(carriers))
    )
    lines.extend(
        [
            "",
            "TOURNAMENT ANALYSIS (stable labelled n=7 partitions)",
            f"carriers={tuple(carriers)}",
            f"separated_pairs={tuple(int(value) for value in separated)}",
            f"retention {fingerprint(retention_tournament)}",
            f"economy {fingerprint(economy_tournament)}",
            f"edge_flips={flips}",
            f"tie_order={tuple(carriers)}",
            f"retention_Hamiltonian_path={tuple(carriers[i] for i in retention_path or ())}",
            f"economy_Hamiltonian_path={tuple(carriers[i] for i in economy_path or ())}",
            "",
            "PRESERVATION VERDICT",
            "for each carrier, stable labelled = its coarsest stable refinement on bounded states n=3..7",
            "for each carrier, stable unordered = coarsest branching-multiset refinement on those states",
            "unordered refinement retains O and unary reflection; it is not a state/reflection quotient",
            "phase and exact induce singleton partitions; n=3,4 use an exact-mask phase fallback",
            "the tested combinatorial line state does not encode metric phase, owners, h_red, or Delta_M",
            "",
            "ALL ASSERTIONS PASSED",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=DEFAULT_SMALL)
    parser.add_argument("--n7-atlas", type=Path, default=DEFAULT_N7)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    report = analyze(load_atlases(args.small_atlas, args.n7_atlas))
    args.output.write_text(report)
    if args.stdout:
        print(report, end="")


if __name__ == "__main__":
    main()
