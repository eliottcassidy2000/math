#!/usr/bin/env python3
"""Exact closed CF lift/face operation kernels for THM-846.

Compose THM-838's centered coordinate copy Phi:X9->X10 with each THM-842
face R_10:X10->X9.  Literal coordinate kernels are proved before any quotient
is minimized.  The three composites are then audited on THM-828's 58 cells
against Q, bare merged-node marginals, coupled bar-P, the apex-relative theta
sheet, and the unique B seam.

Tournament Analysis uses information carriers rather than runners as
vertices.  Its pairwise observable is separated witness pairs, its switch is
raw retention versus retention per logarithmic cost, and the displayed
carrier order is the tie Hamiltonian path.

Preserved: literal coordinate maps, complement/reflection action, the finite
defect bank, Q, locally exact tournament-isomorphism partitions, seam, and
theta.  Destroyed: LRC metric gaps, owners, walls, loneliness, and arbitrary
continued-fraction words.  Challenged assumption: one sidecar is not expected
to repair literal inversion, affine-sheet order, and node projection at once.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from continued_fraction_n9_defect_transport_codex_S13 import (
    RHO_4_5,
    RHO_6_7,
    merged_node_partition,
    phase_schedule,
    recursive_rho,
    rho_5_6,
)
from n9_defect_b3_continuation_purity_codex_S13d import (
    complement,
    face,
    fibre_hist,
    gf2_rank,
    q_cell,
    reflect,
    ta_fingerprint,
    theta,
    tile_index,
    tiles,
)


DEFECT_BASIS = (0x0192486, 0x08C2C0C, 0x11B4600, 0x4483414)
ROLES = ("A", "B", "C")


def build_rho() -> tuple[int, ...]:
    schedules = phase_schedule()
    rhos: dict[int, tuple[int, ...]] = {5: rho_5_6(), 6: RHO_6_7}
    rhos[7] = recursive_rho(7, rhos[5], tuple(schedules[7]["increments"]))
    rhos[8] = recursive_rho(8, rhos[6], tuple(schedules[8]["increments"]))
    rhos[9] = recursive_rho(9, rhos[7], tuple(schedules[9]["increments"]))
    rho = rhos[9]
    assert rho == (
        0, 1, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12, 12, 13, 14, 15,
        16, 17, 17, 18, 19, 20, 15, 21, 22, 23, 24, 21, 25, 26, 24, 27, 26, 27,
    )
    return rho


def role_coordinate_maps(rho: tuple[int, ...]) -> dict[str, tuple[int, ...]]:
    upper = tile_index(10)
    maps = {}
    for role in ROLES:
        source_bits = []
        for a, b in tiles(9):
            source_tile = {"A": (a, b), "B": (a + 1, b), "C": (a + 1, b + 1)}[role]
            source_bits.append(rho[upper[source_tile]])
        maps[role] = tuple(source_bits)
    return maps


def apply_coordinate_map(mask: int, mapping: tuple[int, ...]) -> int:
    return sum(((mask >> source) & 1) << target for target, source in enumerate(mapping))


def seam(mask: int) -> int:
    index = tile_index(9)
    return ((mask >> index[(6, 4)]) ^ (mask >> index[(7, 3)])) & 1


def inverse_b(image: int, lost_seam: int) -> int:
    index = tile_index(9)
    bit64, bit73 = index[(6, 4)], index[(7, 3)]
    recovered = image & ~(1 << bit64)
    recovered |= (((image >> bit73) & 1) ^ lost_seam) << bit64
    return recovered


def read_witnesses(path: Path) -> list[dict[str, int]]:
    rows = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({key: int(row[key], 16) for key in ("D", "u", "v")})
    assert len(rows) == 58
    return rows


def partition(values: list[object]) -> tuple[tuple[int, ...], ...]:
    groups: dict[object, list[int]] = defaultdict(list)
    for index, value in enumerate(values):
        groups[value].append(index)
    return tuple(sorted(tuple(group) for group in groups.values()))


def kernel_audit(source: list[object], target: list[object]) -> dict[str, object]:
    groups: dict[object, list[int]] = defaultdict(list)
    for index, value in enumerate(source):
        groups[value].append(index)
    split_groups = []
    for indices in groups.values():
        targets = {target[index] for index in indices}
        if len(targets) > 1:
            split_groups.append({"indices": indices, "target_cells": len(targets)})
    return {
        "source_fibres": fibre_hist(source),
        "target_fibres": fibre_hist(target),
        "split_source_cells": len(split_groups),
        "split_rows": split_groups,
        "kernel_contained": not split_groups,
    }


def tile_names(bits: set[int]) -> list[str]:
    schema = tiles(9)
    return [f"({schema[bit][0]},{schema[bit][1]})" for bit in sorted(bits)]


def run(witness_path: Path) -> dict[str, object]:
    rows = read_witnesses(witness_path)
    rho = build_rho()
    maps = role_coordinate_maps(rho)
    index9 = tile_index(9)
    full9 = (1 << len(tiles(9))) - 1

    coordinate_rows = {}
    all_source_bits = set(range(len(tiles(9))))
    for role in ROLES:
        mapping = maps[role]
        used = set(mapping)
        multiplicity = Counter(mapping)
        coordinate_rows[role] = {
            "rank": len(used),
            "missing_source_tiles": tile_names(all_source_bits - used),
            "duplicated_source_tiles": {
                f"({tiles(9)[bit][0]},{tiles(9)[bit][1]})": count
                for bit, count in sorted(multiplicity.items()) if count > 1
            },
            "output_to_source_tiles": [
                {
                    "output": list(tiles(9)[target]),
                    "source": list(tiles(9)[source]),
                }
                for target, source in enumerate(mapping)
            ],
        }

    expected_a_missing = {index9[(9, b)] for b in range(1, 8)}
    expected_c_missing = {index9[(a, 1)] for a in range(3, 10)}
    assert set(maps["A"]) == all_source_bits - expected_a_missing
    assert set(maps["C"]) == all_source_bits - expected_c_missing
    apex_bit = index9[(9, 1)]
    assert set(maps["A"]) | set(maps["C"]) == all_source_bits - {apex_bit}
    assert {role: coordinate_rows[role]["rank"] for role in ROLES} == {"A": 21, "B": 27, "C": 21}

    # Prove the linear parts of the full-cube symmetry identities directly
    # as coordinate permutations, rather than only sampling masks below.
    sigma = tuple((reflect(1 << bit, 9)).bit_length() - 1 for bit in range(28))
    assert sorted(sigma) == list(range(28))
    assert all(sigma[maps["A"][target]] == maps["C"][sigma[target]] for target in range(28))
    assert all(sigma[maps["C"][target]] == maps["A"][sigma[target]] for target in range(28))
    assert all(sigma[maps["B"][target]] == maps["B"][sigma[target]] for target in range(28))

    exceptional = {target: source for target, source in enumerate(maps["B"]) if target != source}
    assert exceptional == {index9[(6, 4)]: index9[(7, 3)]}
    assert tuple(maps["B"][maps["B"][i]] for i in range(28)) == maps["B"]
    # Reflection commutation plus equality on the affine constant proves the
    # full affine theta commutation identity for B.
    assert apply_coordinate_map(theta(0, 9), maps["B"]) == theta(0, 9)
    theta_constant_seam = seam(theta(0, 9))
    assert all(
        seam(theta(1 << bit, 9)) ^ theta_constant_seam == seam(1 << bit)
        for bit in range(28)
    )

    # An explicit coordinate left inverse for ordered (F_A,F_C), with the
    # common apex token.  On the canonical apex-zero line slice no extra bit
    # is needed.
    inverse_ac: dict[int, tuple[str, int]] = {}
    for role in ("A", "C"):
        for output, source in enumerate(maps[role]):
            inverse_ac.setdefault(source, (role, output))
    assert set(inverse_ac) == all_source_bits - {apex_bit}

    def reconstruct_ac(a_image: int, c_image: int, apex: int) -> int:
        answer = apex << apex_bit
        for source, (role, output) in inverse_ac.items():
            image = a_image if role == "A" else c_image
            answer |= ((image >> output) & 1) << source
        return answer

    # The pure Phi image has eight equality constraints, one per doubled rho fibre.
    rho_fibres: dict[int, list[int]] = defaultdict(list)
    for target, source in enumerate(rho):
        rho_fibres[source].append(target)
    doubled = {source: targets for source, targets in rho_fibres.items() if len(targets) > 1}
    expected_doubled_tiles = {(3, 1), (4, 2), (7, 3), (5, 3), (6, 4), (7, 5), (8, 6), (9, 7)}
    assert {tiles(9)[source] for source in doubled} == expected_doubled_tiles
    assert len(doubled) == 8 and all(len(targets) == 2 for targets in doubled.values())

    tests = {0, full9, *DEFECT_BASIS}
    for row in rows:
        tests.update((row["u"], row["v"], complement(row["u"], 9), complement(row["v"], 9)))
    for mask in tests:
        fa = apply_coordinate_map(mask, maps["A"])
        fb = apply_coordinate_map(mask, maps["B"])
        fc = apply_coordinate_map(mask, maps["C"])
        assert reconstruct_ac(fa, fc, (mask >> apex_bit) & 1) == mask
        assert inverse_b(fb, seam(mask)) == mask
        assert apply_coordinate_map(fb, maps["B"]) == fb
        assert seam(fb) == 0
        for role, image in (("A", fa), ("B", fb), ("C", fc)):
            assert apply_coordinate_map(complement(mask, 9), maps[role]) == complement(image, 9)
        assert apply_coordinate_map(reflect(mask, 9), maps["A"]) == reflect(fc, 9)
        assert apply_coordinate_map(reflect(mask, 9), maps["C"]) == reflect(fa, 9)
        assert apply_coordinate_map(reflect(mask, 9), maps["B"]) == reflect(fb, 9)
        assert apply_coordinate_map(theta(mask, 9), maps["B"]) == theta(fb, 9)

    f = {role: (lambda mask, mapping=maps[role]: apply_coordinate_map(mask, mapping)) for role in ROLES}
    sectors = sorted({row["D"] for row in rows})
    assert len(sectors) == 11
    assert all(f["B"](word) == word and seam(word) == 0 for word in (*DEFECT_BASIS, *sectors))
    defect_action = {
        role: {
            "basis_rank": gf2_rank([f[role](word) for word in DEFECT_BASIS]),
            "occupied_sector_images": len({f[role](word) for word in sectors}),
        }
        for role in ROLES
    }

    # Classify every source and composite endpoint in one exact local node atlas.
    masks: set[int] = set()
    for row in rows:
        for mask in (row["u"], row["v"]):
            for image in (mask, *(f[role](mask) for role in ROLES)):
                masks.update((image, complement(image, 9)))
    node, classifier = merged_node_partition(masks, 9)
    pbar = lambda mask: tuple(sorted((node[mask], node[complement(mask, 9)])))

    source_u = [row["u"] for row in rows]
    source_v = [row["v"] for row in rows]
    assert all(((mask >> apex_bit) & 1) == 0 for mask in (*source_u, *source_v))
    images_u = {role: [f[role](mask) for mask in source_u] for role in ROLES}
    images_v = {role: [f[role](mask) for mask in source_v] for role in ROLES}
    q_values = {role: [q_cell(mask, 9) for mask in images_u[role]] for role in ROLES}
    p_values = {role: [pbar(mask) for mask in images_u[role]] for role in ROLES}
    q_pair_descent = {
        role: sum(q_cell(a, 9) == q_cell(b, 9) for a, b in zip(images_u[role], images_v[role]))
        for role in ROLES
    }
    p_pair_descent = {
        role: sum(pbar(a) == pbar(b) for a, b in zip(images_u[role], images_v[role]))
        for role in ROLES
    }
    assert q_pair_descent == p_pair_descent == {"A": 0, "B": 58, "C": 0}
    assert all(fibre_hist(q_values[role]) == (58, {1: 58}) for role in ROLES)
    assert all(fibre_hist(p_values[role]) == (58, {1: 58}) for role in ROLES)

    ordered_ac_q = [(q_values["A"][i], q_values["C"][i]) for i in range(58)]
    ordered_ac_p = [(p_values["A"][i], p_values["C"][i]) for i in range(58)]
    ac_q = [tuple(sorted(value)) for value in ordered_ac_q]
    ac_p = [tuple(sorted(value)) for value in ordered_ac_p]
    ac_plus_b_p = [(ac_p[i], p_values["B"][i]) for i in range(58)]
    assert fibre_hist(ordered_ac_q) == fibre_hist(ordered_ac_p) == (58, {1: 58})
    assert fibre_hist(ac_q) == (29, {2: 29})
    assert fibre_hist(ac_p) == (29, {2: 29})
    assert fibre_hist(ac_plus_b_p) == (58, {1: 58})

    # The theta involution acts freely on the chosen 58 representatives.
    source_index = {mask: i for i, mask in enumerate(source_u)}
    theta_mate = [source_index[theta(mask, 9)] for mask in source_u]
    theta_orbits = [(i, mate) for i, mate in enumerate(theta_mate) if i < mate]
    assert len(theta_orbits) == 29 and all(theta_mate[mate] == i for i, mate in theta_orbits)
    assert all(ac_p[i] == ac_p[j] and ac_q[i] == ac_q[j] for i, j in theta_orbits)
    assert all(ordered_ac_q[i] != ordered_ac_q[j] and ordered_ac_p[i] != ordered_ac_p[j]
               for i, j in theta_orbits)
    assert all(seam(source_u[i]) == seam(source_u[j]) for i, j in theta_orbits)
    assert all(p_values["B"][i] != p_values["B"][j] for i, j in theta_orbits)
    theta_sheet_preserved = sum(
        (q_cell(source_u[i], 9) > q_cell(source_u[j], 9))
        == (q_values["B"][i] > q_values["B"][j])
        for i, j in theta_orbits
    )
    assert theta_sheet_preserved == 29
    orientation_triples = []
    for i, j in theta_orbits:
        source_sign = (q_cell(source_u[i], 9) > q_cell(source_u[j], 9)) - (
            q_cell(source_u[i], 9) < q_cell(source_u[j], 9)
        )
        b_sign = (q_values["B"][i] > q_values["B"][j]) - (
            q_values["B"][i] < q_values["B"][j]
        )
        ac_sign = (q_values["A"][i] > q_values["C"][i]) - (
            q_values["A"][i] < q_values["C"][i]
        )
        orientation_triples.append((source_sign, b_sign, ac_sign))
    assert all(source_sign == b_sign == -ac_sign != 0
               for source_sign, b_sign, ac_sign in orientation_triples)
    seam_orbits = Counter(seam(source_u[i]) for i, _ in theta_orbits)
    assert seam_orbits == {0: 13, 1: 16}

    source_direct = [node[mask] for mask in source_u]
    source_partner = [node[complement(mask, 9)] for mask in source_u]
    source_p = [pbar(mask) for mask in source_u]
    source_q = [q_cell(mask, 9) for mask in source_u]
    b_direct = [node[mask] for mask in images_u["B"]]
    b_partner = [node[complement(mask, 9)] for mask in images_u["B"]]
    b_p = p_values["B"]
    source_seam = [seam(mask) for mask in source_u]
    target_seam = [seam(mask) for mask in images_u["B"]]

    assert fibre_hist(source_direct) == (54, {1: 51, 2: 2, 3: 1})
    assert fibre_hist(source_partner) == (53, {1: 48, 2: 5})
    assert fibre_hist(b_direct) == fibre_hist(b_partner) == (55, {1: 52, 2: 3})
    assert fibre_hist(b_p) == (58, {1: 58})
    assert fibre_hist(source_p) == fibre_hist(source_q) == (58, {1: 58})
    assert not any(left == right for left, right in b_p)
    assert not any(target_seam)
    assert fibre_hist(list(zip(b_direct, source_seam))) == (55, {1: 52, 2: 3})
    assert fibre_hist(list(zip(b_partner, source_seam))) == (55, {1: 52, 2: 3})

    source_to_b_p = list(zip(source_p, b_p))
    assert fibre_hist(source_to_b_p) == (58, {1: 58})
    b_image_masks = {
        image
        for mask in (*source_u, *source_v)
        for image in (f["B"](mask), complement(f["B"](mask), 9))
    }
    assert len(b_image_masks) == 232
    assert len({node[mask] for mask in b_image_masks}) == 105

    operation_kernels = {
        "direct_node": kernel_audit(source_direct, b_direct),
        "partner_node": kernel_audit(source_partner, b_partner),
        "direct_node_plus_seam": kernel_audit(
            list(zip(source_direct, source_seam)), list(zip(b_direct, target_seam))
        ),
        "partner_node_plus_seam": kernel_audit(
            list(zip(source_partner, source_seam)), list(zip(b_partner, target_seam))
        ),
        "coupled_bar_P": kernel_audit(source_p, b_p),
        "literal_Q": kernel_audit(source_q, q_values["B"]),
    }
    assert operation_kernels["direct_node"]["split_source_cells"] == 1
    assert operation_kernels["partner_node"]["split_source_cells"] == 2
    assert operation_kernels["direct_node_plus_seam"]["split_source_cells"] == 0
    assert operation_kernels["partner_node_plus_seam"]["split_source_cells"] == 1
    assert operation_kernels["coupled_bar_P"]["split_source_cells"] == 0
    assert operation_kernels["literal_Q"]["split_source_cells"] == 0

    changed_u = sum(f["B"](mask) != mask for mask in source_u)
    assert changed_u == 32 and Counter(source_seam) == {0: 26, 1: 32}

    carriers = [
        ("defect sector", [row["D"] for row in rows]),
        ("B seam", source_seam),
        ("source direct node", source_direct),
        ("F_B direct node", b_direct),
        ("A/C bar-P deck", ac_p),
        ("F_B coupled bar-P", b_p),
        ("A/C literal-Q deck", ac_q),
    ]
    tournament_analysis = ta_fingerprint(carriers)

    return {
        "schema_version": 1,
        "theorem": "THM-846",
        "source_theorems": ["THM-828", "THM-834", "THM-838", "THM-840", "THM-842"],
        "rho_9_10": list(rho),
        "pure_Phi_image": {
            "dimension": 28,
            "target_dimension": 36,
            "codimension": 8,
            "doubled_coordinate_fibres": [
                {
                    "source_tile": list(tiles(9)[source]),
                    "target_tiles": [list(tiles(10)[target]) for target in targets],
                }
                for source, targets in sorted(doubled.items())
            ],
        },
        "closed_coordinate_maps": coordinate_rows,
        "coordinate_theorems": {
            "A_C_used_coordinate_union": 27,
            "A_C_common_missing_tile": [9, 1],
            "ordered_A_C_has_left_inverse_with_apex": True,
            "ordered_A_C_has_left_inverse_on_apex_zero_Q_slice": True,
            "B_idempotent": True,
            "B_identity_except": {"output": [6, 4], "source": [7, 3]},
            "B_lost_seam": "x_(6,4) xor x_(7,3)",
            "B_exact_inverse_with_seam": True,
            "B_commutes_with_complement_reflection_theta": True,
            "A_C_exchange_under_reflection": True,
        },
        "defect_action": {
            "roles": defect_action,
            "B_fixes_all_four_basis_words": True,
            "B_fixes_all_eleven_occupied_sectors": True,
            "B_seam_zero_on_basis_and_sectors": True,
        },
        "finite_bank": {
            "cells": 58,
            "B_changed_representatives": changed_u,
            "B_source_seam_histogram": dict(sorted(Counter(source_seam).items())),
            "rolewise_source_pair_Q_descent": q_pair_descent,
            "rolewise_source_pair_bar_P_descent": p_pair_descent,
            "role_representative_Q_fibres": {role: fibre_hist(q_values[role]) for role in ROLES},
            "role_representative_bar_P_fibres": {role: fibre_hist(p_values[role]) for role in ROLES},
            "ordered_A_C_Q_deck": fibre_hist(ordered_ac_q),
            "ordered_A_C_bar_P_deck": fibre_hist(ordered_ac_p),
            "unordered_A_C_Q_deck": fibre_hist(ac_q),
            "unordered_A_C_bar_P_deck": fibre_hist(ac_p),
            "A_C_bar_P_plus_B_bar_P": fibre_hist(ac_plus_b_p),
            "theta_orbits": len(theta_orbits),
            "theta_seam_orbit_histogram": dict(sorted(seam_orbits.items())),
            "theta_sheet_preserved_by_B": theta_sheet_preserved,
            "theta_source_B_orientation_equals_negative_A_C_role_order": len(
                orientation_triples
            ),
            "B_direct_node_fibres": fibre_hist(b_direct),
            "B_partner_node_fibres": fibre_hist(b_partner),
            "B_direct_node_plus_input_seam_fibres": fibre_hist(
                list(zip(b_direct, source_seam))
            ),
            "B_partner_node_plus_input_seam_fibres": fibre_hist(
                list(zip(b_partner, source_seam))
            ),
            "B_coupled_bar_P_fibres": fibre_hist(b_p),
            "B_coupled_bar_P_loops": sum(left == right for left, right in b_p),
            "B_image_masks": len(b_image_masks),
            "B_image_merged_nodes": len({node[mask] for mask in b_image_masks}),
            "source_to_B_bar_P_graph": fibre_hist(source_to_b_p),
            "source_Q_fibres": fibre_hist(source_q),
            "source_bar_P_fibres": fibre_hist(source_p),
        },
        "operation_kernels": operation_kernels,
        "classifier": {"masks": len(masks), **classifier},
        "tournament_analysis": tournament_analysis,
        "preservation": {
            "preserves": [
                "literal closed coordinate maps", "complement", "reflection", "theta under B",
                "finite defect bank", "Q", "local merged-node equality", "coupled bar-P", "seam",
            ],
            "destroys": [
                "LRC gaps", "owners", "walls", "metric loneliness", "arbitrary CF words",
            ],
            "challenged_assumption": (
                "seam, theta sheet, and coupled node fibre solve three distinct forgetting problems"
            ),
        },
    }


def render(result: dict[str, object]) -> str:
    coordinate = result["coordinate_theorems"]
    bank = result["finite_bank"]
    kernels = result["operation_kernels"]
    ta = result["tournament_analysis"]
    lines = [
        "THM-846 N=9 CLOSED CF LIFT/FACE OPERATION KERNELS",
        "=" * 70,
        f"coordinate ranks A/B/C={tuple(result['closed_coordinate_maps'][r]['rank'] for r in ROLES)}",
        f"A/C used union={coordinate['A_C_used_coordinate_union']} missing={coordinate['A_C_common_missing_tile']} "
        f"left-invertible with apex={coordinate['ordered_A_C_has_left_inverse_with_apex']}",
        f"B idempotent={coordinate['B_idempotent']} exceptional={coordinate['B_identity_except']}",
        f"B lost seam={coordinate['B_lost_seam']} exact inverse={coordinate['B_exact_inverse_with_seam']}",
        f"Phi image codimension={result['pure_Phi_image']['codimension']} doubled fibres={result['pure_Phi_image']['doubled_coordinate_fibres']}",
        f"defect action={result['defect_action']}",
        "",
        "FINITE 58-CELL BANK",
        f"B changed/seam histogram={bank['B_changed_representatives']}/{bank['B_source_seam_histogram']}",
        f"rolewise Q/bar-P descent={bank['rolewise_source_pair_Q_descent']}/{bank['rolewise_source_pair_bar_P_descent']}",
        f"role Q fibres={bank['role_representative_Q_fibres']}",
        f"role bar-P fibres={bank['role_representative_bar_P_fibres']}",
        f"ordered A/C Q/bar-P={bank['ordered_A_C_Q_deck']}/{bank['ordered_A_C_bar_P_deck']}",
        f"unordered A/C Q/bar-P={bank['unordered_A_C_Q_deck']}/{bank['unordered_A_C_bar_P_deck']}",
        f"A/C bar-P + B bar-P={bank['A_C_bar_P_plus_B_bar_P']}",
        f"theta orbits/seam/sheet preserved={bank['theta_orbits']}/{bank['theta_seam_orbit_histogram']}/{bank['theta_sheet_preserved_by_B']}",
        f"theta source=B=-A/C-role orientation={bank['theta_source_B_orientation_equals_negative_A_C_role_order']}",
        f"B node margins/coupled={bank['B_direct_node_fibres']}/{bank['B_partner_node_fibres']}/{bank['B_coupled_bar_P_fibres']}",
        f"B target node + input seam={bank['B_direct_node_plus_input_seam_fibres']}/{bank['B_partner_node_plus_input_seam_fibres']}",
        f"B image masks/nodes/loops={bank['B_image_masks']}/{bank['B_image_merged_nodes']}/{bank['B_coupled_bar_P_loops']}",
        "",
        "THM-840 OPERATION KERNELS (split source cells)",
        f"  { {name: row['split_source_cells'] for name, row in kernels.items()} }",
        f"TOURNAMENT ANALYSIS vertices={len(ta['vertices'])} edge-flips={ta['edge_flips']}",
        f"  retention={ta['retention']}",
        "  both gauges transitive: C3=0, singleton SCCs, one Hamiltonian path",
        "PRESERVATION: literal maps/Q/local nodes/coupled bar-P/seam/theta on the finite bank",
        "DESTROYS: all LRC metric/owner/wall data and arbitrary CF-word continuation",
        "CHALLENGED ASSUMPTION: seam, sheet, and node fibre are not one universal bit",
        "ALL ASSERTIONS PASSED",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--witnesses",
        type=Path,
        default=Path("05-knowledge/results/mobius_cech_n9_exact_join_witnesses_codex_S13.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("05-knowledge/results/n9_closed_cf_lift_face_kernels_codex_S13f.out"),
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("05-knowledge/results/n9_closed_cf_lift_face_kernels_codex_S13f.json"),
    )
    args = parser.parse_args()
    result = run(args.witnesses)
    text = render(result)
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(text, end="")


if __name__ == "__main__":
    main()
