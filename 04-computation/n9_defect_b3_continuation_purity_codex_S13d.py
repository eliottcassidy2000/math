#!/usr/bin/env python3
"""Exact A/B/C and deletion continuation audit for THM-842.

The input is THM-828's 58-witness table plus the certified n=8 merged-node
atlas.  Literal Q means the orbit of a complement line under staircase
reflection.  The script keeps face roles, repaired fixed-path cards, actual
induced tournaments, path position, shortcut seams, defect masks, skew words,
and endpoint node coupling distinct throughout.

Tournament Analysis uses information carriers as vertices.  The pairwise
observable is the number of unordered source-witness pairs separated.  The
two gauges are raw retention and retention per logarithmic cell cost; ties use
the displayed carrier order.  These carrier tournaments preserve only finite
bank separation.  They destroy LRC gaps, owners, wall schedules, and metric
loneliness.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import struct
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path


def tiles(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((a, b) for b in range(1, n - 1) for a in range(n, b + 1, -1) if a - b >= 2)


@lru_cache(maxsize=None)
def tile_index(n: int) -> dict[tuple[int, int], int]:
    return {cell: bit for bit, cell in enumerate(tiles(n))}


def complement(x: int, n: int) -> int:
    return x ^ ((1 << len(tiles(n))) - 1)


def reflect(x: int, n: int) -> int:
    index = tile_index(n)
    result = 0
    for bit, (a, b) in enumerate(tiles(n)):
        result |= ((x >> bit) & 1) << index[(n - b + 1, n - a + 1)]
    return result


def q_cell(x: int, n: int) -> int:
    return min(x, complement(x, n), reflect(x, n), complement(reflect(x, n), n))


def face(x: int, n: int, role: str) -> int:
    upper = tile_index(n)
    result = 0
    for bit, (a, b) in enumerate(tiles(n - 1)):
        source = {"A": (a, b), "B": (a + 1, b), "C": (a + 1, b + 1)}[role]
        result |= ((x >> upper[source]) & 1) << bit
    return result


def repaired_card(x: int, n: int, deleted: int) -> tuple[int, int | None]:
    upper = tile_index(n)
    result = 0
    for bit, (a, b) in enumerate(tiles(n - 1)):
        source = (a + int(a >= deleted), b + int(b >= deleted))
        result |= ((x >> upper[source]) & 1) << bit
    seam = None if deleted in (1, n) else (x >> upper[(deleted + 1, deleted - 1)]) & 1
    return result, seam


def skew_word(x: int, n: int) -> tuple[int, ...]:
    index = tile_index(n)
    position: dict[int, int] = defaultdict(int)
    word = [0] * (n - 3)
    for a, b in tiles(n):
        tau = a + b - 1
        if tau < n:
            p = position[tau]
            position[tau] += 1
            i = index[(a, b)]
            j = index[(n - b + 1, n - a + 1)]
            word[tau - 3] += p * (((x >> i) & 1) - ((x >> j) & 1))
    return tuple(word)


def eta(word: tuple[int, ...]) -> int:
    return next((1 if z > 0 else -1 for z in word if z), 0)


def first_skew_layer(word: tuple[int, ...]) -> int | None:
    return next((i + 3 for i, z in enumerate(word) if z), None)


def theta(x: int, n: int) -> int:
    apex = 1 << tile_index(n)[(n, 1)]
    nonapex = ((1 << len(tiles(n))) - 1) ^ apex
    return reflect(x ^ nonapex, n)


def gf2_rank(rows: list[int]) -> int:
    rows = rows[:]
    rank = 0
    while rows:
        pivot = max(rows)
        if not pivot:
            break
        rows.remove(pivot)
        bit = pivot.bit_length() - 1
        rows = [row ^ pivot if (row >> bit) & 1 else row for row in rows]
        rank += 1
    return rank


def fibre_hist(values: list[object]) -> tuple[int, dict[int, int]]:
    counts = Counter(values)
    return len(counts), dict(sorted(Counter(counts.values()).items()))


def partition(values: list[object]) -> tuple[tuple[int, ...], ...]:
    groups: dict[object, list[int]] = defaultdict(list)
    for i, value in enumerate(values):
        groups[value].append(i)
    return tuple(sorted(tuple(group) for group in groups.values()))


def adjacency(x: int, n: int) -> tuple[int, ...]:
    adj = [0] * n
    for v in range(2, n + 1):
        adj[v - 1] |= 1 << (v - 2)
    for bit, (a, b) in enumerate(tiles(n)):
        if (x >> bit) & 1:
            adj[a - 1] |= 1 << (b - 1)
        else:
            adj[b - 1] |= 1 << (a - 1)
    return tuple(adj)


def delete_adjacency(adj: tuple[int, ...], deleted: int) -> tuple[int, ...]:
    keep = [v for v in range(len(adj)) if v != deleted - 1]
    return tuple(sum(int(bool(adj[v] & (1 << w))) << j for j, w in enumerate(keep)) for v in keep)


def redei_path(adj: tuple[int, ...]) -> tuple[int, ...]:
    path: list[int] = []
    for v in range(len(adj)):
        place = next((i for i, w in enumerate(path) if adj[v] & (1 << w)), len(path))
        path.insert(place, v)
    assert all(adj[path[i]] & (1 << path[i + 1]) for i in range(len(path) - 1))
    return tuple(path)


def encode_along_path(adj: tuple[int, ...], path: tuple[int, ...]) -> int:
    n = len(adj)
    by_label = {n - i: vertex for i, vertex in enumerate(path)}
    result = 0
    for bit, (a, b) in enumerate(tiles(n)):
        result |= int(bool(adj[by_label[a]] & (1 << by_label[b]))) << bit
    return result


def induced_node(x: int, deleted: int, node8: tuple[int, ...]) -> int:
    card_adj = delete_adjacency(adjacency(x, 9), deleted)
    mask = encode_along_path(card_adj, redei_path(card_adj))
    return node8[mask]


def ta_fingerprint(carriers: list[tuple[str, list[object]]]) -> dict[str, object]:
    total_pairs = math.comb(len(carriers[0][1]), 2)
    data = []
    for order, (name, values) in enumerate(carriers):
        counts = Counter(values)
        separated = total_pairs - sum(math.comb(c, 2) for c in counts.values())
        cells = len(counts)
        data.append((name, separated, cells, order))
    retention = sorted(data, key=lambda z: (z[1], -z[3]))
    efficiency = sorted(data, key=lambda z: (z[1] / math.log2(z[2] + 1), -z[3]))
    rp = {z[0]: i for i, z in enumerate(retention)}
    ep = {z[0]: i for i, z in enumerate(efficiency)}
    flips = sum((rp[a[0]] < rp[b[0]]) != (ep[a[0]] < ep[b[0]]) for i, a in enumerate(data) for b in data[i + 1 :])
    return {
        "vertices": [z[0] for z in data],
        "retention": {z[0]: {"separated_pairs": z[1], "cells": z[2]} for z in data},
        "switches": ["raw retention", "retention per logarithmic cell cost"],
        "tie_hamiltonian_path": [z[0] for z in carriers],
        "score_histogram": {str(i): 1 for i in range(len(data))},
        "directed_triangles": 0,
        "scc_sizes": [1] * len(data),
        "hamiltonian_paths": 1,
        "edge_flips": flips,
    }


def run(witness_path: Path, atlas_path: Path, forward_path: Path) -> dict[str, object]:
    rows = []
    with witness_path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append({key: int(row[key], 16) for key in ("D", "u", "v", "xA", "xB", "xC")})
    assert len(rows) == 58
    raw = atlas_path.read_bytes()
    assert len(raw) == 2 * (1 << 21)
    node8 = struct.unpack("<" + str(1 << 21) + "H", raw)
    forward = json.loads(forward_path.read_text())
    upper_nodes = [tuple(forward["pairs"][i]["merged_node_key"]) for i in range(58)]
    assert [int(forward["pairs"][i]["u"], 16) for i in range(58)] == [row["u"] for row in rows]

    full8 = (1 << 21) - 1
    pbar = lambda x: tuple(sorted((node8[x], node8[x ^ full8])))
    parrow = lambda x: (node8[x], node8[x ^ full8])
    roles = ("A", "B", "C")
    role_masks: dict[str, list[tuple[int, int]]] = {r: [] for r in roles}
    d_images: dict[int, dict[str, int]] = {}
    upper_eta, upper_j = [], []
    for row in rows:
        upper_j.append(skew_word(row["u"], 9))
        upper_eta.append(eta(upper_j[-1]))
        d_images.setdefault(row["D"], {})
        for role in roles:
            a, b = face(row["u"], 9, role), face(row["v"], 9, role)
            assert a == row["x" + role] and (a ^ b) == face(row["D"], 9, role)
            role_masks[role].append((a, b))
            d_images[row["D"]][role] = a ^ b

    basis = [0x0192486, 0x08C2C0C, 0x11B4600, 0x4483414]
    for role in roles:
        assert gf2_rank([face(d, 9, role) for d in basis]) == 4
        assert len({data[role] for data in d_images.values()}) == 11
    assert all(data["A"] == reflect(data["C"], 8) and data["B"] == reflect(data["B"], 8) for data in d_images.values())

    q_role = {r: [q_cell(a, 8) for a, _ in role_masks[r]] for r in roles}
    p_role = {r: [pbar(a) for a, _ in role_masks[r]] for r in roles}
    q_pair_descent = {r: sum(q_cell(a, 8) == q_cell(b, 8) for a, b in role_masks[r]) for r in roles}
    p_pair_descent = {r: sum(pbar(a) == pbar(b) for a, b in role_masks[r]) for r in roles}
    assert q_pair_descent == {"A": 0, "B": 58, "C": 0}
    assert p_pair_descent == {"A": 58, "B": 58, "C": 58}
    assert all(q_cell(a, 8) == q_cell(b, 8) for a, b in role_masks["B"])
    assert all(pbar(a) == pbar(b) for a, b in role_masks["B"])
    assert fibre_hist(q_role["B"]) == (58, {1: 58})
    assert fibre_hist(p_role["B"]) == (58, {1: 58})

    endpoint_q = [tuple(sorted((q_role["A"][i], q_role["C"][i]))) for i in range(58)]
    endpoint_p = [tuple(sorted((p_role["A"][i], p_role["C"][i]))) for i in range(58)]
    assert all(parrow(role_masks["A"][i][0]) == parrow(role_masks["C"][i][0]) for i in range(58))
    assert fibre_hist(endpoint_q) == (29, {2: 29}) and all(a != b for a, b in endpoint_q)
    assert fibre_hist(endpoint_p) == (29, {2: 29}) and all(a == b for a, b in endpoint_p)
    joint_q = [(endpoint_q[i], q_role["B"][i]) for i in range(58)]
    joint_p = [(endpoint_p[i], p_role["B"][i]) for i in range(58)]
    assert fibre_hist(joint_q) == (58, {1: 58}) and fibre_hist(joint_p) == (58, {1: 58})

    reps = {row["u"]: i for i, row in enumerate(rows)}
    companions = []
    seen = set()
    for i, row in enumerate(rows):
        mate = theta(row["u"], 9)
        assert mate in reps and mate != row["u"]
        j = reps[mate]
        assert rows[j]["D"] == row["D"] and upper_j[j] == upper_j[i]
        assert endpoint_q[j] == endpoint_q[i] and endpoint_p[j] == endpoint_p[i]
        assert face(mate, 9, "A") == reflect(complement(face(row["u"], 9, "C"), 8), 8)
        assert face(mate, 9, "C") == reflect(complement(face(row["u"], 9, "A"), 8), 8)
        assert face(mate, 9, "B") == theta(face(row["u"], 9, "B"), 8)
        assert q_role["B"][j] != q_role["B"][i]
        assert parrow(role_masks["A"][i][0]) == tuple(reversed(parrow(role_masks["A"][j][0])))
        if i not in seen:
            seen.update((i, j))
            companions.append((min(i, j), max(i, j), row["D"]))
    assert len(companions) == 29 and partition(endpoint_q) == tuple(sorted((a, b) for a, b, _ in companions))

    relative: dict[str, Counter[str]] = {r: Counter() for r in roles}
    eta_values: dict[str, Counter[int]] = {r: Counter() for r in roles}
    first_layers: dict[str, Counter[object]] = {r: Counter() for r in roles}
    role_j: dict[str, list[tuple[int, ...]]] = {r: [] for r in roles}
    for role in roles:
        for i, (a, b) in enumerate(role_masks[role]):
            word = skew_word(a, 8)
            role_j[role].append(word)
            e = eta(word)
            eta_values[role][e] += 1
            first_layers[role][first_skew_layer(word)] += 1
            relative[role]["zero" if e == 0 else "preserve" if e == upper_eta[i] else "reverse"] += 1
            if role == "B":
                assert skew_word(b, 8) == tuple(-z for z in word)
    assert relative["A"] == relative["C"] == Counter({"preserve": 31, "reverse": 23, "zero": 4})
    assert eta_values["A"] == eta_values["C"] == Counter({1: 33, -1: 21, 0: 4})
    assert first_layers["A"] == first_layers["C"] == Counter({5: 25, 6: 18, 7: 11, None: 4})
    assert relative["B"] == Counter({"reverse": 38, "preserve": 20})
    assert eta_values["B"] == Counter({-1: 52, 1: 6})
    assert first_layers["B"] == Counter({5: 44, 6: 12, 7: 2})
    b_words_by_d: dict[int, set[tuple[int, ...]]] = defaultdict(set)
    b_relative_by_d: dict[int, set[str]] = defaultdict(set)
    for i, row in enumerate(rows):
        b_words_by_d[row["D"]].add(role_j["B"][i])
        b_relative_by_d[row["D"]].add("preserve" if eta(role_j["B"][i]) == upper_eta[i] else "reverse")
    assert all(len(words) == 1 and next(iter(words)) != (0,) * 5 for words in b_words_by_d.values())
    preserve_sectors = {0x0192486, 0x08C2C0C, 0x095088A, 0x18E4E8A, 0x1976A0C, 0x5C67A9E}
    assert all(states == {"preserve" if D in preserve_sectors else "reverse"} for D, states in b_relative_by_d.items())
    endpoint_p_eta = [(endpoint_p[i], eta(role_j["A"][i]), eta(role_j["C"][i])) for i in range(58)]
    assert fibre_hist(endpoint_p_eta) == (50, {1: 42, 2: 8})

    repaired_q: dict[int, list[int]] = {}
    repaired_p: dict[int, list[tuple[int, int]]] = {}
    seam_words = []
    repaired_relative: dict[int, Counter[str]] = {}
    b_equals_repaired = {}
    for deleted in range(1, 10):
        repaired_q[deleted], repaired_p[deleted] = [], []
        rel = Counter()
        exact_equal = q_equal = 0
        for i, row in enumerate(rows):
            card, _ = repaired_card(row["u"], 9, deleted)
            repaired_q[deleted].append(q_cell(card, 8))
            repaired_p[deleted].append(pbar(card))
            e = eta(skew_word(card, 8))
            rel["zero" if e == 0 else "preserve" if e == upper_eta[i] else "reverse"] += 1
            if 2 <= deleted <= 8:
                bmask = role_masks["B"][i][0]
                exact_equal += card == bmask
                q_equal += q_cell(card, 8) == q_cell(bmask, 8)
        repaired_relative[deleted] = rel
        if 2 <= deleted <= 8:
            b_equals_repaired[deleted] = (exact_equal, q_equal)
    expected_relative = {
        1: (31, 23, 4), 2: (31, 23, 4), 3: (21, 26, 11), 4: (20, 28, 10),
        5: (52, 6, 0), 6: (20, 28, 10), 7: (21, 26, 11), 8: (31, 23, 4), 9: (31, 23, 4),
    }
    assert all(repaired_relative[j] == Counter({"preserve": p, "reverse": r, "zero": z}) for j, (p, r, z) in expected_relative.items())
    assert b_equals_repaired == {2: (0, 0), 3: (0, 0), 4: (1, 1), 5: (2, 2), 6: (1, 1), 7: (0, 0), 8: (0, 0)}

    orbit_expected = {
        (1, 9): ((29, {2: 29}), (29, {2: 29}), 0, 58),
        (2, 8): ((46, {1: 34, 2: 12}), (46, {1: 34, 2: 12}), 0, 8),
        (3, 7): ((58, {1: 58}), (58, {1: 58}), 0, 0),
        (4, 6): ((58, {1: 58}), (58, {1: 58}), 0, 0),
        (5,): ((58, {1: 58}), (58, {1: 58}), 0, 0),
    }
    repaired_orbits = {}
    for orbit, expected in orbit_expected.items():
        q_values = [tuple(sorted(repaired_q[j][i] for j in orbit)) for i in range(58)]
        p_values = [tuple(sorted(repaired_p[j][i] for j in orbit)) for i in range(58)]
        q_coincide = sum(len(set(value)) < len(value) for value in q_values)
        p_coincide = sum(len(set(value)) < len(value) for value in p_values)
        actual = (fibre_hist(q_values), fibre_hist(p_values), q_coincide, p_coincide)
        assert actual == expected
        q_seam_values = []
        p_seam_values = []
        for i, row in enumerate(rows):
            q_tokens, p_tokens = [], []
            for j in orbit:
                _, seam = repaired_card(row["u"], 9, j)
                token = -1 if seam is None else seam
                q_tokens.append((repaired_q[j][i], token))
                p_tokens.append((repaired_p[j][i], token))
            q_seam_values.append(tuple(sorted(q_tokens)))
            p_seam_values.append(tuple(sorted(p_tokens)))
        assert fibre_hist(q_seam_values) == actual[0] and fibre_hist(p_seam_values) == actual[1]
        repaired_orbits[str(orbit)] = {
            "Q": actual[0], "bar_P": actual[1],
            "Q_role_coincidences": q_coincide, "bar_P_role_coincidences": p_coincide,
            "seam_augmented_same_fibres": True,
        }
    full_repaired_q = [min(tuple(repaired_q[j][i] for j in range(1, 10)), tuple(repaired_q[j][i] for j in range(9, 0, -1))) for i in range(58)]
    assert fibre_hist(full_repaired_q) == (58, {1: 58})

    for row in rows:
        seams_u = tuple(repaired_card(row["u"], 9, j)[1] for j in range(2, 9))
        seams_v = tuple(repaired_card(row["v"], 9, j)[1] for j in range(2, 9))
        assert seams_v == tuple(reversed(seams_u))
        seam_words.append(seams_u)
    seam_ones = [sum(word[j - 2] for word in seam_words) for j in range(2, 9)]
    assert seam_ones == [29, 48, 31, 29, 27, 10, 29]
    assert fibre_hist(seam_words) == (22, {1: 8, 2: 4, 3: 2, 4: 6, 6: 2})
    assert Counter(map(sum, seam_words)) == Counter({1: 5, 2: 9, 3: 15, 4: 15, 5: 9, 6: 5})
    assert Counter(sum(a != b for a, b in zip(word, reversed(word))) for word in seam_words) == Counter({0: 10, 2: 42, 4: 6})

    induced = [[induced_node(row["u"], j, node8) for j in range(1, 10)] for row in rows]
    induced_v = [[induced_node(row["v"], j, node8) for j in range(1, 10)] for row in rows]
    assert all(induced_v[i][j - 1] == induced[i][9 - j] for i in range(58) for j in range(1, 10))
    induced_expected = {
        (1, 9): ((49, {1: 41, 2: 7, 3: 1}), 58),
        (2, 8): ((45, {1: 33, 2: 11, 3: 1}), 30),
        (3, 7): ((55, {1: 52, 2: 3}), 9),
        (4, 6): ((55, {1: 52, 2: 3}), 0),
        (5,): ((53, {1: 48, 2: 5}), 0),
    }
    induced_orbits = {}
    for orbit, expected in induced_expected.items():
        values = [tuple(sorted(induced[i][j - 1] for j in orbit)) for i in range(58)]
        actual = (fibre_hist(values), sum(len(set(value)) < len(value) for value in values))
        assert actual == expected
        induced_orbits[str(orbit)] = {"fibres": actual[0], "role_coincidences": actual[1]}
    ordered_induced = [min(tuple(row), tuple(reversed(row))) for row in induced]
    unpositioned_multiset = [tuple(sorted(row)) for row in induced]
    unpositioned_set = [tuple(sorted(set(row))) for row in induced]
    assert fibre_hist(ordered_induced) == (58, {1: 58})
    assert fibre_hist(unpositioned_multiset) == fibre_hist(unpositioned_set) == (53, {1: 48, 2: 5})
    assert partition(unpositioned_multiset) == partition(unpositioned_set) == partition(upper_nodes)
    assert partition([row[4] for row in induced]) == partition(upper_nodes)

    carriers = [
        ("D sector", [row["D"] for row in rows]),
        ("endpoint Q", endpoint_q),
        ("endpoint bar P", endpoint_p),
        ("upper node", upper_nodes),
        ("unpositioned induced deck", unpositioned_multiset),
        ("B Q", q_role["B"]),
        ("B bar P", p_role["B"]),
        ("endpoint+B Q", joint_q),
        ("ordered induced deck", ordered_induced),
    ]
    ta = ta_fingerprint(carriers)

    return {
        "schema_version": 1,
        "theorem": "THM-842",
        "inputs": {"witnesses": str(witness_path), "node_atlas": str(atlas_path), "forward_map": str(forward_path)},
        "source_cells": 58,
        "role_geometry": {
            "A": "delete fixed-path vertex 9",
            "B": "gap contraction, not one vertex deletion",
            "C": "delete fixed-path vertex 1",
            "reflection": ["A sigma = sigma C", "B sigma = sigma B", "C sigma = sigma A"],
        },
        "role_rank_on_V": {role: 4 for role in roles},
        "role_D_images": {f"0x{D:07x}": {role: f"0x{value:06x}" for role, value in data.items()} for D, data in sorted(d_images.items())},
        "face_continuation": {
            "rolewise_Q_descent": q_pair_descent, "rolewise_bar_P_descent": p_pair_descent,
            "B_Q": fibre_hist(q_role["B"]), "B_bar_P": fibre_hist(p_role["B"]),
            "endpoint_Q": fibre_hist(endpoint_q), "endpoint_bar_P": fibre_hist(endpoint_p),
            "endpoint_plus_B_Q": fibre_hist(joint_q), "endpoint_plus_B_bar_P": fibre_hist(joint_p),
            "endpoint_Q_role_coincidences": 0, "endpoint_bar_P_role_coincidences": 58,
        },
        "apex_relative_antipode": {
            "formula": "theta_n(u)=sigma_n(u xor (FULL_n xor apex))",
            "fixed_points": 0, "orbits": 29,
            "preserves": ["D", "full upper skew word J", "upper eta"],
            "endpoint_identity": "A theta = sigma kappa C; C theta = sigma kappa A",
            "gap_identity": "B theta = theta B",
            "golden_companions": [[i, j, f"0x{D:07x}"] for i, j, D in sorted(companions)],
        },
        "chirality": {
            role: {
                "eta_histogram": dict(sorted(eta_values[role].items())),
                "first_layer_histogram": {"none" if k is None else str(k): v for k, v in sorted(first_layers[role].items(), key=lambda z: -1 if z[0] is None else z[0])},
                "relative_to_upper": dict(relative[role]),
                "sector_constant": role == "B",
            } for role in roles
        },
        "endpoint_bar_P_plus_A_C_eta": fibre_hist(endpoint_p_eta),
        "repaired_cards": {
            "position_orbits": repaired_orbits,
            "full_ordered_mod_reversal": fibre_hist(full_repaired_q),
            "B_exact_and_Q_equal_by_internal_position": {str(j): list(v) for j, v in b_equals_repaired.items()},
            "eta_relative_by_position": {str(j): dict(c) for j, c in repaired_relative.items()},
            "shortcut_seams": {
                "ones_by_j2_to_j8": seam_ones,
                "word_fibres": fibre_hist(seam_words),
                "weight_histogram": dict(sorted(Counter(map(sum, seam_words)).items())),
            },
        },
        "induced_tournament_decks": {
            "position_orbits": induced_orbits,
            "ordered_nine_cards_mod_reversal": fibre_hist(ordered_induced),
            "unpositioned_multiset": fibre_hist(unpositioned_multiset),
            "unpositioned_set": fibre_hist(unpositioned_set),
            "unpositioned_partition_equals_upper_node": True,
            "central_node_partition_equals_upper_node": True,
        },
        "bare_sector_guardrail": "all four missing sectors have B footprints before A/C or raw-S2 death",
        "tournament_analysis": ta,
        "methodology": {
            "alternate_vertices": ["face roles", "proof obligations", "information carriers", "path positions"],
            "challenged_assumption": "tournament vertices and continuation states need not be runners or unmarked classes",
            "preserves": ["finite-bank equality", "D action", "Q/bar-P fibres", "skew", "path position when retained"],
            "destroys": ["LRC gaps", "owners", "walls", "metric loneliness"],
        },
    }


def render(result: dict[str, object]) -> str:
    face_data = result["face_continuation"]
    anti = result["apex_relative_antipode"]
    repaired = result["repaired_cards"]
    induced = result["induced_tournament_decks"]
    ta = result["tournament_analysis"]
    lines = [
        "THM-842 N=9 DEFECT B3 CONTINUATION PURITY",
        "=" * 68,
        f"source cells={result['source_cells']} role ranks A/B/C=4/4/4",
        "geometry: A=delete vertex9, C=delete vertex1, B=gap contraction (not vertex deletion)",
        f"B Q/bar-P={face_data['B_Q']}/{face_data['B_bar_P']}",
        f"rolewise source-pair descent Q={face_data['rolewise_Q_descent']} bar-P={face_data['rolewise_bar_P_descent']}",
        f"endpoint Q/bar-P={face_data['endpoint_Q']}/{face_data['endpoint_bar_P']}",
        f"endpoint+B Q/bar-P={face_data['endpoint_plus_B_Q']}/{face_data['endpoint_plus_B_bar_P']}",
        f"endpoint role coincidences Q/bar-P={face_data['endpoint_Q_role_coincidences']}/{face_data['endpoint_bar_P_role_coincidences']}",
        f"theta fixed/orbits={anti['fixed_points']}/{anti['orbits']} preserves D+full-J+eta",
        "theta identities: A theta=sigma kappa C; C theta=sigma kappa A; B theta=theta B",
        f"chirality A={result['chirality']['A']} C={result['chirality']['C']} B={result['chirality']['B']}",
        f"endpoint bar-P+(etaA,etaC)={result['endpoint_bar_P_plus_A_C_eta']}",
        f"repaired position orbits={repaired['position_orbits']}",
        f"B equals repaired internal (exact,Q)={repaired['B_exact_and_Q_equal_by_internal_position']}",
        f"shortcut seams={repaired['shortcut_seams']}",
        f"induced position orbits={induced['position_orbits']}",
        f"induced ordered/unpositioned={induced['ordered_nine_cards_mod_reversal']}/{induced['unpositioned_multiset']}",
        "unpositioned and central induced-node partitions equal the upper bare-node partition: YES",
        f"TOURNAMENT ANALYSIS vertices={len(ta['vertices'])} edge-flips={ta['edge_flips']}",
        f"  retention={ta['retention']}",
        "  both gauges transitive: C3=0, singleton SCCs, one Hamiltonian path",
        "PRESERVATION: finite bank, D/Q/bar-P/skew/path-position only",
        "DESTROYS: LRC gaps, owners, walls, metric loneliness",
        "CHALLENGED ASSUMPTION: use face roles/information carriers, not only runners or unmarked nodes",
        "ALL ASSERTIONS PASSED",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--witnesses", type=Path, default=Path("05-knowledge/results/mobius_cech_n9_exact_join_witnesses_codex_S13.tsv"))
    parser.add_argument("--node-atlas", type=Path, default=Path("/tmp/n8_merged_node_rank_u16.bin"))
    parser.add_argument("--forward-map", type=Path, default=Path("05-knowledge/results/n9_false_palindrome_node_forward_map_codex_S13b.json"))
    parser.add_argument("--output", type=Path, default=Path("05-knowledge/results/n9_defect_b3_continuation_purity_codex_S13d.out"))
    parser.add_argument("--json", type=Path, default=Path("05-knowledge/results/n9_defect_b3_continuation_purity_codex_S13d.json"))
    args = parser.parse_args()
    result = run(args.witnesses, args.node_atlas, args.forward_map)
    text = render(result)
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(text, end="")


if __name__ == "__main__":
    main()
