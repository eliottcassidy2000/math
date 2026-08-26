#!/usr/bin/env python3
"""Exact minimal-anchor hierarchy for zero-original q=50 bodies.

The load-bearing blocker enumeration in this audit is deliberately complete:
if a repair hypergraph has a cover of size below the target budget, every
extension through that budget is also enumerated.  Stopping at first-hit or
inclusion-minimal covers is not sufficient for blocker descent.
"""

from collections import Counter, defaultdict
import hashlib
import importlib.util
from itertools import combinations
from math import comb, gcd
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")

HERE = Path(__file__).resolve().parent
FOUR_PATH = HERE / "lrc14_q50_zero_original_four_anchor_transfer.py"
FOUR_SHA256 = "1e1ec137b27b8555b2c9759a3e283746e09f445a8bd6bdba25ccd94187b05958"

EXPECTED_ANCHOR_COUNTS = {1: 0, 2: 0, 3: 0, 4: 32, 5: 297, 6: 24}
EXPECTED_ANCHOR_SHA256 = {
    1: "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    2: "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    3: "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    4: "359881969de938d942ddf69501c0d66fee783fe12910e348b915a23c82e3b560",
    5: "fdde1bd87afcb906da126ee5a8cd43c7f3eb0cabc1a67f844c150fb3f7cc96d6",
    6: "34564111d2efd367ca6a23f02eef4099760d055b08e5cd9e4e0b30f981ce07d5",
}
EXPECTED_SIZE6 = tuple(sorted(
    tuple(sorted((42, 63, 132, 286, left, right)))
    for left in (8, 16, 88, 176)
    for right in (10, 20, 30, 170, 190, 290)
))
EXPECTED_LAYER_DATA = {
    5: (3_017, "2c01c096640ff8690496a424fce357cdb38d0a29cafc81bacfe94eea64729f06"),
    6: (85_324, "503ffc5e72825f275a9b351c1ceab1b13c452c5082ec4009020bb858c44aafdd"),
    7: (821_737, "2a2b26ffd0c26000260ead13f507955b332d5c61927219418a49225b1a004b7d"),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_four():
    require(hashlib.sha256(FOUR_PATH.read_bytes()).hexdigest() == FOUR_SHA256,
            "THM-4182 primary dependency hash changed")
    spec = importlib.util.spec_from_file_location("thm4182_base", FOUR_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-4182 audit")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def valid_anchor(labels):
    return gcd(*labels) == 1 and all(
        any(label % modulus == 0 for label in labels)
        for modulus in range(2, 15)
    )


def enumerate_minimal_anchors(base, ground):
    layers = {}
    previous = []
    for size in range(1, 7):
        anchors = tuple(
            labels
            for labels in combinations(ground, size)
            if valid_anchor(labels)
            and not any(set(old).issubset(labels) for old in previous)
        )
        require(len(anchors) == EXPECTED_ANCHOR_COUNTS[size],
                ("minimal-anchor count changed", size, len(anchors)))
        digest = hashlib.sha256()
        for anchor in anchors:
            digest.update(base.mask_of(anchor).to_bytes(4, "little"))
        require(digest.hexdigest() == EXPECTED_ANCHOR_SHA256[size],
                ("minimal-anchor labels changed", size, digest.hexdigest()))
        layers[size] = anchors
        previous.extend(anchors)
    require(layers[6] == EXPECTED_SIZE6, "size-six product classification changed")
    require(all(286 in anchor for anchor in layers[4] + layers[5]),
            "286 is not forced at sizes four and five")
    return layers


def edge_incidence(edges, vertices):
    incidence = {vertex: 0 for vertex in vertices}
    for edge_index, edge in enumerate(edges):
        bits = edge
        while bits:
            bit = bits & -bits
            vertex = bit.bit_length() - 1
            if vertex in incidence:
                incidence[vertex] |= 1 << edge_index
            bits ^= bit
    return incidence, (1 << len(edges)) - 1


def complete_covers(edges, vertices, budget):
    """Return every labelled cover of cardinality at most ``budget``."""
    incidence, all_edges = edge_incidence(edges, vertices)
    covers = []
    for size in range(budget + 1):
        for chosen_vertices in combinations(vertices, size):
            covered = 0
            chosen = 0
            for vertex in chosen_vertices:
                covered |= incidence[vertex]
                chosen |= 1 << vertex
            if covered == all_edges:
                covers.append(chosen)
    return tuple(sorted(covers))


def filter_covers(covers, edges, vertices):
    incidence, all_edges = edge_incidence(edges, vertices)
    survivors = []
    for chosen in covers:
        covered = 0
        bits = chosen
        while bits:
            bit = bits & -bits
            covered |= incidence[bit.bit_length() - 1]
            bits ^= bit
        if covered == all_edges:
            survivors.append(chosen)
    return tuple(survivors)


def fnv1a64_add_u64(value, word):
    for shift in range(0, 64, 8):
        value ^= (word >> shift) & 0xFF
        value = (value * 0x100000001B3) & ((1 << 64) - 1)
    return value


def update_row_digest(value, anchor_mask, active_count, covers, survivors):
    for word in (anchor_mask, active_count, len(covers)):
        value = fnv1a64_add_u64(value, word)
    for cover in covers:
        value = fnv1a64_add_u64(value, cover)
    value = fnv1a64_add_u64(value, len(survivors))
    for survivor in survivors:
        value = fnv1a64_add_u64(value, survivor)
    return value


def audit_size_five(base, anchors, edges5, edges6, edges7):
    total_histogram = Counter()
    no_cover_rows = 0
    rows_with_covers = 0
    residual_presentations = []
    digest = 0xCBF29CE484222325

    for anchor in anchors:
        anchor_mask = base.mask_of(anchor)
        vertices = tuple(index for index in range(30) if not anchor_mask >> index & 1)
        active5 = tuple(edge for edge in edges5 if not edge & anchor_mask)
        covers5 = complete_covers(active5, vertices, 5)
        if covers5:
            rows_with_covers += 1
        else:
            no_cover_rows += 1
        total_histogram.update(cover.bit_count() for cover in covers5)

        active6 = tuple(edge for edge in edges6 if not edge & anchor_mask)
        covers6 = filter_covers(covers5, active6, vertices)
        residual_presentations.extend((anchor, anchor_mask, cover) for cover in covers6)
        digest = update_row_digest(digest, anchor_mask, len(active5), covers5, covers6)
        print(
            "SIZE5_ROW", "ANCHOR", anchor, "E5_EDGES", len(active5),
            "E5_COVERS_THROUGH_5", len(covers5),
            "E5_COVER_SIZE_HIST", tuple(sorted(Counter(
                cover.bit_count() for cover in covers5
            ).items())),
            "E6_COVERS_THROUGH_5", len(covers6),
        )

    require((no_cover_rows, rows_with_covers) == (41, 256),
            ("size-five 41/256 split changed", no_cover_rows, rows_with_covers))
    require(sum(total_histogram.values()) == 21_506,
            ("complete size-five cover count changed", total_histogram))
    require(dict(total_histogram) == {3: 33, 4: 1_229, 5: 20_244},
            ("complete size-five cover histogram changed", total_histogram))
    require(len(residual_presentations) == 15,
            ("size-five E6 residual count changed", len(residual_presentations)))

    bodies = defaultdict(list)
    for anchor, anchor_mask, cover in residual_presentations:
        body = anchor_mask | cover
        bodies[body].append((anchor, cover))
    expected_bodies = {
        base.mask_of((16, 88, 95, 145, 168, 170, 193, 240, 252, 286)): (9, 736),
        base.mask_of((80, 88, 95, 145, 168, 170, 193, 240, 252, 286)): (6, 514),
    }
    require(set(bodies) == set(expected_bodies), "size-five residual bodies changed")
    for anchor, anchor_mask, cover in residual_presentations:
        print("SIZE5_E6_BLOCKER", "ANCHOR", anchor, "K", base.labels_of(cover))
        require(all(edge & cover for edge in edges6 if not edge & anchor_mask),
                ("declared E6 blocker failed", anchor, base.labels_of(cover)))
    for body, presentations in bodies.items():
        expected_multiplicity, expected_witnesses = expected_bodies[body]
        witnesses = tuple(edge for edge in edges7 if not edge & body)
        require((len(presentations), len(witnesses)) ==
                (expected_multiplicity, expected_witnesses),
                ("size-five body crossing changed", base.labels_of(body)))
        print(
            "SIZE5_RESIDUAL_BODY", base.labels_of(body),
            "PRESENTATIONS", len(presentations),
            "E7_WITNESSES", len(witnesses),
            "FIRST_WITNESS", base.labels_of(witnesses[0]),
        )

    print(
        "SIZE5_SUMMARY", "NO_E5_COVER_ROWS", no_cover_rows,
        "E5_COVER_ROWS", rows_with_covers,
        "COMPLETE_E5_COVERS", sum(total_histogram.values()),
        "COVER_SIZE_HIST", tuple(sorted(total_histogram.items())),
        "E6_BLOCKER_PRESENTATIONS", len(residual_presentations),
        "E6_BLOCKER_BODIES", len(bodies),
        "ROW_LEDGER_FNV1A64_LE", f"{digest:016x}",
    )
    return residual_presentations, digest


def audit_size_six(base, anchors, edges5, edges6):
    total_histogram = Counter()
    no_cover_rows = 0
    rows_with_covers = 0
    all_presentations = []
    digest = 0xCBF29CE484222325

    for anchor in anchors:
        anchor_mask = base.mask_of(anchor)
        vertices = tuple(index for index in range(30) if not anchor_mask >> index & 1)
        active5 = tuple(edge for edge in edges5 if not edge & anchor_mask)
        covers5 = complete_covers(active5, vertices, 4)
        total_histogram.update(cover.bit_count() for cover in covers5)
        if covers5:
            rows_with_covers += 1
        else:
            no_cover_rows += 1
        active6 = tuple(edge for edge in edges6 if not edge & anchor_mask)
        covers6 = filter_covers(covers5, active6, vertices)
        require(not covers6, ("size-six E6 blocker survived", anchor, covers6))
        all_presentations.extend((anchor, cover) for cover in covers5)
        digest = update_row_digest(digest, anchor_mask, len(active5), covers5, covers6)
        print(
            "SIZE6_ROW", "ANCHOR", anchor, "E5_EDGES", len(active5),
            "E5_COVERS_THROUGH_4", len(covers5),
            "COVERS", tuple(base.labels_of(cover) for cover in covers5),
            "E6_COVERS_THROUGH_4", len(covers6),
        )

    require((no_cover_rows, rows_with_covers) == (22, 2),
            ("size-six 22/2 split changed", no_cover_rows, rows_with_covers))
    require(dict(total_histogram) == {4: 4},
            ("size-six cover histogram changed", total_histogram))
    expected_rows = {
        (42, 63, 88, 132, 286, 290): 3,
        (42, 63, 132, 176, 286, 290): 1,
    }
    require(Counter(anchor for anchor, _ in all_presentations) == expected_rows,
            ("exceptional size-six rows changed", all_presentations))
    for anchor, cover in all_presentations:
        body = base.mask_of(anchor) | cover
        witnesses = tuple(edge for edge in edges6 if not edge & body)
        require(witnesses, ("size-six cover did not cross E6", anchor, cover))
        print(
            "SIZE6_E5_BLOCKER", "ANCHOR", anchor, "K", base.labels_of(cover),
            "E6_WITNESSES", len(witnesses),
            "FIRST_WITNESS", base.labels_of(witnesses[0]),
        )
    print(
        "SIZE6_SUMMARY", "NO_E5_COVER_ROWS", no_cover_rows,
        "E5_COVER_ROWS", rows_with_covers,
        "COMPLETE_E5_COVERS", sum(total_histogram.values()),
        "COVER_SIZE_HIST", tuple(sorted(total_histogram.items())),
        "E6_BLOCKERS", 0,
        "ROW_LEDGER_FNV1A64_LE", f"{digest:016x}",
    )
    return all_presentations, digest


def body_hierarchy(base, ground, anchor_layers):
    label_bit = {label: 1 << index for index, label in enumerate(base.POOL)}
    covered = set()
    union_sizes = {}
    presentation_counts = {}
    for size in (4, 5, 6):
        presentations = 0
        for anchor in anchor_layers[size]:
            anchor_mask = base.mask_of(anchor)
            optional = tuple(label for label in ground if label not in anchor)
            for choice in combinations(optional, 10 - size):
                body = anchor_mask
                for label in choice:
                    body |= label_bit[label]
                covered.add(body)
                presentations += 1
        presentation_counts[size] = presentations
        union_sizes[size] = len(covered)
    require(presentation_counts == {4: 3_230_304, 5: 7_821_198, 6: 143_640},
            ("hierarchy presentation counts changed", presentation_counts))
    require(union_sizes == {4: 1_138_494, 5: 1_487_206, 6: 1_491_665},
            ("hierarchy union sizes changed", union_sizes))

    modulus_bits = {
        label: sum(1 << (modulus - 2) for modulus in range(2, 15)
                   if label % modulus == 0)
        for label in ground
    }
    full_modulus_mask = (1 << 13) - 1
    valid_count = 0
    for body in combinations(ground, 10):
        divisor_mask = 0
        common_gcd = 0
        for label in body:
            divisor_mask |= modulus_bits[label]
            common_gcd = gcd(common_gcd, label)
        valid_count += divisor_mask == full_modulus_mask and common_gcd == 1
    require(valid_count == len(covered) == 1_491_665,
            ("ambient valid-body count or hierarchy coverage changed", valid_count, len(covered)))
    print(
        "BODY_HIERARCHY", "AMBIENT", comb(27, 10),
        "VALID_PRIMITIVE_DIVISOR_COMPLETE", valid_count,
        "SIZE4_UNION", union_sizes[4], "AFTER_SIZE5", union_sizes[5],
        "AFTER_SIZE6", union_sizes[6],
        "REMAINDER_AFTER_SIZE4", valid_count - union_sizes[4],
        "REMAINDER_AFTER_SIZE5", valid_count - union_sizes[5],
        "REMAINDER_AFTER_SIZE6", valid_count - union_sizes[6],
        "PRESENTATIONS_BY_ANCHOR_SIZE", tuple(sorted(presentation_counts.items())),
    )
    return valid_count


def main():
    four = load_four()
    base = four.load_base()
    _, atom_mass, _ = base.build_full_failure_atoms()
    ground = tuple(label for label in base.POOL if label not in four.ORIGINAL)
    anchor_layers = enumerate_minimal_anchors(base, ground)

    print("AUDIT q50_zero_original_minimal_anchor_hierarchy_python")
    print("FOUR_ANCHOR_SCRIPT_SHA256", FOUR_SHA256)
    for size in range(1, 7):
        print(
            "MINIMAL_ANCHOR_LAYER", size, "COUNT", len(anchor_layers[size]),
            "MASK_SHA256", EXPECTED_ANCHOR_SHA256[size],
        )
    print("SIZE6_PRODUCT", "CORE", (42, 63, 132, 286),
          "LEFT", (8, 16, 88, 176), "RIGHT", (10, 20, 30, 170, 190, 290))

    layers = {}
    for arity in (5, 6, 7):
        layer = four.global_layer(base, atom_mass, arity)
        edges, equalities, min_edge, max_nonedge, mask_hash = layer
        expected_count, expected_hash = EXPECTED_LAYER_DATA[arity]
        require((len(edges), equalities, mask_hash) ==
                (expected_count, 0, expected_hash),
                ("global layer changed", arity, len(edges), equalities, mask_hash))
        layers[arity] = edges
        print(
            "GLOBAL_LAYER", arity, "CANDIDATES", comb(30, arity),
            "EDGES", len(edges), "EQUALITIES", equalities,
            "MIN_EDGE_DELTA", min_edge, "MAX_NONEDGE_DELTA", max_nonedge,
            "EDGE_MASK_SHA256", mask_hash,
        )

    audit_size_five(base, anchor_layers[5], layers[5], layers[6], layers[7])
    audit_size_six(base, anchor_layers[6], layers[5], layers[6])
    valid_count = body_hierarchy(base, ground, anchor_layers)
    cumulative = 888_030 + 6_660_225 + 1_071_961 + valid_count
    require(cumulative == 10_111_881, "complete q50 four-slice count changed")
    print("Q50_FOUR_DISJOINT_ORIGINAL_ANCHOR_SLICES", cumulative)
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
