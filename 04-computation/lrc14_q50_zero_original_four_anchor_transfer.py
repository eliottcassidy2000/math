#!/usr/bin/env python3
"""Exact q=50 four-anchor repair-incidence audit.

This audit starts from the frozen THM-4178 labelled pool-wall atoms, but its
target and combinatorics are new.  It classifies primitive divisor-complete
four-anchors avoiding the three original anchors, searches their depth-six
repair hypergraphs through transversal budget six, carries every surviving
blocker across depth seven, and performs a literal zero-original body-union
census.
"""

from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import importlib.util
from itertools import combinations
from math import comb, gcd
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178.py"
BASE_SHA256 = "ff9748f07cbc1bab20832860ff5f52bca59573262725a4e69c5b42ab2d078251"

ORIGINAL = frozenset((120, 126, 143))
X = frozenset((85, 95, 145, 193))

EXPECTED_ANCHORS = (
    (10, 63, 168, 286),
    (15, 40, 252, 286),
    (15, 80, 252, 286),
    (15, 240, 252, 286),
    (20, 63, 168, 286),
    (30, 63, 168, 286),
    (40, 63, 84, 286),
    (40, 63, 168, 286),
    (40, 63, 252, 286),
    (40, 85, 252, 286),
    (40, 95, 252, 286),
    (40, 145, 252, 286),
    (40, 193, 252, 286),
    (42, 63, 240, 286),
    (60, 63, 168, 286),
    (63, 80, 84, 286),
    (63, 80, 168, 286),
    (63, 80, 252, 286),
    (63, 84, 240, 286),
    (63, 168, 170, 286),
    (63, 168, 190, 286),
    (63, 168, 240, 286),
    (63, 168, 286, 290),
    (63, 240, 252, 286),
    (80, 85, 252, 286),
    (80, 95, 252, 286),
    (80, 145, 252, 286),
    (80, 193, 252, 286),
    (85, 240, 252, 286),
    (95, 240, 252, 286),
    (145, 240, 252, 286),
    (193, 240, 252, 286),
)

EXPECTED_BLOCKER_COUNTS = {
    (80, 85, 252, 286): 1,
    (80, 95, 252, 286): 3,
    (80, 145, 252, 286): 3,
    (80, 193, 252, 286): 3,
    (85, 240, 252, 286): 2,
    (95, 240, 252, 286): 5,
    (145, 240, 252, 286): 5,
    (193, 240, 252, 286): 5,
}

EXPECTED_D7_COVERS = {
    (80, 85, 252, 286): (8, 88, 95, 145, 168, 193, 240),
    (80, 95, 252, 286): (8, 85, 88, 145, 168, 193, 240),
    (80, 145, 252, 286): (8, 85, 88, 95, 168, 193, 240),
    (80, 193, 252, 286): (8, 85, 88, 95, 145, 168, 240),
    (85, 240, 252, 286): (8, 80, 88, 95, 145, 168, 193),
    (95, 240, 252, 286): (8, 80, 85, 88, 145, 168, 193),
    (145, 240, 252, 286): (8, 80, 85, 88, 95, 168, 193),
    (193, 240, 252, 286): (8, 80, 85, 88, 95, 145, 168),
}

BODY_SPINE = frozenset((88, 95, 145, 168, 193, 240, 252, 286))
BODY_SELECTOR_PAIRS = ((8, 80), (16, 85), (16, 170), (80, 85), (80, 170))
EXPECTED_BODY_PRESENTATIONS = (6, 4, 3, 8, 6)
EXPECTED_BODY_D7_WITNESSES = (243, 439, 736, 278, 514)
EXPECTED_PATTERN_COUNTS = {
    0b11110: 27,
    0b10101: 5,
    0b00110: 284,
    0b01010: 11,
    0b10100: 91,
    0b11000: 172,
    0b00001: 238,
    0b00010: 117,
    0b00100: 329,
    0b01000: 68,
    0b10000: 219,
}
EXPECTED_UNION_HISTOGRAM = {
    1: 262_761,
    2: 332_137,
    3: 204_900,
    4: 178_178,
    5: 59_560,
    6: 66_399,
    7: 13_188,
    8: 11_229,
    9: 6_149,
    10: 2_666,
    11: 658,
    12: 493,
    13: 140,
    15: 36,
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_base():
    require(
        hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
        "THM-4178 pool-wall dependency hash changed",
    )
    spec = importlib.util.spec_from_file_location("thm4178_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-4178")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def divisor_complete(anchor):
    return all(any(label % modulus == 0 for label in anchor) for modulus in range(2, 15))


def enumerate_four_anchors(base):
    zero_original_pool = tuple(label for label in base.POOL if label not in ORIGINAL)
    anchors = tuple(
        anchor
        for anchor in combinations(zero_original_pool, 4)
        if gcd(*anchor) == 1 and divisor_complete(anchor)
    )
    require(anchors == EXPECTED_ANCHORS, ("four-anchor classification changed", anchors))
    require(all(286 in anchor for anchor in anchors), "13-pin 286 was not forced")
    return anchors, zero_original_pool


def global_layer(base, atom_mass, arity):
    edges = []
    equalities = 0
    minimum_edge_delta = None
    maximum_nonedge_delta = None
    edge_hash = hashlib.sha256()
    for vertices in combinations(range(len(base.POOL)), arity):
        mask = sum(1 << vertex for vertex in vertices)
        numerator = base.deletion_numerator(mask, atom_mass)
        delta = 9 * numerator - 8 * base.Q * base.COMMON
        equalities += delta == 0
        if delta >= 0:
            edges.append(mask)
            minimum_edge_delta = (
                delta if minimum_edge_delta is None else min(minimum_edge_delta, delta)
            )
            edge_hash.update(mask.to_bytes(4, "little"))
        else:
            maximum_nonedge_delta = (
                delta if maximum_nonedge_delta is None else max(maximum_nonedge_delta, delta)
            )
    require(edges, ("empty global repair layer", arity))
    return (
        tuple(edges),
        equalities,
        minimum_edge_delta,
        maximum_nonedge_delta,
        edge_hash.hexdigest(),
    )


def enumerate_covers_through_budget(edges, vertex_count, budget):
    """Enumerate every inclusion-first cover through ``budget`` exactly."""
    incidence = [0] * vertex_count
    for edge_index, edge in enumerate(edges):
        vertices = edge
        while vertices:
            bit = vertices & -vertices
            incidence[bit.bit_length() - 1] |= 1 << edge_index
            vertices ^= bit

    all_edges = (1 << len(edges)) - 1
    seen = [set() for _ in range(budget + 1)]
    seen[0].add(0)
    covers = []
    nodes = 0

    def search(chosen, covered):
        nonlocal nodes
        nodes += 1
        depth = chosen.bit_count()
        if covered == all_edges:
            covers.append(chosen)
            return
        if depth == budget:
            return
        uncovered = all_edges ^ covered
        first_index = (uncovered & -uncovered).bit_length() - 1
        branch = edges[first_index]
        while branch:
            bit = branch & -branch
            child = chosen | bit
            if child not in seen[depth + 1]:
                seen[depth + 1].add(child)
                search(child, covered | incidence[bit.bit_length() - 1])
            branch ^= bit

    search(0, 0)
    require(
        all(all(edge & cover for edge in edges) for cover in covers),
        "cover enumeration returned a non-cover",
    )
    return tuple(covers), nodes, tuple(len(level) for level in seen)


def zero_original_body_union(base, anchors, zero_original_pool):
    multiplicity = Counter()
    for anchor in anchors:
        anchor_mask = base.mask_of(anchor)
        optional = tuple(label for label in zero_original_pool if label not in anchor)
        require(len(optional) == 23, ("zero-original optional size", anchor))
        for choice in combinations(optional, 6):
            multiplicity[anchor_mask | base.mask_of(choice)] += 1
    histogram = Counter(multiplicity.values())
    require(dict(sorted(histogram.items())) == EXPECTED_UNION_HISTOGRAM,
            ("body multiplicity histogram changed", histogram))
    require(len(multiplicity) == 1_138_494, "zero-original body union changed")
    require(sum(multiplicity.values()) == 32 * comb(23, 6) == 3_230_304,
            "zero-original presentation count changed")
    return histogram, len(multiplicity), sum(multiplicity.values())


def main():
    base = load_base()
    _, atom_mass, _ = base.build_full_failure_atoms()
    anchors, zero_original_pool = enumerate_four_anchors(base)

    print("AUDIT q50_zero_original_four_anchor_repair_incidence_python")
    print("BASE_SCRIPT_SHA256", BASE_SHA256)
    print("FOUR_ANCHORS", len(anchors), "ALL_CONTAIN_286", True)
    for anchor in anchors:
        print("ANCHOR", anchor)

    layers = {}
    for arity, expected_edges in ((6, 85_324), (7, 821_737)):
        layer = global_layer(base, atom_mass, arity)
        edges, equalities, min_edge, max_nonedge, edge_hash = layer
        require(len(edges) == expected_edges and equalities == 0,
                ("global repair layer changed", arity, len(edges), equalities))
        layers[arity] = layer
        print(
            "GLOBAL_LAYER", arity, "CANDIDATES", comb(30, arity),
            "EDGES", len(edges), "EQUALITIES", equalities,
            "MIN_EDGE_DELTA", min_edge, "MAX_NONEDGE_DELTA", max_nonedge,
            "EDGE_MASK_SHA256", edge_hash,
        )

    edges6 = layers[6][0]
    edges7 = layers[7][0]
    edge7_set = frozenset(edges7)
    blocker_presentations = []
    positive_d6 = 0

    for anchor in anchors:
        anchor_mask = base.mask_of(anchor)
        active6 = tuple(edge for edge in edges6 if not edge & anchor_mask)
        blockers, nodes, levels = enumerate_covers_through_budget(active6, 30, 6)
        require(all(not blocker & anchor_mask for blocker in blockers),
                ("blocker used anchor", anchor))
        if anchor in EXPECTED_BLOCKER_COUNTS:
            require(len(blockers) == EXPECTED_BLOCKER_COUNTS[anchor],
                    ("depth-six blocker count changed", anchor, len(blockers)))
            require(all(blocker.bit_count() == 6 for blocker in blockers),
                    ("depth-six cover below six", anchor))
        else:
            require(not blockers, ("unexpected depth-six blocker", anchor, blockers))
            positive_d6 += 1
        print(
            "D6_ROW", "ANCHOR", anchor, "EDGES", len(active6),
            "BLOCKERS_THROUGH_6", len(blockers), "ENUM_NODES", nodes,
            "ENUM_LEVELS", levels,
        )
        for blocker in blockers:
            print("D6_BLOCKER", "ANCHOR", anchor, "K", base.labels_of(blocker))
            blocker_presentations.append((anchor, anchor_mask, blocker))

    require(positive_d6 == 24, ("depth-six positive row count", positive_d6))
    require(len(blocker_presentations) == 27, "depth-six blocker presentation count")

    for anchor in EXPECTED_BLOCKER_COUNTS:
        anchor_mask = base.mask_of(anchor)
        active7 = tuple(edge for edge in edges7 if not edge & anchor_mask)
        cover = base.mask_of(EXPECTED_D7_COVERS[anchor])
        require(not cover & anchor_mask and cover.bit_count() == 7,
                ("invalid declared depth-seven cover", anchor))
        require(all(edge & cover for edge in active7),
                ("declared depth-seven cover misses a repair", anchor))
        print(
            "D7_ROW", "ANCHOR", anchor, "EDGES", len(active7),
            "TAU", 7, "COVER", base.labels_of(cover),
            "NO_COVER_THROUGH_6_BY_BLOCKER_DESCENT", True,
        )

    body_presentations = defaultdict(list)
    for anchor, anchor_mask, blocker in blocker_presentations:
        body = anchor_mask | blocker
        witness = next((edge for edge in edges7 if not edge & body), None)
        require(witness is not None, ("uncrossed depth-six blocker", anchor, blocker))
        body_presentations[body].append((anchor, blocker))
    require(len(body_presentations) == 5, "blocker presentations did not collapse to five bodies")

    expected_bodies = tuple(base.mask_of(BODY_SPINE | frozenset(pair))
                            for pair in BODY_SELECTOR_PAIRS)
    require(set(body_presentations) == set(expected_bodies),
            ("five-state obstruction graph changed", tuple(map(base.labels_of, body_presentations))))
    for index, body in enumerate(expected_bodies):
        presentations = body_presentations[body]
        witness_count = sum(not edge & body for edge in edges7)
        require(len(presentations) == EXPECTED_BODY_PRESENTATIONS[index],
                ("body presentation multiplicity", index, len(presentations)))
        require(witness_count == EXPECTED_BODY_D7_WITNESSES[index],
                ("body witness count", index, witness_count))
        first_witness = next(edge for edge in edges7 if not edge & body)
        print(
            "BLOCKER_BODY", index, "SELECTOR_PAIR", BODY_SELECTOR_PAIRS[index],
            "LABELS", base.labels_of(body), "PRESENTATIONS", len(presentations),
            "D7_WITNESSES", witness_count,
            "FIRST_WITNESS", base.labels_of(first_witness),
        )

    pattern_counts = Counter()
    pattern_examples = {}
    for edge in edges7:
        pattern = sum(1 << index for index, body in enumerate(expected_bodies)
                      if not edge & body)
        if pattern:
            pattern_counts[pattern] += 1
            pattern_examples.setdefault(pattern, edge)
    require(dict(pattern_counts) == EXPECTED_PATTERN_COUNTS,
            ("body-witness incidence pattern changed", pattern_counts))
    require(0b11111 not in pattern_counts, "a single repair unexpectedly covers all five bodies")
    for pattern in sorted(pattern_counts, key=lambda value: (-value.bit_count(), value)):
        states = tuple(index for index in range(5) if pattern >> index & 1)
        print(
            "WITNESS_PATTERN", states, "COUNT", pattern_counts[pattern],
            "EXAMPLE", base.labels_of(pattern_examples[pattern]),
        )

    bank = (
        (8, 10, 15, 20, 143, 176, 290),
        (10, 15, 20, 85, 120, 143, 176),
    )
    bank_patterns = []
    for repair in bank:
        mask = base.mask_of(repair)
        require(mask in edge7_set, ("bank member is not a repair", repair))
        pattern = sum(1 << index for index, body in enumerate(expected_bodies)
                      if not mask & body)
        bank_patterns.append(pattern)
        numerator = base.deletion_numerator(mask, atom_mass)
        mass = Fraction(numerator, 14 * base.Q * base.COMMON)
        require(mass > Fraction(4, 63), ("bank member not strict", repair, mass))
        print(
            "BANK_REPAIR", repair, "BODY_PATTERN",
            tuple(index for index in range(5) if pattern >> index & 1),
            "MASS", mass, "MARGIN", mass - Fraction(4, 63),
        )
    require(bank_patterns[0] | bank_patterns[1] == 0b11111,
            "two-repair bank does not cover all obstruction bodies")
    print("COMMON_ONE_REPAIR_BANK", False, "EXPLICIT_TWO_REPAIR_BANK", True)

    histogram, distinct_bodies, presentations = zero_original_body_union(
        base, anchors, zero_original_pool
    )
    print(
        "ZERO_ORIGINAL_BODY_UNION", "PRESENTATIONS", presentations,
        "DISTINCT_BODIES", distinct_bodies,
        "MULTIPLICITY_HISTOGRAM", tuple(sorted(histogram.items())),
    )
    cumulative = 888_030 + 6_660_225 + 1_071_961 + distinct_bodies
    require(cumulative == 9_758_710, "four-slice cumulative count changed")
    print("Q50_FOUR_DISJOINT_ORIGINAL_ANCHOR_SLICES", cumulative)
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
