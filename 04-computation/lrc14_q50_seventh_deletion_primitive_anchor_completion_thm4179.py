#!/usr/bin/env python3
"""Exact blocker-descent audit for THM-4179.

This path reuses THM-4178's frozen pool-wall failure atoms, but changes the
combinatorics: it enumerates every seven-cover of each hostile depth-six
repair hypergraph with Python big-integer incidence bitsets, then exhibits a
depth-seven repair missed by every such blocker.  An explicit eight-cover
proves the new transversal number is exactly eight.
"""

from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
from itertools import combinations
from math import comb
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178.py"
BASE_SHA256 = "ff9748f07cbc1bab20832860ff5f52bca59573262725a4e69c5b42ab2d078251"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_base():
    require(hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
            "THM-4178 pool-wall dependency hash changed")
    spec = importlib.util.spec_from_file_location("thm4178_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-4178 dependency")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def local_layer(base, atom_mass, anchor, arity):
    anchor_mask = base.mask_of(anchor)
    local_to_pool = tuple(
        index for index in range(len(base.POOL)) if not (anchor_mask >> index & 1)
    )
    require(len(local_to_pool) == 27, "optional ground set changed")

    edges = []
    equalities = 0
    minimum_edge_delta = None
    maximum_nonedge_delta = None
    for local_vertices in combinations(range(27), arity):
        local_mask = sum(1 << vertex for vertex in local_vertices)
        global_mask = sum(1 << local_to_pool[vertex] for vertex in local_vertices)
        numerator = base.deletion_numerator(global_mask, atom_mass)
        delta = 9 * numerator - 8 * base.Q * base.COMMON
        equalities += delta == 0
        if delta >= 0:
            edges.append(local_mask)
            minimum_edge_delta = (
                delta if minimum_edge_delta is None else min(minimum_edge_delta, delta)
            )
        else:
            maximum_nonedge_delta = (
                delta if maximum_nonedge_delta is None else max(maximum_nonedge_delta, delta)
            )
    require(len(edges) > 0, ("empty repair layer", anchor, arity))
    return local_to_pool, tuple(edges), equalities, minimum_edge_delta, maximum_nonedge_delta


def label_tuple(base, local_to_pool, mask):
    return tuple(
        base.POOL[local_to_pool[vertex]] for vertex in range(27) if mask >> vertex & 1
    )


def mask_from_labels(base, local_to_pool, labels):
    wanted = set(labels)
    mask = sum(
        1 << vertex
        for vertex, pool_index in enumerate(local_to_pool)
        if base.POOL[pool_index] in wanted
    )
    require(mask.bit_count() == len(wanted), ("labels outside optional ground set", labels))
    return mask


def enumerate_covers_through_budget(edges, vertex_count, budget):
    """Enumerate every first-hit cover of size at most ``budget`` exactly.

    Each edge index is one bit of a Python integer.  At a non-cover state the
    first uncovered edge must contribute a vertex to every completion, so
    branching on all its vertices is exhaustive.  Covered incidence is a
    deterministic function of the chosen vertex mask; per-depth deduplication
    therefore loses no completion.
    """
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
    require(all(all(edge & cover for edge in edges) for cover in covers),
            "enumerated an invalid cover")
    return tuple(covers), nodes, tuple(len(level) for level in seen)


def blocker_body_mass(base, atom_mass, anchor, local_to_pool, blocker):
    kept_pool_mask = base.mask_of(anchor)
    for vertex in range(27):
        if blocker >> vertex & 1:
            kept_pool_mask |= 1 << local_to_pool[vertex]
    numerator = sum(
        weight for failure, weight in atom_mass.items() if not failure & kept_pool_mask
    )
    denominator = 14 * base.Q * base.COMMON
    return numerator, denominator, Fraction(numerator, denominator)


def audit_body_union(base):
    original = set(base.ORIGINAL_ANCHORS)
    nonoriginal = tuple(label for label in base.POOL if label not in original)
    require(len(nonoriginal) == 27 and 252 in nonoriginal, "body-count ground set changed")
    alternatives = {40, 80, 240}
    others = tuple(label for label in nonoriginal if label != 252)
    histogram = Counter()
    for tail in combinations(others, 8):
        multiplicity = len(alternatives.intersection(tail))
        if multiplicity:
            histogram[multiplicity] += 1

    expected = {1: 3 * comb(23, 7), 2: 3 * comb(23, 6), 3: comb(23, 5)}
    require(dict(histogram) == expected, ("body presentation histogram changed", histogram))
    union_count = sum(histogram.values())
    presentation_count = sum(key * value for key, value in histogram.items())
    require(union_count == comb(26, 8) - comb(23, 8) == 1_071_961,
            "body union formula failed")
    require(presentation_count == 3 * comb(25, 7) == 1_442_100,
            "anchor-presentation count failed")
    require(union_count - comb(25, 7) == 591_261, "increment count failed")
    return histogram, union_count, presentation_count


def main():
    base = load_base()
    _, atom_mass, _ = base.build_full_failure_atoms()
    blocked_rows = (
        ((80, 143, 252), (8, 85, 88, 95, 145, 168, 193, 240), 298_279),
        ((143, 240, 252), (8, 80, 85, 88, 95, 145, 168, 193), 286_291),
    )

    print("AUDIT thm4179_q50_blocker_descent_python_bigint")
    print("BASE_SCRIPT_SHA256", BASE_SHA256)
    for anchor, eight_cover_labels, expected_d7_edges in blocked_rows:
        local_to_pool, edges6, equalities6, min6, max6 = local_layer(
            base, atom_mass, anchor, 6
        )
        blockers, nodes, levels = enumerate_covers_through_budget(edges6, 27, 7)
        require(equalities6 == 0, ("depth-six equality", anchor))
        require(min(blocker.bit_count() for blocker in blockers) == 7,
                ("depth-six transversal below seven", anchor))
        expected_blocker_count = 1 if anchor == (80, 143, 252) else 3
        require(len(blockers) == expected_blocker_count,
                ("depth-six blocker count changed", anchor, len(blockers)))

        local_to_pool7, edges7, equalities7, min7, max7 = local_layer(
            base, atom_mass, anchor, 7
        )
        require(local_to_pool7 == local_to_pool, "optional coordinate changed across arities")
        require(len(edges7) == expected_d7_edges and equalities7 == 0,
                ("depth-seven ledger changed", anchor, len(edges7), equalities7))

        eight_cover = mask_from_labels(base, local_to_pool, eight_cover_labels)
        require(eight_cover.bit_count() == 8 and all(edge & eight_cover for edge in edges7),
                ("invalid depth-seven eight-cover", anchor))

        print(
            "ROW", "ANCHOR", anchor,
            "D6_EDGES", len(edges6), "D6_EQUALITIES", equalities6,
            "D6_MIN_EDGE_DELTA", min6, "D6_MAX_NONEDGE_DELTA", max6,
            "D6_EXACT_SEVEN_BLOCKERS", len(blockers),
            "ENUM_NODES", nodes, "ENUM_LEVELS", levels,
        )
        for blocker in blockers:
            witness = next((edge for edge in edges7 if not edge & blocker), None)
            require(witness is not None, ("depth-six blocker survived depth seven", anchor))
            numerator, denominator, body_mass = blocker_body_mass(
                base, atom_mass, anchor, local_to_pool, blocker
            )
            require(body_mass > Fraction(4, 63), ("blocker body below Haar threshold", anchor))
            print(
                "BLOCKER", "ANCHOR", anchor,
                "K", label_tuple(base, local_to_pool, blocker),
                "D7_MISSED_REPAIR", label_tuple(base, local_to_pool, witness),
                "BODY_MASS_NUM_DEN", (numerator, denominator),
                "BODY_MASS_REDUCED", body_mass,
                "BODY_MARGIN", body_mass - Fraction(4, 63),
            )
        print(
            "D7", "ANCHOR", anchor, "CANDIDATES", comb(27, 7),
            "EDGES", len(edges7), "EQUALITIES", equalities7,
            "MIN_EDGE_DELTA", min7, "MAX_NONEDGE_DELTA", max7,
            "TAU", 8, "EIGHT_COVER", label_tuple(base, local_to_pool, eight_cover),
        )

    histogram, union_count, presentation_count = audit_body_union(base)
    print(
        "EXACTLY_ONE_ORIGINAL_ANCHOR_BODY_UNION",
        "PRESENTATION_MULTIPLICITY_HIST", tuple(sorted(histogram.items())),
        "DISTINCT_BODIES", union_count,
        "ANCHOR_PRESENTATIONS", presentation_count,
        "PAIR_INTERSECTION", comb(24, 6),
        "TRIPLE_INTERSECTION", comb(23, 5),
        "INCREMENT_BEYOND_40_ANCHOR", union_count - comb(25, 7),
    )
    print("MAXIMAL_DELETION_DUALITY_CONTROLS", 4, "ALL_STRICTLY_ABOVE_4_OVER_63", True)
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
