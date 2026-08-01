#!/usr/bin/env python3
r"""Uniform closure of the complete 21-body robust-edge-8 block.

The three sum-threshold creation codes are ``DIIDI`` (13 bodies), ``IDDII``
(4 bodies), and ``IDIDD`` (4 bodies).  The last two have connected low graphs
for every realized high subgraph.  A ``DIIDI`` low graph disconnects exactly
when its five star edges are high; the star center is then a free projective
scale and the other five vertices induce one of ``K5``, ``K5-e``,
``K5-e`` (the other labelled edge), or ``K5-P3``.

The edge-9 clique-anchor chart is complete whenever every vertex outside a
maximum low clique is directly low-adjacent to that clique.  Edge 8 contains
the first connected graph outside this chart: its last vertex attaches at
low distance two.  The repair is general.  Normalize a maximum low clique and
grow through the connected low graph.  At every stage, choose an unassigned
vertex with a low edge to the assigned set.  Its value is its chosen parent's
value times one of the finite symmetric low ratios.  Exact pruning against all
assigned vertices then enumerates every labelled projective word.  This staged
low-spanning-tree chart contributes 2,624 words for the unique missed graph.

For each disconnected case, choose a pair in the five-level component and
bound the free center's debt by its level-one value.  Component dilation grows
the pair gcd and decreases component debt.  Thus all four cases are uniform
cylinders, not finite center checks.

The exact census has 652,688 rows.  Positive exact margins close all 21 bodies
for arbitrary positive reflected levels.  The cumulative reflected THM-2941
body count becomes 2,442/3,003, leaving 561 bodies with at most seven robust
edges.  This does not prove physical LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EDGE9 = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_edge9_threshold_block_uniform_closure_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_edge8_threshold_block_uniform_closure_thm2941.out"
)
EXPECTED_EDGE9_SHA256 = "8540d1689ca002f284e9dcfb1730c382bb8eb6a79b62d09084f53707629ab621"
EXPECTED_SEMANTIC_SHA256 = "f197d1a7f237526b03c3ca97a47fb48900d560d9a6bd75838140c74d0c59a04b"
EXPECTED_BODY_DIGEST = "17c120a9fb01b2e30cf42debad83dc80433f65d4eeba496414b537543c1b1641"
EXPECTED_ASSIGNMENT_DIGEST = "76c30b33786c99b75b931ad688d9db2840c7ec7466550a727fca39fa65a41adb"

BODY_COUNT = 3003
PRIOR_CLOSED_COUNT = 2421
ROBUST_EDGE8_BODY_COUNT = 21
CUMULATIVE_CLOSED_COUNT = 2442
REMAINING_BODY_COUNT = 561
EXPECTED_CODE_HIST = (("DIIDI", 13), ("IDDII", 4), ("IDIDD", 4))
EXPECTED_CODE_ROW_COUNTS = (
    ("DIIDI", 430352),
    ("IDDII", 119808),
    ("IDIDD", 102528),
)
CONNECTED_ASSIGNMENT_COUNT = 641144
CYLINDER_TEMPLATE_COUNT = 11544
CERTIFICATE_ROW_COUNT = 652688
K5_COMPONENT_WORD_COUNT = 240
K5E_COMPONENT_WORD_COUNT = 168
K5P3_COMPONENT_WORD_COUNT = 312
STAGED_FALLBACK_WORD_COUNT = 2624
DIRECT_CONTROL_SCALES = (1, 2, 5)
DIRECT_CONTROL_COUNT = 63
EXPECTED_WEAKEST = (
    F(424512216668439000719269, 7282605102107128237009800),
    "IDDII",
    "(5, (3, 3, 2, 1, 1, 0), 1)",
    (1, 3, 4, 9, 12, 13),
    ((1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)),
    ((1, 4), (2, 4), (2, 5), (3, 5), (4, 5)),
    (12, 3, 2, 4, 16, 9),
    F(409, 6552),
    2,
    3,
    1,
    2,
    2,
    5966,
    1,
    586,
    "low",
)

DIIDI_FULL = (
    (0, 5), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5),
)
DIIDI_STAR5 = ((0, 5), (1, 5), (2, 5), (3, 5), (4, 5))
IDIDD_FULL = (
    (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5),
)
IDDII_FULL = (
    (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5),
)
STAGED_FALLBACK_HIGH = (
    (0, 5), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


require(sha256(EDGE9) == EXPECTED_EDGE9_SHA256,
        ("edge-9 theorem changed", sha256(EDGE9), EXPECTED_EDGE9_SHA256))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E9 = import_module("robust_edge8_edge9", EDGE9)
E10 = E9.E10
R = E9.R
LOW = E9.LOW
U = E9.U
E11 = E9.E11
E12 = E9.E12
SYMMETRIC_LOW_RATIOS = E9.SYMMETRIC_LOW_RATIOS


def maximum_low_clique(vertex_count: int, high_edges) -> tuple[int, ...]:
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    for size in range(vertex_count, 0, -1):
        for core in combinations(range(vertex_count), size):
            if all(edge not in high_set for edge in combinations(core, 2)):
                return core
    raise RuntimeError(("no low clique", vertex_count, high_set))


def staged_classified_words(vertex_count: int, high_edges) -> tuple[tuple[int, ...], ...]:
    """Enumerate a connected labelled low graph by a low spanning tree."""
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    core = maximum_low_clique(vertex_count, high_set)
    words = set()
    for clique in E9.clique_types(len(core)):
        for clique_values in permutations(clique):
            partial = [None] * vertex_count
            for index, value in zip(core, clique_values):
                partial[index] = value

            def extend() -> None:
                unassigned = tuple(
                    index for index, value in enumerate(partial) if value is None
                )
                if not unassigned:
                    minimum = min(partial)
                    normalized = tuple(value / minimum for value in partial)
                    words.add(E9.primitive_word(normalized))
                    return
                eligible = tuple(
                    index
                    for index in unassigned
                    if any(
                        tuple(sorted((index, assigned))) not in high_set
                        for assigned, value in enumerate(partial)
                        if value is not None
                    )
                )
                require(eligible, ("disconnected low graph", vertex_count, high_set, partial))
                index = eligible[0]
                anchor = next(
                    assigned
                    for assigned, value in enumerate(partial)
                    if value is not None
                    and tuple(sorted((index, assigned))) not in high_set
                )
                candidates = tuple(sorted({
                    partial[anchor] * ratio for ratio in SYMMETRIC_LOW_RATIOS
                }))
                for candidate in candidates:
                    if candidate in partial:
                        continue
                    if not all(
                        (not LOW.low_adjacent(candidate, value))
                        == (tuple(sorted((index, assigned))) in high_set)
                        for assigned, value in enumerate(partial)
                        if value is not None
                    ):
                        continue
                    partial[index] = candidate
                    extend()
                    partial[index] = None

            extend()

    for word in words:
        actual_high = frozenset(
            edge
            for edge in combinations(range(vertex_count), 2)
            if not LOW.low_adjacent(F(word[edge[0]]), F(word[edge[1]]))
        )
        require(actual_high == high_set, (word, actual_high, high_set))
    return tuple(sorted(words))


def classified_words(vertex_count: int, high_edges) -> tuple[tuple[int, ...], ...]:
    """Use the fast edge-9 chart, with the complete staged chart as fallback."""
    normalized = tuple(sorted(tuple(sorted(edge)) for edge in high_edges))
    try:
        return E9.classified_words(vertex_count, normalized)
    except RuntimeError:
        require(normalized == STAGED_FALLBACK_HIGH,
                ("unexpected clique-anchor miss", vertex_count, normalized))
        words = staged_classified_words(vertex_count, normalized)
        require(len(words) == STAGED_FALLBACK_WORD_COUNT,
                (normalized, len(words), STAGED_FALLBACK_WORD_COUNT))
        return words


def low_graph_is_connected(vertex_count: int, high_edges) -> bool:
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    reached = {0}
    frontier = [0]
    while frontier:
        first = frontier.pop()
        for second in range(vertex_count):
            if second in reached or second == first:
                continue
            if tuple(sorted((first, second))) not in high_set:
                reached.add(second)
                frontier.append(second)
    return len(reached) == vertex_count


def orbit_count(words, vertex_count: int, high_edges) -> int:
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    automorphisms = tuple(
        permutation
        for permutation in permutations(range(vertex_count))
        if frozenset(
            tuple(sorted((permutation[i], permutation[j]))) for i, j in high_set
        ) == high_set
    )
    representatives = {
        min(
            tuple(word[permutation[index]] for index in range(vertex_count))
            for permutation in automorphisms
        )
        for word in words
    }
    return len(representatives)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    component_specs = (
        ("k5_cylinder", (), K5_COMPONENT_WORD_COUNT, 2),
        ("k5e_cylinder_a", ((2, 4),), K5E_COMPONENT_WORD_COUNT, 14),
        ("k5e_cylinder_b", ((3, 4),), K5E_COMPONENT_WORD_COUNT, 14),
        ("k5p3_cylinder", ((2, 4), (3, 4)), K5P3_COMPONENT_WORD_COUNT, 78),
    )
    components = []
    for category, leaf_high, expected_words, expected_orbits in component_specs:
        if category == "k5_cylinder":
            # Preserve the established edge-10 cylinder ordering in the
            # assignment digest: primitive type first, then permutations.
            words = tuple(
                word
                for primitive in E10.LOW_K5_PRIMITIVES
                for word in permutations(primitive)
            )
        else:
            words = E9.classified_words(5, leaf_high)
        require(len(words) == expected_words, (category, len(words), expected_words))
        require(orbit_count(words, 5, leaf_high) == expected_orbits,
                (category, orbit_count(words, 5, leaf_high), expected_orbits))
        components.append((category, tuple(leaf_high), words))

    bodies = []
    code_hist: Counter[str] = Counter()
    body_digest = hashlib.sha256()
    universal_exceptions = {row[0] for row in U.EXPECTED_EXCEPTIONS}
    expected_complements = {
        "DIIDI": DIIDI_FULL,
        "IDIDD": IDIDD_FULL,
        "IDDII": IDDII_FULL,
    }
    for body in combinations(range(1, 15), 6):
        ruler, universal_debt, robust_edges = LOW.robust_edges(body)
        if len(robust_edges) != 8:
            continue
        require(body not in universal_exceptions,
                ("edge-8 body lacks complete same-level graph", body))
        threshold = E12.nonrobust_threshold(body, ruler, universal_debt)
        code = E12.threshold_creation_code(body, threshold)
        nonedges = tuple(
            edge for edge in combinations(range(6), 2) if edge not in set(robust_edges)
        )
        require(code in expected_complements, (body, code, nonedges))
        require(nonedges == expected_complements[code],
                ("edge-8 complement changed", body, code, nonedges))
        code_hist[code] += 1
        bodies.append((body, ruler, universal_debt, robust_edges, nonedges, threshold, code))
        body_digest.update(
            f"{body}|{ruler}|{universal_debt}|{threshold}|{robust_edges}|{nonedges}|{code}\n".encode()
        )
    require(tuple(sorted(code_hist.items())) == EXPECTED_CODE_HIST, code_hist)
    require(len(bodies) == ROBUST_EDGE8_BODY_COUNT, len(bodies))
    require(body_digest.hexdigest() == EXPECTED_BODY_DIGEST,
            (body_digest.hexdigest(), EXPECTED_BODY_DIGEST))

    word_cache = {}
    rows = []
    code_row_counts: Counter[str] = Counter()
    category_row_counts: Counter[str] = Counter()
    certificate_kinds: Counter[str] = Counter()
    assignment_digest = hashlib.sha256()
    weakest_by_body = {}
    connected_rows = 0
    cylinder_rows = 0

    for body, ruler, _, robust_edges, nonedges, _, code in bodies:
        safe_ruler, safe_ranges = R.safe_cell_ranges(body)
        require(safe_ruler == ruler, (body, safe_ruler, ruler))
        profile_cache = {}
        robust_set = set(robust_edges)

        for size in range(1, len(nonedges) + 1):
            for high_edges in combinations(nonedges, size):
                high_edges = tuple(sorted(high_edges))
                if not low_graph_is_connected(6, high_edges):
                    require(code == "DIIDI" and set(DIIDI_STAR5) <= set(high_edges),
                            ("unexpected disconnected low graph", code, high_edges))
                    continue
                if high_edges not in word_cache:
                    word_cache[high_edges] = classified_words(6, high_edges)
                category = str(E11.graph_key(high_edges))
                for levels in word_cache[high_edges]:
                    induced_high = tuple(
                        edge for edge in combinations(range(6), 2)
                        if not LOW.low_adjacent(F(levels[edge[0]]), F(levels[edge[1]]))
                    )
                    require(induced_high == high_edges,
                            (body, high_edges, levels, induced_high))
                    require(all(
                        LOW.low_adjacent(F(levels[i]), F(levels[j]))
                        for i, j in robust_set
                    ), (body, levels, robust_edges))
                    debt = E10.singleton_debt(body, ruler, levels)
                    best = E10.best_certificate(
                        body, ruler, safe_ranges, levels, debt, profile_cache
                    )
                    require(best[0] > 0,
                            ("finite edge-8 assignment failed", code, category,
                             body, nonedges, high_edges, levels, best))
                    row = (
                        best[0], code, category, body, nonedges, high_edges, levels,
                    ) + best[1:]
                    rows.append(row)
                    connected_rows += 1
                    code_row_counts[code] += 1
                    category_row_counts[category] += 1
                    certificate_kinds[best[-1]] += 1
                    assignment_digest.update(
                        f"{code}|{category}|{body}|{ruler}|{nonedges}|{high_edges}|{levels}|{debt}|{best}\n".encode()
                    )
                    if body not in weakest_by_body or row < weakest_by_body[body]:
                        weakest_by_body[body] = row

        if code == "DIIDI":
            leaf_pairs = tuple(combinations(range(5), 2))
            for category, leaf_high, component_words in components:
                high_edges = tuple(sorted(DIIDI_STAR5 + leaf_high))
                require(not low_graph_is_connected(6, high_edges),
                        (category, high_edges))
                for leaf_levels in component_words:
                    debt = E10.free_center_debt_bound(body, ruler, leaf_levels)
                    best = E10.best_certificate(
                        body, ruler, safe_ranges, leaf_levels, debt,
                        profile_cache, pair_indices=leaf_pairs,
                    )
                    require(best[0] > 0,
                            ("edge-8 cylinder failed", category, body,
                             leaf_levels, best, debt))
                    row = (
                        best[0], code, category, body, nonedges,
                        high_edges, leaf_levels,
                    ) + best[1:]
                    rows.append(row)
                    cylinder_rows += 1
                    code_row_counts[code] += 1
                    category_row_counts[category] += 1
                    certificate_kinds[best[-1]] += 1
                    assignment_digest.update(
                        f"{code}|{category}|{body}|{ruler}|{nonedges}|{high_edges}|{leaf_levels}|{debt}|{best}\n".encode()
                    )
                    if body not in weakest_by_body or row < weakest_by_body[body]:
                        weakest_by_body[body] = row

    require(connected_rows == CONNECTED_ASSIGNMENT_COUNT, connected_rows)
    require(cylinder_rows == CYLINDER_TEMPLATE_COUNT, cylinder_rows)
    require(len(rows) == CERTIFICATE_ROW_COUNT, len(rows))
    require(tuple(sorted(code_row_counts.items())) == EXPECTED_CODE_ROW_COUNTS,
            code_row_counts)
    require(certificate_kinds == Counter({"low": CERTIFICATE_ROW_COUNT}),
            certificate_kinds)
    require(len(weakest_by_body) == ROBUST_EDGE8_BODY_COUNT, len(weakest_by_body))
    require(assignment_digest.hexdigest() == EXPECTED_ASSIGNMENT_DIGEST,
            (assignment_digest.hexdigest(), EXPECTED_ASSIGNMENT_DIGEST))
    weakest = min(rows)
    if EXPECTED_WEAKEST is not None:
        require(weakest == EXPECTED_WEAKEST, (weakest, EXPECTED_WEAKEST))

    direct_controls = 0
    ruler_by_body = {body: ruler for body, ruler, *_ in bodies}
    for body, row in sorted(weakest_by_body.items()):
        category = row[2]
        levels = row[6]
        floor, i, j, divisor, cell = row[7], row[8], row[9], row[12], row[13]
        ruler = ruler_by_body[body]
        a, b = body[i], body[j]
        p, q = levels[i], levels[j]
        require(R.body_cell_is_safe(ruler, body, cell),
                ("certificate cell unsafe", body, levels, (i, j), cell))
        for scale in DIRECT_CONTROL_SCALES:
            actual = E10.intersection_mass(
                R.reflected_level_arcs(ruler, a, scale * p, cell),
                R.reflected_level_arcs(ruler, b, scale * q, cell),
            )
            transported = floor - F(4 * (a + b), scale * divisor * ruler)
            if category.endswith("_cylinder") or "_cylinder_" in category:
                scaled_debt = E10.free_center_debt_bound(
                    body, ruler, tuple(scale * value for value in levels)
                )
            else:
                scaled_debt = E10.singleton_debt(
                    body, ruler, tuple(scale * value for value in levels)
                )
            require(actual >= transported > scaled_debt,
                    (body, category, levels, (i, j), scale,
                     actual, transported, scaled_debt, row))
            direct_controls += 1
    require(direct_controls == DIRECT_CONTROL_COUNT, direct_controls)
    require(PRIOR_CLOSED_COUNT + len(bodies) == CUMULATIVE_CLOSED_COUNT,
            (PRIOR_CLOSED_COUNT, len(bodies), CUMULATIVE_CLOSED_COUNT))
    require(BODY_COUNT - CUMULATIVE_CLOSED_COUNT == REMAINING_BODY_COUNT,
            (BODY_COUNT, CUMULATIVE_CLOSED_COUNT, REMAINING_BODY_COUNT))

    semantic_payload = (
        sha256(EDGE9),
        tuple(sorted(code_hist.items())),
        tuple(sorted(code_row_counts.items())),
        connected_rows,
        cylinder_rows,
        len(rows),
        (
            K5_COMPONENT_WORD_COUNT,
            K5E_COMPONENT_WORD_COUNT,
            K5E_COMPONENT_WORD_COUNT,
            K5P3_COMPONENT_WORD_COUNT,
        ),
        STAGED_FALLBACK_HIGH,
        STAGED_FALLBACK_WORD_COUNT,
        weakest,
        tuple(sorted(certificate_kinds.items())),
        DIRECT_CONTROL_SCALES,
        direct_controls,
        body_digest.hexdigest(),
        assignment_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 reflected robust-edge-8 threshold-block uniform closure exact referee",
        f"edge8_code_hist={tuple(sorted(code_hist.items()))};"
        f"code_row_counts={tuple(sorted(code_row_counts.items()))}",
        f"connected_projective_rows={connected_rows};"
        f"cylinder_templates={cylinder_rows};"
        f"category_row_counts={tuple(sorted(category_row_counts.items()))};"
        f"total_certificate_rows={len(rows)};"
        f"certificate_kinds={tuple(sorted(certificate_kinds.items()))}",
        f"component_word_counts={(K5_COMPONENT_WORD_COUNT, K5E_COMPONENT_WORD_COUNT, K5E_COMPONENT_WORD_COUNT, K5P3_COMPONENT_WORD_COUNT)};"
        f"staged_fallback_words={STAGED_FALLBACK_WORD_COUNT};"
        f"staged_fallback_high={STAGED_FALLBACK_HIGH}",
        f"weakest_margin={qtext(weakest[0])};code={weakest[1]};"
        f"category={weakest[2]};body={weakest[3]};"
        f"nonrobust_edges={weakest[4]};high_edges={weakest[5]};"
        f"levels={weakest[6]};floor={qtext(weakest[7])};"
        f"certificate_pair={(weakest[8], weakest[9])};"
        f"reduced_channel={(weakest[10], weakest[11])};gcd={weakest[12]};"
        f"j={weakest[13]};lift={weakest[14]};endpoint_radius={weakest[15]};"
        f"kind={weakest[16]}",
        f"direct_control_scales={DIRECT_CONTROL_SCALES};"
        f"per_body_weakest_controls={direct_controls}",
        "classifier_law=a maximum low clique plus a low spanning-tree growth order enumerates every connected projective word",
        "component_law=every disconnected edge-8 low graph is a finite five-level component plus one free star center",
        "scale_law=dilating the connected component preserves its reduced channel and decreases transport loss and component debt",
        "conclusion=all 21 robust-edge-8 bodies close for every assignment of positive reflected levels",
        f"corollary=arbitrary-level body closure rises from {PRIOR_CLOSED_COUNT} to {CUMULATIVE_CLOSED_COUNT};remaining_bodies={REMAINING_BODY_COUNT}",
        "scope=reflected THM-2941 residual family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"edge9_sha256={sha256(EDGE9)}",
        f"body_digest={body_digest.hexdigest()}",
        f"assignment_digest={assignment_digest.hexdigest()}",
        f"source_sha256={source_sha}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
