#!/usr/bin/env python3
r"""Uniform closure of the complete 35-body robust-edge-9 block.

The two sum-threshold creation codes are

    DIIID (18 bodies),       IDIDI (17 bodies).

The full ``DIIID`` complement is ``K1,5`` plus one edge between two leaves.
It has two disconnected low-graph subcases: the five star edges alone leave a
low ``K5`` component, while all six high edges leave a low ``K5-e`` component.
In both cases the star center is a free projective scale.  The low ``K5`` bank
has the two familiar primitive types.  A new exact five-vertex clique-anchor
classification proves that a fixed labelled low ``K5-e`` component has 168
primitive words, or 14 projective orbits.

For either component, choose a pair inside the five-level component and bound
the free center's debt by its level-one value.  Dilating the component
multiplies the pair gcd and decreases component debt, while every positive
center level lies below the same debt bound.  Thus these are uniform cylinders,
not finite center checks.

Every other nonempty high subgraph has an anchorable maximum low clique.  The
enhanced classifier tries all maximum independent sets until every outside
vertex has a low anchor, repairing the first-choice limitation exposed by one
five-edge ``DIIID`` subgraph.

The resulting exact census contains 440,352 certificate rows:

    DIIID 259,200,       IDIDI 181,152.

All margins are positive and every selected certificate is low-phase.  Exact
reflected-arc controls at scales ``1,2,5`` replay the weakest row on each body,
giving 105 transport checks.  Therefore all 35 robust-edge-9 bodies close for
arbitrary positive reflected levels.  The cumulative body count is
2,421/3,003, leaving 582 bodies with at most eight robust edges.  This is
scoped to the THM-2941 reflected residual, not physical LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EDGE10 = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_edge10_threshold_block_uniform_closure_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_edge9_threshold_block_uniform_closure_thm2941.out"
)
EXPECTED_EDGE10_SHA256 = "cd2a4e84b2527f3fc7bf79980d2816ad207d770639fae9b42a1b9202d54cd2cd"
EXPECTED_SEMANTIC_SHA256 = "d8e5718d66156a8ddb51227d7bbb9f0474e2fcaea687e2662e029b0db580974a"

BODY_COUNT = 3003
PRIOR_CLOSED_COUNT = 2386
ROBUST_EDGE9_BODY_COUNT = 35
CUMULATIVE_CLOSED_COUNT = 2421
REMAINING_BODY_COUNT = 582
EXPECTED_CODE_HIST = (("DIIID", 18), ("IDIDI", 17))
EXPECTED_CODE_ROW_COUNTS = (("DIIID", 259200), ("IDIDI", 181152))
CONNECTED_ASSIGNMENT_COUNT = 433008
STAR5_CYLINDER_TEMPLATE_COUNT = 4320
K5E_CYLINDER_TEMPLATE_COUNT = 3024
CERTIFICATE_ROW_COUNT = 440352
K5E_COMPONENT_WORD_COUNT = 168
K5E_COMPONENT_ORBIT_COUNT = 14
DIRECT_CONTROL_SCALES = (1, 2, 5)
DIRECT_CONTROL_COUNT = 105
EXPECTED_WEAKEST = (
    F(15922411415672011132257703, 273113614278530404150800600),
    "IDIDI",
    "(5, (3, 3, 2, 1, 1, 0), 1)",
    (1, 2, 4, 9, 12, 13),
    ((1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)),
    ((1, 5), (2, 4), (2, 5), (3, 4), (4, 5)),
    (12, 3, 2, 4, 9, 16),
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
EXPECTED_STAR5_CYLINDER_WEAKEST = (
    F(131561040198654363844061, 1976502746201064220713840),
    "DIIID",
    "star5_cylinder",
    (1, 4, 5, 6, 7, 14),
    ((0, 5), (1, 5), (2, 5), (3, 4), (3, 5), (4, 5)),
    ((0, 5), (1, 5), (2, 5), (3, 5), (4, 5)),
    (1, 3, 4, 2, 6),
    F(17, 240),
    2,
    3,
    2,
    1,
    2,
    5100,
    6,
    427,
    "low",
)
EXPECTED_K5E_CYLINDER_WEAKEST = (
    F(621859723357427846213739, 9331732975464874113195364),
    "DIIID",
    "k5e_cylinder",
    (1, 2, 3, 7, 11, 14),
    ((0, 5), (1, 5), (2, 5), (3, 4), (3, 5), (4, 5)),
    ((0, 5), (1, 5), (2, 5), (3, 4), (3, 5), (4, 5)),
    (6, 2, 4, 8, 1),
    F(445, 6468),
    2,
    3,
    1,
    2,
    4,
    5972,
    1,
    496,
    "low",
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


require(sha256(EDGE10) == EXPECTED_EDGE10_SHA256,
        ("edge-10 theorem changed", sha256(EDGE10), EXPECTED_EDGE10_SHA256))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E10 = import_module("robust_edge9_edge10", EDGE10)
R = E10.R
LOW = E10.LOW
X = E10.X
U = E10.U
E11 = E10.E11
E12 = E10.E12
SYMMETRIC_LOW_RATIOS = E11.E13.SYMMETRIC_LOW_RATIOS


def anchorable_maximum_independent_set(vertex_count: int, high_edges):
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    for size in range(vertex_count, 0, -1):
        cores = tuple(
            vertices
            for vertices in combinations(range(vertex_count), size)
            if all(edge not in high_set for edge in combinations(vertices, 2))
        )
        if not cores:
            continue
        for core in cores:
            outside = tuple(index for index in range(vertex_count) if index not in core)
            if all(
                any(tuple(sorted((outside_index, core_index))) not in high_set
                    for core_index in core)
                for outside_index in outside
            ):
                return core
        raise RuntimeError(("no anchorable maximum low clique", vertex_count, high_set, cores))
    raise RuntimeError(("no low clique", vertex_count, high_set))


def clique_types(size: int):
    return tuple(
        row
        for row in combinations(LOW.NORMALIZED_VERTICES, size)
        if row[0] == 1 and LOW.is_low_clique(row)
    )


def primitive_word(normalized: tuple[F, ...]) -> tuple[int, ...]:
    denominator = lcm(*(value.denominator for value in normalized))
    integers = tuple(int(value * denominator) for value in normalized)
    common = gcd(*integers)
    primitive = tuple(value // common for value in integers)
    require(gcd(*primitive) == 1, (normalized, primitive))
    return primitive


def classified_words(vertex_count: int, high_edges) -> tuple[tuple[int, ...], ...]:
    """All primitive words with this labelled high graph and an anchorable core."""
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    core = anchorable_maximum_independent_set(vertex_count, high_set)
    outside = tuple(index for index in range(vertex_count) if index not in core)
    words = set()
    for clique in clique_types(len(core)):
        for clique_values in permutations(clique):
            partial = [None] * vertex_count
            for index, value in zip(core, clique_values):
                partial[index] = value
            candidate_banks = []
            for outside_index in outside:
                anchor = next(
                    core_index
                    for core_index in core
                    if tuple(sorted((outside_index, core_index))) not in high_set
                )
                candidate_banks.append(
                    tuple(sorted({
                        partial[anchor] * ratio
                        for ratio in SYMMETRIC_LOW_RATIOS
                    }))
                )
            for outside_values in product(*candidate_banks):
                if len(set(clique_values + outside_values)) != vertex_count:
                    continue
                word = partial[:]
                for index, value in zip(outside, outside_values):
                    word[index] = value
                if not all(
                    (not LOW.low_adjacent(word[i], word[j]))
                    == (tuple(sorted((i, j))) in high_set)
                    for i, j in combinations(range(vertex_count), 2)
                ):
                    continue
                minimum = min(word)
                words.add(primitive_word(tuple(value / minimum for value in word)))
    for word in words:
        actual_high = frozenset(
            edge for edge in combinations(range(vertex_count), 2)
            if not LOW.low_adjacent(F(word[edge[0]]), F(word[edge[1]]))
        )
        require(actual_high == high_set, (word, actual_high, high_set))
    return tuple(sorted(words))


def orbit_count(words, vertex_count: int, high_edges) -> int:
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    automorphisms = tuple(
        permutation
        for permutation in permutations(range(vertex_count))
        if frozenset(
            tuple(sorted((permutation[i], permutation[j])))
            for i, j in high_set
        ) == high_set
    )
    representatives = {
        min(tuple(word[permutation[index]] for index in range(vertex_count))
            for permutation in automorphisms)
        for word in words
    }
    return len(representatives)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    k5e_high_edge = ((3, 4),)
    k5e_words = classified_words(5, k5e_high_edge)
    require(len(k5e_words) == K5E_COMPONENT_WORD_COUNT, len(k5e_words))
    require(orbit_count(k5e_words, 5, k5e_high_edge) == K5E_COMPONENT_ORBIT_COUNT,
            orbit_count(k5e_words, 5, k5e_high_edge))

    bodies = []
    code_hist: Counter[str] = Counter()
    body_digest = hashlib.sha256()
    universal_exceptions = {row[0] for row in U.EXPECTED_EXCEPTIONS}
    for body in combinations(range(1, 15), 6):
        ruler, universal_debt, robust_edges = LOW.robust_edges(body)
        if len(robust_edges) != 9:
            continue
        require(body not in universal_exceptions,
                ("edge-9 body lacks complete same-level graph", body))
        threshold = E12.nonrobust_threshold(body, ruler, universal_debt)
        code = E12.threshold_creation_code(body, threshold)
        nonedges = tuple(
            edge for edge in combinations(range(6), 2) if edge not in set(robust_edges)
        )
        require(code in {"DIIID", "IDIDI"}, (body, code, nonedges))
        if code == "DIIID":
            require(
                nonedges == ((0, 5), (1, 5), (2, 5), (3, 4), (3, 5), (4, 5)),
                ("DIIID complement changed", body, nonedges),
            )
        code_hist[code] += 1
        bodies.append((body, ruler, universal_debt, robust_edges, nonedges, threshold, code))
        body_digest.update(
            f"{body}|{ruler}|{universal_debt}|{threshold}|{robust_edges}|{nonedges}|{code}\n".encode()
        )
    require(tuple(sorted(code_hist.items())) == EXPECTED_CODE_HIST, code_hist)
    require(len(bodies) == ROBUST_EDGE9_BODY_COUNT, len(bodies))

    word_cache = {}
    rows = []
    code_row_counts: Counter[str] = Counter()
    certificate_kinds: Counter[str] = Counter()
    assignment_digest = hashlib.sha256()
    weakest_by_body = {}
    connected_rows = 0
    star5_rows = 0
    k5e_rows = 0

    diiid_full = ((0, 5), (1, 5), (2, 5), (3, 4), (3, 5), (4, 5))
    diiid_star5 = tuple(edge for edge in diiid_full if edge != (3, 4))
    for body, ruler, _, robust_edges, nonedges, _, code in bodies:
        safe_ruler, safe_ranges = R.safe_cell_ranges(body)
        require(safe_ruler == ruler, (body, safe_ruler, ruler))
        profile_cache = {}
        robust_set = set(robust_edges)

        for size in range(1, len(nonedges) + 1):
            for high_edges in combinations(nonedges, size):
                high_edges = tuple(sorted(high_edges))
                if code == "DIIID" and high_edges in (diiid_star5, diiid_full):
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
                    require(all(LOW.low_adjacent(F(levels[i]), F(levels[j]))
                                for i, j in robust_set),
                            (body, levels, robust_edges))
                    debt = E10.singleton_debt(body, ruler, levels)
                    best = E10.best_certificate(
                        body, ruler, safe_ranges, levels, debt, profile_cache
                    )
                    require(best[0] > 0,
                            ("finite edge-9 assignment failed", code, category,
                             body, nonedges, high_edges, levels, best))
                    row = (
                        best[0], code, category, body, nonedges,
                        high_edges, levels,
                    ) + best[1:]
                    rows.append(row)
                    connected_rows += 1
                    code_row_counts[code] += 1
                    certificate_kinds[best[-1]] += 1
                    assignment_digest.update(
                        f"{code}|{category}|{body}|{ruler}|{nonedges}|{high_edges}|{levels}|{debt}|{best}\n".encode()
                    )
                    if body not in weakest_by_body or row < weakest_by_body[body]:
                        weakest_by_body[body] = row

        if code == "DIIID":
            leaf_pairs = tuple(combinations(range(5), 2))
            for leaf_type in E10.LOW_K5_PRIMITIVES:
                for leaf_levels in permutations(leaf_type):
                    debt = E10.free_center_debt_bound(body, ruler, leaf_levels)
                    best = E10.best_certificate(
                        body, ruler, safe_ranges, leaf_levels, debt,
                        profile_cache, pair_indices=leaf_pairs,
                    )
                    require(best[0] > 0,
                            ("star5 edge-9 cylinder failed", body,
                             leaf_type, leaf_levels, best, debt))
                    row = (
                        best[0], code, "star5_cylinder", body,
                        nonedges, diiid_star5, leaf_levels,
                    ) + best[1:]
                    rows.append(row)
                    star5_rows += 1
                    code_row_counts[code] += 1
                    certificate_kinds[best[-1]] += 1
                    assignment_digest.update(
                        f"{code}|star5_cylinder|{body}|{ruler}|{nonedges}|{leaf_levels}|{debt}|{best}\n".encode()
                    )
                    if body not in weakest_by_body or row < weakest_by_body[body]:
                        weakest_by_body[body] = row

            for leaf_levels in k5e_words:
                debt = E10.free_center_debt_bound(body, ruler, leaf_levels)
                best = E10.best_certificate(
                    body, ruler, safe_ranges, leaf_levels, debt,
                    profile_cache, pair_indices=leaf_pairs,
                )
                require(best[0] > 0,
                        ("K5-e edge-9 cylinder failed", body,
                         leaf_levels, best, debt))
                row = (
                    best[0], code, "k5e_cylinder", body,
                    nonedges, diiid_full, leaf_levels,
                ) + best[1:]
                rows.append(row)
                k5e_rows += 1
                code_row_counts[code] += 1
                certificate_kinds[best[-1]] += 1
                assignment_digest.update(
                    f"{code}|k5e_cylinder|{body}|{ruler}|{nonedges}|{leaf_levels}|{debt}|{best}\n".encode()
                )
                if body not in weakest_by_body or row < weakest_by_body[body]:
                    weakest_by_body[body] = row

    require(connected_rows == CONNECTED_ASSIGNMENT_COUNT, connected_rows)
    require(star5_rows == STAR5_CYLINDER_TEMPLATE_COUNT, star5_rows)
    require(k5e_rows == K5E_CYLINDER_TEMPLATE_COUNT, k5e_rows)
    require(len(rows) == CERTIFICATE_ROW_COUNT, len(rows))
    require(tuple(sorted(code_row_counts.items())) == EXPECTED_CODE_ROW_COUNTS,
            code_row_counts)
    require(min(rows) == EXPECTED_WEAKEST, min(rows))
    star5_weakest = min(row for row in rows if row[2] == "star5_cylinder")
    k5e_weakest = min(row for row in rows if row[2] == "k5e_cylinder")
    require(star5_weakest == EXPECTED_STAR5_CYLINDER_WEAKEST, star5_weakest)
    require(k5e_weakest == EXPECTED_K5E_CYLINDER_WEAKEST, k5e_weakest)
    require(certificate_kinds == Counter({"low": CERTIFICATE_ROW_COUNT}),
            certificate_kinds)
    require(len(weakest_by_body) == ROBUST_EDGE9_BODY_COUNT, len(weakest_by_body))

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
            if category in {"star5_cylinder", "k5e_cylinder"}:
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
        k5e_words,
        tuple(bodies),
        tuple(sorted((edges, words) for edges, words in word_cache.items())),
        min(rows),
        star5_weakest,
        k5e_weakest,
        tuple(sorted(code_row_counts.items())),
        tuple(sorted(certificate_kinds.items())),
        connected_rows,
        star5_rows,
        k5e_rows,
        tuple(sorted(weakest_by_body.items())),
        DIRECT_CONTROL_SCALES,
        direct_controls,
        body_digest.hexdigest(),
        assignment_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    weakest = min(rows)
    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 reflected robust-edge-9 threshold-block uniform closure exact referee",
        f"edge9_code_hist={tuple(sorted(code_hist.items()))};"
        f"code_row_counts={tuple(sorted(code_row_counts.items()))}",
        f"connected_projective_rows={connected_rows};"
        f"star5_cylinder_templates={star5_rows};"
        f"k5e_component_words={len(k5e_words)};"
        f"k5e_component_orbits={K5E_COMPONENT_ORBIT_COUNT};"
        f"k5e_cylinder_templates={k5e_rows};"
        f"total_certificate_rows={len(rows)};"
        f"certificate_kinds={tuple(sorted(certificate_kinds.items()))}",
        f"weakest_margin={qtext(weakest[0])};code={weakest[1]};"
        f"category={weakest[2]};body={weakest[3]};"
        f"nonrobust_edges={weakest[4]};high_edges={weakest[5]};"
        f"levels={weakest[6]};floor={qtext(weakest[7])};"
        f"certificate_pair={(weakest[8],weakest[9])};"
        f"reduced_channel={(weakest[10],weakest[11])};gcd={weakest[12]};"
        f"j={weakest[13]};lift={weakest[14]};endpoint_radius={weakest[15]};"
        f"kind={weakest[16]}",
        f"star5_cylinder_weakest={qtext(star5_weakest[0])};"
        f"k5e_cylinder_weakest={qtext(k5e_weakest[0])}",
        f"direct_control_scales={DIRECT_CONTROL_SCALES};"
        f"per_body_weakest_controls={direct_controls}",
        "component_law=each disconnected low component has a finite primitive bank while the isolated center contributes only worst-case level-one debt",
        "scale_law=dilating a connected component preserves its reduced channel and cell and decreases loss and component debt",
        "conclusion=all 35 robust-edge-9 bodies close for every assignment of positive reflected levels",
        f"corollary=arbitrary-level body closure rises from {PRIOR_CLOSED_COUNT} to {CUMULATIVE_CLOSED_COUNT};remaining_bodies={REMAINING_BODY_COUNT}",
        "scope=reflected THM-2941 residual family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"edge10_sha256={sha256(EDGE10)}",
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
