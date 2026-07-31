#!/usr/bin/env python3
r"""Clique-anchor closure of all 15 robust-edge-11 bodies.

The sum-threshold complement has one of two creation codes at eleven robust
edges:

* ``IDIII`` on 7 bodies, whose four nonrobust edges form ``K1,4``;
* ``IIDID`` on 8 bodies, whose four nonrobust edges form a paw.

An uncertified word has six distinct levels and its nonempty high-phase graph
is a subgraph of that complement.  Besides the earlier one-edge, path, and
three-star cases, this introduces ``2K2``, ``P4``, ``K1,4``, and the paw.  The
triangle case remains projectively empty.

This referee uses a uniform finite classifier.  For a prescribed labelled high
graph ``H``, choose a maximum independent set ``C``.  Its levels form a
low-phase clique, one of the exact normalized clique bank.  Every vertex
outside ``C`` has a low neighbour in ``C`` for the graphs at hand, so its ratio
to the first such neighbour lies in the finite symmetric low-ratio bank.  All
permutations of the clique values and all resulting outside candidates are
checked against the complete high/low adjacency matrix.  The exact fixed-label
word counts (and projective orbit counts) are

    edge   96 (2),       P3     96 (8),
    2K2    16 (1),       K1,3  312 (26),
    K3      0 (0),       P4    128 (32),
    K1,4  672 (28),      paw   184 (46).

Enumerating every nonempty high subgraph on all 15 bodies gives 33,216
primitive labelled assignments.  Every assignment has a positive exact
body-safe one-pair margin over its actual singleton debt, and every selected
certificate is low-phase.  Integer dilation preserves the reduced channel and
cell while decreasing transport loss and debt.  Direct reflected-arc controls
at scales ``1,2,5`` verify all 99,648 chosen certificates.

Hence all 15 robust-edge-11 bodies close for arbitrary positive reflected
levels.  The cumulative body count is 2,354/3,003, leaving 649 bodies with at
most ten robust edges.  This is scoped to the THM-2941 reflected residual, not
physical LRC(14).
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
BASE = ROOT / "04-computation" / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
LOW_PHASE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
)
EXCEPTIONAL = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_exceptional_low_channel_uniform_closure_thm2941.py"
)
UNIVERSAL = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
)
EDGE13 = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_k6_minus_adjacent_pair_uniform_closure_thm2941.py"
)
EDGE13_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_k6_minus_adjacent_pair_uniform_closure_thm2941.out"
)
EDGE12 = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_edge12_threshold_shapes_uniform_closure_thm2941.py"
)
EDGE12_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_edge12_threshold_shapes_uniform_closure_thm2941.out"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_edge11_threshold_shapes_uniform_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_LOW_PHASE_SHA256 = "416c36f16f7c821feb8d260882711d2717069147b8604a93ba60432785cf1d1c"
EXPECTED_EXCEPTIONAL_SHA256 = "bde992db1edbd9dd744ff22744a1afef79cf4bcc54a2f918793c2603f062df7c"
EXPECTED_UNIVERSAL_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"
EXPECTED_EDGE13_SHA256 = "2fdde052252ea6067fd587dcf1d41db2a4ce381da4b4ab21f8cd5c0efb7d8e3e"
EXPECTED_EDGE13_OUTPUT_SHA256 = "c64d3fbc77d0e917b6bd3b8d498cfe4b568f2cd6fcc68890475654f1be1527d9"
EXPECTED_EDGE12_SHA256 = "a47791c4cbb211e401b080f7b56380bae5c52b1ce865d596be594210cdb6eee6"
EXPECTED_EDGE12_OUTPUT_SHA256 = "e7dc1bba56ae281685e76b8b855fb6b702f1116c9797246b959249d85c63d576"
EXPECTED_SEMANTIC_SHA256 = "3e3ed7b68c58233aef3b4f2450f60b945ec59cd7c7826783abece85eb16a739f"

BODY_COUNT = 3003
PRIOR_CLOSED_COUNT = 2339
ROBUST_EDGE11_BODY_COUNT = 15
CUMULATIVE_CLOSED_COUNT = 2354
REMAINING_BODY_COUNT = 649
EXPECTED_EDGE11_CODE_HIST = (("IDIII", 7), ("IIDID", 8))
EXPECTED_SHAPE_BANK = (
    ("edge", 96, 2),
    ("p3", 96, 8),
    ("2k2", 16, 1),
    ("star3", 312, 26),
    ("triangle", 0, 0),
    ("p4", 128, 32),
    ("star4", 672, 28),
    ("paw", 184, 46),
)
EXPECTED_ASSIGNMENT_CATEGORIES = (
    ("2k2", 128),
    ("edge", 5760),
    ("p3", 7872),
    ("p4", 2048),
    ("paw", 1472),
    ("star3", 11232),
    ("star4", 4704),
)
ASSIGNMENT_COUNT = 33216
DIRECT_CONTROL_SCALES = (1, 2, 5)
DIRECT_CONTROL_COUNT = 99648
EXPECTED_WEAKEST = (
    F(2703015168079847245426187, 43021744871641066038616290),
    "star3",
    (1, 2, 3, 7, 13, 14),
    ((2, 5), (3, 4), (3, 5), (4, 5)),
    ((2, 5), (3, 5), (4, 5)),
    (6, 3, 1, 2, 4, 9),
    F(81, 1274),
    0,
    1,
    2,
    1,
    3,
    4874,
    2,
    666,
    "low",
)

SHAPE_KEYS = {
    (1, (1, 1, 0, 0, 0, 0), 0): "edge",
    (2, (2, 1, 1, 0, 0, 0), 0): "p3",
    (2, (1, 1, 1, 1, 0, 0), 0): "2k2",
    (3, (3, 1, 1, 1, 0, 0), 0): "star3",
    (3, (2, 2, 2, 0, 0, 0), 1): "triangle",
    (3, (2, 2, 1, 1, 0, 0), 0): "p4",
    (4, (4, 1, 1, 1, 1, 0), 0): "star4",
    (4, (3, 2, 2, 1, 0, 0), 1): "paw",
}
CANONICAL_SHAPES = (
    ("edge", ((0, 1),)),
    ("p3", ((0, 1), (1, 2))),
    ("2k2", ((0, 1), (2, 3))),
    ("star3", ((0, 1), (0, 2), (0, 3))),
    ("triangle", ((0, 1), (0, 2), (1, 2))),
    ("p4", ((0, 1), (1, 2), (2, 3))),
    ("star4", ((0, 1), (0, 2), (0, 3), (0, 4))),
    ("paw", ((0, 1), (0, 2), (0, 3), (1, 2))),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


for path, expected in (
    (BASE, EXPECTED_BASE_SHA256),
    (LOW_PHASE, EXPECTED_LOW_PHASE_SHA256),
    (EXCEPTIONAL, EXPECTED_EXCEPTIONAL_SHA256),
    (UNIVERSAL, EXPECTED_UNIVERSAL_SHA256),
    (EDGE13, EXPECTED_EDGE13_SHA256),
    (EDGE13_OUTPUT, EXPECTED_EDGE13_OUTPUT_SHA256),
    (EDGE12, EXPECTED_EDGE12_SHA256),
    (EDGE12_OUTPUT, EXPECTED_EDGE12_OUTPUT_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("robust_edge11_base", BASE)
LOW = import_module("robust_edge11_low_phase", LOW_PHASE)
X = import_module("robust_edge11_optimizer", EXCEPTIONAL)
U = import_module("robust_edge11_universal", UNIVERSAL)
E13 = import_module("robust_edge11_edge13", EDGE13)
E12 = import_module("robust_edge11_edge12", EDGE12)


def singleton_debt(body: tuple[int, ...], ruler: int, levels: tuple[int, ...]) -> F:
    return sum(
        (F(e, 7 * (level * ruler - e)) for e, level in zip(body, levels)),
        F(0),
    )


def intersection_mass(first, second) -> F:
    i = 0
    j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(
            F(0),
            min(first[i][1], second[j][1]) - max(first[i][0], second[j][0]),
        )
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def graph_key(edges) -> tuple[int, tuple[int, ...], int]:
    edge_set = frozenset(tuple(sorted(edge)) for edge in edges)
    degrees = tuple(
        sorted((sum(vertex in edge for edge in edge_set) for vertex in range(6)), reverse=True)
    )
    triangles = sum(
        all(edge in edge_set for edge in combinations(vertices, 2))
        for vertices in combinations(range(6), 3)
    )
    return len(edge_set), degrees, triangles


def shape_name(edges) -> str:
    key = graph_key(edges)
    require(key in SHAPE_KEYS, ("unknown high-graph shape", edges, key))
    return SHAPE_KEYS[key]


def clique_types(size: int):
    return tuple(
        row
        for row in combinations(LOW.NORMALIZED_VERTICES, size)
        if row[0] == 1 and LOW.is_low_clique(row)
    )


def maximum_independent_set(high_edges) -> tuple[int, ...]:
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    for size in range(6, 0, -1):
        rows = tuple(
            vertices
            for vertices in combinations(range(6), size)
            if all(edge not in high_set for edge in combinations(vertices, 2))
        )
        if rows:
            return rows[0]
    raise RuntimeError(("no independent set", high_edges))


def primitive_word(normalized: tuple[F, ...]) -> tuple[int, ...]:
    denominator = lcm(*(value.denominator for value in normalized))
    integers = tuple(int(value * denominator) for value in normalized)
    common = gcd(*integers)
    primitive = tuple(value // common for value in integers)
    require(gcd(*primitive) == 1, (normalized, primitive))
    return primitive


def classified_words(high_edges) -> tuple[tuple[int, ...], ...]:
    """All primitive words with exactly this labelled high graph."""
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    core = maximum_independent_set(high_set)
    outside = tuple(index for index in range(6) if index not in core)
    words = set()
    for clique in clique_types(len(core)):
        for clique_values in permutations(clique):
            partial = [None] * 6
            for index, value in zip(core, clique_values):
                partial[index] = value
            candidate_banks = []
            for outside_index in outside:
                low_anchors = tuple(
                    core_index
                    for core_index in core
                    if tuple(sorted((outside_index, core_index))) not in high_set
                )
                require(low_anchors,
                        ("finite clique-anchor condition failed", high_set, core, outside_index))
                anchor = low_anchors[0]
                candidate_banks.append(
                    tuple(sorted({
                        partial[anchor] * ratio
                        for ratio in E13.SYMMETRIC_LOW_RATIOS
                    }))
                )
            for outside_values in product(*candidate_banks):
                if len(set(clique_values + outside_values)) != 6:
                    continue
                word = partial[:]
                for index, value in zip(outside, outside_values):
                    word[index] = value
                if not all(
                    (not LOW.low_adjacent(word[i], word[j]))
                    == (tuple(sorted((i, j))) in high_set)
                    for i, j in combinations(range(6), 2)
                ):
                    continue
                minimum = min(word)
                normalized = tuple(value / minimum for value in word)
                words.add(primitive_word(normalized))
    for word in words:
        actual_high = frozenset(
            edge for edge in combinations(range(6), 2)
            if not LOW.low_adjacent(F(word[edge[0]]), F(word[edge[1]]))
        )
        require(actual_high == high_set, (word, actual_high, high_set))
    return tuple(sorted(words))


def graph_automorphisms(high_edges):
    high_set = frozenset(tuple(sorted(edge)) for edge in high_edges)
    return tuple(
        permutation
        for permutation in permutations(range(6))
        if frozenset(
            tuple(sorted((permutation[i], permutation[j])))
            for i, j in high_set
        ) == high_set
    )


def orbit_count(words, high_edges) -> int:
    automorphisms = graph_automorphisms(high_edges)
    representatives = {
        min(tuple(word[permutation[index]] for index in range(6))
            for permutation in automorphisms)
        for word in words
    }
    return len(representatives)


def best_certificate(
    body: tuple[int, ...],
    ruler: int,
    safe_ranges: tuple[tuple[int, int], ...],
    levels: tuple[int, ...],
    debt: F,
    profile_cache: dict,
):
    certificates = []
    for i, j in combinations(range(6), 2):
        p, q = levels[i], levels[j]
        divisor = gcd(p, q)
        P, Q = p // divisor, q // divisor
        require(P != Q and gcd(P, Q) == 1, (levels, i, j, P, Q))
        if P + Q <= 7:
            key = (i, j, P, Q)
            if key not in profile_cache:
                profile_cache[key] = X.best_skeleton_profile(
                    safe_ranges,
                    ruler,
                    P * body[j] - Q * body[i],
                    P,
                    Q,
                )
            floor, cell, lift, endpoint_radius = profile_cache[key]
            kind = "low"
        else:
            floor = LOW.FIBER_FLOOR
            cell = safe_ranges[0][0]
            lift = -1
            endpoint_radius = -1
            kind = "high"
        margin = floor - F(4 * (body[i] + body[j]), divisor * ruler) - debt
        certificates.append(
            (
                margin,
                floor,
                i,
                j,
                P,
                Q,
                divisor,
                cell,
                lift,
                endpoint_radius,
                kind,
            )
        )
    return max(certificates)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    shape_bank = []
    for name, canonical_edges in CANONICAL_SHAPES:
        require(shape_name(canonical_edges) == name, (name, canonical_edges))
        words = classified_words(canonical_edges)
        shape_bank.append((name, len(words), orbit_count(words, canonical_edges)))
    require(tuple(shape_bank) == EXPECTED_SHAPE_BANK, shape_bank)

    bodies = []
    code_hist: Counter[str] = Counter()
    body_digest = hashlib.sha256()
    universal_exceptions = {row[0] for row in U.EXPECTED_EXCEPTIONS}
    for body in combinations(range(1, 15), 6):
        ruler, debt, robust_edges = LOW.robust_edges(body)
        if len(robust_edges) != 11:
            continue
        require(body not in universal_exceptions,
                ("edge-11 body lacks complete same-level graph", body))
        threshold = E12.nonrobust_threshold(body, ruler, debt)
        code = E12.threshold_creation_code(body, threshold)
        nonedges = tuple(
            edge for edge in combinations(range(6), 2) if edge not in set(robust_edges)
        )
        complement_shape = shape_name(nonedges)
        require((code, complement_shape) in (("IDIII", "star4"), ("IIDID", "paw")),
                (body, code, nonedges, complement_shape))
        code_hist[code] += 1
        bodies.append((body, ruler, debt, robust_edges, nonedges, threshold, code, complement_shape))
        body_digest.update(
            f"{body}|{ruler}|{debt}|{threshold}|{robust_edges}|{nonedges}|{code}|{complement_shape}\n".encode()
        )
    require(tuple(sorted(code_hist.items())) == EXPECTED_EDGE11_CODE_HIST, code_hist)
    require(len(bodies) == ROBUST_EDGE11_BODY_COUNT, len(bodies))

    word_cache = {}
    assignment_rows = []
    category_counts: Counter[str] = Counter()
    certificate_kinds: Counter[str] = Counter()
    certificate_floors: Counter[F] = Counter()
    direct_controls = 0
    assignment_digest = hashlib.sha256()

    for body, ruler, _, robust_edges, nonedges, _, _, _ in bodies:
        safe_ruler, safe_ranges = R.safe_cell_ranges(body)
        require(safe_ruler == ruler, (body, safe_ruler, ruler))
        profile_cache = {}
        robust_set = set(robust_edges)
        for size in range(1, len(nonedges) + 1):
            for high_edges in combinations(nonedges, size):
                high_edges = tuple(sorted(high_edges))
                if high_edges not in word_cache:
                    word_cache[high_edges] = classified_words(high_edges)
                category = shape_name(high_edges)
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
                    debt = singleton_debt(body, ruler, levels)
                    best = best_certificate(
                        body, ruler, safe_ranges, levels, debt, profile_cache
                    )
                    require(best[0] > 0,
                            ("primitive edge-11 assignment failed", category,
                             body, nonedges, high_edges, levels, best))
                    (
                        margin,
                        floor,
                        i,
                        j,
                        P,
                        Q,
                        divisor,
                        cell,
                        lift,
                        endpoint_radius,
                        kind,
                    ) = best
                    a, b = body[i], body[j]
                    p, q = levels[i], levels[j]
                    require(R.body_cell_is_safe(ruler, body, cell),
                            ("certificate cell unsafe", body, levels, (i, j), cell))
                    for scale in DIRECT_CONTROL_SCALES:
                        actual = intersection_mass(
                            R.reflected_level_arcs(ruler, a, scale * p, cell),
                            R.reflected_level_arcs(ruler, b, scale * q, cell),
                        )
                        transported = floor - F(
                            4 * (a + b), scale * divisor * ruler
                        )
                        scaled_debt = singleton_debt(
                            body, ruler, tuple(scale * value for value in levels)
                        )
                        require(actual >= transported > scaled_debt,
                                (body, levels, (i, j), scale, actual,
                                 transported, scaled_debt, best))
                        direct_controls += 1
                    row = (
                        margin,
                        category,
                        body,
                        nonedges,
                        high_edges,
                        levels,
                        floor,
                        i,
                        j,
                        P,
                        Q,
                        divisor,
                        cell,
                        lift,
                        endpoint_radius,
                        kind,
                    )
                    assignment_rows.append(row)
                    category_counts[category] += 1
                    certificate_kinds[kind] += 1
                    certificate_floors[floor] += 1
                    assignment_digest.update(
                        f"{category}|{body}|{ruler}|{nonedges}|{high_edges}|{levels}|{debt}|{best}\n".encode()
                    )

    require(tuple(sorted(category_counts.items())) == EXPECTED_ASSIGNMENT_CATEGORIES,
            category_counts)
    require(len(assignment_rows) == ASSIGNMENT_COUNT, len(assignment_rows))
    require(direct_controls == DIRECT_CONTROL_COUNT, direct_controls)
    require(min(assignment_rows) == EXPECTED_WEAKEST, min(assignment_rows))
    require(certificate_kinds == Counter({"low": ASSIGNMENT_COUNT}), certificate_kinds)
    require(PRIOR_CLOSED_COUNT + len(bodies) == CUMULATIVE_CLOSED_COUNT,
            (PRIOR_CLOSED_COUNT, len(bodies), CUMULATIVE_CLOSED_COUNT))
    require(BODY_COUNT - CUMULATIVE_CLOSED_COUNT == REMAINING_BODY_COUNT,
            (BODY_COUNT, CUMULATIVE_CLOSED_COUNT, REMAINING_BODY_COUNT))

    semantic_payload = (
        tuple(shape_bank),
        tuple(bodies),
        tuple(sorted((edges, words) for edges, words in word_cache.items())),
        min(assignment_rows),
        tuple(sorted(category_counts.items())),
        tuple(sorted(certificate_kinds.items())),
        tuple(sorted(certificate_floors.items())),
        len(assignment_rows),
        DIRECT_CONTROL_SCALES,
        direct_controls,
        body_digest.hexdigest(),
        assignment_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    (
        weakest_margin,
        weakest_category,
        weakest_body,
        weakest_nonedges,
        weakest_high_edges,
        weakest_levels,
        weakest_floor,
        weakest_i,
        weakest_j,
        weakest_P,
        weakest_Q,
        weakest_divisor,
        weakest_cell,
        weakest_lift,
        weakest_radius,
        weakest_kind,
    ) = min(assignment_rows)
    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 reflected robust-edge-11 clique-anchor uniform closure exact referee",
        f"edge11_code_hist={tuple(sorted(code_hist.items()))};shape_bank={tuple(shape_bank)}",
        f"assignment_categories={tuple(sorted(category_counts.items()))};"
        f"total_primitive_assignments={len(assignment_rows)};"
        f"certificate_kinds={tuple(sorted(certificate_kinds.items()))};"
        f"certificate_floor_values={len(certificate_floors)}",
        f"weakest_margin={qtext(weakest_margin)};category={weakest_category};"
        f"body={weakest_body};nonrobust_edges={weakest_nonedges};"
        f"high_edges={weakest_high_edges};levels={weakest_levels};"
        f"floor={qtext(weakest_floor)};certificate_pair={(weakest_i,weakest_j)};"
        f"reduced_channel={(weakest_P,weakest_Q)};gcd={weakest_divisor};"
        f"j={weakest_cell};lift={weakest_lift};endpoint_radius={weakest_radius};"
        f"kind={weakest_kind}",
        f"direct_control_scales={DIRECT_CONTROL_SCALES};direct_controls={direct_controls}",
        "classification_law=maximum low clique plus one finite symmetric-ratio bank per outside vertex",
        "scale_law=multiplying all levels preserves the reduced channel and cell, decreases transport loss and singleton debt",
        "conclusion=all 15 robust-edge-11 bodies close for every assignment of positive reflected levels",
        f"corollary=arbitrary-level body closure rises from {PRIOR_CLOSED_COUNT} to {CUMULATIVE_CLOSED_COUNT};remaining_bodies={REMAINING_BODY_COUNT}",
        "scope=reflected THM-2941 residual family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"low_phase_sha256={sha256(LOW_PHASE)}",
        f"exceptional_sha256={sha256(EXCEPTIONAL)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"edge13_sha256={sha256(EDGE13)}",
        f"edge13_output_sha256={sha256(EDGE13_OUTPUT)}",
        f"edge12_sha256={sha256(EDGE12)}",
        f"edge12_output_sha256={sha256(EDGE12_OUTPUT)}",
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
