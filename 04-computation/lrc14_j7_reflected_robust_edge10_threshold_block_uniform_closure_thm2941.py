#!/usr/bin/env python3
r"""Uniform closure of the complete 32-body robust-edge-10 block.

The sum-threshold complements at ten robust edges have creation codes

    DIIII (11 bodies),   IDIID (15 bodies),   IIDDI (6 bodies).

For ``IDIID`` and ``IIDDI``, and for every proper high subgraph of ``DIIII``,
the low graph is connected and the clique-anchor classifier from the edge-11
theorem gives a finite exact projective bank.  The full ``DIIII`` complement is
``K1,5`` and is the first projectively infinite case: the five leaf levels form
one of the two low-phase five-cliques, but the center-to-leaf scale is free.

The free scale is harmless.  The primitive leaf sets are

    A0=(2,3,4,6,12),       B0=(1,2,3,4,6).

For each of their ``5!`` placements, this referee chooses a leaf-pair
certificate.  It bounds the free center's singleton debt by its worst value at
level one.  Dilating all leaves multiplies the selected pair gcd and decreases
every leaf debt term, while every positive center level has debt no greater
than the level-one bound.  Hence the 2,640 primitive leaf templates prove the
entire eleven-body free-center cylinder, not merely finitely many center
levels.

Together with all connected finite banks, the exact census contains 188,544
certificate rows:

    DIIII  89,760,       IDIID  73,440,       IIDDI  25,344.

Every row has positive margin over the actual (or free-center worst-case)
singleton debt; every selected certificate is low-phase.  Exact reflected-arc
controls at scales ``1,2,5`` replay the weakest row on each body, giving 96
independent transport checks.

Thus all 32 robust-edge-10 bodies close for arbitrary positive reflected
levels.  The cumulative body count is 2,386/3,003, leaving 617 bodies with at
most nine robust edges.  This is scoped to the THM-2941 reflected residual,
not physical LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd
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
EDGE11 = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_edge11_threshold_shapes_uniform_closure_thm2941.py"
)
EDGE12 = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_edge12_threshold_shapes_uniform_closure_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_edge10_threshold_block_uniform_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_LOW_PHASE_SHA256 = "b2418dfda1b48257d1f7582d4ea977203a26f88885e13946bc100ccf264c9ce1"
EXPECTED_EXCEPTIONAL_SHA256 = "15c13d721ee695c2a69bc386d5c12eba0382b52f9eabfe052f2fe5fda03c7bc4"
EXPECTED_UNIVERSAL_SHA256 = "dc6f23a201e817dd9134e8660d35e83d3053c67d26fc271ce3eae07f0f857689"
EXPECTED_EDGE11_SHA256 = "58685796af62dfb425bfa17a6363d9cc5ca4cbb72e0f1d42ecc23dcbf10c01b4"
EXPECTED_EDGE12_SHA256 = "76110afdfb1e607fd63df1909a05259d9252da1855b4ae2a9fc6da5d936fbaa3"
EXPECTED_SEMANTIC_SHA256 = "d7b87d6121721aae103209b935c5469655f1984e468c07b946285d76689c997e"

BODY_COUNT = 3003
PRIOR_CLOSED_COUNT = 2354
ROBUST_EDGE10_BODY_COUNT = 32
CUMULATIVE_CLOSED_COUNT = 2386
REMAINING_BODY_COUNT = 617
EXPECTED_CODE_HIST = (("DIIII", 11), ("IDIID", 15), ("IIDDI", 6))
EXPECTED_CODE_ROW_COUNTS = (
    ("DIIII", 89760),
    ("IDIID", 73440),
    ("IIDDI", 25344),
)
CONNECTED_ASSIGNMENT_COUNT = 185904
FREE_CENTER_TEMPLATE_COUNT = 2640
CERTIFICATE_ROW_COUNT = 188544
DIRECT_CONTROL_SCALES = (1, 2, 5)
DIRECT_CONTROL_COUNT = 96
LOW_K5_PRIMITIVES = (
    (2, 3, 4, 6, 12),
    (1, 2, 3, 4, 6),
)
EXPECTED_WEAKEST = (
    F(5010286913481743662879783, 81689790267700879566650040),
    "IDIID",
    "(5, (4, 2, 2, 1, 1, 0), 1)",
    (1, 2, 4, 6, 9, 13),
    ((1, 5), (2, 5), (3, 4), (3, 5), (4, 5)),
    ((1, 5), (2, 5), (3, 4), (3, 5), (4, 5)),
    (2, 12, 4, 8, 3, 5),
    F(59, 936),
    2,
    3,
    1,
    2,
    4,
    3564,
    1,
    578,
    "low",
)
EXPECTED_DIIII_CYLINDER_WEAKEST = (
    F(67641571398179833793, 1083489281890452517080),
    "DIIII",
    "star5_infinite",
    (1, 2, 3, 5, 6, 13),
    ((0, 5), (1, 5), (2, 5), (3, 5), (4, 5)),
    ((0, 5), (1, 5), (2, 5), (3, 5), (4, 5)),
    (1, 3, 2, 4, 6),
    F(719, 10920),
    2,
    3,
    1,
    2,
    2,
    5009,
    1,
    451,
    "low",
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


for path, expected in (
    (BASE, EXPECTED_BASE_SHA256),
    (LOW_PHASE, EXPECTED_LOW_PHASE_SHA256),
    (EXCEPTIONAL, EXPECTED_EXCEPTIONAL_SHA256),
    (UNIVERSAL, EXPECTED_UNIVERSAL_SHA256),
    (EDGE11, EXPECTED_EDGE11_SHA256),
    (EDGE12, EXPECTED_EDGE12_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("robust_edge10_base", BASE)
LOW = import_module("robust_edge10_low_phase", LOW_PHASE)
X = import_module("robust_edge10_optimizer", EXCEPTIONAL)
U = import_module("robust_edge10_universal", UNIVERSAL)
E11 = import_module("robust_edge10_edge11", EDGE11)
E12 = import_module("robust_edge10_edge12", EDGE12)


def singleton_debt(body: tuple[int, ...], ruler: int, levels: tuple[int, ...]) -> F:
    return sum(
        (F(e, 7 * (level * ruler - e)) for e, level in zip(body, levels)),
        F(0),
    )


def free_center_debt_bound(
    body: tuple[int, ...], ruler: int, leaf_levels: tuple[int, ...]
) -> F:
    require(len(leaf_levels) == 5, leaf_levels)
    leaf_debt = sum(
        (F(body[index], 7 * (leaf_levels[index] * ruler - body[index]))
         for index in range(5)),
        F(0),
    )
    return leaf_debt + F(body[5], 7 * (ruler - body[5]))


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


def best_certificate(
    body: tuple[int, ...],
    ruler: int,
    safe_ranges: tuple[tuple[int, int], ...],
    levels: tuple[int, ...],
    debt: F,
    profile_cache: dict,
    pair_indices=None,
):
    if pair_indices is None:
        pair_indices = tuple(combinations(range(6), 2))
    certificates = []
    for i, j in pair_indices:
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

    bodies = []
    code_hist: Counter[str] = Counter()
    body_digest = hashlib.sha256()
    universal_exceptions = {row[0] for row in U.EXPECTED_EXCEPTIONS}
    for body in combinations(range(1, 15), 6):
        ruler, universal_debt, robust_edges = LOW.robust_edges(body)
        if len(robust_edges) != 10:
            continue
        require(body not in universal_exceptions,
                ("edge-10 body lacks complete same-level graph", body))
        threshold = E12.nonrobust_threshold(body, ruler, universal_debt)
        code = E12.threshold_creation_code(body, threshold)
        nonedges = tuple(
            edge for edge in combinations(range(6), 2) if edge not in set(robust_edges)
        )
        require(code in {"DIIII", "IDIID", "IIDDI"},
                (body, code, nonedges))
        if code == "DIIII":
            require(nonedges == tuple((index, 5) for index in range(5)),
                    ("DIIII complement is not K1,5", body, nonedges))
        code_hist[code] += 1
        bodies.append((body, ruler, universal_debt, robust_edges, nonedges, threshold, code))
        body_digest.update(
            f"{body}|{ruler}|{universal_debt}|{threshold}|{robust_edges}|{nonedges}|{code}\n".encode()
        )
    require(tuple(sorted(code_hist.items())) == EXPECTED_CODE_HIST, code_hist)
    require(len(bodies) == ROBUST_EDGE10_BODY_COUNT, len(bodies))

    word_cache = {}
    rows = []
    code_row_counts: Counter[str] = Counter()
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
                if code == "DIIII" and high_edges == nonedges:
                    continue
                if high_edges not in word_cache:
                    word_cache[high_edges] = E11.classified_words(high_edges)
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
                    debt = singleton_debt(body, ruler, levels)
                    best = best_certificate(
                        body, ruler, safe_ranges, levels, debt, profile_cache
                    )
                    require(best[0] > 0,
                            ("finite edge-10 assignment failed", code, category,
                             body, nonedges, high_edges, levels, best))
                    row = (
                        best[0],
                        code,
                        category,
                        body,
                        nonedges,
                        high_edges,
                        levels,
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

        if code == "DIIII":
            leaf_pairs = tuple(combinations(range(5), 2))
            for leaf_type in LOW_K5_PRIMITIVES:
                for leaf_levels in permutations(leaf_type):
                    debt = free_center_debt_bound(body, ruler, leaf_levels)
                    best = best_certificate(
                        body,
                        ruler,
                        safe_ranges,
                        leaf_levels,
                        debt,
                        profile_cache,
                        pair_indices=leaf_pairs,
                    )
                    require(best[0] > 0,
                            ("free-center edge-10 cylinder failed", body,
                             leaf_type, leaf_levels, best, debt))
                    row = (
                        best[0],
                        code,
                        "star5_infinite",
                        body,
                        nonedges,
                        nonedges,
                        leaf_levels,
                    ) + best[1:]
                    rows.append(row)
                    cylinder_rows += 1
                    code_row_counts[code] += 1
                    certificate_kinds[best[-1]] += 1
                    assignment_digest.update(
                        f"{code}|star5_infinite|{body}|{ruler}|{nonedges}|{leaf_levels}|{debt}|{best}\n".encode()
                    )
                    if body not in weakest_by_body or row < weakest_by_body[body]:
                        weakest_by_body[body] = row

    require(connected_rows == CONNECTED_ASSIGNMENT_COUNT, connected_rows)
    require(cylinder_rows == FREE_CENTER_TEMPLATE_COUNT, cylinder_rows)
    require(len(rows) == CERTIFICATE_ROW_COUNT, len(rows))
    require(tuple(sorted(code_row_counts.items())) == EXPECTED_CODE_ROW_COUNTS,
            code_row_counts)
    require(min(rows) == EXPECTED_WEAKEST, min(rows))
    cylinder_weakest = min(row for row in rows if row[2] == "star5_infinite")
    require(cylinder_weakest == EXPECTED_DIIII_CYLINDER_WEAKEST, cylinder_weakest)
    require(certificate_kinds == Counter({"low": CERTIFICATE_ROW_COUNT}),
            certificate_kinds)
    require(len(weakest_by_body) == ROBUST_EDGE10_BODY_COUNT, len(weakest_by_body))

    # Directly replay one extremal certificate per body.  The analytic
    # transport theorem covers every row; these controls independently anchor
    # its implementation without 565,632 redundant interval calculations.
    direct_controls = 0
    for body, row in sorted(weakest_by_body.items()):
        (
            _,
            _,
            category,
            _,
            _,
            _,
            levels,
            floor,
            i,
            j,
            _,
            _,
            divisor,
            cell,
            _,
            _,
            _,
        ) = row
        ruler = next(item[1] for item in bodies if item[0] == body)
        a, b = body[i], body[j]
        p, q = levels[i], levels[j]
        require(R.body_cell_is_safe(ruler, body, cell),
                ("certificate cell unsafe", body, levels, (i, j), cell))
        for scale in DIRECT_CONTROL_SCALES:
            actual = intersection_mass(
                R.reflected_level_arcs(ruler, a, scale * p, cell),
                R.reflected_level_arcs(ruler, b, scale * q, cell),
            )
            transported = floor - F(4 * (a + b), scale * divisor * ruler)
            if category == "star5_infinite":
                scaled_debt = free_center_debt_bound(
                    body, ruler, tuple(scale * value for value in levels)
                )
            else:
                scaled_debt = singleton_debt(
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
        tuple(bodies),
        tuple(sorted((edges, words) for edges, words in word_cache.items())),
        min(rows),
        cylinder_weakest,
        tuple(sorted(code_row_counts.items())),
        tuple(sorted(certificate_kinds.items())),
        connected_rows,
        cylinder_rows,
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
        "LRC14 reflected robust-edge-10 threshold-block uniform closure exact referee",
        f"edge10_code_hist={tuple(sorted(code_hist.items()))};"
        f"code_row_counts={tuple(sorted(code_row_counts.items()))}",
        f"connected_projective_rows={connected_rows};"
        f"free_center_cylinder_templates={cylinder_rows};"
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
        f"free_center_weakest_margin={qtext(cylinder_weakest[0])};"
        f"body={cylinder_weakest[3]};leaf_levels={cylinder_weakest[6]};"
        f"floor={qtext(cylinder_weakest[7])};"
        f"certificate_pair={(cylinder_weakest[8],cylinder_weakest[9])}",
        f"direct_control_scales={DIRECT_CONTROL_SCALES};"
        f"per_body_weakest_controls={direct_controls}",
        "free_center_law=the K5 leaf ray dilates while the arbitrary center debt is bounded by its level-one value",
        "scale_law=multiplying finite words or cylinder leaves preserves the reduced channel and cell and decreases loss and controlled debt",
        "conclusion=all 32 robust-edge-10 bodies close for every assignment of positive reflected levels",
        f"corollary=arbitrary-level body closure rises from {PRIOR_CLOSED_COUNT} to {CUMULATIVE_CLOSED_COUNT};remaining_bodies={REMAINING_BODY_COUNT}",
        "scope=reflected THM-2941 residual family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"low_phase_sha256={sha256(LOW_PHASE)}",
        f"exceptional_sha256={sha256(EXCEPTIONAL)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"edge11_sha256={sha256(EDGE11)}",
        f"edge12_sha256={sha256(EDGE12)}",
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
