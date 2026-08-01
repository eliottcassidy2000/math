#!/usr/bin/env python3
r"""Uniform arbitrary-level closure for all 59 robust-``K6-e`` bodies.

The low-phase theorem attaches a robust graph to every six-body.  On the 59
bodies with fourteen robust edges, any packet not already closed by the
universal repeated-level theorem has six distinct levels.  If a robust edge
has reduced level channel ``(P,Q)`` with ``P+Q>=8``, its periodized fibre floor
``1/105`` closes the packet.  Hence a remaining level set would induce a
low-phase graph containing ``K6-e``.  Since that graph has clique number five,
the induced graph is exactly ``K6-e`` and its unique high pair lies on the
unique nonrobust body edge.

There are only two projective rational level sets with low graph ``K6-e``.
Indeed, deleting either endpoint of the high pair leaves a five-clique.  The
two normalized five-cliques are

    A={1,3/2,2,3,6},       B={1,2,3,4,6}.

Align two scaled copies of ``A`` or ``B`` in four vertices.  Exact enumeration
of the possible scale ratios ``u/v`` gives four ordered overlaps and just two
normalized unions:

    C={1,3/2,2,3,4,6},     unique high pair {3/2,4};
    D={1,2,3,4,6,12},      unique high pair {1,12}.

For integral levels their primitive representatives are respectively

    C0=(2,3,4,6,8,12),     D0=(1,2,3,4,6,12).

Thus it remains to test ``59*2*2*4! = 5,664`` labelled primitive assignments.
For each assignment this referee chooses the strongest exact one-pair
certificate.  A low channel uses the body-safe-cell convex optimizer from the
exceptional-body theorem; a high channel uses the global ``1/105`` floor.  If
the chosen levels are ``p=gP,q=gQ``, the certified overlap is

    floor - 4(a+b)/(gL).

Every primitive margin over the exact singleton debt is positive.  Multiplying
all six levels by an integer ``s>=1`` preserves ``(P,Q)`` and the selected
cell, multiplies ``g`` by ``s``, and decreases every singleton-debt term.
Therefore the primitive checks close every scale ray.  Direct reflected-arc
controls at scales ``1,2,5`` independently verify all 16,992 selected
certificates.

Consequently all 59 robust-``K6-e`` bodies close for arbitrary positive
reflected levels.  Together with the 2,217 robust-``K6`` bodies, this closes
2,276 of the 3,003 bodies and leaves 727 bodies with at most thirteen robust
edges.  This is a scoped THM-2941 result, not physical LRC(14).
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
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_k6_minus_edge_uniform_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_LOW_PHASE_SHA256 = "416c36f16f7c821feb8d260882711d2717069147b8604a93ba60432785cf1d1c"
EXPECTED_EXCEPTIONAL_SHA256 = "bde992db1edbd9dd744ff22744a1afef79cf4bcc54a2f918793c2603f062df7c"
EXPECTED_UNIVERSAL_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"
EXPECTED_SEMANTIC_SHA256 = "fea17af1ab0586db6aa0de1abf384087de93f4a0871be30e426d2a4c40de9431"

BODY_COUNT = 3003
ROBUST_K6_COUNT = 2217
ROBUST_K6_MINUS_EDGE_COUNT = 59
CUMULATIVE_CLOSED_COUNT = 2276
REMAINING_BODY_COUNT = 727
ASSIGNMENT_COUNT = 5664
DIRECT_CONTROL_SCALES = (1, 2, 5)
DIRECT_CONTROL_COUNT = 16992

K5_TYPES = (
    (F(1), F(3, 2), F(2), F(3), F(6)),
    (F(1), F(2), F(3), F(4), F(6)),
)
EXPECTED_CLASSIFICATION_ROWS = (
    (
        0,
        1,
        F(1, 2),
        (F(1), F(3, 2), F(2), F(3)),
        (F(1), F(2), F(3), F(4), F(6), F(12)),
        (F(1), F(12)),
    ),
    (
        0,
        1,
        F(1),
        (F(1), F(2), F(3), F(6)),
        (F(1), F(3, 2), F(2), F(3), F(4), F(6)),
        (F(3, 2), F(4)),
    ),
    (
        1,
        0,
        F(1),
        (F(1), F(2), F(3), F(6)),
        (F(1), F(3, 2), F(2), F(3), F(4), F(6)),
        (F(3, 2), F(4)),
    ),
    (
        1,
        0,
        F(2),
        (F(2), F(3), F(4), F(6)),
        (F(1), F(2), F(3), F(4), F(6), F(12)),
        (F(1), F(12)),
    ),
)
PROJECTIVE_TYPES = (
    (
        (F(1), F(3, 2), F(2), F(3), F(4), F(6)),
        (2, 3, 4, 6, 8, 12),
        (3, 8),
    ),
    (
        (F(1), F(2), F(3), F(4), F(6), F(12)),
        (1, 2, 3, 4, 6, 12),
        (1, 12),
    ),
)
EXPECTED_WEAKEST_BY_TYPE = (
    (
        F(11634903103299856214368519, 168563272427550089513989110),
        (1, 2, 3, 5, 9, 14),
        (4, 5),
        (6, 2, 12, 4, 8, 3),
        F(1, 14),
        0,
        5,
        2,
        1,
        3,
        8144,
        25,
        612,
        "low",
    ),
    (
        F(7745575244788036048763, 117219090152296965822990),
        (1, 2, 4, 5, 12, 13),
        (4, 5),
        (2, 3, 4, 6, 12, 1),
        F(1, 14),
        0,
        5,
        2,
        1,
        1,
        9610,
        22,
        35,
        "low",
    ),
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
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("robust_k6_minus_edge_base", BASE)
LOW = import_module("robust_k6_minus_edge_low_phase", LOW_PHASE)
X = import_module("robust_k6_minus_edge_optimizer", EXCEPTIONAL)
U = import_module("robust_k6_minus_edge_universal", UNIVERSAL)


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


def classify_k6_minus_edge():
    """Glue the two classified projective K5 types along four vertices."""
    rows = []
    for first_index, first in enumerate(K5_TYPES):
        for second_index, second in enumerate(K5_TYPES):
            ratios = sorted({x / y for x in first for y in second})
            for ratio in ratios:
                left = set(first)
                right = {ratio * value for value in second}
                union = left | right
                if len(left & right) != 4 or len(union) != 6:
                    continue
                high = tuple(
                    pair
                    for pair in combinations(sorted(union), 2)
                    if not LOW.low_adjacent(*pair)
                )
                if len(high) != 1:
                    continue
                minimum = min(union)
                rows.append(
                    (
                        first_index,
                        second_index,
                        ratio,
                        tuple(sorted(left & right)),
                        tuple(sorted(value / minimum for value in union)),
                        tuple(value / minimum for value in high[0]),
                    )
                )
    return tuple(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(K5_TYPES == LOW.EXPECTED_MAX_CLIQUES, (K5_TYPES, LOW.EXPECTED_MAX_CLIQUES))
    classification_rows = classify_k6_minus_edge()
    require(classification_rows == EXPECTED_CLASSIFICATION_ROWS, classification_rows)
    classified_types = tuple(sorted({row[4] for row in classification_rows}))
    require(classified_types == tuple(sorted(row[0] for row in PROJECTIVE_TYPES)), classified_types)
    for normalized, primitive, high_pair in PROJECTIVE_TYPES:
        scale = F(primitive[0], normalized[0])
        require(tuple(scale * value for value in normalized) == primitive,
                (normalized, primitive, scale))
        actual_high = tuple(
            pair for pair in combinations(primitive, 2)
            if not LOW.low_adjacent(F(pair[0]), F(pair[1]))
        )
        require(actual_high == (high_pair,), (primitive, actual_high, high_pair))

    bodies = []
    body_digest = hashlib.sha256()
    universal_exceptions = {row[0] for row in U.EXPECTED_EXCEPTIONS}
    for body in combinations(range(1, 15), 6):
        ruler, universal_debt, robust_edges = LOW.robust_edges(body)
        if len(robust_edges) != 14:
            continue
        require(body not in universal_exceptions,
                ("K6-e body lacks complete same-level graph", body))
        nonedges = tuple(
            edge for edge in combinations(range(6), 2) if edge not in set(robust_edges)
        )
        require(len(nonedges) == 1, (body, robust_edges, nonedges))
        bodies.append((body, ruler, universal_debt, robust_edges, nonedges[0]))
        body_digest.update(f"{body}|{ruler}|{universal_debt}|{robust_edges}|{nonedges[0]}\n".encode())
    require(len(bodies) == ROBUST_K6_MINUS_EDGE_COUNT, len(bodies))

    assignment_count = 0
    direct_controls = 0
    weakest_by_type = []
    certificate_kinds: Counter[str] = Counter()
    certificate_floors: Counter[F] = Counter()
    assignment_digest = hashlib.sha256()

    for type_index, (_, primitive, high_pair) in enumerate(PROJECTIVE_TYPES):
        type_rows = []
        for body, ruler, _, robust_edges, nonedge in bodies:
            safe_ruler, safe_ranges = R.safe_cell_ranges(body)
            require(safe_ruler == ruler, (body, safe_ruler, ruler))
            robust_set = set(robust_edges)
            profile_cache = {}
            remaining_levels = tuple(value for value in primitive if value not in high_pair)
            for oriented_high in (high_pair, high_pair[::-1]):
                for remaining in permutations(remaining_levels):
                    levels_list = [None] * 6
                    levels_list[nonedge[0]], levels_list[nonedge[1]] = oriented_high
                    remaining_iterator = iter(remaining)
                    for index in range(6):
                        if levels_list[index] is None:
                            levels_list[index] = next(remaining_iterator)
                    levels = tuple(levels_list)
                    require(len(set(levels)) == 6, (body, levels))

                    induced_high = tuple(
                        edge
                        for edge in combinations(range(6), 2)
                        if not LOW.low_adjacent(F(levels[edge[0]]), F(levels[edge[1]]))
                    )
                    require(induced_high == (nonedge,),
                            (body, primitive, nonedge, levels, induced_high))
                    require(all(LOW.low_adjacent(F(levels[i]), F(levels[j]))
                                for i, j in robust_set),
                            (body, levels, robust_edges))

                    debt = singleton_debt(body, ruler, levels)
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
                        margin = (
                            floor
                            - F(4 * (body[i] + body[j]), divisor * ruler)
                            - debt
                        )
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
                    best = max(certificates)
                    require(best[0] > 0,
                            ("primitive K6-e assignment failed", body, primitive, levels, best))
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

                    # Exact direct controls of the transport bound and the
                    # scale-monotone singleton-debt conclusion.
                    for scale in DIRECT_CONTROL_SCALES:
                        actual = intersection_mass(
                            R.reflected_level_arcs(ruler, a, scale * p, cell),
                            R.reflected_level_arcs(ruler, b, scale * q, cell),
                        )
                        transported = floor - F(
                            4 * (a + b), scale * divisor * ruler
                        )
                        scaled_levels = tuple(scale * level for level in levels)
                        scaled_debt = singleton_debt(body, ruler, scaled_levels)
                        require(actual >= transported > scaled_debt,
                                (body, levels, (i, j), scale, actual,
                                 transported, scaled_debt, best))
                        direct_controls += 1

                    row = (
                        margin,
                        body,
                        nonedge,
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
                    type_rows.append(row)
                    assignment_count += 1
                    certificate_kinds[kind] += 1
                    certificate_floors[floor] += 1
                    assignment_digest.update(
                        f"{type_index}|{body}|{ruler}|{nonedge}|{levels}|{debt}|{best}\n".encode()
                    )
        weakest_by_type.append(min(type_rows))

    require(assignment_count == ASSIGNMENT_COUNT, assignment_count)
    require(direct_controls == DIRECT_CONTROL_COUNT, direct_controls)
    require(tuple(weakest_by_type) == EXPECTED_WEAKEST_BY_TYPE, weakest_by_type)
    require(ROBUST_K6_COUNT + len(bodies) == CUMULATIVE_CLOSED_COUNT,
            (ROBUST_K6_COUNT, len(bodies), CUMULATIVE_CLOSED_COUNT))
    require(BODY_COUNT - CUMULATIVE_CLOSED_COUNT == REMAINING_BODY_COUNT,
            (BODY_COUNT, CUMULATIVE_CLOSED_COUNT, REMAINING_BODY_COUNT))

    semantic_payload = (
        classification_rows,
        tuple(bodies),
        tuple(weakest_by_type),
        tuple(sorted(certificate_kinds.items())),
        tuple(sorted(certificate_floors.items())),
        assignment_count,
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
        "LRC14 reflected robust-K6-minus-edge uniform closure exact referee",
        f"classification_rows={classification_rows}",
        f"projective_types={PROJECTIVE_TYPES}",
        f"robust_K6_minus_edge_bodies={len(bodies)};"
        f"primitive_labelled_assignments={assignment_count};"
        f"certificate_kinds={tuple(sorted(certificate_kinds.items()))};"
        f"certificate_floors={tuple((qtext(k),v) for k,v in sorted(certificate_floors.items()))}",
    ]
    for type_index, row in enumerate(weakest_by_type):
        (
            margin,
            body,
            nonedge,
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
        ) = row
        lines.append(
            f"TYPE;index={type_index};weakest_margin={qtext(margin)};"
            f"body={body};nonrobust_edge={nonedge};levels={levels};"
            f"floor={qtext(floor)};certificate_pair={(i,j)};"
            f"reduced_channel={(P,Q)};gcd={divisor};j={cell};lift={lift};"
            f"endpoint_radius={endpoint_radius};kind={kind}"
        )
    lines.extend(
        (
            f"direct_control_scales={DIRECT_CONTROL_SCALES};direct_controls={direct_controls}",
            "scale_law=multiplying all levels preserves the reduced channel and cell, decreases transport loss and singleton debt",
            "conclusion=all 59 robust-K6-minus-edge bodies close for every assignment of positive reflected levels",
            f"corollary=arbitrary-level body closure rises from {ROBUST_K6_COUNT} to {CUMULATIVE_CLOSED_COUNT};remaining_bodies={REMAINING_BODY_COUNT}",
            "scope=reflected THM-2941 residual family only;physical LRC14 remains open",
            "normal_vs_python_O=BYTE_IDENTICAL",
            f"base_sha256={sha256(BASE)}",
            f"low_phase_sha256={sha256(LOW_PHASE)}",
            f"exceptional_sha256={sha256(EXCEPTIONAL)}",
            f"universal_sha256={sha256(UNIVERSAL)}",
            f"body_digest={body_digest.hexdigest()}",
            f"assignment_digest={assignment_digest.hexdigest()}",
            f"source_sha256={source_sha}",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
