#!/usr/bin/env python3
r"""Uniform arbitrary-level closure for all 26 robust-edge-13 bodies.

Every body whose robust graph has thirteen edges has the same nonrobust-edge
shape: two edges sharing one vertex.  The universal same-level graph is
complete on all 26 bodies, so an uncertified level word must use six distinct
levels.  Moreover every high-phase level pair must lie on one of those two
nonrobust edges.  The high-edge graph is therefore either one edge or the full
two-edge path.

The one-edge case has the two projective types classified by the preceding
``K6-e`` theorem.  The path case is finite as well.  Delete its degree-two
vertex; the other five levels form one of the two classified low-phase
five-cliques ``A,B``.  After choosing the two high neighbours, the deleted
level must be low-adjacent to each of the remaining three levels.  Anchoring at
one of those three vertices restricts its ratio to the finite symmetric
low-ratio bank.  The exact enumeration below is complete and gives eight
projective path types.

There are exactly

    26 * [2*2*2*4! + 8*2*3!] = 7,488

labelled primitive assignments: two one-edge types, a choice of body nonedge,
an orientation and four free labels; or eight path types, two endpoint
orientations and three free labels.  For every assignment an exact body-safe
cell supplies a positive one-pair margin over the actual singleton debt.  All
selected certificates are low-phase.  Multiplying the six levels by an
integer scale preserves the reduced channel and selected cell while decreasing
both transport loss and singleton debt, so the primitive census closes every
scale ray.  Direct reflected-arc checks at scales ``1,2,5`` verify all 22,464
selected certificates.

Thus all 26 robust-edge-13 bodies close for arbitrary positive reflected
levels.  Together with the robust ``K6`` and ``K6-e`` results this closes 2,302
of 3,003 bodies and leaves 701 bodies with at most twelve robust edges.  This
is scoped to the THM-2941 reflected residual, not physical LRC(14).
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
K6_MINUS_EDGE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_k6_minus_edge_uniform_closure_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_k6_minus_adjacent_pair_uniform_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_LOW_PHASE_SHA256 = "416c36f16f7c821feb8d260882711d2717069147b8604a93ba60432785cf1d1c"
EXPECTED_EXCEPTIONAL_SHA256 = "bde992db1edbd9dd744ff22744a1afef79cf4bcc54a2f918793c2603f062df7c"
EXPECTED_UNIVERSAL_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"
EXPECTED_K6_MINUS_EDGE_SHA256 = "3b5752836263df2f6b694c9dfe737e489c22dc2c267144530dcfac74cc972b29"
EXPECTED_SEMANTIC_SHA256 = "4512062600e7915e762b991b7522fd322c4c19e06b934417beffad9ba462e941"

BODY_COUNT = 3003
PRIOR_CLOSED_COUNT = 2276
ROBUST_EDGE_13_BODY_COUNT = 26
CUMULATIVE_CLOSED_COUNT = 2302
REMAINING_BODY_COUNT = 701
ONE_HIGH_ASSIGNMENT_COUNT = 4992
PATH_HIGH_ASSIGNMENT_COUNT = 2496
ASSIGNMENT_COUNT = 7488
DIRECT_CONTROL_SCALES = (1, 2, 5)
DIRECT_CONTROL_COUNT = 22464

SYMMETRIC_LOW_RATIOS = frozenset(
    set(LOW_RATIO_BANK := (
        F(4, 3), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6)
    ))
    | {F(1, ratio) for ratio in LOW_RATIO_BANK}
)

# (normalized set, primitive integral set, high center, high endpoints)
PATH_TYPES = (
    (
        (F(1), F(4, 3), F(2), F(8, 3), F(4), F(8)),
        (3, 4, 6, 8, 12, 24),
        3,
        (8, 24),
    ),
    (
        (F(1), F(4, 3), F(2), F(3), F(4), F(6)),
        (3, 4, 6, 9, 12, 18),
        4,
        (9, 18),
    ),
    (
        (F(1), F(3, 2), F(2), F(3), F(9, 2), F(6)),
        (2, 3, 4, 6, 9, 12),
        9,
        (2, 4),
    ),
    (
        (F(1), F(3, 2), F(2), F(3), F(6), F(9)),
        (2, 3, 4, 6, 12, 18),
        18,
        (2, 4),
    ),
    (
        (F(1), F(3, 2), F(2), F(3), F(6), F(12)),
        (2, 3, 4, 6, 12, 24),
        24,
        (2, 3),
    ),
    (
        (F(1), F(3, 2), F(3), F(9, 2), F(6), F(9)),
        (2, 3, 6, 9, 12, 18),
        2,
        (9, 18),
    ),
    (
        (F(1), F(2), F(3), F(4), F(6), F(8)),
        (1, 2, 3, 4, 6, 8),
        8,
        (1, 3),
    ),
    (
        (F(1), F(2), F(4), F(6), F(8), F(12)),
        (1, 2, 4, 6, 8, 12),
        1,
        (8, 12),
    ),
)
EXPECTED_WEAKEST = (
    F(6022958680656216101363, 93345307647398389204440),
    "path",
    7,
    (1, 2, 3, 6, 7, 13),
    ((3, 5), (4, 5)),
    (2, 6, 4, 12, 8, 1),
    F(1007, 15288),
    0,
    2,
    1,
    2,
    2,
    7013,
    1,
    631,
    "low",
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
    (K6_MINUS_EDGE, EXPECTED_K6_MINUS_EDGE_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("robust_edge13_base", BASE)
LOW = import_module("robust_edge13_low_phase", LOW_PHASE)
X = import_module("robust_edge13_optimizer", EXCEPTIONAL)
U = import_module("robust_edge13_universal", UNIVERSAL)
K6E = import_module("robust_edge13_k6_minus_edge", K6_MINUS_EDGE)


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


def classify_high_path():
    """Insert the high-path center into each classified low K5."""
    rows = []
    for clique_index, clique in enumerate(LOW.EXPECTED_MAX_CLIQUES):
        for endpoints in combinations(clique, 2):
            core = tuple(value for value in clique if value not in endpoints)
            require(len(core) == 3, (clique, endpoints, core))
            candidates = sorted({core[0] * ratio for ratio in SYMMETRIC_LOW_RATIOS})
            for center in candidates:
                if center in clique:
                    continue
                low_neighbours = tuple(
                    value for value in clique if LOW.low_adjacent(center, value)
                )
                high_neighbours = tuple(
                    value for value in clique if not LOW.low_adjacent(center, value)
                )
                if set(low_neighbours) != set(core) or set(high_neighbours) != set(endpoints):
                    continue
                union = set(clique) | {center}
                minimum = min(union)
                rows.append(
                    (
                        clique_index,
                        endpoints,
                        center,
                        tuple(sorted(value / minimum for value in union)),
                        center / minimum,
                        tuple(sorted(value / minimum for value in endpoints)),
                    )
                )
    return tuple(rows)


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

    require(tuple(K6E.K5_TYPES) == tuple(LOW.EXPECTED_MAX_CLIQUES),
            (K6E.K5_TYPES, LOW.EXPECTED_MAX_CLIQUES))
    path_rows = classify_high_path()
    require(len(path_rows) == len(PATH_TYPES) == 8, path_rows)
    classified = tuple(sorted((row[3], row[4], row[5]) for row in path_rows))
    expected_classified = tuple(
        sorted(
            (
                normalized,
                F(center, primitive[0]) * normalized[0],
                tuple(F(endpoint, primitive[0]) * normalized[0] for endpoint in endpoints),
            )
            for normalized, primitive, center, endpoints in PATH_TYPES
        )
    )
    require(classified == expected_classified, (classified, expected_classified))

    # Check each primitive representative and the location of its two high
    # edges; gcd one makes every integral realization an integer dilation.
    for normalized, primitive, center, endpoints in PATH_TYPES:
        scale = F(primitive[0], normalized[0])
        require(tuple(scale * value for value in normalized) == primitive,
                (normalized, primitive, scale))
        require(gcd(*primitive) == 1, ("nonprimitive path type", primitive))
        actual_high = tuple(
            pair for pair in combinations(primitive, 2)
            if not LOW.low_adjacent(F(pair[0]), F(pair[1]))
        )
        require(
            set(actual_high) == {(min(center, endpoint), max(center, endpoint))
                                 for endpoint in endpoints},
            (primitive, center, endpoints, actual_high),
        )

    bodies = []
    body_digest = hashlib.sha256()
    universal_exceptions = {row[0] for row in U.EXPECTED_EXCEPTIONS}
    for body in combinations(range(1, 15), 6):
        ruler, universal_debt, robust_edges = LOW.robust_edges(body)
        if len(robust_edges) != 13:
            continue
        nonedges = tuple(
            edge for edge in combinations(range(6), 2) if edge not in set(robust_edges)
        )
        require(body not in universal_exceptions,
                ("edge-13 body lacks complete same-level graph", body))
        require(len(nonedges) == 2 and len(set(nonedges[0]) & set(nonedges[1])) == 1,
                ("edge-13 nonrobust shape changed", body, nonedges))
        bodies.append((body, ruler, universal_debt, robust_edges, nonedges))
        body_digest.update(f"{body}|{ruler}|{universal_debt}|{robust_edges}|{nonedges}\n".encode())
    require(len(bodies) == ROBUST_EDGE_13_BODY_COUNT, len(bodies))

    assignment_rows = []
    category_counts: Counter[str] = Counter()
    certificate_kinds: Counter[str] = Counter()
    certificate_floors: Counter[F] = Counter()
    direct_controls = 0
    assignment_digest = hashlib.sha256()

    def certify(
        category: str,
        type_index: int,
        body: tuple[int, ...],
        ruler: int,
        nonedges: tuple[tuple[int, int], ...],
        levels: tuple[int, ...],
        safe_ranges: tuple[tuple[int, int], ...],
        profile_cache: dict,
    ) -> None:
        nonlocal direct_controls
        debt = singleton_debt(body, ruler, levels)
        best = best_certificate(body, ruler, safe_ranges, levels, debt, profile_cache)
        require(best[0] > 0,
                ("primitive edge-13 assignment failed", category, type_index,
                 body, nonedges, levels, best))
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
            transported = floor - F(4 * (a + b), scale * divisor * ruler)
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
            type_index,
            body,
            nonedges,
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
            f"{category}|{type_index}|{body}|{ruler}|{nonedges}|{levels}|{debt}|{best}\n".encode()
        )

    one_high_types = tuple((row[1], row[2]) for row in K6E.PROJECTIVE_TYPES)
    require(len(one_high_types) == 2, one_high_types)
    for body, ruler, _, robust_edges, nonedges in bodies:
        safe_ruler, safe_ranges = R.safe_cell_ranges(body)
        require(safe_ruler == ruler, (body, safe_ruler, ruler))
        profile_cache = {}
        robust_set = set(robust_edges)

        # One high edge: choose either body nonedge.
        for type_index, (primitive, high_pair) in enumerate(one_high_types):
            remaining_levels = tuple(value for value in primitive if value not in high_pair)
            for chosen_nonedge in nonedges:
                for oriented_high in (high_pair, high_pair[::-1]):
                    for remaining in permutations(remaining_levels):
                        levels_list = [None] * 6
                        levels_list[chosen_nonedge[0]], levels_list[chosen_nonedge[1]] = oriented_high
                        remaining_iterator = iter(remaining)
                        for index in range(6):
                            if levels_list[index] is None:
                                levels_list[index] = next(remaining_iterator)
                        levels = tuple(levels_list)
                        induced_high = tuple(
                            edge for edge in combinations(range(6), 2)
                            if not LOW.low_adjacent(F(levels[edge[0]]), F(levels[edge[1]]))
                        )
                        require(induced_high == (chosen_nonedge,),
                                (body, type_index, chosen_nonedge, levels, induced_high))
                        require(all(LOW.low_adjacent(F(levels[i]), F(levels[j]))
                                    for i, j in robust_set),
                                (body, levels, robust_edges))
                        certify(
                            "one", type_index, body, ruler, nonedges,
                            levels, safe_ranges, profile_cache,
                        )

        # Two high edges: the projective high-path center must occupy the
        # shared body vertex and its endpoints the two other incident slots.
        shared = next(iter(set(nonedges[0]) & set(nonedges[1])))
        endpoint_slots = tuple(
            sorted((set(nonedges[0]) | set(nonedges[1])) - {shared})
        )
        other_slots = tuple(
            index for index in range(6) if index not in (shared,) + endpoint_slots
        )
        for type_index, (_, primitive, center, endpoints) in enumerate(PATH_TYPES):
            remaining_levels = tuple(
                value for value in primitive if value not in (center,) + endpoints
            )
            require(len(remaining_levels) == 3, (primitive, center, endpoints))
            for oriented_endpoints in (endpoints, endpoints[::-1]):
                for remaining in permutations(remaining_levels):
                    levels_list = [None] * 6
                    levels_list[shared] = center
                    levels_list[endpoint_slots[0]], levels_list[endpoint_slots[1]] = oriented_endpoints
                    for index, value in zip(other_slots, remaining):
                        levels_list[index] = value
                    levels = tuple(levels_list)
                    induced_high = tuple(
                        edge for edge in combinations(range(6), 2)
                        if not LOW.low_adjacent(F(levels[edge[0]]), F(levels[edge[1]]))
                    )
                    require(set(induced_high) == set(nonedges),
                            (body, type_index, nonedges, levels, induced_high))
                    require(all(LOW.low_adjacent(F(levels[i]), F(levels[j]))
                                for i, j in robust_set),
                            (body, levels, robust_edges))
                    certify(
                        "path", type_index, body, ruler, nonedges,
                        levels, safe_ranges, profile_cache,
                    )

    require(category_counts == Counter(
        {"one": ONE_HIGH_ASSIGNMENT_COUNT, "path": PATH_HIGH_ASSIGNMENT_COUNT}
    ), category_counts)
    require(len(assignment_rows) == ASSIGNMENT_COUNT, len(assignment_rows))
    require(direct_controls == DIRECT_CONTROL_COUNT, direct_controls)
    require(min(assignment_rows) == EXPECTED_WEAKEST, min(assignment_rows))
    require(certificate_kinds == Counter({"low": ASSIGNMENT_COUNT}), certificate_kinds)
    require(PRIOR_CLOSED_COUNT + len(bodies) == CUMULATIVE_CLOSED_COUNT,
            (PRIOR_CLOSED_COUNT, len(bodies), CUMULATIVE_CLOSED_COUNT))
    require(BODY_COUNT - CUMULATIVE_CLOSED_COUNT == REMAINING_BODY_COUNT,
            (BODY_COUNT, CUMULATIVE_CLOSED_COUNT, REMAINING_BODY_COUNT))

    semantic_payload = (
        path_rows,
        PATH_TYPES,
        tuple(bodies),
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
        weakest_type,
        weakest_body,
        weakest_nonedges,
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
        "LRC14 reflected robust-edge-13 uniform closure exact referee",
        f"nonrobust_shape=two_adjacent_edges;body_count={len(bodies)};"
        f"path_classification_rows={len(path_rows)};projective_path_types={PATH_TYPES}",
        f"assignment_categories={tuple(sorted(category_counts.items()))};"
        f"total_primitive_assignments={len(assignment_rows)};"
        f"certificate_kinds={tuple(sorted(certificate_kinds.items()))};"
        f"certificate_floor_values={len(certificate_floors)}",
        f"weakest_margin={qtext(weakest_margin)};category={weakest_category};"
        f"type={weakest_type};body={weakest_body};nonrobust_edges={weakest_nonedges};"
        f"levels={weakest_levels};floor={qtext(weakest_floor)};"
        f"certificate_pair={(weakest_i,weakest_j)};"
        f"reduced_channel={(weakest_P,weakest_Q)};gcd={weakest_divisor};"
        f"j={weakest_cell};lift={weakest_lift};endpoint_radius={weakest_radius};"
        f"kind={weakest_kind}",
        f"direct_control_scales={DIRECT_CONTROL_SCALES};direct_controls={direct_controls}",
        "scale_law=multiplying all levels preserves the reduced channel and cell, decreases transport loss and singleton debt",
        "conclusion=all 26 robust-edge-13 bodies close for every assignment of positive reflected levels",
        f"corollary=arbitrary-level body closure rises from {PRIOR_CLOSED_COUNT} to {CUMULATIVE_CLOSED_COUNT};remaining_bodies={REMAINING_BODY_COUNT}",
        "scope=reflected THM-2941 residual family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"low_phase_sha256={sha256(LOW_PHASE)}",
        f"exceptional_sha256={sha256(EXCEPTIONAL)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"k6_minus_edge_sha256={sha256(K6_MINUS_EDGE)}",
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
