#!/usr/bin/env python3
r"""Threshold-shape closure of all 37 robust-edge-12 bodies.

The robust graph is not an arbitrary graph.  For a body ``E``, ruler ``L``,
and universal debt ``eps(E)``, put

    theta(E) = L/4 * (1/105-eps(E)).

The robust condition is exactly ``a+b<theta(E)``.  Thus the nonrobust
complement is the upper sum-threshold graph ``a+b>=theta(E)``.  On every
induced set of sorted labels, either the least label is isolated (when its sum
with the greatest is below theta) or the greatest label is dominating.  This
gives an exact extreme peel by letters ``I,D`` and explains the severe graph
shape restriction without a general graph atlas.

At thirteen robust edges the complement has creation code ``IIIDI``.  At
twelve robust edges it has only two shapes:

* 17 bodies have code ``IIDII`` and complement ``K1,3``;
* 20 bodies have code ``IIIDD`` and complement ``K3``.

An uncertified packet has six distinct levels and every high-phase level pair
lies in the nonrobust complement.  One-high-edge and two-high-edge-path level
sets were classified by the preceding two theorems.  For the full three-edge
cases this referee proves:

* no projective level set has high graph exactly ``K3``;
* exactly 26 projective level sets have high graph exactly ``K1,3``.

For the first claim, normalize the three low core levels to one of the 13
low-phase three-cliques.  Each proposed high vertex must belong to the finite
symmetric low-neighbour bank of the least core vertex; no pairwise-high triple
exists.  For the second, delete the high-star center, choose its three high
neighbours in one of the two classified five-cliques, and enumerate the center
from the symmetric low-neighbour bank of either remaining core vertex.

The remaining exact census has 26,616 primitive labelled assignments:

    one high edge:       10,656,
    two-edge high path:  10,656,
    full high star:       5,304.

Every assignment has a positive body-safe one-pair margin over its actual
singleton debt; every selected certificate is low-phase.  Integer dilation
preserves the reduced channel and cell while decreasing transport loss and
singleton debt.  Direct reflected-arc controls at scales ``1,2,5`` verify all
79,848 chosen certificates.  Hence all 37 edge-12 bodies close for arbitrary
positive reflected levels.  The cumulative body count becomes 2,339/3,003,
leaving 664 bodies with at most eleven robust edges.  This is scoped to the
THM-2941 reflected residual, not physical LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
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
K6_MINUS_EDGE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_k6_minus_edge_uniform_closure_thm2941.py"
)
EDGE13 = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_k6_minus_adjacent_pair_uniform_closure_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_robust_edge12_threshold_shapes_uniform_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_LOW_PHASE_SHA256 = "b2418dfda1b48257d1f7582d4ea977203a26f88885e13946bc100ccf264c9ce1"
EXPECTED_EXCEPTIONAL_SHA256 = "15c13d721ee695c2a69bc386d5c12eba0382b52f9eabfe052f2fe5fda03c7bc4"
EXPECTED_UNIVERSAL_SHA256 = "dc6f23a201e817dd9134e8660d35e83d3053c67d26fc271ce3eae07f0f857689"
EXPECTED_K6_MINUS_EDGE_SHA256 = "1a18c91803185033bd1faf8005a88ba817cc93edb0ff76a27287427e2299e97a"
EXPECTED_EDGE13_SHA256 = "027cf455babba5eb7647ffabb6e678e091b3a9c8426fa51703bec05cd5ab1586"
EXPECTED_SEMANTIC_SHA256 = "cbb549b7c81f5cd62b86b0b99d44555d5bb798ca19e34a225cf339a0ca6ca1ca"

BODY_COUNT = 3003
PRIOR_CLOSED_COUNT = 2302
ROBUST_EDGE12_BODY_COUNT = 37
CUMULATIVE_CLOSED_COUNT = 2339
REMAINING_BODY_COUNT = 664
EXPECTED_EDGE12_CODE_HIST = (("IIDII", 17), ("IIIDD", 20))
EXPECTED_THRESHOLD_CODE_COUNT = 29
EXPECTED_NORMALIZED_K3_COUNT = 13
EXPECTED_HIGH_TRIANGLE_TYPE_COUNT = 0
EXPECTED_HIGH_STAR_TYPE_COUNT = 26
ONE_ASSIGNMENT_COUNT = 10656
PATH_ASSIGNMENT_COUNT = 10656
STAR_ASSIGNMENT_COUNT = 5304
ASSIGNMENT_COUNT = 26616
DIRECT_CONTROL_SCALES = (1, 2, 5)
DIRECT_CONTROL_COUNT = 79848
EXPECTED_WEAKEST = (
    F(1610277377987793354881, 25313591961935770016160),
    "path",
    7,
    (1, 2, 4, 5, 8, 13),
    ((2, 5), (3, 5), (4, 5)),
    (6, 2, 12, 4, 8, 1),
    F(137, 2080),
    1,
    3,
    1,
    2,
    2,
    6679,
    1,
    601,
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
    (K6_MINUS_EDGE, EXPECTED_K6_MINUS_EDGE_SHA256),
    (EDGE13, EXPECTED_EDGE13_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("robust_edge12_base", BASE)
LOW = import_module("robust_edge12_low_phase", LOW_PHASE)
X = import_module("robust_edge12_optimizer", EXCEPTIONAL)
U = import_module("robust_edge12_universal", UNIVERSAL)
K6E = import_module("robust_edge12_k6_minus_edge", K6_MINUS_EDGE)
E13 = import_module("robust_edge12_edge13", EDGE13)


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


def nonrobust_threshold(body: tuple[int, ...], ruler: int, debt: F) -> F:
    return F(ruler, 4) * (LOW.FIBER_FLOOR - debt)


def threshold_creation_code(body: tuple[int, ...], threshold: F) -> str:
    """Delete an isolated minimum or a dominating maximum, five times."""
    vertices = list(range(6))
    code = []
    while len(vertices) > 1:
        least = vertices[0]
        greatest = vertices[-1]
        if body[least] + body[greatest] < threshold:
            require(all(body[least] + body[index] < threshold for index in vertices[1:]),
                    ("least is not isolated", body, threshold, vertices))
            code.append("I")
            vertices.pop(0)
        else:
            require(all(body[greatest] + body[index] >= threshold for index in vertices[:-1]),
                    ("greatest is not dominating", body, threshold, vertices))
            code.append("D")
            vertices.pop()
    return "".join(code)


def normalized_low_k3s():
    return tuple(
        row
        for row in combinations(LOW.NORMALIZED_VERTICES, 3)
        if row[0] == 1 and LOW.is_low_clique(row)
    )


def classify_high_triangle(k3_types):
    """Attach three pairwise-high vertices, each low to a low K3 core."""
    rows = []
    for core_index, core in enumerate(k3_types):
        candidates = tuple(
            sorted(
                {
                    core[0] * ratio
                    for ratio in E13.SYMMETRIC_LOW_RATIOS
                    if core[0] * ratio not in core
                    and all(LOW.low_adjacent(core[0] * ratio, value) for value in core)
                }
            )
        )
        for high_vertices in combinations(candidates, 3):
            if not all(
                not LOW.low_adjacent(*pair)
                for pair in combinations(high_vertices, 2)
            ):
                continue
            union = set(core) | set(high_vertices)
            if len(union) != 6:
                continue
            minimum = min(union)
            rows.append(
                (
                    core_index,
                    core,
                    high_vertices,
                    tuple(sorted(value / minimum for value in union)),
                )
            )
    return tuple(rows)


def classify_high_star():
    """Insert a center high to three vertices and low to the other two of a K5."""
    rows = []
    for clique_index, clique in enumerate(LOW.EXPECTED_MAX_CLIQUES):
        for endpoints in combinations(clique, 3):
            core = tuple(value for value in clique if value not in endpoints)
            require(len(core) == 2, (clique, endpoints, core))
            candidates = sorted({
                core[0] * ratio for ratio in E13.SYMMETRIC_LOW_RATIOS
            })
            for center in candidates:
                if center in clique:
                    continue
                if not all(LOW.low_adjacent(center, value) for value in core):
                    continue
                if not all(not LOW.low_adjacent(center, value) for value in endpoints):
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


def primitive_star_type(row):
    normalized, center, endpoints = row[3], row[4], row[5]
    denominator = lcm(*(value.denominator for value in normalized))
    integers = tuple(int(value * denominator) for value in normalized)
    common = gcd(*integers)
    primitive = tuple(value // common for value in integers)
    scale = F(denominator, common)
    primitive_center = int(center * scale)
    primitive_endpoints = tuple(int(value * scale) for value in endpoints)
    require(gcd(*primitive) == 1, ("nonprimitive star type", row, primitive))
    require(tuple(F(value, scale) for value in primitive) == normalized,
            (row, primitive, scale))
    return normalized, primitive, primitive_center, primitive_endpoints


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

    k3_types = normalized_low_k3s()
    require(len(k3_types) == EXPECTED_NORMALIZED_K3_COUNT, k3_types)
    triangle_rows = classify_high_triangle(k3_types)
    require(len(triangle_rows) == EXPECTED_HIGH_TRIANGLE_TYPE_COUNT, triangle_rows)
    star_rows = classify_high_star()
    require(len(star_rows) == EXPECTED_HIGH_STAR_TYPE_COUNT, star_rows)
    star_types = tuple(sorted({primitive_star_type(row) for row in star_rows}))
    require(len(star_types) == len(star_rows) == EXPECTED_HIGH_STAR_TYPE_COUNT,
            (len(star_types), len(star_rows), star_rows))
    for normalized, primitive, center, endpoints in star_types:
        actual_high = tuple(
            pair for pair in combinations(primitive, 2)
            if not LOW.low_adjacent(F(pair[0]), F(pair[1]))
        )
        require(
            set(actual_high) == {
                (min(center, endpoint), max(center, endpoint))
                for endpoint in endpoints
            },
            (normalized, primitive, center, endpoints, actual_high),
        )

    bodies = []
    code_hist: Counter[str] = Counter()
    edge12_code_hist: Counter[str] = Counter()
    body_digest = hashlib.sha256()
    universal_exceptions = {row[0] for row in U.EXPECTED_EXCEPTIONS}
    for body in combinations(range(1, 15), 6):
        ruler, debt, robust_edges = LOW.robust_edges(body)
        threshold = nonrobust_threshold(body, ruler, debt)
        expected_robust = tuple(
            edge for edge in combinations(range(6), 2)
            if body[edge[0]] + body[edge[1]] < threshold
        )
        require(robust_edges == expected_robust,
                ("scalar robust threshold failed", body, threshold,
                 robust_edges, expected_robust))
        code = threshold_creation_code(body, threshold)
        code_hist[code] += 1
        if len(robust_edges) != 12:
            continue
        require(body not in universal_exceptions,
                ("edge-12 body lacks complete same-level graph", body))
        nonedges = tuple(
            edge for edge in combinations(range(6), 2) if edge not in set(robust_edges)
        )
        degrees = tuple(sum(index in edge for edge in nonedges) for index in range(6))
        if sorted(degrees, reverse=True) == [3, 1, 1, 1, 0, 0]:
            shape = "star"
            require(code == "IIDII", (body, nonedges, code))
        elif sorted(degrees, reverse=True) == [2, 2, 2, 0, 0, 0]:
            shape = "triangle"
            require(code == "IIIDD", (body, nonedges, code))
        else:
            raise RuntimeError(("unexpected edge-12 threshold shape", body, nonedges, degrees, code))
        edge12_code_hist[code] += 1
        bodies.append((body, ruler, debt, robust_edges, nonedges, threshold, code, shape))
        body_digest.update(
            f"{body}|{ruler}|{debt}|{threshold}|{robust_edges}|{nonedges}|{code}|{shape}\n".encode()
        )
    require(len(code_hist) == EXPECTED_THRESHOLD_CODE_COUNT, len(code_hist))
    require(tuple(sorted(edge12_code_hist.items())) == EXPECTED_EDGE12_CODE_HIST,
            edge12_code_hist)
    require(len(bodies) == ROBUST_EDGE12_BODY_COUNT, len(bodies))

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
                ("primitive edge-12 assignment failed", category, type_index,
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
    path_types = tuple((row[1], row[2], row[3]) for row in E13.PATH_TYPES)
    require((len(one_high_types), len(path_types)) == (2, 8),
            (one_high_types, path_types))
    for body, ruler, _, robust_edges, nonedges, _, _, shape in bodies:
        safe_ruler, safe_ranges = R.safe_cell_ranges(body)
        require(safe_ruler == ruler, (body, safe_ruler, ruler))
        profile_cache = {}
        robust_set = set(robust_edges)

        # One high edge: choose one of the three nonrobust edges.
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

        # Two high edges.  Every pair of complement edges in K3 or K1,3 is
        # adjacent, so each choice realizes one of the eight path types.
        for chosen_edges in combinations(nonedges, 2):
            common = set(chosen_edges[0]) & set(chosen_edges[1])
            require(len(common) == 1, (body, nonedges, chosen_edges))
            shared = next(iter(common))
            endpoint_slots = tuple(
                sorted((set(chosen_edges[0]) | set(chosen_edges[1])) - {shared})
            )
            other_slots = tuple(
                index for index in range(6) if index not in (shared,) + endpoint_slots
            )
            for type_index, (primitive, center, endpoints) in enumerate(path_types):
                remaining_levels = tuple(
                    value for value in primitive if value not in (center,) + endpoints
                )
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
                        require(set(induced_high) == set(chosen_edges),
                                (body, type_index, chosen_edges, levels, induced_high))
                        require(all(LOW.low_adjacent(F(levels[i]), F(levels[j]))
                                    for i, j in robust_set),
                                (body, levels, robust_edges))
                        certify(
                            "path", type_index, body, ruler, nonedges,
                            levels, safe_ranges, profile_cache,
                        )

        # A full high triangle is impossible by the projective classification.
        # On a star complement, enumerate the 26 full high-star types.
        if shape == "star":
            degrees = {index: sum(index in edge for edge in nonedges) for index in range(6)}
            shared = next(index for index, degree in degrees.items() if degree == 3)
            endpoint_slots = tuple(sorted(index for index, degree in degrees.items() if degree == 1))
            core_slots = tuple(sorted(index for index, degree in degrees.items() if degree == 0))
            for type_index, (_, primitive, center, endpoints) in enumerate(star_types):
                remaining_levels = tuple(
                    value for value in primitive if value not in (center,) + endpoints
                )
                require(len(remaining_levels) == 2,
                        (primitive, center, endpoints, remaining_levels))
                for endpoint_values in permutations(endpoints):
                    for core_values in permutations(remaining_levels):
                        levels_list = [None] * 6
                        levels_list[shared] = center
                        for index, value in zip(endpoint_slots, endpoint_values):
                            levels_list[index] = value
                        for index, value in zip(core_slots, core_values):
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
                            "star", type_index, body, ruler, nonedges,
                            levels, safe_ranges, profile_cache,
                        )

    require(category_counts == Counter({
        "one": ONE_ASSIGNMENT_COUNT,
        "path": PATH_ASSIGNMENT_COUNT,
        "star": STAR_ASSIGNMENT_COUNT,
    }), category_counts)
    require(len(assignment_rows) == ASSIGNMENT_COUNT, len(assignment_rows))
    require(direct_controls == DIRECT_CONTROL_COUNT, direct_controls)
    require(min(assignment_rows) == EXPECTED_WEAKEST, min(assignment_rows))
    require(certificate_kinds == Counter({"low": ASSIGNMENT_COUNT}), certificate_kinds)
    require(PRIOR_CLOSED_COUNT + len(bodies) == CUMULATIVE_CLOSED_COUNT,
            (PRIOR_CLOSED_COUNT, len(bodies), CUMULATIVE_CLOSED_COUNT))
    require(BODY_COUNT - CUMULATIVE_CLOSED_COUNT == REMAINING_BODY_COUNT,
            (BODY_COUNT, CUMULATIVE_CLOSED_COUNT, REMAINING_BODY_COUNT))

    semantic_payload = (
        tuple(sorted(code_hist.items())),
        tuple(sorted(edge12_code_hist.items())),
        k3_types,
        triangle_rows,
        star_rows,
        star_types,
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
        "LRC14 reflected robust-edge-12 threshold-shape uniform closure exact referee",
        f"threshold_creation_code_hist={tuple(sorted(code_hist.items()))};"
        f"edge12_code_hist={tuple(sorted(edge12_code_hist.items()))}",
        f"normalized_low_K3_types={len(k3_types)};high_triangle_types={len(triangle_rows)};"
        f"high_star_types={len(star_types)};projective_star_types={star_types}",
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
        "threshold_law=robust iff a+b<theta(E);each induced complement has isolated minimum or dominating maximum",
        "scale_law=multiplying all levels preserves the reduced channel and cell, decreases transport loss and singleton debt",
        "conclusion=all 37 robust-edge-12 bodies close for every assignment of positive reflected levels",
        f"corollary=arbitrary-level body closure rises from {PRIOR_CLOSED_COUNT} to {CUMULATIVE_CLOSED_COUNT};remaining_bodies={REMAINING_BODY_COUNT}",
        "scope=reflected THM-2941 residual family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"low_phase_sha256={sha256(LOW_PHASE)}",
        f"exceptional_sha256={sha256(EXCEPTIONAL)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"k6_minus_edge_sha256={sha256(K6_MINUS_EDGE)}",
        f"edge13_sha256={sha256(EDGE13)}",
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
