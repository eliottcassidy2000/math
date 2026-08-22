#!/usr/bin/env python3
"""Exact first-open six-complement relation on the THM-3366 raw rows.

For k=2,3, enumerate every literal body F in C([14],6), every inherited
divisor row passing the exact support cutoff, and every six-subset C of the
pool [14] that covers the unsupported open-cell target U_D(F).  Rows with a
cover of size at most five are kept out: they are already THM-3366 terminals.
The reported minimum histogram is deliberately truncated at depth six:
``None`` means no cover by at most six pool clocks, not no pool-14 cover.

The resulting full relation F -> C is a finite-exact boundary atlas.  It does
not prove LRC(14): six complement clocks plus the seven inherited clocks give
thirteen clocks, exactly the first open count, and mutation forgets the next
row's aligned sector, divisor, tail labels, and distinctness sidecars.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/lrc14_k1_body_complement_clock_scan_kps_s172.py"
EXPECTED_BASE_SHA256 = "bdb2001cf22f7e92884e895b0095021e42e8f1febd9adbf779b250a2f6c53507"
EXPECTED = {
    2: {
        "input": 27_163,
        "hist": ((1, 12_662), (2, 2_764), (3, 998), (4, 1_106),
                 (5, 1_668), (6, 4_814), (None, 3_151)),
        "rows": 4_814,
        "incidences": 4_918,
        "completion_hist": ((1, 4_727), (2, 70), (3, 17)),
        "edges": 4_908,
        "outdegree": ((0, 1_539), (1, 1_099), (2, 299), (3, 57),
                      (4, 8), (5, 1)),
        "sinks": 1_539,
    },
    3: {
        "input": 26_970,
        "hist": ((1, 12_659), (2, 2_764), (3, 976), (4, 1_052),
                 (5, 1_602), (6, 4_778), (None, 3_139)),
        "rows": 4_778,
        "incidences": 4_881,
        "completion_hist": ((1, 4_692), (2, 69), (3, 17)),
        "edges": 4_872,
        "outdegree": ((0, 1_559), (1, 1_088), (2, 295), (3, 54),
                      (4, 6), (5, 1)),
        "sinks": 1_559,
    },
}
EXPECTED_SEMANTIC_SHA256 = "ce82dbe6054cea1a250d6632da5b0e2d9c695445564a026b4851b3c7e5bb1fc0"
EXPECTED_SENTINEL_HOSTILE = (
    (1, 2, 4, 6, 9, 10),
    2_520,
    1_260,
    646,
    (1, 2, 3, 5, 8, 9, 10),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_base():
    require(lf_sha256(BASE_PATH) == EXPECTED_BASE_SHA256, "base dependency changed")
    spec = spec_from_file_location("exact_six_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "base import")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def strongly_connected_components(vertices, adjacency):
    """Deterministic Kosaraju decomposition."""
    reverse = {vertex: set() for vertex in vertices}
    for vertex, targets in adjacency.items():
        for target in targets:
            reverse[target].add(vertex)

    sys.setrecursionlimit(10_000)
    seen = set()
    order = []

    def forward(vertex):
        seen.add(vertex)
        for target in sorted(adjacency[vertex]):
            if target not in seen:
                forward(target)
        order.append(vertex)

    for vertex in sorted(vertices):
        if vertex not in seen:
            forward(vertex)

    seen.clear()
    components = []

    def backward(vertex, component):
        seen.add(vertex)
        component.append(vertex)
        for source in sorted(reverse[vertex]):
            if source not in seen:
                backward(source, component)

    for vertex in reversed(order):
        if vertex not in seen:
            component = []
            backward(vertex, component)
            components.append(tuple(sorted(component)))
    return tuple(sorted(components, key=lambda component: (len(component), component)))


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert found")

    base = load_base()
    support = base.load_support()
    pool = tuple(range(1, 15))
    bodies = tuple(combinations(pool, 6))
    points = base.arrangement_points(pool)
    atoms = tuple(zip(points, points[1:]))
    samples = points + tuple((left + right) / 2 for left, right in atoms)
    masks = tuple(
        sum(
            (1 << index for index, x_value in enumerate(samples)
             if base.danger(clock, x_value)),
            0,
        )
        for clock in pool
    )
    union_mask = 0
    for mask in masks:
        union_mask |= mask
    covers_by_bit = tuple(
        tuple(index for index, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(len(samples))
    )

    six_unions = []
    for completion in bodies:
        mask = 0
        for clock in completion:
            mask |= masks[clock - 1]
        six_unions.append((completion, mask))
    six_unions = tuple(six_unions)

    @lru_cache(maxsize=None)
    def solve_exact(remaining: int, depth: int):
        if remaining == 0:
            return ()
        if depth == 0:
            return None
        bits = tuple(
            bit for bit in range(len(samples)) if remaining & (1 << bit)
        )
        pivot = min(bits, key=lambda bit: len(covers_by_bit[bit]))
        for index in covers_by_bit[pivot]:
            reduced = remaining & ~masks[index]
            if reduced == remaining:
                continue
            suffix = solve_exact(reduced, depth - 1)
            if suffix is not None:
                return (pool[index],) + suffix
        return None

    def target_mask(divisor, gaps):
        target = 0
        for index, x_value in enumerate(points):
            if ((divisor * x_value).denominator != 1
                    and any(left < x_value < right for left, right in gaps)):
                target |= 1 << index
        offset = len(points)
        for index, (left, right) in enumerate(atoms):
            if any(
                max(left, gap_left) < min(right, gap_right)
                for gap_left, gap_right in gaps
            ):
                target |= 1 << (offset + index)
        return target

    cutoffs = {2: Q(887, 990), 3: Q(125, 143)}
    reports = []
    semantic = sha256()
    sentinel_hostile = None
    for sector in (2, 3):
        input_rows = 0
        minimum_histogram = Counter()
        exact_rows = []
        relation = defaultdict(set)
        completion_count_histogram = Counter()
        completion_incidences = 0
        full_period_unique = 0
        lower_self_rows = []

        for body in bodies:
            ruler, ranges = support.safe_cell_ranges(body)
            for divisor in support.divisors(ruler):
                count = support.support_size_bitset(divisor, ranges)
                if Q(count, divisor) > cutoffs[sector]:
                    continue
                input_rows += 1
                arcs = base.residue_arcs(divisor, ranges)
                gaps = base.unsupported_gaps(divisor, arcs)
                target = target_mask(divisor, gaps)
                if (
                    sector == 2
                    and body == EXPECTED_SENTINEL_HOSTILE[0]
                    and divisor == EXPECTED_SENTINEL_HOSTILE[2]
                ):
                    witness = EXPECTED_SENTINEL_HOSTILE[-1]
                    witness_mask = 0
                    for clock in witness:
                        witness_mask |= masks[clock - 1]
                    require(
                        solve_exact(target, 6) is None,
                        ("sentinel unexpectedly has a six-cover", body, divisor),
                    )
                    require(
                        target & ~witness_mask == 0,
                        ("sentinel seven-cover failed", body, divisor, witness),
                    )
                    sentinel_hostile = (body, ruler, divisor, count, witness)
                if target & ~union_mask:
                    minimum_histogram[None] += 1
                    continue
                minimum = None
                for depth in range(7):
                    if solve_exact(target, depth) is not None:
                        minimum = depth
                        break
                minimum_histogram[minimum] += 1
                if minimum != 6:
                    continue

                completions = tuple(
                    completion
                    for completion, mask in six_unions
                    if target & ~mask == 0
                )
                require(completions, (sector, body, divisor))
                exact_rows.append((body, ruler, divisor, count, completions))
                completion_count_histogram[len(completions)] += 1
                completion_incidences += len(completions)
                relation[body].update(completions)
                if divisor == ruler:
                    require(completions == (body,),
                            ("full-period uniqueness", sector, body, completions))
                    full_period_unique += 1
                elif body in completions:
                    lower_self_rows.append((body, ruler, divisor, count))

        vertices = set(bodies)
        adjacency = {
            body: set(relation.get(body, ())) - {body}
            for body in vertices
        }
        components = strongly_connected_components(vertices, adjacency)
        cycles = tuple(component for component in components if len(component) > 1)
        outdegree = tuple(sorted(Counter(map(len, adjacency.values())).items()))
        sinks = sum(not targets for targets in adjacency.values())
        edge_count = sum(len(targets) for targets in relation.values())
        report = (
            sector,
            input_rows,
            tuple(sorted(minimum_histogram.items(), key=lambda item: (item[0] is None, item[0]))),
            len(exact_rows),
            completion_incidences,
            tuple(sorted(completion_count_histogram.items())),
            full_period_unique,
            tuple(lower_self_rows),
            edge_count,
            outdegree,
            sinks,
            cycles,
        )
        reports.append(report)
        for row in exact_rows:
            semantic.update((repr((sector, row)) + "\n").encode("ascii"))

        expected = EXPECTED[sector]
        require(input_rows == expected["input"], (sector, input_rows))
        require(report[2] == expected["hist"], (sector, report[2]))
        require(len(exact_rows) == expected["rows"], (sector, len(exact_rows)))
        require(completion_incidences == expected["incidences"],
                (sector, completion_incidences))
        require(report[5] == expected["completion_hist"], (sector, report[5]))
        require(full_period_unique == 3_003, (sector, full_period_unique))
        require(len(lower_self_rows) == 6, (sector, lower_self_rows))
        require(edge_count == expected["edges"], (sector, edge_count))
        require(outdegree == expected["outdegree"], (sector, outdegree))
        require(sinks == expected["sinks"], (sector, sinks))
        require(tuple(map(len, cycles)) == (3, 3), (sector, cycles))

    require(
        sentinel_hostile == EXPECTED_SENTINEL_HOSTILE,
        ("truncated-minimum sentinel", sentinel_hostile),
    )

    semantic_digest = semantic.hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256,
                (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("LRC14 exact-six complement mutation relation")
    print("status=FINITE-EXACT first-open boundary atlas;not a terminal or physical iteration theorem")
    print(f"dependency_sha256_lf={lf_sha256(BASE_PATH)}")
    for report in reports:
        (sector, input_rows, minimum_histogram, exact_rows, incidences,
         completion_histogram, full_period_unique, lower_self_rows, edges,
         outdegree, sinks, cycles) = report
        print(
            f"k{sector}_input_rows={input_rows};"
            f"minimum_through_6_histogram={minimum_histogram}"
        )
        print(f"k{sector}_exact6_rows={exact_rows};completion_incidences={incidences};completion_count_histogram={completion_histogram}")
        print(f"k{sector}_full_period_unique={full_period_unique};lower_self_rows={lower_self_rows}")
        print(f"k{sector}_body_relation_edges={edges};nonself_outdegree={outdegree};mutation_sinks={sinks}")
        print(f"k{sector}_nontrivial_sccs={cycles}")
    print(
        "truncated_minimum_sentinel="
        f"{sentinel_hostile};no_cover_by_at_most_6=True;seven_cover_verified=True"
    )
    print("typed_stopping_boundary=six_complements_plus_seven_inherited_clocks_equals_13_open_clocks;mutation_forgets_next_sector_divisor_tail_and_distinctness")
    print(f"semantic_sha256={semantic_digest}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=stdlib_exact_Fraction_and_bitsets;no_float;no_assert;normal_and_O_truth_gates")


if __name__ == "__main__":
    main()
