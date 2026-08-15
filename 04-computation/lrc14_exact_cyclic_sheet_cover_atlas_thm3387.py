#!/usr/bin/env python3
"""Exact Boolean sheet-cover atlas for THM-3387.

THM-3385 gives a sufficient sum-capacity criterion.  Here the true finite
object is used: at base y, each transverse speed blocks a subset of Z/q, and
the core quotient is exact iff every full transverse sheet cover occurs
inside the descended core danger set.  The literal six-body universe is
classified by two independent routes:

* exact aligned-cell interval projection from the pinned THM-3385 companion;
* an integer event sweep of the Boolean sheet hypergraph.

The q=2 hypergraph collapses further to an undirected gcd graph.  Runtime
gates use RuntimeError, so python -O retains every decision.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/lrc14_q_fibre_complement_clock_thm3385.py"
PINNED = (
    (
        "THM-3385",
        ROOT / "01-canon/theorems/THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks.md",
        "f22f36c898ffd508e65210ccae720ce185c9bd5c6bff606d7caa018c1caa8db7",
    ),
    (
        "THM-3385-script",
        BASE_PATH,
        "8d3add2104ca91e89017ae643f573c8c5f6b45a4fdfecd9d65c16f48fd51b356",
    ),
    (
        "THM-3385-output",
        ROOT / "05-knowledge/results/lrc14_q_fibre_complement_clock_thm3385.out",
        "5b7349e13e36e8424346aa2c50f26f760f766a161f046ca09ad7bcd6e4d665b4",
    ),
    (
        "THM-3366-allk-audit",
        ROOT / "05-knowledge/results/lrc14_allk_complement_clock_independent_audit_20260814.out",
        "31c096e528a148cd717f4d9c7485c93cdec6175b33b83bae9d2be627249d743c",
    ),
    (
        "THM-3366-S172-output",
        ROOT / "05-knowledge/results/lrc14_k1_body_complement_clock_scan_kps_s172.out",
        "94b0f25f0c61bbfc3d8acf6abfbda5856475c69f90aed24ba2de8738cc7a59bf",
    ),
)

EXPECTED_CANDIDATES = 23_569
EXPECTED_EXACT_ROWS = 15_393
EXPECTED_GLOBAL_TRANSVERSE_ROWS = 15_381
EXPECTED_CAPACITY_ROWS = 6_420
EXPECTED_Q_HISTOGRAM = (
    (2, 252),
    (3, 588),
    (4, 619),
    (5, 1619),
    (6, 1478),
    (7, 2079),
    (8, 1152),
    (9, 1205),
    (10, 1269),
    (11, 1287),
    (12, 1271),
    (13, 1287),
    (14, 1287),
)
EXPECTED_STRUCTURAL_ROWS = (15_246, 15_393, 15_393, 15_393, 2735, 1104, 291)
EXPECTED_STRUCTURAL_OCCURRENCES = (
    93_732_978_513_930,
    3_659_255_462_265,
    133_947_094_813,
    4_445_930_697,
    15_952_450,
    154_344,
    291,
)
EXPECTED_S172_SIZE_HISTOGRAM = ((1, 12658), (2, 2029), (3, 409), (4, 150), (5, 147))
EXPECTED_S172_OCCURRENCE_HISTOGRAM = (
    (1, 3_427_329_787_389),
    (2, 182_339_009_252),
    (3, 35_772_924_216),
    (4, 1_643_011_511),
    (5, 12_170_729_897),
)
EXPECTED_Q2_INDEPENDENT_SETS = (
    (1,),
    (3,),
    (5,),
    (7,),
    (9,),
    (11,),
    (13,),
    (1, 3),
    (1, 5),
    (3, 9),
)
EXPECTED_CORE_RESCUES = (
    ((1, 3, 6, 8, 11, 13), 48048, 8008, 6, (1,), (1, 3, 8, 11, 13)),
    ((3, 5, 6, 8, 11, 13), 240240, 40040, 6, (1,), (3, 5, 8, 11, 13)),
    ((3, 6, 7, 8, 11, 13), 336336, 56056, 6, (1,), (3, 7, 8, 11, 13)),
    ((3, 6, 8, 9, 11, 13), 144144, 48048, 3, (1, 2, 3), (8, 11, 13)),
    ((3, 6, 8, 9, 11, 13), 144144, 24024, 6, (1,), (3, 8, 9, 11, 13)),
    ((3, 6, 8, 10, 11, 13), 240240, 40040, 6, (1,), (3, 8, 10, 11, 13)),
    ((3, 6, 8, 11, 12, 13), 48048, 16016, 3, (1, 2, 4), (8, 11, 13)),
    ((3, 6, 8, 11, 12, 13), 48048, 8008, 6, (1, 2), (3, 8, 11, 13)),
    ((3, 6, 8, 11, 13, 14), 336336, 56056, 6, (1,), (3, 8, 11, 13, 14)),
    ((5, 6, 11, 12, 13, 14), 840840, 168168, 5, (1,), (6, 11, 12, 13, 14)),
    ((6, 8, 9, 11, 12, 13), 144144, 48048, 3, (2, 3, 4), (8, 11, 13)),
    ((7, 9, 10, 11, 12, 13), 2522520, 504504, 5, (2,), (7, 9, 11, 12, 13)),
)
EXPECTED_SEMANTIC_DIGEST = "0a38ff30719647e20b6eecd9e2ea13e700426989ee209791ef62bdeeedacd496"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_base():
    spec = spec_from_file_location("thm3385_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "base import spec")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(b"\n")

    def hexdigest(self):
        return self._hash.hexdigest()


def strict_danger_numerator(speed, point_numerator, point_denominator):
    residue = (speed * point_numerator) % point_denominator
    distance = min(residue, point_denominator - residue)
    return 14 * distance < point_denominator


def source_sheet_blocked(speed, degree, sheet, sample, scale):
    denominator = 2 * scale * degree
    numerator = sample + 2 * scale * sheet
    return strict_danger_numerator(speed, numerator, denominator)


def transverse_event_positions(degree, transverse, scale):
    multiplier = scale // 14
    events = {0, scale}
    for speed in transverse:
        require(multiplier % speed == 0, ("event scale", degree, transverse, speed))
        for sheet in range(degree):
            for tooth in range(speed + 1):
                for sign in (-1, 1):
                    point = multiplier * degree * (14 * tooth + sign) // speed - sheet * scale
                    if 0 <= point <= scale:
                        events.add(point)
    return tuple(sorted(events))


def event_samples(events):
    return tuple(2 * point for point in events) + tuple(
        left + right for left, right in zip(events, events[1:])
    )


def transverse_full_cover(degree, transverse, sample, scale):
    return all(
        any(source_sheet_blocked(speed, degree, sheet, sample, scale) for speed in transverse)
        for sheet in range(degree)
    )


@lru_cache(maxsize=None)
def global_transverse_survives(degree, transverse):
    scale = 14 * lcm(*transverse)
    events = transverse_event_positions(degree, transverse, scale)
    for sample in event_samples(events):
        if transverse_full_cover(degree, transverse, sample, scale):
            return False
    return True


def full_cover_is_core_contained(degree, core, transverse):
    source_speeds = tuple(degree * clock for clock in core) + transverse
    scale = 14 * lcm(*source_speeds)
    multiplier = scale // 14
    events = set(transverse_event_positions(degree, transverse, scale))
    for clock in core:
        require(multiplier % clock == 0, ("core event scale", degree, core, clock))
        for tooth in range(clock + 1):
            for sign in (-1, 1):
                point = multiplier * (14 * tooth + sign) // clock
                if 0 <= point <= scale:
                    events.add(point)

    full_samples = 0
    for sample in event_samples(tuple(sorted(events))):
        if not transverse_full_cover(degree, transverse, sample, scale):
            continue
        full_samples += 1
        core_safe = all(
            not strict_danger_numerator(clock, sample, 2 * scale)
            for clock in core
        )
        if core_safe:
            return False, full_samples
    return True, full_samples


def q2_edge(left, right):
    require(left % 2 == 1 and right % 2 == 1 and left != right, (left, right))
    return left + right > 7 * gcd(left, right)


def q2_independent(transverse):
    return all(not q2_edge(left, right) for left, right in combinations(transverse, 2))


def all_divisors(number):
    low = []
    high = []
    divisor = 1
    while divisor * divisor <= number:
        if number % divisor == 0:
            low.append(divisor)
            if divisor * divisor != number:
                high.append(number // divisor)
        divisor += 1
    return tuple(low + high[::-1])


@lru_cache(maxsize=None)
def denominator_shapes_downward(modulus, arity):
    total = comb(len(all_divisors(modulus)) + arity - 2, arity)
    return total - sum(
        denominator_shapes_downward(divisor, arity)
        for divisor in all_divisors(modulus)
        if divisor < modulus
    )


def main():
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))
    base = load_base()

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(isinstance(node, ast.Constant) and isinstance(node.value, float) for node in ast.walk(tree)),
        "floating literal",
    )

    # The q=2 pair graph is checked independently far beyond the literal pool.
    pair_checks = 0
    for left, right in combinations(range(1, 200, 2), 2):
        predicted_survival = not q2_edge(left, right)
        exact_survival = global_transverse_survives(2, (left, right))
        require(predicted_survival == exact_survival, ("q2 edge", left, right))
        pair_checks += 1

    literal_odds = tuple(range(1, 15, 2))
    independent_sets = tuple(
        subset
        for size in range(1, 7)
        for subset in combinations(literal_odds, size)
        if q2_independent(subset)
    )
    require(independent_sets == EXPECTED_Q2_INDEPENDENT_SETS, independent_sets)

    semantic = ExactDigest()
    candidates = 0
    exact_records = []
    global_records = 0
    capacity_records = 0
    q_histogram = Counter()
    rescued = []

    for body in combinations(range(1, 15), 6):
        modulus, safe_ranges = base.safe_cell_ranges(body)
        for degree in range(2, 15):
            core = tuple(speed // degree for speed in body if speed % degree == 0)
            transverse = tuple(speed for speed in body if speed % degree != 0)
            if not core or not transverse:
                continue
            candidates += 1
            quotient = modulus // degree

            projected = base.projected_support_ranges(quotient, safe_ranges)
            unsupported = base.complementary_ranges(quotient, projected)
            descended = base.danger_cell_ranges(quotient, core)
            cell_exact = unsupported == descended

            globally_safe = global_transverse_survives(degree, transverse)
            if globally_safe:
                require(cell_exact, ("global survival must descend", body, degree))
            if not cell_exact:
                continue

            support_count = sum(right - left for left, right in projected)
            capacity = sum(base.fibre_capacity(speed, degree) for speed in transverse)
            record = (
                body,
                modulus,
                quotient,
                degree,
                core,
                transverse,
                capacity,
                support_count,
            )
            exact_records.append(record)
            q_histogram[degree] += 1
            semantic.add(record)

            if capacity < degree:
                capacity_records += 1
            if globally_safe:
                global_records += 1
            else:
                contained, full_samples = full_cover_is_core_contained(
                    degree, core, transverse
                )
                require(contained and full_samples > 0, ("core rescue", record, full_samples))
                rescued.append((body, modulus, quotient, degree, core, transverse))

    require(candidates == EXPECTED_CANDIDATES, candidates)
    require(len(exact_records) == EXPECTED_EXACT_ROWS, len(exact_records))
    require(global_records == EXPECTED_GLOBAL_TRANSVERSE_ROWS, global_records)
    require(capacity_records == EXPECTED_CAPACITY_ROWS, capacity_records)
    require(tuple(sorted(q_histogram.items())) == EXPECTED_Q_HISTOGRAM, q_histogram)
    require(tuple(rescued) == EXPECTED_CORE_RESCUES, rescued)

    # Core danger unions are irredundant on every exact row.  Since S172 uses
    # exactly the descended-body candidate tuple, the injection from this
    # atlas into S172's closed rows preserves its least-completion size.
    s172_size_histogram = Counter()
    s172_occurrence_histogram = Counter()
    for record in exact_records:
        core = record[4]
        full_union = base.danger_cell_ranges(record[2], core)
        for removed in core:
            reduced = tuple(clock for clock in core if clock != removed)
            require(
                base.danger_cell_ranges(record[2], reduced) != full_union,
                ("redundant descended core", record, removed),
            )
        size = len(core)
        weight = denominator_shapes_downward(record[2], 6)
        s172_size_histogram[size] += 1
        s172_occurrence_histogram[size] += weight
    require(
        tuple(sorted(s172_size_histogram.items())) == EXPECTED_S172_SIZE_HISTOGRAM,
        s172_size_histogram,
    )
    require(
        tuple(sorted(s172_occurrence_histogram.items()))
        == EXPECTED_S172_OCCURRENCE_HISTOGRAM,
        s172_occurrence_histogram,
    )

    q2_records = tuple(record for record in exact_records if record[3] == 2)
    require(len(q2_records) == 252, len(q2_records))
    require(all(q2_independent(record[5]) for record in q2_records), "q2 graph mismatch")
    require(Counter(len(record[5]) for record in q2_records) == {1: 147, 2: 105}, "q2 sizes")

    structural_rows = []
    structural_occurrences = []
    for k in range(7):
        selected = tuple(
            record
            for record in exact_records
            if Q(record[7], record[2]) <= base.CUTOFFS[k]
            and (k >= 1 or len(record[4]) <= 4)
        )
        structural_rows.append(len(selected))
        arity = 7 - k
        total = 0
        for record in selected:
            downward = denominator_shapes_downward(record[2], arity)
            require(
                downward == base.denominator_shapes(record[2], arity),
                ("occurrence route", record[2], arity, downward),
            )
            total += downward
        structural_occurrences.append(total)
        semantic.add((k, tuple((record[0], record[2], record[3]) for record in selected)))

    require(tuple(structural_rows) == EXPECTED_STRUCTURAL_ROWS, structural_rows)
    require(
        tuple(structural_occurrences) == EXPECTED_STRUCTURAL_OCCURRENCES,
        structural_occurrences,
    )

    # Strict/open and graph hostiles.
    strict_q7 = sum(
        base.danger(1, (Q(1, 2) + sheet) / 7)
        for sheet in range(7)
    )
    closed_q7 = sum(
        base.circular_distance((Q(1, 2) + sheet) / 7) <= Q(1, 14)
        for sheet in range(7)
    )
    require((strict_q7, closed_q7) == (0, 2), (strict_q7, closed_q7))
    require(global_transverse_survives(2, (1, 3)), "q2 nonedge positive")
    require(not global_transverse_survives(2, (1, 9)), "q2 edge hostile")
    require(not global_transverse_survives(3, (8, 11, 13)), "core rescue setup")

    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("THM-3387 EXACT CYCLIC SHEET-COVER ATLAS")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=FINITE-EXACT Boolean-fibre iff atlas plus PROVED q2 gcd graph")
    print("exact_iff=every_full_transverse_sheet_cover_lies_inside_descended_core_danger")
    print(f"candidate_body_degree_rows={candidates};exact_rows={len(exact_records)};failed_rows={candidates-len(exact_records)}")
    print(f"global_transverse_rows={global_records};core_rescued_rows={len(rescued)};capacity_subclass={capacity_records}")
    print(f"q_histogram={tuple(sorted(q_histogram.items()))}")
    print(f"structural_rows_k0_to_k6={tuple(structural_rows)}")
    print(f"structural_occurrences_k0_to_k6={tuple(structural_occurrences)}")
    print(f"s172_exact_identification=size_histogram:{tuple(sorted(s172_size_histogram.items()))};occurrence_histogram:{tuple(sorted(s172_occurrence_histogram.items()))}")
    print("q2_pair_edge=left+right>7*gcd(left,right);relation=symmetric_not_tournament")
    print(f"q2_pair_checks={pair_checks};literal_independent_sets={independent_sets}")
    print("q2_exact_rows=252=147_single_transverse+105_two_transverse;new_beyond_capacity=105")
    print(f"core_rescue_records={tuple(rescued)}")
    print("strictness_hostile=q7,u1,y1/2:strict_blocks_0,closed_blocks_2")
    print("capacity_hostile=q2:{1,3}_survives,{1,9}_covers;sum_capacity_equal")
    print("typing=Boolean_sheet_hypergraph;grid_owner_retained;counts_are_subclass_not_additive")
    print("scope=no_refined_ledger_decrement,no_physical_tail_realization,no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
