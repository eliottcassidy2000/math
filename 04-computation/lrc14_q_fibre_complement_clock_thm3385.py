#!/usr/bin/env python3
"""Exact companion for THM-3385's q-fibre complement-clock theorem.

For rho=1/14 and pi_q(x)=q*x, speeds divisible by q are constant on a
q-sheet fibre after division by q.  A transverse speed u blocks at most

    gcd(u,q) * ceil((q/gcd(u,q))/7)

sheets.  This companion checks the finite-fibre bound, the exact aligned
cell projection for every literal six-speed body in {1,...,14}, the induced
row/occurrence census, the seven refined k=3 instances, and hostile boundary
cases.  All gates are RuntimeError checks and remain active under python -O.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import comb, gcd, lcm
from pathlib import Path


RHO = Q(1, 14)
BODY_POOL = tuple(range(1, 15))
BODY_SIZE = 6

PINNED = (
    (
        "THM-3366",
        Path("01-canon/theorems/THM-3366-all-sector-complement-clock-completion.md"),
        "9d6a356cf7440f0eae271ec8163b53895df9ac70332be27a1b0fe62da5172e43",
    ),
    (
        "THM-3381",
        Path("01-canon/theorems/THM-3381-reflected-residue-affine-phase-transport-and-frozen-tree-stability.md"),
        "2aae4c0d7112c2d55e347016fb6fe194b921df63c33a7d7e429ec4d320676137",
    ),
    (
        "THM-2928-support",
        Path("04-computation/lrc14_two_drift_body_projection_support_thm2928.py"),
        "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a",
    ),
    (
        "THM-3366-k3-audit",
        Path("05-knowledge/results/lrc14_k3_refined_complement_clock_independent_audit_20260814.out"),
        "0bf0d3b3c6b530ed3b37b3cf3d6c9dd38edebf331e2f6c7fcbcca1741a5e15c0",
    ),
)

SAFE_FLOORS = {
    0: Q(1),
    1: Q(6, 7),
    2: Q(66, 91),
    3: Q(55, 91),
    4: Q(558, 1183),
    5: Q(478, 1365),
    6: Q(61, 273),
    7: Q(15, 154),
}
CUTOFFS = {
    k: (Q(1) - SAFE_FLOORS[7 - k]) / SAFE_FLOORS[k]
    for k in range(8)
}
EXPECTED_CUTOFFS = {
    0: Q(139, 154),
    1: Q(106, 117),
    2: Q(887, 990),
    3: Q(125, 143),
    4: Q(26, 31),
    5: Q(375, 478),
    6: Q(39, 61),
    7: Q(0),
}

EXPECTED_Q_HISTOGRAM = (
    (2, 147),
    (3, 45),
    (4, 35),
    (5, 495),
    (6, 57),
    (7, 2079),
    (11, 1287),
    (12, 196),
    (13, 1287),
    (14, 792),
)
EXPECTED_STRUCTURAL_ROWS = (6273, 6420, 6420, 6420, 1272, 227, 192)
EXPECTED_STRUCTURAL_OCCURRENCES = (
    83_942_791_283_821,
    3_115_391_844_730,
    106_719_668_256,
    3_253_046_409,
    8_230_092,
    35_818,
    192,
)
EXPECTED_SEMANTIC_DIGEST = "76345873bbfbe37fde097e9a23e92744fb376ad0eb94d7f0eb1ff33c5f7223de"

EXPECTED_SEVEN = (
    ((1, 2, 6, 8, 10, 14), 11760, 5880, 388),
    ((2, 3, 6, 8, 10, 14), 11760, 5880, 388),
    ((2, 5, 6, 8, 10, 14), 11760, 5880, 388),
    ((2, 6, 7, 8, 10, 14), 11760, 5880, 388),
    ((2, 6, 8, 9, 10, 14), 35280, 17640, 1008),
    ((2, 6, 8, 10, 11, 14), 129360, 64680, 2544),
    ((2, 6, 8, 10, 13, 14), 152880, 76440, 2544),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(b"\n")

    def hexdigest(self):
        return self._hash.hexdigest()


def floor_q(value):
    return value.numerator // value.denominator


def mod_one(value):
    return value - floor_q(value)


def circular_distance(value):
    residue = mod_one(value)
    return min(residue, 1 - residue)


def danger(speed, point, rho=RHO):
    return circular_distance(speed * point) < rho


def ceil_q(value):
    return -((-value.numerator) // value.denominator)


def fibre_capacity(speed, degree, rho=RHO):
    g = gcd(speed, degree)
    orbit = degree // g
    return g * ceil_q(2 * rho * orbit)


def merge_ranges(ranges):
    merged = []
    for left, right in sorted(ranges):
        if left >= right:
            continue
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def safe_cell_ranges(body):
    """Exact THM-2928 body-safe L-cell ranges, independently rebuilt."""
    modulus = lcm(*(14 * speed for speed in body))
    blocked = []
    for speed in body:
        half = modulus // (14 * speed)
        period = modulus // speed
        for tooth in range(speed + 1):
            center = tooth * period
            blocked.append((max(0, center - half), min(modulus, center + half)))
    blocked = merge_ranges(blocked)
    safe = []
    cursor = 0
    for left, right in blocked:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < modulus:
        safe.append((cursor, modulus))
    return modulus, tuple(safe)


def projected_support_ranges(modulus, safe_ranges):
    pieces = []
    for left, right in safe_ranges:
        length = right - left
        if length >= modulus:
            return ((0, modulus),)
        start = left % modulus
        end = start + length
        if end <= modulus:
            pieces.append((start, end))
        else:
            pieces.append((start, modulus))
            pieces.append((0, end - modulus))
    return merge_ranges(pieces)


def complementary_ranges(modulus, ranges):
    gaps = []
    cursor = 0
    for left, right in ranges:
        if cursor < left:
            gaps.append((cursor, left))
        cursor = right
    if cursor < modulus:
        gaps.append((cursor, modulus))
    return tuple(gaps)


def danger_cell_ranges(modulus, clocks):
    """Open-cell labels contained in the union of D_c, on a D-grid."""
    blocked = []
    for clock in clocks:
        require(modulus % (14 * clock) == 0, ("unaligned clock", modulus, clock))
        half = modulus // (14 * clock)
        period = modulus // clock
        for tooth in range(clock + 1):
            center = tooth * period
            blocked.append((max(0, center - half), min(modulus, center + half)))
    return merge_ranges(blocked)


def event_points(degree, core, transverse, rho=RHO):
    points = {Q(0), Q(1)}

    # Boundaries of the descended target combs D_c(y).
    for clock in core:
        for tooth in range(clock + 1):
            for sign in (-1, 1):
                point = (Q(tooth) + sign * rho) / clock
                if Q(0) <= point <= Q(1):
                    points.add(point)

    # Boundaries of every source clause D_s((y+k)/q).
    source_speeds = tuple(degree * clock for clock in core) + tuple(transverse)
    for speed in source_speeds:
        for sheet in range(degree):
            for tooth in range(speed + 1):
                for sign in (-1, 1):
                    point = degree * (Q(tooth) + sign * rho) / speed - sheet
                    if Q(0) <= point <= Q(1):
                        points.add(point)
    return tuple(sorted(points))


def source_has_safe_sheet(degree, core, transverse, base, rho=RHO):
    for sheet in range(degree):
        point = (base + sheet) / degree
        speeds = tuple(degree * clock for clock in core) + tuple(transverse)
        if all(not danger(speed, point, rho) for speed in speeds):
            return True
    return False


def exact_event_sweep(degree, core, transverse, rho=RHO):
    points = event_points(degree, core, transverse, rho)
    samples = points + tuple((left + right) / 2 for left, right in zip(points, points[1:]))
    for base in samples:
        target_safe = all(not danger(clock, base, rho) for clock in core)
        require(
            source_has_safe_sheet(degree, core, transverse, base, rho) == target_safe,
            ("event identity", degree, core, transverse, rho, base),
        )
    return len(samples)


def max_grid_hits(orbit, rho=RHO):
    events = {Q(0), Q(1)}
    for residue in range(orbit):
        for sign in (-1, 1):
            events.add(mod_one(sign * rho - Q(residue, orbit)))
    points = tuple(sorted(events))
    samples = points + tuple((left + right) / 2 for left, right in zip(points, points[1:]))
    return max(
        sum(danger(1, phase + Q(residue, orbit), rho) for residue in range(orbit))
        for phase in samples
    )


def arrangement_points(clocks):
    points = {Q(0), Q(1)}
    for clock in clocks:
        for tooth in range(clock + 1):
            for sign in (-1, 1):
                point = (Q(tooth) + sign * RHO) / clock
                if Q(0) <= point <= Q(1):
                    points.add(point)
    return tuple(sorted(points))


@lru_cache(maxsize=None)
def candidate_safe_structure(clocks):
    points = arrangement_points(clocks)
    safe_points = tuple(
        point for point in points if all(not danger(clock, point) for clock in clocks)
    )
    safe_atoms = tuple(
        (left, right)
        for left, right in zip(points, points[1:])
        if all(not danger(clock, (left + right) / 2) for clock in clocks)
    )
    return safe_points, safe_atoms


def covers_unsupported_cells(clocks, modulus, integer_gaps):
    gaps = tuple((Q(left, modulus), Q(right, modulus)) for left, right in integer_gaps)
    safe_points, safe_atoms = candidate_safe_structure(clocks)
    for point in safe_points:
        if (modulus * point).denominator == 1:
            continue
        if any(left < point < right for left, right in gaps):
            return False
    for left, right in safe_atoms:
        if any(max(left, gap_left) < min(right, gap_right) for gap_left, gap_right in gaps):
            return False
    return True


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


def mobius(number):
    result = 1
    remaining = number
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime:
            prime += 1
            continue
        remaining //= prime
        if remaining % prime == 0:
            return 0
        result = -result
        while remaining % prime == 0:
            remaining //= prime
        prime += 1
    return -result if remaining > 1 else result


@lru_cache(maxsize=None)
def denominator_shapes(modulus, arity):
    return sum(
        mobius(modulus // divisor)
        * comb(len(all_divisors(divisor)) + arity - 2, arity)
        for divisor in all_divisors(modulus)
    )


def capacity_positive_control(transverse, degree):
    return sum(fibre_capacity(speed, degree) for speed in transverse) < degree


def main():
    require(CUTOFFS == EXPECTED_CUTOFFS, ("cutoffs", CUTOFFS))
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(isinstance(node, ast.Constant) and isinstance(node.value, float) for node in ast.walk(tree)),
        "floating literal",
    )

    # The open-arc fibre capacity is exact, not merely a union bound per speed.
    capacity_checks = 0
    for orbit in range(1, 101):
        exact = max_grid_hits(orbit)
        predicted = (orbit + 6) // 7
        require(exact == predicted, ("grid capacity", orbit, exact, predicted))
        capacity_checks += 1

    records = []
    q_histogram = Counter()
    row_keys = set()
    semantic = ExactDigest()
    first_by_q = {}

    for body in combinations(BODY_POOL, BODY_SIZE):
        modulus, safe_ranges = safe_cell_ranges(body)
        for degree in range(2, 15):
            core = tuple(speed // degree for speed in body if speed % degree == 0)
            transverse = tuple(speed for speed in body if speed % degree != 0)
            if not core or not transverse or modulus % degree:
                continue
            capacity = sum(fibre_capacity(speed, degree) for speed in transverse)
            if capacity >= degree:
                continue

            quotient = modulus // degree
            projected = projected_support_ranges(quotient, safe_ranges)
            unsupported = complementary_ranges(quotient, projected)
            expected_unsupported = danger_cell_ranges(quotient, core)
            require(
                unsupported == expected_unsupported,
                ("discrete q-fold identity", body, modulus, degree, unsupported, expected_unsupported),
            )
            support_count = sum(right - left for left, right in projected)
            require(
                support_count + sum(right - left for left, right in unsupported) == quotient,
                ("partition", body, quotient),
            )

            key = (body, quotient)
            require(key not in row_keys, ("duplicate body/divisor certificate", key, degree))
            row_keys.add(key)
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
            records.append(record)
            q_histogram[degree] += 1
            first_by_q.setdefault(degree, record)
            semantic.add(record)

    require(len(records) == 6420, ("certificate rows", len(records)))
    require(tuple(sorted(q_histogram.items())) == EXPECTED_Q_HISTOGRAM, q_histogram)

    # Independent continuous event sweeps for every occurring degree and all
    # seven refined q=2 rows retain strict boundary samples.
    event_sweep_checks = 0
    for record in first_by_q.values():
        event_sweep_checks += exact_event_sweep(record[3], record[4], record[5])

    structural_rows = []
    structural_occurrences = []
    for k in range(7):
        arity = 7 - k
        selected = [
            record
            for record in records
            if Q(record[7], record[2]) <= CUTOFFS[k]
            and (k >= 1 or len(record[4]) <= 4)
        ]
        structural_rows.append(len(selected))
        structural_occurrences.append(
            sum(denominator_shapes(record[2], arity) for record in selected)
        )
        semantic.add((k, tuple((record[0], record[2], record[3]) for record in selected)))

    require(tuple(structural_rows) == EXPECTED_STRUCTURAL_ROWS, structural_rows)
    require(
        tuple(structural_occurrences) == EXPECTED_STRUCTURAL_OCCURRENCES,
        structural_occurrences,
    )

    # The seven current refined k=3 hits are exactly the fixed half-even core.
    seven = []
    certificate_by_key = {(record[0], record[2]): record for record in records}
    for body, modulus, quotient, occurrences in EXPECTED_SEVEN:
        record = certificate_by_key[(body, quotient)]
        require(record[1] == modulus and record[3] == 2, ("seven modulus", record))
        require(record[4] == (1, 3, 4, 5, 7), ("seven core", record))
        require(len(record[5]) == 1 and record[5][0] % 2 == 1, ("seven odd", record))
        require(record[6] == 1, ("seven capacity", record))
        event_sweep_checks += exact_event_sweep(2, record[4], record[5])
        seven.append((body, modulus, quotient, record[4], occurrences))
    require(sum(item[-1] for item in seven) == 7648, seven)
    semantic.add(tuple(seven))

    # q=4 and q=6 are finite Boolean sheet-cover problems, with hyperedge
    # capacities determined by gcd type.  No orientation is introduced.
    q4_capacity = tuple((residue, fibre_capacity(residue, 4)) for residue in (1, 2, 3))
    q6_capacity = tuple((residue, fibre_capacity(residue, 6)) for residue in (1, 2, 3, 4, 5))
    require(q4_capacity == ((1, 1), (2, 2), (3, 1)), q4_capacity)
    require(q6_capacity == ((1, 1), (2, 2), (3, 3), (4, 2), (5, 1)), q6_capacity)

    # Sharp and hostile controls.
    require(exact_event_sweep(2, (), (1,), Q(1, 4)) > 0, "rho=1/4 strict boundary")
    x = Q(1, 4)
    require(
        danger(1, x, Q(2, 7)) and danger(1, x + Q(1, 2), Q(2, 7)),
        "rho>1/4 hostile",
    )

    # Capacity equality is not an iff: (1,3) survives, while (1,9) blocks
    # both q=2 sheets at y=1/9.
    require(not capacity_positive_control((1, 3), 2), "equality positive setup")
    require(exact_event_sweep(2, (), (1, 3)) > 0, "capacity equality positive")
    hostile_y = Q(1, 9)
    hostile_sheets = tuple((hostile_y + sheet) / 2 for sheet in range(2))
    require(
        all(any(danger(speed, point) for speed in (1, 9)) for point in hostile_sheets),
        ("multiple transverse hostile", hostile_y, hostile_sheets),
    )

    # A q-divisible speed is core data, not a transverse speed.
    require(
        all(danger(2, Q(sheet, 2)) for sheet in range(2)),
        "divisible speed omitted from core",
    )

    # The quotient divisor and endpoint convention are load-bearing.
    hostile_body = (1, 2, 6, 8, 10, 14)
    hostile_modulus, hostile_safe = safe_cell_ranges(hostile_body)
    wrong_quotient = hostile_modulus // 4
    wrong_gaps = complementary_ranges(
        wrong_quotient,
        projected_support_ranges(wrong_quotient, hostile_safe),
    )
    wrong_first_gap = (Q(wrong_gaps[0][0], wrong_quotient), Q(wrong_gaps[0][1], wrong_quotient))
    expected_half_even_first_gap = (Q(0), Q(1, 14))
    require(wrong_first_gap == (Q(0), Q(1, 28)), wrong_first_gap)
    require(wrong_first_gap != expected_half_even_first_gap, "wrong divisor accidentally exact")
    require(danger(1, Q(0)), "grid point belongs to open danger comb")
    # U_D contains no D-grid point, hence the theorem removes the grid from
    # the exact discrete equality and delegates it to THM-3366's owner clock.

    # The q=2 five-core certificate is genuinely outside the k=0 displayed
    # budget: all 147 body/divisor rows need five clocks within pool 1..14.
    # Check the actual THM-3366 target convention, including strict candidate
    # endpoints and removal of D-grid points, against every subset through 4.
    pool = tuple(range(1, 15))
    q2_five_core_records = tuple(
        record for record in records if record[3] == 2 and len(record[4]) == 5
    )
    require(len(q2_five_core_records) == 147, len(q2_five_core_records))
    four_clock_tests = 0
    for record in q2_five_core_records:
        integer_gaps = danger_cell_ranges(record[2], record[4])
        for size in range(5):
            for candidate in combinations(pool, size):
                four_clock_tests += 1
                require(
                    not covers_unsupported_cells(candidate, record[2], integer_gaps),
                    ("unexpected four-clock cover", record[0], record[2], candidate),
                )

    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("THM-3385 Q-FIBRE COMPLEMENT-CLOCK COMPANION")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=PROVED analytic q-fibre identity plus FINITE-EXACT literal-body census")
    print("fibre_capacity=gcd(u,q)*ceil((q/gcd(u,q))/7);sufficient_if=sum_capacity<q")
    print(f"capacity_exact_orbit_checks={capacity_checks};orbit_range=1..100")
    print(f"literal_body_certificates={len(records)};q_histogram={tuple(sorted(q_histogram.items()))}")
    print(f"structural_rows_k0_to_k6={tuple(structural_rows)}")
    print(f"structural_occurrences_k0_to_k6={tuple(structural_occurrences)}")
    print(f"continuous_event_samples={event_sweep_checks};strict_boundaries_included=True")
    print(f"refined_k3_seven={tuple(seven)};occurrences=7648;completion=(1,3,4,5,7)")
    print(f"boolean_fibre_q4={q4_capacity};boolean_fibre_q6={q6_capacity};object=hypergraph_not_tournament")
    print("hostile_radius=q2,rho>1/4,x=1/4 blocks both sheets;rho=1/4 strict survives")
    print("hostile_capacity=sum=q is inconclusive:(1,3) survives but (1,9) blocks y=1/9")
    print("hostile_divisibility=q2,u=2 blocks both sheets at y=0 and must descend as core clock 1")
    print(f"hostile_wrong_divisor=L/4 first_gap={wrong_first_gap} expected_half_even_gap={expected_half_even_first_gap}")
    print(f"hostile_k0_pool14=q2_rows:{len(q2_five_core_records)};strict_grid_aware_subset_tests:{four_clock_tests};all_need_at_least_5")
    print("grid_convention=U_D equals descended danger union minus the D-grid;owner clock remains required")
    print("scope=no arbitrary reflected phase,no necessity of capacity bound,no physical drift realization,no LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
