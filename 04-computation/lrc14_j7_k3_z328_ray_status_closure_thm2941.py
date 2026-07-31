#!/usr/bin/env python3
"""Exact all-label closure of the nine projected ``k=3, z1=328`` rows.

The scalar supplier reconstructs the nine bodies globally from all 3,003
six-body carriers.  The residue-ray/reversal law then computes the exact
infinite-label maximum for every denominator class.  Every remaining state is
rejected by either a crude fibre overload or the common 16-cell Hunter status
table.  Each status rejection is replayed below by a second exact rational
matrix checker independent of the status solver.  It is then normalized by a
basis-independent minimal split-Smith barcode circuit and a canonical exact
affine tariff.  Solver-selected Farkas bases are never frozen into partitions,
representatives, or semantic hashes.  No finite label horizon remains after
the scalar handoff.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_ENGINE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
)
SCALAR_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z328_scalar_slice_thm2941.py"
)
SCALAR_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z328_scalar_slice_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z328_ray_status_closure_thm2941.out"
)
EXPECTED_RAY_ENGINE_SHA256 = (
    "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2"
)
EXPECTED_SCALAR_SHA256 = (
    "75fd846d9070b267003c70e175a9af1225815c3567a5a2558891f5af7a61f8f3"
)
EXPECTED_SCALAR_OUTPUT_SHA256 = (
    "f5c11f364a626141af181d84f39d48030ef91a8ddf7d74b9602bb15cd7eb626e"
)
EXPECTED_SEMANTIC_SHA256 = (
    "70164f3126f95ecfd6c38f88fc1d3c65ba90264e9fe8e1ecd61d394f5d91c5ab"
)
EXPECTED_INSTANCE_TARIFF_SHA256 = (
    "3959d0bf62626f84758da687687d83156f7f3cee47e63ac2aff722a3a27bdd8b"
)
FIRST = 328
EXPECTED_COUNTS = {
    (1, 2, 6, 8, 12, 14): (1, 0, 1),
    (1, 2, 8, 10, 12, 14): (3, 0, 3),
    (1, 4, 6, 8, 12, 14): (1, 0, 1),
    (1, 4, 8, 10, 12, 14): (35, 16, 19),
    (1, 6, 8, 10, 12, 14): (2, 0, 2),
    (2, 4, 6, 8, 12, 14): (1, 0, 1),
    (2, 4, 8, 10, 12, 14): (8, 2, 6),
    (2, 6, 8, 10, 12, 14): (8, 0, 8),
    (2, 8, 10, 11, 12, 14): (26, 18, 8),
}
EXPECTED_TOTALS = (85, 36, 49, 0)
EXPECTED_M_HISTOGRAM = ((4, 2), (5, 22), (6, 4), (7, 21))
EXPECTED_CIRCUIT_CENSUS = (
    ((2,), 7),
    ((2, 4), 3),
    ((2, 5, 6), 1),
    ((3,), 29),
    ((4,), 9),
)
EXPECTED_DISTINCT_TARIFFS = 25


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(RAY_ENGINE_PATH) == EXPECTED_RAY_ENGINE_SHA256, "ray engine changed")
require(file_sha256(SCALAR_PATH) == EXPECTED_SCALAR_SHA256, "scalar supplier changed")
if EXPECTED_SCALAR_OUTPUT_SHA256 is not None:
    require(
        file_sha256(SCALAR_OUTPUT_PATH) == EXPECTED_SCALAR_OUTPUT_SHA256,
        "scalar transcript changed",
    )
ray = load_module("z328_ray_engine", RAY_ENGINE_PATH)
scalar = load_module("z328_scalar_supplier", SCALAR_PATH)
ray.FIRST = FIRST
BODIES = scalar.EXPECTED_BODIES
require(tuple(EXPECTED_COUNTS) == BODIES, "body/count order mismatch")


def independent_farkas_check(q, marginals, capacities, histogram, certificate):
    """Rebuild and verify the 16-cell Farkas inequalities exactly."""
    thresholds, alpha, z = certificate
    tail_rows = []
    tail_rhs = []
    rebuilt_thresholds = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(count for load, count in histogram if load >= threshold)
        good = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(good):
            continue
        rebuilt_thresholds.append(threshold)
        tail_rows.append(good)
        tail_rhs.append(demand)
    equality_rows = [
        (1,) * 16,
        *[
            tuple((pattern >> index) & 1 for pattern in range(16))
            for index in range(4)
        ],
    ]
    equality_rhs = (q, *marginals)
    require(tuple(rebuilt_thresholds) == tuple(thresholds), "threshold mismatch")
    require(len(alpha) == len(tail_rows) and len(z) == 5, "certificate shape")
    require(all(value >= 0 for value in alpha), "negative Farkas alpha")
    slacks = []
    for pattern in range(16):
        value = sum(z[row] * equality_rows[row][pattern] for row in range(5))
        value -= sum(
            alpha[row] * tail_rows[row][pattern] for row in range(len(alpha))
        )
        slacks.append(value)
    contradiction = sum(z[row] * equality_rhs[row] for row in range(5))
    contradiction -= sum(
        alpha[row] * tail_rhs[row] for row in range(len(alpha))
    )
    require(all(value >= 0 for value in slacks), "negative Farkas column")
    require(contradiction < 0, "nonnegative Farkas contradiction")
    return tuple(slacks), contradiction < 0


def solve_square(rows):
    """Solve one square augmented rational system, or return ``None``."""
    matrix = [[Q(value) for value in row] for row in rows]
    n = len(matrix)
    require(n and all(len(row) == n + 1 for row in matrix), "nonsquare system")
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if matrix[row][column]),
            None,
        )
        if pivot is None:
            return None
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = matrix[column][column]
        matrix[column] = [value / scale for value in matrix[column]]
        for row in range(n):
            if row == column or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                left - scale * right
                for left, right in zip(matrix[row], matrix[column])
            ]
    return tuple(matrix[row][-1] for row in range(n))


STATUS_COLUMNS = tuple(
    (Q(1), *(Q((pattern >> bit) & 1) for bit in range(4)))
    for pattern in range(16)
)


@lru_cache(maxsize=None)
def tariff_vertices(good):
    """Enumerate all exact affine-majorant vertices for a 16-cell good mask."""
    good = tuple(map(Q, good))
    vertices = set()
    for active in combinations(range(16), 5):
        answer = solve_square(
            ((*STATUS_COLUMNS[pattern], good[pattern]) for pattern in active)
        )
        if answer is None:
            continue
        if all(
            sum(
                coefficient * value
                for coefficient, value in zip(answer, column)
            )
            >= good[pattern]
            for pattern, column in enumerate(STATUS_COLUMNS)
        ):
            vertices.add(answer)
    require(vertices, ("no tariff vertex", good))
    return tuple(sorted(vertices))


def canonical_tariff(q, marginals, good):
    """Minimum objective, then lexicographically first exact affine tariff."""
    objective = (Q(q), *(Q(value) for value in marginals))
    return min(
        (
            sum(
                coefficient * value
                for coefficient, value in zip(alpha, objective)
            ),
            alpha,
        )
        for alpha in tariff_vertices(tuple(good))
    )


def canonical_primal(q, marginals, good):
    """Exact maximum and lexicographically first basic common status table."""
    rhs = (Q(q), *(Q(value) for value in marginals))
    candidates = []
    for basis in combinations(range(16), 5):
        rows = [
            tuple(STATUS_COLUMNS[pattern][coordinate] for pattern in basis)
            + (rhs[coordinate],)
            for coordinate in range(5)
        ]
        values = solve_square(rows)
        if values is None or any(value < 0 for value in values):
            continue
        table = tuple(
            next(
                (
                    values[index]
                    for index, pattern in enumerate(basis)
                    if pattern == cell
                ),
                Q(0),
            )
            for cell in range(16)
        )
        require(
            tuple(
                sum(
                    table[pattern] * STATUS_COLUMNS[pattern][coordinate]
                    for pattern in range(16)
                )
                for coordinate in range(5)
            )
            == rhs,
            ("bad canonical primal", basis),
        )
        value = sum(table[pattern] * good[pattern] for pattern in range(16))
        candidates.append((value, table))
    require(candidates, ("empty marginal polytope", q, marginals))
    optimum = max(value for value, _table in candidates)
    table = min(table for value, table in candidates if value == optimum)
    return optimum, table


def demand_at(histogram, grade):
    return sum(count for load, count in histogram if load >= grade)


def histogram_for_selected(q, selected):
    """Realize selected nested tail demands as a literal load histogram."""
    selected = tuple(sorted(selected))
    rows = []
    if selected and q > selected[0][1]:
        rows.append((0, q - selected[0][1]))
    for index, (grade, demand) in enumerate(selected):
        next_demand = (
            selected[index + 1][1] if index + 1 < len(selected) else 0
        )
        count = demand - next_demand
        require(count >= 0, ("nonnested tails", selected))
        if count:
            rows.append((grade, count))
    require(sum(count for _load, count in rows) == q, (q, selected, rows))
    return tuple(rows)


def minimal_circuit(q, marginals, capacities, histogram):
    """Cardinality-first, then lexicographic split-Smith tail circuit."""
    grades = tuple(
        load
        for load, _count in histogram
        if load > 0 and not all(capacity >= load for capacity in capacities)
    )
    demands = {grade: demand_at(histogram, grade) for grade in grades}
    for size in range(1, len(grades) + 1):
        for support in combinations(grades, size):
            if size == 1:
                good = tuple(Q(capacity >= support[0]) for capacity in capacities)
                upper, _alpha = canonical_tariff(q, marginals, good)
                feasible = upper >= demands[support[0]]
            else:
                selected = tuple((grade, demands[grade]) for grade in support)
                synthetic = histogram_for_selected(q, selected)
                feasible, exact = ray.local.common_status_feasible(
                    q, marginals, capacities, synthetic
                )
                require(exact is not None, ("missing exact circuit check", support))
            if not feasible:
                return support, demands
    raise RuntimeError(("status kill has no barcode circuit", grades))


def canonical_barcode_tariff(q, marginals, capacities, histogram):
    """Return the basis-independent circuit and canonical exact tariff."""
    support, demand_map = minimal_circuit(
        q, marginals, capacities, histogram
    )
    goods = tuple(
        tuple(Q(capacity >= grade) for capacity in capacities)
        for grade in support
    )
    demands = tuple(Q(demand_map[grade]) for grade in support)
    if len(support) == 1:
        lambdas = (Q(1),)
        upper, alpha = canonical_tariff(q, marginals, goods[0])
        primal, table = canonical_primal(q, marginals, goods[0])
        require(upper == primal, ("primal/dual mismatch", support, upper, primal))
    else:
        # The normalized uniform circuit is exact on the complete z=328
        # universe.  This is a deterministic certificate, not a solver basis.
        lambdas = tuple(Q(1, len(support)) for _grade in support)
        combined_good = tuple(
            sum(
                weight * good[pattern]
                for weight, good in zip(lambdas, goods)
            )
            for pattern in range(16)
        )
        upper, alpha = canonical_tariff(q, marginals, combined_good)
        table = ()
    weighted_demand = sum(
        weight * demand for weight, demand in zip(lambdas, demands)
    )
    deficit = weighted_demand - upper
    require(deficit > 0, ("nonpositive canonical deficit", support, deficit))
    return support, demands, upper, deficit, lambdas, alpha, table


def evaluate_body(body):
    stream = ray.Stream(body)
    trials, states, checks, signs = ray.ray_quotient_states(stream)
    crude, status, survivors = ray.common_status_screen(stream, states)
    require(
        (len(states), len(crude), len(status)) == EXPECTED_COUNTS[body],
        ("body counts changed", body, len(states), len(crude), len(status)),
    )
    require(set(states) == set(crude) | set(status) | set(survivors), "bad partition")
    require(not (set(crude) & set(status)), "kill stages overlap")

    for ds, witness in crude.items():
        gap, q, M, target, capacity = witness
        D = lcm(*ds)
        require(M == D // q, "crude modulus changed")
        require(target == dict(stream.target_data(D))[q], "crude target changed")
        require(
            capacity == sum(ray.local.fibre_cap(D, d, q) for d in ds),
            "crude capacity changed",
        )
        require(gap == target - capacity and gap > 0, "invalid crude witness")

    verified_farkas_checks = 0
    all_negative = True
    canonical_records = []
    for ds, witness in sorted(status.items()):
        q, M, marginals, cap_set, histogram, certificate = witness
        D = lcm(*ds)
        require(M == D // q, "status modulus changed")
        rebuilt_marginals, capacities = ray.local.hunter_status_data(D, ds, q)
        require(rebuilt_marginals == marginals, "status marginals changed")
        require(tuple(sorted(set(capacities))) == cap_set, "status caps changed")
        arcs = ray.fibre.projected_support_arcs(D, stream.ranges)
        require(
            ray.fibre.residue_load_histogram(arcs, q) == histogram,
            "status load histogram changed",
        )
        slacks, certificate_negative = independent_farkas_check(
            q, marginals, capacities, histogram, certificate
        )
        verified_farkas_checks += 1
        all_negative = all_negative and certificate_negative
        require(all(value >= 0 for value in slacks), "internal Farkas replay failed")
        (
            support,
            demands,
            upper,
            deficit,
            lambdas,
            alpha,
            table,
        ) = canonical_barcode_tariff(q, marginals, capacities, histogram)
        canonical_records.append(
            (
                body,
                ds,
                witness[:-1],
                support,
                demands,
                upper,
                deficit,
                lambdas,
                alpha,
                table,
            )
        )

    sign_totals = {
        sign: sum(
            count
            for (_d, candidate_sign), count in signs.items()
            if candidate_sign == sign
        )
        for sign in (-1, 0, 1)
    }
    require(sign_totals[-1] == sign_totals[1], "ray sign imbalance")
    deterministic_status = tuple(
        (ds, witness[:-1]) for ds, witness in sorted(status.items())
    )
    partition = (
        tuple(sorted(states.items())),
        tuple(sorted(crude.items())),
        deterministic_status,
        tuple(survivors),
    )
    return (
        body,
        stream.h,
        len(stream.carrier),
        stream.L,
        stream.high_floor,
        stream.first_d,
        checks,
        tuple(sorted(sign_totals.items())),
        len(ray.support.divisors(stream.L)) - 1,
        trials,
        partition,
        tuple(sorted(Counter(w[1] for w in status.values()).items())),
        max(states.items(), key=lambda item: (item[1]["excess"], item[0]), default=None),
        hashlib.sha256(repr(partition).encode()).hexdigest(),
        verified_farkas_checks,
        all_negative,
        tuple(canonical_records),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    pair_rows, control_marginals, control_caps, control_certificate = ray.local.controls()
    records = tuple(evaluate_body(body) for body in BODIES)
    totals = tuple(
        sum(len(record[10][index]) for record in records) for index in range(4)
    )
    require(totals == EXPECTED_TOTALS, ("global counts changed", totals))
    require(
        sum(record[14] for record in records) == totals[2],
        "exact Farkas check count changed",
    )
    require(all(record[15] for record in records), "nonnegative Farkas replay")
    global_m = Counter()
    for record in records:
        global_m.update(dict(record[11]))
    require(tuple(sorted(global_m.items())) == EXPECTED_M_HISTOGRAM, "M histogram changed")
    circuit_census = Counter()
    tariff_census = Counter()
    instance_digest = hashlib.sha256()
    all_canonical_records = []
    for record in records:
        for canonical in record[16]:
            circuit_census[canonical[3]] += 1
            tariff_census[
                (
                    canonical[3],
                    canonical[4],
                    canonical[5],
                    canonical[6],
                    canonical[7],
                    canonical[8],
                )
            ] += 1
            instance_digest.update(f"{canonical}\n".encode())
            all_canonical_records.append(canonical)
    require(
        tuple(sorted(circuit_census.items())) == EXPECTED_CIRCUIT_CENSUS,
        ("canonical circuit census changed", circuit_census),
    )
    require(
        len(tariff_census) == EXPECTED_DISTINCT_TARIFFS,
        ("canonical tariff census changed", len(tariff_census)),
    )
    require(
        instance_digest.hexdigest() == EXPECTED_INSTANCE_TARIFF_SHA256,
        ("deterministic instance/tariff digest changed", instance_digest.hexdigest()),
    )
    representative = all_canonical_records[0]
    semantic_payload = (
        FIRST,
        BODIES,
        records,
        totals,
        tuple(sorted(global_m.items())),
        pair_rows,
        control_marginals,
        control_caps,
        tuple(sorted(circuit_census.items())),
        tuple(sorted(tariff_census.items(), key=repr)),
        instance_digest.hexdigest(),
        representative,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=328 all-label ray/status closure",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        "independent_audit=basis-independent minimal barcode/canonical tariff plus internal exact matrix replay",
        f"scalar_source_sha256={file_sha256(SCALAR_PATH)}",
        f"scalar_output_sha256={file_sha256(SCALAR_OUTPUT_PATH)}",
        "scope=9 globally reconstructed scalar bodies;no finite label horizon",
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "A(L-b)=-A(b);all denominator maxima attained"
        ),
        f"pair_overlap_exhaustive_controls={pair_rows}",
        (
            f"totals=states:{totals[0]};crude_kills:{totals[1]};"
            f"status_kills:{totals[2]};survivors:{totals[3]}"
        ),
        f"status_M_histogram={dict(sorted(global_m.items()))}",
        f"minimal_barcode_circuit_census={dict(sorted(circuit_census.items()))}",
        f"distinct_canonical_tariffs={len(tariff_census)}",
        f"deterministic_instance_tariff_sha256={instance_digest.hexdigest()}",
    ]
    for record in records:
        (
            body,
            h,
            components,
            L,
            high_floor,
            first_d,
            checks,
            signs,
            denominator_classes,
            trials,
            partition,
            status_histogram,
            maximum,
            partition_digest,
            verified_checks,
            all_negative,
            canonical_records,
        ) = record
        states, crude, status, survivors = partition
        local_circuits = Counter(canonical[3] for canonical in canonical_records)
        lines.append(
            f"E={body};h={h};r={components};L={L};high={high_floor};d1={first_d};"
            f"ray_checks={checks};ray_signs={dict(signs)};"
            f"denominator_classes={denominator_classes};denominator_trials={trials};"
            f"states={len(states)};crude_kills={len(crude)};"
            f"status_kills={len(status)};survivors={len(survivors)};"
            f"status_M={dict(status_histogram)};max_state={maximum};"
            f"partition_sha256={partition_digest};"
            f"minimal_circuits={dict(sorted(local_circuits.items()))};"
            f"verified_farkas_checks={verified_checks};"
            f"all_farkas_contradictions_negative={int(all_negative)}"
        )
    lines.extend(
        (
            f"representative_canonical_barcode_tariff={representative}",
            "internal_solver_farkas_replay=49/49:PASS;solver_bases_not_frozen",
            "conclusion=all 9 projected k=3 z1=328 scalar rows are empty",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
