#!/usr/bin/env python3
"""Bounded exact q-profile hostile for two aligned and five drifting combs.

This is the next relaxation after the small-denominator common-u blocker.
For each q|D it keeps the complete two-level fibre profile of every enlarged
needle and applies the upward-event/fractional-cover theorem to every target
load threshold.  The q-fibre status law is a discrete counting law, so unlike
the compact/open common-u test, equality with the real cover bound is allowed.

The census is intentionally bounded to D<=100 and enumerates only 3,680 exact
denominator shapes.  It does not claim a global k=2 census.
"""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[2]
QCAP_PATH = (
    ROOT
    / ".scratch"
    / "lrc14_k2_five_drift_pattern_compression_20260729"
    / "qfiber_hostile.py"
)
spec = spec_from_file_location("lrc14_k2_qcap_hostile", QCAP_PATH)
qcap = module_from_spec(spec)
spec.loader.exec_module(qcap)
base = qcap.base
support = base.support
combined = base.combined

EXPECTED_D100 = {
    "positive_checks": 12,
    "positive_union_size": 49,
    "stages": {
        "support_status": (
            6089,
            1097,
            53,
            38,
            5,
            "cb318bc8bba79a3e1bc493dbfe2d4d7bfddf8e010a45b2f5cb6858a2bf056113",
        ),
        "status_and_q_profile": (
            500,
            120,
            21,
            15,
            3,
            "8a9df03e38f100dcc6f833bd8d325277b161fc1ea4ce3c51ce7746c4a4d48b55",
        ),
    },
    "final_divisor_occurrences": {28: 3, 56: 401, 84: 96},
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def q_type(D, d, q):
    ell = (d + 6) // 7
    common = gcd(d, q)
    quotient, remainder = divmod(ell, common)
    height = D // lcm(d, q)
    low = height * quotient
    increment = height if remainder else 0
    high_count = (q // common) * remainder
    require(0 <= high_count <= q, "high-count marginal changed")
    return low, increment, high_count


def target_thresholds(histogram):
    result = []
    running = 0
    for load, count in reversed(histogram):
        running += count
        if load:
            result.append((load, running))
    return tuple(result)


@lru_cache(maxsize=None)
def profile_passes(q, baseline, profile, thresholds):
    increments = tuple(item[0] for item in profile)
    marginals = tuple(Q(item[1], q) for item in profile)
    for target_load, target_count in thresholds:
        need = target_load - baseline
        if need <= 0:
            continue
        maximum_mass = base.maximum_upward_mass(
            increments, marginals, need
        )
        if Q(target_count) > q * maximum_mass:
            return False
    return True


def profile_slack(q, baseline, profile, thresholds):
    """Minimum slack over nontrivial thresholds needing exceptional load."""
    increments = tuple(item[0] for item in profile)
    marginals = tuple(Q(item[1], q) for item in profile)
    slacks = []
    for target_load, target_count in thresholds:
        need = target_load - baseline
        if need <= 0:
            continue
        maximum_mass = base.maximum_upward_mass(increments, marginals, need)
        slacks.append(q * maximum_mass - target_count)
    return min(slacks) if slacks else None


@lru_cache(maxsize=None)
def integer_status_law(q, baseline, profile, thresholds):
    """Exact integer q-column status law, or None when infeasible."""
    increments = tuple(item[0] for item in profile)
    high_counts = tuple(item[1] for item in profile)
    count = len(increments)
    require(all(0 <= value <= q for value in high_counts), "bad high count")
    model = cp_model.CpModel()
    state_counts = [
        model.NewIntVar(0, q, f"n_{state}")
        for state in range(1 << count)
    ]
    model.Add(sum(state_counts) == q)
    for index, high_count in enumerate(high_counts):
        model.Add(
            sum(
                state_counts[state]
                for state in range(1 << count)
                if (state >> index) & 1
            )
            == high_count
        )
    state_loads = tuple(
        sum(
            increments[index]
            for index in range(count)
            if (state >> index) & 1
        )
        for state in range(1 << count)
    )
    for target_load, target_count in thresholds:
        need = target_load - baseline
        if need <= 0:
            continue
        heavy_states = [
            state_counts[state]
            for state, load in enumerate(state_loads)
            if load >= need
        ]
        if not heavy_states:
            return None
        model.Add(sum(heavy_states) >= target_count)

    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    status = solver.Solve(model)
    if status == cp_model.INFEASIBLE:
        return None
    require(
        status in (cp_model.FEASIBLE, cp_model.OPTIMAL),
        ("integer status solver returned unknown", status),
    )
    law = tuple(solver.Value(variable) for variable in state_counts)
    require(sum(law) == q, "integer law total changed")
    for index, high_count in enumerate(high_counts):
        require(
            sum(
                law[state]
                for state in range(1 << count)
                if (state >> index) & 1
            )
            == high_count,
            ("integer law marginal changed", index),
        )
    for target_load, target_count in thresholds:
        need = target_load - baseline
        if need <= 0:
            continue
        require(
            sum(
                law[state]
                for state, load in enumerate(state_loads)
                if load >= need
            )
            >= target_count,
            ("integer law target threshold changed", target_load),
        )
    return law


def literal_profile_control():
    D = 84
    denominators = (6, 7, 12, 14, 21)
    phases = (1, 2, 3, 4, 5)
    steps = (1, 1, 5, 3, 2)
    require(lcm(*denominators) == D, "positive-control lcm changed")
    union = set()
    for d, phase, step in zip(denominators, phases, steps):
        union.update(qcap.exact_enlarged_needle(D, d, phase, step))
    checks = 0
    for q in support.divisors(D):
        baseline = 0
        profile = []
        for d in denominators:
            low, increment, high_count = q_type(D, d, q)
            baseline += low
            if increment:
                profile.append((increment, high_count))
        profile.sort()
        histogram = Counter(Counter(point % q for point in union).values())
        histogram.setdefault(0, q - sum(histogram.values()))
        thresholds = target_thresholds(tuple(sorted(histogram.items())))
        require(
            profile_passes(q, baseline, tuple(profile), thresholds),
            ("literal positive control rejected", q),
        )
        require(
            integer_status_law(
                q, baseline, tuple(profile), thresholds
            )
            is not None,
            ("literal integer positive control rejected", q),
        )
        checks += 1
    return checks, len(union)


def main(max_D):
    positive_checks, positive_union_size = literal_profile_control()
    rows_by_D = defaultdict(list)
    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > base.SUPPORT_CUTOFF or D > max_D:
                continue
            arcs = combined.projected_support_arcs(D, ranges)
            small_loads = {}
            for divisor in base.SMALL_DIVISORS:
                if D % divisor == 0:
                    histogram = combined.residue_load_histogram(arcs, divisor)
                    small_loads[divisor] = combined.top_class_load(histogram, 1)
            q_thresholds = tuple(
                (
                    q,
                    target_thresholds(
                        combined.residue_load_histogram(arcs, q)
                    ),
                )
                for q in support.divisors(D)
            )
            rows_by_D[D].append(
                (support_count, body, L, small_loads, q_thresholds)
            )

    stage_names = (
        "support_status",
        "status_and_q_profile",
        "status_q_profile_integer",
    )
    occurrences = Counter()
    shape_sets = {stage: set() for stage in stage_names}
    row_sets = {stage: set() for stage in stage_names}
    body_sets = {stage: set() for stage in stage_names}
    divisor_sets = {stage: set() for stage in stage_names}
    divisor_occurrences = {stage: Counter() for stage in stage_names}
    pattern_occurrences = {stage: Counter() for stage in stage_names}
    semantic = {stage: sha256() for stage in stage_names}
    first_failed_q = Counter()
    final_rows = Counter()
    final_shapes = Counter()
    final_slacks = Counter()
    final_tight_q = Counter()
    integer_first_failure = Counter()
    integer_survivor_rows = Counter()
    integer_survivor_laws = Counter()

    for D in sorted(rows_by_D):
        divisors = tuple(support.divisors(D))
        features = []
        for values in qcap.denominator_shapes(D):
            pattern = tuple(values.count(d) for d in base.SMALL_DIVISORS)
            large_capacity = sum(
                (D // d) * ((d + 6) // 7)
                for d in values
                if d >= 6
            )
            profiles = []
            for q in divisors:
                baseline = 0
                profile = []
                for d in values:
                    low, increment, high_count = q_type(D, d, q)
                    baseline += low
                    if increment:
                        profile.append((increment, high_count))
                profiles.append((baseline, tuple(sorted(profile))))
            features.append(
                (values, pattern, large_capacity, tuple(profiles))
            )

        for support_count, body, L, small_loads, q_thresholds in rows_by_D[D]:
            require(
                tuple(q for q, _thresholds in q_thresholds) == divisors,
                "q-threshold alphabet changed",
            )
            for values, pattern, large_capacity, profiles in features:
                local_weights, marginals = base.small_vectors(
                    pattern, D, small_loads
                )
                local_limit = base.status_load_limit(local_weights, marginals)
                support_status = (
                    large_capacity + local_limit >= support_count
                )
                failures = []
                if support_status:
                    tests = list(zip(q_thresholds, profiles))
                    tests.sort(
                        key=lambda test: (
                            test[1][0]
                            + sum(increment for increment, _count in test[1][1])
                            - test[0][1][0][0],
                            test[0][0],
                        )
                    )
                    for (q, thresholds), (baseline, profile) in tests:
                        if not profile_passes(q, baseline, profile, thresholds):
                            failures.append(q)
                            break
                stage_pass = {
                    "support_status": support_status,
                    "status_and_q_profile": support_status and not failures,
                }
                integer_failures = []
                integer_laws = []
                if support_status and not failures:
                    for (q, thresholds), (baseline, profile) in zip(
                        q_thresholds, profiles
                    ):
                        law = integer_status_law(
                            q, baseline, profile, thresholds
                        )
                        if law is None:
                            integer_failures.append(q)
                            break
                        integer_laws.append((q, profile, law))
                stage_pass["status_q_profile_integer"] = (
                    support_status and not failures and not integer_failures
                )
                if support_status and failures:
                    first_failed_q[(D, min(failures))] += 1
                if support_status and not failures:
                    slacks = []
                    for (q, thresholds), (baseline, profile) in zip(
                        q_thresholds, profiles
                    ):
                        slack = profile_slack(
                            q, baseline, profile, thresholds
                        )
                        if slack is not None:
                            slacks.append((q, slack))
                    require(slacks, "real survivor has no nontrivial q threshold")
                    minimum_slack = min(slack for _q, slack in slacks)
                    final_slacks[minimum_slack] += 1
                    for q, slack in slacks:
                        if slack == minimum_slack:
                            final_tight_q[(D, q)] += 1
                    final_rows[(D, body, L, support_count)] += 1
                    final_shapes[(D, values)] += 1
                    if integer_failures:
                        integer_first_failure[(D, integer_failures[0])] += 1
                    else:
                        integer_survivor_rows[
                            (D, body, L, support_count)
                        ] += 1
                        for q, profile, law in integer_laws:
                            integer_survivor_laws[
                                (
                                    D,
                                    q,
                                    profile,
                                    tuple(
                                        (state, count)
                                        for state, count in enumerate(law)
                                        if count
                                    ),
                                )
                            ] += 1
                for stage in stage_names:
                    if not stage_pass[stage]:
                        continue
                    occurrences[stage] += 1
                    shape_sets[stage].add((D, values))
                    row_sets[stage].add((body, D))
                    body_sets[stage].add(body)
                    divisor_sets[stage].add(D)
                    divisor_occurrences[stage][D] += 1
                    pattern_occurrences[stage][pattern] += 1
                    semantic[stage].update(
                        f"{body}|{L}|{D}|{support_count}|{values}\n".encode()
                    )

    if max_D == 100:
        require(
            (positive_checks, positive_union_size)
            == (
                EXPECTED_D100["positive_checks"],
                EXPECTED_D100["positive_union_size"],
            ),
            "literal positive profile control changed",
        )
        for stage in ("support_status", "status_and_q_profile"):
            actual = (
                occurrences[stage],
                len(shape_sets[stage]),
                len(row_sets[stage]),
                len(body_sets[stage]),
                len(divisor_sets[stage]),
                semantic[stage].hexdigest(),
            )
            require(
                actual == EXPECTED_D100["stages"][stage],
                ("D100 q-profile ledger changed", stage, actual),
            )
        require(
            dict(divisor_occurrences["status_and_q_profile"])
            == EXPECTED_D100["final_divisor_occurrences"],
            "D100 final divisor ledger changed",
        )
        require(
            final_slacks == Counter({Q(0): 500}),
            ("D100 equality wall changed", final_slacks),
        )

    print("LRC14 k=2/p=5 bounded exact q-profile hostile")
    print(f"max_D={max_D}")
    print(f"literal_positive_profile_checks={positive_checks}")
    print(f"literal_positive_union_size={positive_union_size}")
    print(
        "support_rows_by_D="
        f"{Counter({D: len(rows) for D, rows in rows_by_D.items()})}"
    )
    for stage in stage_names:
        print(
            f"stage={stage},occurrences={occurrences[stage]},"
            f"shapes={len(shape_sets[stage])},"
            f"rows={len(row_sets[stage])},"
            f"bodies={len(body_sets[stage])},"
            f"divisors={len(divisor_sets[stage])},"
            f"semantic_sha256={semantic[stage].hexdigest()}"
        )
        print(
            f"stage={stage},divisor_occurrences="
            f"{Counter(dict(sorted(divisor_occurrences[stage].items())))}"
        )
        print(
            f"stage={stage},pattern_occurrences="
            f"{Counter(dict(sorted(pattern_occurrences[stage].items())))}"
        )
    print(f"first_failed_q_histogram={first_failed_q}")
    print(f"final_slack_histogram={final_slacks}")
    print(f"final_tight_q_histogram={final_tight_q}")
    print(f"integer_first_failure_histogram={integer_first_failure}")
    print(
        "integer_survivor_rows="
        f"{[(*key, count) for key, count in sorted(integer_survivor_rows.items())]}"
    )
    print(
        "integer_survivor_law_types_top30="
        f"{[(*key, count) for key, count in integer_survivor_laws.most_common(30)]}"
    )
    print(
        "final_rows=(D,body,L,support_count,occurrences)="
        f"{[(*key, count) for key, count in sorted(final_rows.items())]}"
    )
    print(
        "final_shape_occurrence_histogram="
        f"{Counter(final_shapes.values())}"
    )
    print(
        "final_shapes_top30=(D,denominators,rows)="
        f"{[(*key, count) for key, count in sorted(final_shapes.items(), key=lambda item: (-item[1], item[0]))[:30]]}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--max-D", type=int, default=100)
    arguments = parser.parse_args()
    main(arguments.max_D)
