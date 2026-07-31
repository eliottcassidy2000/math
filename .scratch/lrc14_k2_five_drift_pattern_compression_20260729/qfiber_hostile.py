#!/usr/bin/env python3
"""Bounded all-divisor fibre-cap hostile for the k=2/p=5 compression.

Unlike ``scout.py``'s scalable divisor-Mobius aggregation, this sidecar
materializes only the 3,680 denominator shapes occurring in the hostile
universe D<=100.  It intersects the exact support-status blocker with every
phase-free q-fibre cap, q|D.

For d,q|D, ell=ceil(d/7), g=gcd(d,q), H=D/lcm(d,q), an enlarged d-needle has
at most

    H ceil(ell/g)

points in any q-fibre.  Therefore a target row with largest q-fibre load T_q
requires

    T_q <= sum_i H_i ceil(ell_i/g_i)

for every divisor q.  This is an upper relaxation: it grants each needle its
largest fibre independently and does not retain which fibre is exceptional.
"""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE_PATH = (
    ROOT
    / ".scratch"
    / "lrc14_k2_five_drift_pattern_compression_20260729"
    / "scout.py"
)
spec = spec_from_file_location("lrc14_k2_five_drift_scout", BASE_PATH)
base = module_from_spec(spec)
spec.loader.exec_module(base)
support = base.support
combined = base.combined

EXPECTED_D100 = {
    "support_rows_by_D": {14: 1, 28: 6, 42: 14, 56: 15, 70: 23, 84: 65, 98: 23},
    "shape_counts_by_D": {
        14: 19,
        28: 100,
        42: 402,
        56: 321,
        70: 402,
        84: 2336,
        98: 100,
    },
    "stages": {
        "raw": (
            174448,
            3680,
            147,
            78,
            7,
            "3c7062cd459ee8044361bf3d6138aa292a4b9a72266e83ddaf0864e511b3a0b5",
        ),
        "scalar": (
            117415,
            3369,
            147,
            78,
            7,
            "f89b6aab1a702aa0b2f578fcc182b0a7c8788339a013f924585c14de24943eb9",
        ),
        "support_status": (
            6089,
            1097,
            53,
            38,
            5,
            "cb318bc8bba79a3e1bc493dbfe2d4d7bfddf8e010a45b2f5cb6858a2bf056113",
        ),
        "all_q_cap": (
            116338,
            3146,
            147,
            78,
            7,
            "44511f3639202ddc844c27a7d1c12db69d1f138fec813360136294a9313799c1",
        ),
        "status_and_all_q": (
            5059,
            874,
            53,
            38,
            5,
            "540daeee9f92478a5e01d9eac59225931870f571ff95119bc21572eaefd749dc",
        ),
    },
    "decisive_q": {(42, 7): 3, (56, 4): 35, (84, 7): 53, (84, 12): 905, (84, 14): 34},
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def q_fibre_cap(D, d, q):
    require(D % d == 0 and D % q == 0, "nondivisor q-cap request")
    ell = (d + 6) // 7
    common = gcd(d, q)
    height = D // lcm(d, q)
    return height * ((ell + common - 1) // common)


@lru_cache(maxsize=None)
def denominator_shapes(D):
    alphabet = tuple(d for d in support.divisors(D) if d > 1)
    result = tuple(
        values
        for values in combinations_with_replacement(alphabet, base.ARITY)
        if lcm(*values) == D
    )
    require(
        len(result) == base.lcm_multiset_shapes(D, base.ARITY),
        ("bounded shape enumeration changed", D),
    )
    return result


def exact_enlarged_needle(D, d, phase, step):
    require(D % d == 0, "control denominator does not divide D")
    require(gcd(step, d) == 1, "control step is not a unit")
    width = (d + 6) // 7
    classes = {
        (phase + index * step) % d
        for index in range(width)
    }
    return tuple(x for x in range(D) if x % d in classes)


def fibre_cap_controls():
    """Brute the sharp cap for every D<=60, d|D, q|D, and unit step."""
    cases = 0
    for D in range(1, 61):
        divisors = support.divisors(D)
        for d in (value for value in divisors if value > 1):
            for step in range(1, d):
                if gcd(step, d) != 1:
                    continue
                needle = exact_enlarged_needle(D, d, 0, step)
                for q in divisors:
                    loads = Counter(x % q for x in needle)
                    actual = max(loads.values())
                    expected = q_fibre_cap(D, d, q)
                    require(
                        actual == expected,
                        ("sharp q-fibre cap failed", D, d, q, step, actual, expected),
                    )
                    cases += 1
    return cases


def stage_summary(max_D):
    rows_by_D = defaultdict(list)
    body_count = 0
    body_divisor_rows = 0
    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            body_divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > base.SUPPORT_CUTOFF or D > max_D:
                continue
            arcs = combined.projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support projection changed", body, D),
            )
            small_loads = {}
            for divisor in base.SMALL_DIVISORS:
                if D % divisor == 0:
                    histogram = combined.residue_load_histogram(arcs, divisor)
                    small_loads[divisor] = combined.top_class_load(histogram, 1)
            q_maxima = tuple(
                (
                    q,
                    max(
                        load
                        for load, _count in combined.residue_load_histogram(
                            arcs, q
                        )
                    ),
                )
                for q in support.divisors(D)
            )
            rows_by_D[D].append(
                (support_count, body, L, small_loads, q_maxima)
            )

    require(body_count == 3003, "body universe changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")

    stage_names = (
        "raw",
        "scalar",
        "support_status",
        "all_q_cap",
        "status_and_all_q",
    )
    occurrences = Counter()
    shapes = Counter()
    row_sets = {stage: set() for stage in stage_names}
    body_sets = {stage: set() for stage in stage_names}
    divisor_sets = {stage: set() for stage in stage_names}
    divisor_occurrences = {stage: Counter() for stage in stage_names}
    pattern_occurrences = {stage: Counter() for stage in stage_names}
    semantic = {stage: sha256() for stage in stage_names}
    decisive_q = Counter()
    weighted_distribution_controls = 0

    for D in sorted(rows_by_D):
        divisors = support.divisors(D)
        shape_features = []
        literal_distribution = Counter()
        for values in denominator_shapes(D):
            pattern = tuple(values.count(d) for d in base.SMALL_DIVISORS)
            ambient_weights, marginals = base.small_vectors(pattern, D)
            large_capacity = sum(
                (D // d) * ((d + 6) // 7)
                for d in values
                if d >= 6
            )
            total_capacity = large_capacity + sum(ambient_weights)
            literal_distribution[(pattern, large_capacity)] += 1
            q_caps = tuple(
                sum(q_fibre_cap(D, d, q) for d in values)
                for q in divisors
            )
            shape_features.append(
                (values, pattern, marginals, large_capacity, total_capacity, q_caps)
            )
        require(
            literal_distribution == base.shape_distribution(D),
            ("independent weighted exact-lcm distribution failed", D),
        )
        weighted_distribution_controls += 1

        shape_survives = {stage: set() for stage in stage_names}
        for support_count, body, L, small_loads, q_maxima in rows_by_D[D]:
            require(
                tuple(q for q, _load in q_maxima) == tuple(divisors),
                "q alphabet changed",
            )
            for (
                values,
                pattern,
                marginals,
                large_capacity,
                total_capacity,
                q_caps,
            ) in shape_features:
                local_weights, local_marginals = base.small_vectors(
                    pattern, D, small_loads
                )
                require(local_marginals == marginals, "status marginals changed")
                local_limit = base.status_load_limit(local_weights, marginals)
                passes = {
                    "raw": True,
                    "scalar": total_capacity >= support_count,
                    "support_status": (
                        large_capacity + local_limit >= support_count
                    ),
                    "all_q_cap": all(
                        target_load <= cap
                        for (_q, target_load), cap in zip(q_maxima, q_caps)
                    ),
                }
                passes["status_and_all_q"] = (
                    passes["support_status"] and passes["all_q_cap"]
                )
                if passes["support_status"] and not passes["all_q_cap"]:
                    bad = [
                        (q, target_load - cap)
                        for (q, target_load), cap in zip(q_maxima, q_caps)
                        if target_load > cap
                    ]
                    maximum_defect = max(defect for _q, defect in bad)
                    first_q = min(
                        q for q, defect in bad if defect == maximum_defect
                    )
                    decisive_q[(D, first_q)] += 1
                for stage in stage_names:
                    if not passes[stage]:
                        continue
                    occurrences[stage] += 1
                    shape_survives[stage].add(values)
                    row_key = (body, D)
                    row_sets[stage].add(row_key)
                    body_sets[stage].add(body)
                    divisor_sets[stage].add(D)
                    divisor_occurrences[stage][D] += 1
                    pattern_occurrences[stage][pattern] += 1
                    semantic[stage].update(
                        (
                            f"{body}|{L}|{D}|{support_count}|{values}\n"
                        ).encode()
                    )
        for stage in stage_names:
            shapes[stage] += len(shape_survives[stage])

    return {
        "rows_by_D": rows_by_D,
        "stages": stage_names,
        "occurrences": occurrences,
        "shapes": shapes,
        "rows": row_sets,
        "bodies": body_sets,
        "divisors": divisor_sets,
        "divisor_occurrences": divisor_occurrences,
        "pattern_occurrences": pattern_occurrences,
        "semantic": semantic,
        "decisive_q": decisive_q,
        "weighted_distribution_controls": weighted_distribution_controls,
    }


def main(max_D):
    control_cases = fibre_cap_controls()
    result = stage_summary(max_D)
    if max_D == 100:
        actual_rows = {
            D: len(rows) for D, rows in result["rows_by_D"].items()
        }
        actual_shapes = {
            D: len(denominator_shapes(D)) for D in result["rows_by_D"]
        }
        require(
            actual_rows == EXPECTED_D100["support_rows_by_D"],
            ("D100 support-row ledger changed", actual_rows),
        )
        require(
            actual_shapes == EXPECTED_D100["shape_counts_by_D"],
            ("D100 shape ledger changed", actual_shapes),
        )
        for stage in result["stages"]:
            actual = (
                result["occurrences"][stage],
                result["shapes"][stage],
                len(result["rows"][stage]),
                len(result["bodies"][stage]),
                len(result["divisors"][stage]),
                result["semantic"][stage].hexdigest(),
            )
            require(
                actual == EXPECTED_D100["stages"][stage],
                ("D100 q-stage ledger changed", stage, actual),
            )
        require(
            dict(result["decisive_q"]) == EXPECTED_D100["decisive_q"],
            ("D100 decisive-q ledger changed", result["decisive_q"]),
        )
    print("LRC14 k=2/p=5 bounded all-divisor q-fibre hostile")
    print(f"max_D={max_D}")
    print(f"sharp_fibre_cap_control_cases={control_cases}")
    print(
        "independent_weighted_distribution_controls="
        f"{result['weighted_distribution_controls']}"
    )
    print(
        "support_rows_by_D="
        f"{Counter({D: len(rows) for D, rows in result['rows_by_D'].items()})}"
    )
    print(
        "shape_counts_by_D="
        f"{Counter({D: len(denominator_shapes(D)) for D in result['rows_by_D']})}"
    )
    for stage in result["stages"]:
        print(
            f"stage={stage},occurrences={result['occurrences'][stage]},"
            f"shapes={result['shapes'][stage]},"
            f"rows={len(result['rows'][stage])},"
            f"bodies={len(result['bodies'][stage])},"
            f"divisors={len(result['divisors'][stage])},"
            f"semantic_sha256={result['semantic'][stage].hexdigest()}"
        )
        print(
            f"stage={stage},divisor_occurrences="
            f"{Counter(dict(sorted(result['divisor_occurrences'][stage].items())))}"
        )
        print(
            f"stage={stage},pattern_occurrences="
            f"{Counter(dict(sorted(result['pattern_occurrences'][stage].items())))}"
        )
    print(f"decisive_q_histogram={result['decisive_q']}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--max-D", type=int, default=100)
    arguments = parser.parse_args()
    main(arguments.max_D)
