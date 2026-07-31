#!/usr/bin/env python3
"""Compose common-u activity with all-q fibre profiles on the 500-row ledger.

The earlier screens were intentionally independent:

* common-u status bounded which small denominators can be active;
* q-profile status granted every denominator its enlarged trace.

Here they are joined.  For a fixed body row and denominator shape, and for
each subset E of the d<=5 masks declared active, test every q|D after removing
the complete q-profile of each inactive small mask.  The feasible subsets form
an upward Boolean event A.  Coverage forces the actual common-u activity
status into A throughout the compact aligned carrier.  The exact
fractional-cover value of A, with marginals d/7, must therefore be strictly
larger than 66/91.

This remains an upper relaxation: large masks keep their full enlarged
profiles, every q gets an independently optimized real status table, and the
locations of exceptional fibres and lifted unit APs are forgotten.
"""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
QPROFILE_PATH = (
    ROOT
    / ".scratch"
    / "lrc14_k2_five_drift_pattern_compression_20260729"
    / "qprofile_hostile.py"
)
spec = spec_from_file_location("lrc14_k2_qprofile_hostile", QPROFILE_PATH)
qp = module_from_spec(spec)
spec.loader.exec_module(qp)
base = qp.base
qcap = qp.qcap
support = base.support
combined = base.combined

EXPECTED_D100 = {
    "control_count": 17,
    "incoming": (500, 120, 21),
    "surviving": (
        60,
        4,
        15,
        15,
        1,
        "be9d41336a4be0a62ef980305b6260c8e68184fa025198c5f0251013251fdc03",
    ),
    "killed_by_D": {28: 3, 56: 341, 84: 96},
    "survived_by_D": {56: 60},
    "surviving_shapes": {
        (56, (7, 8, 8, 8, 8)),
        (56, (8, 8, 8, 8, 14)),
        (56, (8, 8, 8, 8, 28)),
        (56, (8, 8, 8, 8, 56)),
    },
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def maximum_event_mass(minimal_masks, marginals):
    """Exact real maximum of an upward event from its minimal true masks."""
    minimal_masks = tuple(minimal_masks)
    marginals = tuple(marginals)
    count = len(marginals)
    if not minimal_masks:
        return Q(0)
    if 0 in minimal_masks:
        return Q(1)
    rows = tuple(
        tuple((mask >> index) & 1 for index in range(count))
        for mask in minimal_masks
    )
    constraints = [(row, Q(1)) for row in rows]
    constraints.extend(
        (
            tuple(int(index == coordinate) for index in range(count)),
            Q(0),
        )
        for coordinate in range(count)
    )
    optimum = None
    for chosen in combinations(range(len(constraints)), count):
        matrix = tuple(constraints[index][0] for index in chosen)
        rhs = tuple(constraints[index][1] for index in chosen)
        point = base.solve_square(matrix, rhs)
        if point is None or any(value < 0 for value in point):
            continue
        if any(
            sum(value * coefficient for value, coefficient in zip(point, row))
            < 1
            for row in rows
        ):
            continue
        objective = sum(
            marginal * value
            for marginal, value in zip(marginals, point)
        )
        if optimum is None or objective < optimum:
            optimum = objective
    require(optimum is not None, "upward event cover has no vertex")
    return min(Q(1), optimum)


def minimal_true_masks(feasible_masks, count):
    feasible = set(feasible_masks)
    result = []
    for mask in sorted(feasible):
        if all(
            (mask ^ (1 << index)) not in feasible
            for index in range(count)
            if (mask >> index) & 1
        ):
            result.append(mask)
    return tuple(result)


def event_controls():
    checks = 0
    # Complete h-of-m events reproduce the known symmetric value.
    for count in range(1, 6):
        marginals = (Q(4, 7),) * count
        for needed in range(1, count + 1):
            feasible = [
                mask
                for mask in range(1 << count)
                if mask.bit_count() >= needed
            ]
            minimal = minimal_true_masks(feasible, count)
            expected = min(Q(1), Q(4 * count, 7 * needed))
            require(
                maximum_event_mass(minimal, marginals) == expected,
                ("symmetric event control changed", count, needed),
            )
            checks += 1
    # A heterogeneous union and intersection.
    marginals = (Q(2, 7), Q(5, 7))
    require(maximum_event_mass((1, 2), marginals) == 1, "union control changed")
    require(maximum_event_mass((3,), marginals) == Q(2, 7), "intersection control changed")
    return checks + 2


def row_ledger(max_D):
    rows_by_D = defaultdict(list)
    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            if D > max_D or D not in (28, 56, 84):
                continue
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > base.SUPPORT_CUTOFF:
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
                    qp.target_thresholds(
                        combined.residue_load_histogram(arcs, q)
                    ),
                )
                for q in support.divisors(D)
            )
            rows_by_D[D].append(
                (support_count, body, L, small_loads, q_thresholds)
            )
    return rows_by_D


def denominator_q_components(D, values):
    """Per q: always-large profile and labelled optional-small components."""
    result = []
    for q in support.divisors(D):
        large_baseline = 0
        large_profile = []
        small_components = []
        small_denominators = []
        for d in values:
            low, increment, high_count = qp.q_type(D, d, q)
            component = (low, increment, high_count)
            if d in base.SMALL_DIVISORS:
                small_components.append(component)
                small_denominators.append(d)
            else:
                large_baseline += low
                if increment:
                    large_profile.append((increment, high_count))
        result.append(
            (
                q,
                large_baseline,
                tuple(sorted(large_profile)),
                tuple(small_components),
                tuple(small_denominators),
            )
        )
    return tuple(result)


@lru_cache(maxsize=None)
def subset_profile_passes(q, large_baseline, large_profile, small_components, mask, thresholds):
    baseline = large_baseline
    profile = list(large_profile)
    for index, (low, increment, high_count) in enumerate(small_components):
        if not (mask >> index) & 1:
            continue
        baseline += low
        if increment:
            profile.append((increment, high_count))
    profile.sort()
    return qp.profile_passes(q, baseline, tuple(profile), thresholds)


def main(max_D):
    control_count = event_controls()
    rows_by_D = row_ledger(max_D)
    original_occurrences = 0
    original_shapes = set()
    original_rows = set()
    composed_occurrences = 0
    composed_shapes = set()
    composed_rows = set()
    composed_bodies = set()
    composed_divisors = set()
    killed_by_D = Counter()
    survived_by_D = Counter()
    killed_by_pattern = Counter()
    survived_by_pattern = Counter()
    event_types = Counter()
    semantic = sha256()
    minimal_survivor = None
    minimal_kill = None

    for D in sorted(rows_by_D):
        q_values = tuple(support.divisors(D))
        for values in qcap.denominator_shapes(D):
            pattern = tuple(values.count(d) for d in base.SMALL_DIVISORS)
            small_denominators = tuple(
                d for d in values if d in base.SMALL_DIVISORS
            )
            large_capacity = sum(
                (D // d) * ((d + 6) // 7)
                for d in values
                if d >= 6
            )
            components = denominator_q_components(D, values)
            require(
                tuple(item[0] for item in components) == q_values,
                "q component alphabet changed",
            )
            for support_count, body, L, small_loads, q_thresholds in rows_by_D[D]:
                local_weights, local_marginals = base.small_vectors(
                    pattern, D, small_loads
                )
                if (
                    large_capacity
                    + base.status_load_limit(local_weights, local_marginals)
                    < support_count
                ):
                    continue

                # Reproduce the 500 occurrence-level real q-profile ledger.
                original_pass = True
                for (q, thresholds), (
                    _q,
                    large_baseline,
                    large_profile,
                    small_components,
                    _small_denominators,
                ) in zip(q_thresholds, components):
                    full_mask = (1 << len(small_components)) - 1
                    if not subset_profile_passes(
                        q,
                        large_baseline,
                        large_profile,
                        small_components,
                        full_mask,
                        thresholds,
                    ):
                        original_pass = False
                        break
                if not original_pass:
                    continue
                original_occurrences += 1
                original_shapes.add((D, values))
                original_rows.add((D, body))

                count = len(small_denominators)
                feasible_masks = []
                for mask in range(1 << count):
                    feasible = True
                    for (q, thresholds), (
                        _q,
                        large_baseline,
                        large_profile,
                        small_components,
                        _small_denominators,
                    ) in zip(q_thresholds, components):
                        if not subset_profile_passes(
                            q,
                            large_baseline,
                            large_profile,
                            small_components,
                            mask,
                            thresholds,
                        ):
                            feasible = False
                            break
                    if feasible:
                        feasible_masks.append(mask)
                for mask in feasible_masks:
                    for index in range(count):
                        require(
                            (mask | (1 << index)) in feasible_masks,
                            ("composed event is not upward", D, values, mask),
                        )
                minimal = minimal_true_masks(feasible_masks, count)
                marginals = tuple(Q(d, 7) for d in small_denominators)
                maximum_mass = maximum_event_mass(minimal, marginals)
                event_key = (
                    tuple(small_denominators),
                    minimal,
                    maximum_mass,
                )
                event_types[event_key] += 1
                survives = maximum_mass > base.TWO_SAFE_FLOOR
                record = (
                    D,
                    body,
                    L,
                    support_count,
                    values,
                    minimal,
                    maximum_mass,
                )
                if survives:
                    composed_occurrences += 1
                    composed_shapes.add((D, values))
                    composed_rows.add((D, body))
                    composed_bodies.add(body)
                    composed_divisors.add(D)
                    survived_by_D[D] += 1
                    survived_by_pattern[pattern] += 1
                    semantic.update(f"{record}\n".encode())
                    if minimal_survivor is None or record < minimal_survivor:
                        minimal_survivor = record
                else:
                    killed_by_D[D] += 1
                    killed_by_pattern[pattern] += 1
                    if minimal_kill is None or record < minimal_kill:
                        minimal_kill = record

    require(
        (original_occurrences, len(original_shapes), len(original_rows))
        == EXPECTED_D100["incoming"],
        (
            "incoming 500-ledger changed",
            original_occurrences,
            len(original_shapes),
            len(original_rows),
        ),
    )
    if max_D == 100:
        require(
            control_count == EXPECTED_D100["control_count"],
            "event control count changed",
        )
        actual_surviving = (
            composed_occurrences,
            len(composed_shapes),
            len(composed_rows),
            len(composed_bodies),
            len(composed_divisors),
            semantic.hexdigest(),
        )
        require(
            actual_surviving == EXPECTED_D100["surviving"],
            ("composed D100 ledger changed", actual_surviving),
        )
        require(
            dict(killed_by_D) == EXPECTED_D100["killed_by_D"],
            ("composed kill divisors changed", killed_by_D),
        )
        require(
            dict(survived_by_D) == EXPECTED_D100["survived_by_D"],
            ("composed survivor divisors changed", survived_by_D),
        )
        require(
            composed_shapes == EXPECTED_D100["surviving_shapes"],
            ("composed survivor shapes changed", composed_shapes),
        )

    print("LRC14 k=2/p=5 composed common-u x all-q profile scout")
    print(f"max_D={max_D}")
    print(f"event_cover_control_count={control_count}")
    print(
        f"incoming_occurrences={original_occurrences},"
        f"incoming_shapes={len(original_shapes)},"
        f"incoming_rows={len(original_rows)}"
    )
    print(
        f"surviving_occurrences={composed_occurrences},"
        f"surviving_shapes={len(composed_shapes)},"
        f"surviving_rows={len(composed_rows)},"
        f"surviving_bodies={len(composed_bodies)},"
        f"surviving_divisors={len(composed_divisors)},"
        f"semantic_sha256={semantic.hexdigest()}"
    )
    print(f"killed_by_D={killed_by_D}")
    print(f"survived_by_D={survived_by_D}")
    print(f"killed_by_pattern={killed_by_pattern}")
    print(f"survived_by_pattern={survived_by_pattern}")
    print(f"event_types={event_types}")
    print(f"minimal_kill={minimal_kill}")
    print(f"minimal_survivor={minimal_survivor}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--max-D", type=int, default=100)
    arguments = parser.parse_args()
    main(arguments.max_D)
