#!/usr/bin/env python3
"""Exact bounded compression scout for two aligned combs and five drifts.

This is a scratch research instrument, not a canonical theorem referee.
It never enumerates the 50,874,159,718 denominator five-multisets, much less
their 951,545,890,235 body-row occurrences.

For a common resolving denominator D and d|D, an enlarged drift needle has
ambient capacity

    C_D(d) = (D/d) ceil(d/7).

When d<7, its exact common-phase activity marginal is d/7.  The safe carrier
of the two aligned combs has mass at least 66/91.  Since 5/7 < 66/91 < 6/7,
denominators 2,3,4,5 are "small": no one of them can be forced active on the
whole aligned carrier.  Denominators at least six are initially granted their
full ambient capacity.

For a multiplicity pattern m=(m_2,m_3,m_4,m_5), let the labelled small-mask
loads be a_i and marginals r_i=d_i/7.  Coverage forces the upward event

    sum_i a_i X_i >= |S_D| - C_large

throughout the aligned carrier.  The exact real transport theorem says

    max Pr(A) = min(1, tau_min(A)(r)),

where tau is the weighted fractional-cover value of the minimal true clutter.
The carrier is compact and the proper activity event is open, so equality at
66/91 is fatal: a surviving threshold must have maximum mass strictly larger.

The divisor-poset generating function is grouped by only

    (m_2,m_3,m_4,m_5, total capacity from d>=6).

For every E|D, divisors in the E-alphabet with a common feature are selected
with repetition, and Mobius inversion extracts lcm exactly D.  Thus every raw
denominator multiset is counted exactly without materializing it.

Stages:
  raw             exact lcm-D denominator multiset universe;
  scalar          total ambient capacity at least |S_D|;
  ambient_status  exact fractional blocker on ambient small-class loads;
  support_status  same blocker with each small weight reduced to the largest
                  exact S_D residue-class load modulo d.

All stages are necessary relaxations.  They grant all d>=6 masks their full
ambient capacity and allow arbitrary joint laws with the correct one-marginals.
They do not assert one common arithmetic status law across thresholds, located
unit-needle phases, or simultaneous q-fibre realizability.
"""

from argparse import ArgumentParser
from bisect import bisect_left
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
COMBINED_PATH = (
    ROOT / "04-computation" / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
EXPECTED_ROWS = 27163
EXPECTED_DIVISORS = 219
EXPECTED_SHAPES = 50874159718
EXPECTED_OCCURRENCES = 951545890235
EXPECTED_D100 = {
    "support_rows": 147,
    "support_divisors": 7,
    "features": 200,
    "stages": {
        "raw": (
            174448,
            3680,
            147,
            78,
            7,
            "b82780f84965bb8b0a3bea24d27a97e1be3eba39cb4c8658a76393127bafb498",
        ),
        "scalar": (
            117415,
            3369,
            147,
            78,
            7,
            "33ba125aed732a69c374b19088c888ee332e89ecaf4a423a7fb55fd0370ec04a",
        ),
        "ambient_status": (
            9972,
            1618,
            110,
            73,
            5,
            "35c2880326404f7579fac1b51d324a7bc1693ae4ab47babc10b2aa4ce1bc3210",
        ),
        "support_status": (
            6089,
            1097,
            53,
            38,
            5,
            "d76475f312471bf542a082cf951565d2c6efed83a65484f1673be589ebb7de09",
        ),
    },
}

TWO_SAFE_FLOOR = Q(66, 91)
FIVE_SAFE_FLOOR = Q(478, 1365)
SUPPORT_CUTOFF = (Q(1) - FIVE_SAFE_FLOOR) / TWO_SAFE_FLOOR
SMALL_DIVISORS = (2, 3, 4, 5)
ARITY = 5


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(COMBINED_PATH) == EXPECTED_COMBINED_SHA256,
    "frozen body-projection dependency changed",
)
require(SUPPORT_CUTOFF == Q(887, 990), "support cutoff changed")
spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


def solve_square(rows, rhs):
    """Solve a square rational system, returning None when singular."""
    size = len(rows)
    matrix = [
        [Q(value) for value in row] + [Q(rhs[index])]
        for index, row in enumerate(rows)
    ]
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            return None
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = matrix[column][column]
        matrix[column] = [value / scale for value in matrix[column]]
        for row in range(size):
            if row == column or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(
                    matrix[row], matrix[column]
                )
            ]
    return tuple(matrix[index][-1] for index in range(size))


@lru_cache(maxsize=None)
def maximum_upward_mass(weights, marginals, need):
    """Exact maximum mass of ``sum a_i X_i >= need``."""
    weights = tuple(weights)
    marginals = tuple(marginals)
    count = len(weights)
    require(len(marginals) == count, "weight/marginal arity changed")
    require(all(weight > 0 for weight in weights), "nonpositive status weight")
    require(all(0 <= marginal <= 1 for marginal in marginals), "bad marginal")
    if need <= 0:
        return Q(1)
    if need > sum(weights):
        return Q(0)
    common_scale = gcd(need, *weights)
    if common_scale > 1:
        return maximum_upward_mass(
            tuple(weight // common_scale for weight in weights),
            marginals,
            need // common_scale,
        )

    minimal = []
    for mask in range(1 << count):
        weight = sum(
            weights[index]
            for index in range(count)
            if (mask >> index) & 1
        )
        if weight < need:
            continue
        if any(
            weight - weights[index] >= need
            for index in range(count)
            if (mask >> index) & 1
        ):
            continue
        minimal.append(
            tuple((mask >> index) & 1 for index in range(count))
        )
    require(minimal, "nonempty threshold event has no minimal state")

    constraints = [(row, Q(1)) for row in minimal]
    constraints.extend(
        (
            tuple(int(index == coordinate) for index in range(count)),
            Q(0),
        )
        for coordinate in range(count)
    )
    optimum = None
    for chosen in combinations(range(len(constraints)), count):
        rows = tuple(constraints[index][0] for index in chosen)
        rhs = tuple(constraints[index][1] for index in chosen)
        point = solve_square(rows, rhs)
        if point is None or any(value < 0 for value in point):
            continue
        if any(
            sum(value * coefficient for value, coefficient in zip(point, row))
            < 1
            for row in minimal
        ):
            continue
        objective = sum(
            marginal * value
            for marginal, value in zip(marginals, point)
        )
        if optimum is None or objective < optimum:
            optimum = objective
    require(optimum is not None, "fractional cover LP has no vertex")
    return min(Q(1), optimum)


@lru_cache(maxsize=None)
def status_load_limit(weights, marginals):
    """Largest attainable active load with relaxed mass strictly above 66/91."""
    weights = tuple(weights)
    marginals = tuple(marginals)
    if not weights:
        return 0
    attainable = sorted(
        {
            sum(
                weights[index]
                for index in range(len(weights))
                if (mask >> index) & 1
            )
            for mask in range(1 << len(weights))
        }
    )
    require(
        maximum_upward_mass(weights, marginals, attainable[0]) > TWO_SAFE_FLOOR,
        "zero threshold must pass",
    )
    low = 0
    high = len(attainable)
    while high - low > 1:
        middle = (low + high) // 2
        if (
            maximum_upward_mass(weights, marginals, attainable[middle])
            > TWO_SAFE_FLOOR
        ):
            low = middle
        else:
            high = middle
    result = attainable[low]
    larger = [load for load in attainable if load > result]
    if larger:
        require(
            maximum_upward_mass(weights, marginals, min(larger))
            <= TWO_SAFE_FLOOR,
            "status threshold is not maximal",
        )
    return result


def activity_controls():
    """Exact direction, equality, and heterogeneous blocker controls."""
    require(
        maximum_upward_mass((1,), (Q(5, 7),), 1) == Q(5, 7),
        "single d5 activity mass changed",
    )
    require(
        maximum_upward_mass((1, 1), (Q(5, 7), Q(5, 7)), 1) == 1,
        "two-d5 union control changed",
    )
    require(
        maximum_upward_mass((1, 1), (Q(5, 7), Q(5, 7)), 2) == Q(5, 7),
        "two-d5 intersection control changed",
    )
    require(
        status_load_limit((1, 1), (Q(5, 7), Q(5, 7))) == 1,
        "strict d5 status boundary changed",
    )
    require(
        maximum_upward_mass((3, 2), (Q(2, 7), Q(3, 7)), 2) == Q(5, 7),
        "heterogeneous union control changed",
    )
    require(
        maximum_upward_mass((3, 2), (Q(2, 7), Q(3, 7)), 5) == Q(2, 7),
        "heterogeneous intersection control changed",
    )
    # This deliberately checks the compact/open strict boundary: 66/91 is
    # not enough to contain a compact carrier of measure at least 66/91.
    require(
        not Q(66, 91) > TWO_SAFE_FLOOR,
        "strict equality boundary changed",
    )
    return 7


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
    if remaining > 1:
        result = -result
    return result


def lcm_multiset_shapes(D, arity):
    """Divisor-lattice Mobius check for the unweighted coefficient."""
    total = 0
    for divisor in support.divisors(D):
        alphabet = len(support.divisors(divisor)) - 1
        total += mobius(D // divisor) * comb(alphabet + arity - 1, arity)
    return total


@lru_cache(maxsize=None)
def shape_distribution(D):
    """Map (small multiplicity pattern, large capacity) to shape count."""
    result = Counter()
    for E in support.divisors(D):
        sign = mobius(D // E)
        if not sign:
            continue
        groups = Counter()
        for divisor in (d for d in support.divisors(E) if d > 1):
            capacity = (D // divisor) * ((divisor + 6) // 7)
            feature = tuple(int(divisor == small) for small in SMALL_DIVISORS)
            feature += (capacity if divisor >= 6 else 0,)
            groups[feature] += 1
        # State: used, m2,m3,m4,m5,large_capacity.
        states = {(0, 0, 0, 0, 0, 0): 1}
        for feature, alphabet_size in groups.items():
            units = feature[:-1]
            unit_large = feature[-1]
            additions = []
            for state, multiplicity in states.items():
                used = state[0]
                pattern = state[1:-1]
                large = state[-1]
                for copies in range(1, ARITY + 1 - used):
                    coefficient = comb(alphabet_size + copies - 1, copies)
                    additions.append(
                        (
                            (
                                used + copies,
                                *(
                                    old + copies * unit
                                    for old, unit in zip(pattern, units)
                                ),
                                large + copies * unit_large,
                            ),
                            multiplicity * coefficient,
                        )
                    )
            for state, multiplicity in additions:
                states[state] = states.get(state, 0) + multiplicity
        for state, multiplicity in states.items():
            used = state[0]
            pattern = state[1:-1]
            large = state[-1]
            if used == ARITY:
                result[(pattern, large)] += sign * multiplicity
    require(
        all(multiplicity >= 0 for multiplicity in result.values()),
        ("negative weighted Mobius coefficient", D),
    )
    result += Counter()
    expected = lcm_multiset_shapes(D, ARITY)
    require(sum(result.values()) == expected, ("shape GF changed", D))
    return result


def brute_shape_controls():
    cases = 0
    for D in range(1, 61):
        alphabet = tuple(d for d in support.divisors(D) if d > 1)
        brute = Counter()
        for values in combinations_with_replacement(alphabet, ARITY):
            if lcm(*values) != D:
                continue
            pattern = tuple(values.count(d) for d in SMALL_DIVISORS)
            large = sum(
                (D // d) * ((d + 6) // 7)
                for d in values
                if d >= 6
            )
            brute[(pattern, large)] += 1
        require(brute == shape_distribution(D), ("weighted GF failed", D))
        cases += 1
    return cases


def small_vectors(pattern, D, row_loads=None):
    """Return labelled small-mask weights and exact activity marginals."""
    weights = []
    marginals = []
    for divisor, copies in zip(SMALL_DIVISORS, pattern):
        if not copies:
            continue
        require(D % divisor == 0, ("small divisor does not divide D", D, divisor))
        if row_loads is None:
            weight = D // divisor
        else:
            require(
                divisor in row_loads,
                ("missing row small load", D, divisor, pattern, row_loads),
            )
            weight = row_loads[divisor]
            require(0 < weight <= D // divisor, ("small load out of range", D, divisor))
        weights.extend([weight] * copies)
        marginals.extend([Q(divisor, 7)] * copies)
    return tuple(weights), tuple(marginals)


def pattern_threshold_table():
    """All exact normalized ambient thresholds for at most five small masks."""
    table = {}
    for count2 in range(ARITY + 1):
        for count3 in range(ARITY + 1 - count2):
            for count4 in range(ARITY + 1 - count2 - count3):
                for count5 in range(ARITY + 1 - count2 - count3 - count4):
                    pattern = (count2, count3, count4, count5)
                    present = [
                        divisor
                        for divisor, copies in zip(SMALL_DIVISORS, pattern)
                        if copies
                    ]
                    scale = lcm(*present) if present else 1
                    weights, marginals = small_vectors(pattern, scale)
                    limit = status_load_limit(weights, marginals)
                    table[pattern] = Q(limit, scale)
                    # Scaling cannot change the status event or its strict
                    # blocker value; check a hostile second scale exactly.
                    weights2, marginals2 = small_vectors(pattern, 2 * scale)
                    require(marginals2 == marginals, "scaled marginals changed")
                    require(
                        status_load_limit(weights2, marginals2) == 2 * limit,
                        ("normalized pattern threshold changed", pattern),
                    )
    require(len(table) == 126, "pattern threshold universe changed")
    return table


def suffix_counter(counter):
    keys = sorted(counter)
    suffix = [0] * (len(keys) + 1)
    for index in range(len(keys) - 1, -1, -1):
        suffix[index] = suffix[index + 1] + counter[keys[index]]
    return keys, suffix


def count_at_least(keys, suffix, threshold):
    return suffix[bisect_left(keys, threshold)]


def main(max_D, progress=False):
    activity_control_count = activity_controls()
    shape_control_cases = brute_shape_controls()
    threshold_table = pattern_threshold_table()

    by_divisor = defaultdict(list)
    body_count = 0
    body_divisor_rows = 0
    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            body_divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > SUPPORT_CUTOFF:
                continue
            if max_D is not None and D > max_D:
                continue
            arcs = combined.projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("projected support mismatch", body, D),
            )
            small_loads = {}
            for divisor in SMALL_DIVISORS:
                if D % divisor == 0:
                    histogram = combined.residue_load_histogram(arcs, divisor)
                    small_loads[divisor] = combined.top_class_load(histogram, 1)
            by_divisor[D].append((support_count, body, L, small_loads))

    require(body_count == 3003, "body universe changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")
    if max_D is None:
        require(sum(map(len, by_divisor.values())) == EXPECTED_ROWS, "k2 row count changed")
        require(len(by_divisor) == EXPECTED_DIVISORS, "k2 divisor alphabet changed")

    stage_names = ("raw", "scalar", "ambient_status", "support_status")
    occurrences = Counter()
    shapes = Counter()
    rows = {stage: set() for stage in stage_names}
    bodies = {stage: set() for stage in stage_names}
    divisors = {stage: set() for stage in stage_names}
    pattern_occurrences = {stage: Counter() for stage in stage_names}
    pattern_shapes = {stage: Counter() for stage in stage_names}
    divisor_occurrences = {stage: Counter() for stage in stage_names}
    semantic = {stage: sha256() for stage in stage_names}
    feature_count = 0

    for D in sorted(by_divisor):
        records = sorted(by_divisor[D], key=lambda row: (row[0], row[1], row[2]))
        distribution = shape_distribution(D)
        by_pattern = defaultdict(Counter)
        for (pattern, large_capacity), multiplicity in distribution.items():
            by_pattern[pattern][large_capacity] += multiplicity
        feature_count += len(distribution)

        for pattern, capacity_counter in by_pattern.items():
            keys, suffix = suffix_counter(capacity_counter)
            ambient_weights, marginals = small_vectors(pattern, D)
            ambient_small_capacity = sum(ambient_weights)
            ambient_limit = status_load_limit(ambient_weights, marginals)
            require(
                Q(ambient_limit, D) == threshold_table[pattern],
                ("ambient pattern threshold mismatch", D, pattern),
            )

            minimum_requirements = {
                "raw": min(keys),
                "scalar": min(
                    support_count - ambient_small_capacity
                    for support_count, _body, _L, _loads in records
                ),
                "ambient_status": min(
                    support_count - ambient_limit
                    for support_count, _body, _L, _loads in records
                ),
                "support_status": None,
            }
            support_requirements = []
            for support_count, _body, _L, small_loads in records:
                local_weights, local_marginals = small_vectors(
                    pattern, D, small_loads
                )
                require(local_marginals == marginals, "activity marginals changed")
                local_limit = status_load_limit(local_weights, local_marginals)
                support_requirements.append(support_count - local_limit)
            minimum_requirements["support_status"] = min(support_requirements)

            for stage in stage_names:
                stage_shape_count = count_at_least(
                    keys, suffix, minimum_requirements[stage]
                )
                shapes[stage] += stage_shape_count
                if stage_shape_count:
                    pattern_shapes[stage][pattern] += stage_shape_count

            for support_count, body, L, small_loads in records:
                local_weights, local_marginals = small_vectors(
                    pattern, D, small_loads
                )
                local_limit = status_load_limit(local_weights, local_marginals)
                thresholds = {
                    "raw": min(keys),
                    "scalar": support_count - ambient_small_capacity,
                    "ambient_status": support_count - ambient_limit,
                    "support_status": support_count - local_limit,
                }
                for stage in stage_names:
                    count = count_at_least(keys, suffix, thresholds[stage])
                    if not count:
                        continue
                    occurrences[stage] += count
                    pattern_occurrences[stage][pattern] += count
                    divisor_occurrences[stage][D] += count
                    row_key = (body, D)
                    rows[stage].add(row_key)
                    bodies[stage].add(body)
                    divisors[stage].add(D)
                    semantic[stage].update(
                        (
                            f"{body}|{L}|{D}|{support_count}|{pattern}|"
                            f"{count}\n"
                        ).encode()
                    )
        if progress:
            print(
                f"progress D={D},rows={len(records)},"
                f"features={len(distribution)},"
                f"raw_shapes={sum(distribution.values())}"
            )

    if max_D is None:
        require(shapes["raw"] == EXPECTED_SHAPES, "raw shape universe changed")
        require(
            occurrences["raw"] == EXPECTED_OCCURRENCES,
            "raw occurrence universe changed",
        )
    for left, right in zip(stage_names, stage_names[1:]):
        require(occurrences[right] <= occurrences[left], ("occurrence monotonicity", left, right))
        require(shapes[right] <= shapes[left], ("shape monotonicity", left, right))
        require(rows[right] <= rows[left], ("row monotonicity", left, right))
    if max_D == 100:
        require(
            sum(map(len, by_divisor.values())) == EXPECTED_D100["support_rows"],
            "D100 support-row ledger changed",
        )
        require(
            len(by_divisor) == EXPECTED_D100["support_divisors"],
            "D100 support-divisor ledger changed",
        )
        require(
            feature_count == EXPECTED_D100["features"],
            "D100 feature ledger changed",
        )
        for stage in stage_names:
            actual = (
                occurrences[stage],
                shapes[stage],
                len(rows[stage]),
                len(bodies[stage]),
                len(divisors[stage]),
                semantic[stage].hexdigest(),
            )
            require(
                actual == EXPECTED_D100["stages"][stage],
                ("D100 stage ledger changed", stage, actual),
            )

    threshold_histogram = Counter(threshold_table.values())
    positive_thresholds = {
        pattern: threshold
        for pattern, threshold in threshold_table.items()
        if threshold
    }

    print("LRC14 k=2 aligned / five-drift divisor-status GF scout")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"support_cutoff={SUPPORT_CUTOFF}")
    print(f"two_aligned_safe_floor={TWO_SAFE_FLOOR}")
    print(f"small_divisors={SMALL_DIVISORS}")
    print(f"activity_control_count={activity_control_count}")
    print(f"weighted_shape_control_cases={shape_control_cases}")
    print(f"pattern_threshold_count={len(threshold_table)}")
    print(f"pattern_threshold_histogram={threshold_histogram}")
    print(f"positive_pattern_thresholds={positive_thresholds}")
    print(f"max_D={max_D}")
    print(f"body_count={body_count}")
    print(f"body_divisor_rows={body_divisor_rows}")
    print(f"support_rows={sum(map(len, by_divisor.values()))}")
    print(f"support_divisors={len(by_divisor)}")
    print(f"aggregated_features={feature_count}")
    for stage in stage_names:
        print(
            f"stage={stage},occurrences={occurrences[stage]},"
            f"shapes={shapes[stage]},rows={len(rows[stage])},"
            f"bodies={len(bodies[stage])},divisors={len(divisors[stage])},"
            f"semantic_sha256={semantic[stage].hexdigest()}"
        )
        print(
            f"stage={stage},pattern_occurrences="
            f"{Counter(dict(sorted(pattern_occurrences[stage].items())))}"
        )
        print(
            f"stage={stage},pattern_shapes="
            f"{Counter(dict(sorted(pattern_shapes[stage].items())))}"
        )
        print(
            f"stage={stage},divisor_occurrences="
            f"{Counter(dict(sorted(divisor_occurrences[stage].items())))}"
        )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--max-D", type=int, default=100)
    parser.add_argument("--progress", action="store_true")
    arguments = parser.parse_args()
    main(max_D=arguments.max_D, progress=arguments.progress)
