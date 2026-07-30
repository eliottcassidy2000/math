#!/usr/bin/env python3
"""Exact aggregate scout for three aligned combs and four drifts.

This is a scratch research instrument, not a canonical theorem referee.
It never enumerates the 694,921,995 denominator four-multisets row by row.

Let ``D`` be the common resolving denominator and, for ``d|D``,

    C_D(d) = (D/d) ceil(d/7).

The truncated multiset generating function is

    G_D(y,x,z) =
      product_{d|D,d>1} (1-y*x^C_D(d)*z_d)^(-1)  mod y^5.

The implementation extracts the lcm-exact part by divisor-lattice Mobius
inversion.  Divisors with the same proof feature are grouped inside each
downward divisor alphabet before the truncated product is multiplied.  After
extracting ``y^4`` at exact lcm ``D``, it retains only

    (multiplicity of d=2,3,4; total capacity of d>4).

Thus every one of the raw denominator shapes is counted, but shapes with the
same proof-relevant data are aggregated.

There is a second exact relaxation.  For a fixed common aligned phase ``u``,
a denominator ``d<7`` mask is nonempty iff

    ||c*u|| < d/14,

and then occupies one lifted residue class.  Its activity set has Haar mass
``d/7``.  The safe set of three aligned masks has mass at least 55/91, so a
cover forces the small-mask status to lie in the upward threshold event

    sum(active small-mask loads) >= support - large-mask capacity

on a set of mass at least 55/91.

For an upward event A with exact one-marginals r_i, the general transport
theorem gives

    max Pr(A) = min(1, min_w sum_i r_i*w_i),

where w is a fractional cover of the inclusion-minimal true status sets.
The cover LP has at most four variables here and is solved exactly over
Fractions.  Equality is also impossible here: the aligned safe carrier is
compact, the heavy-status event is open, and a compact subset of a proper
open subset of the circle has strictly smaller measure.  Thus a status
survives this relaxation only when its maximum mass is *strictly greater*
than 55/91.

Two status ledgers are reported:

* ``ambient_status`` grants an active small mask its full ambient capacity;
* ``support_status`` grants only the largest exact S_D residue-class load.

Both grant every d>4 mask its full ambient capacity, so both are necessary
conditions and neither asserts realizability.  The second is stronger while
still depending only on a tiny per-row sidecar for d=2,3,4.
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


ROOT = Path(__file__).resolve().parents[1]
COMBINED_PATH = (
    ROOT / "04-computation" / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
EXPECTED_ROWS = 26970
EXPECTED_DIVISORS = 217
EXPECTED_SHAPES = 694921995
EXPECTED_OCCURRENCES = 21357714101

THREE_SAFE_FLOOR = Q(55, 91)
FOUR_SAFE_FLOOR = Q(558, 1183)
SUPPORT_CUTOFF = (Q(1) - FOUR_SAFE_FLOOR) / THREE_SAFE_FLOOR


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(COMBINED_PATH) == EXPECTED_COMBINED_SHA256,
    "frozen body-projection dependency changed",
)
require(SUPPORT_CUTOFF == Q(125, 143), "support cutoff changed")
spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


@lru_cache(maxsize=None)
def divisors_of(number):
    return tuple(support.divisors(number))


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
    """Exact maximum mass of ``sum_{i in status} weights[i] >= need``."""
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

    # Vertices of {w>=0, h.w>=1 for h minimal true}.  Nonnegativity rows
    # are included as equations w_i=0 among candidate active constraints.
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
def status_need_limit(weights, marginals):
    """Largest integral threshold whose relaxed mass is strictly above 55/91."""
    weights = tuple(weights)
    marginals = tuple(marginals)
    if not weights:
        return 0
    # The upward event changes only at a subset-sum load.  On the integer
    # interval between consecutive subset sums it is constant, so evaluating
    # the upper endpoint of each interval is exact and costs at most 2^4-1
    # cover LPs, independent of the magnitude of D.
    thresholds = sorted(
        {
            sum(
                weights[index]
                for index in range(len(weights))
                if (mask >> index) & 1
            )
            for mask in range(1, 1 << len(weights))
        }
    )
    low = 0
    for threshold in thresholds:
        if (
            maximum_upward_mass(weights, marginals, threshold)
            > THREE_SAFE_FLOOR
        ):
            low = threshold
        else:
            break
    require(
        maximum_upward_mass(weights, marginals, low) > THREE_SAFE_FLOOR,
        "status limit is not admissible",
    )
    if low < sum(weights):
        require(
            maximum_upward_mass(weights, marginals, low + 1)
            <= THREE_SAFE_FLOOR,
            "status limit is not maximal",
        )
    return low


def activity_controls():
    """Exact direction/equality controls for the status theorem use."""
    require(
        maximum_upward_mass((1,), (Q(4, 7),), 1) == Q(4, 7),
        "single d=4 activity mass changed",
    )
    require(
        maximum_upward_mass((1, 1), (Q(4, 7), Q(4, 7)), 1) == 1,
        "two-d4 union control changed",
    )
    require(
        maximum_upward_mass((1, 1), (Q(4, 7), Q(4, 7)), 2) == Q(4, 7),
        "two-d4 intersection control changed",
    )
    require(
        status_need_limit((1, 1), (Q(4, 7), Q(4, 7))) == 1,
        "strict status boundary changed",
    )
    require(
        maximum_upward_mass((2, 3), (Q(2, 7), Q(3, 7)), 2) == Q(5, 7),
        "weighted union control changed",
    )
    require(
        maximum_upward_mass((2, 3), (Q(2, 7), Q(3, 7)), 3) == Q(3, 7),
        "weighted singleton control changed",
    )
    # Exact ambient (d=2,d=3) breakpoint at residual D/3: equality is
    # admissible because either label can pay.  One integer step above it
    # only d=2 can pay.
    require(
        maximum_upward_mass((3, 2), (Q(2, 7), Q(3, 7)), 2) == Q(5, 7),
        "(2,3) equality endpoint was falsely killed",
    )
    require(
        maximum_upward_mass((3, 2), (Q(2, 7), Q(3, 7)), 3) == Q(2, 7),
        "(2,3) post-endpoint status changed",
    )
    # Closest ambient non-killing cover value: 13/21 exceeds 55/91 by
    # exactly 4/273.
    require(
        maximum_upward_mass(
            (3, 2, 2, 2),
            (Q(2, 7), Q(3, 7), Q(3, 7), Q(3, 7)),
            5,
        )
        == Q(13, 21),
        "closest non-killing activity control changed",
    )
    require(Q(13, 21) - THREE_SAFE_FLOOR == Q(4, 273), "upper margin changed")
    require(THREE_SAFE_FLOOR - Q(4, 7) == Q(3, 91), "lower margin changed")
    return 11


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
    for divisor in divisors_of(D):
        alphabet = len(divisors_of(divisor)) - 1
        total += mobius(D // divisor) * comb(alphabet + arity - 1, arity)
    return total


@lru_cache(maxsize=None)
def shape_distribution(D):
    """Map (small multiplicity pattern, large capacity) to shape count."""
    # Compute the lcm-exact coefficient by Mobius inversion.  Within each
    # divisor alphabet E|D, divisors with the same proof feature are grouped.
    # Choosing ``copies`` symbols from a group of size g with repetition has
    # coefficient C(g+copies-1,copies).
    result = Counter()
    for E in divisors_of(D):
        sign = mobius(D // E)
        if not sign:
            continue
        groups = Counter()
        for divisor in (d for d in divisors_of(E) if d > 1):
            capacity = (D // divisor) * ((divisor + 6) // 7)
            feature = (
                int(divisor == 2),
                int(divisor == 3),
                int(divisor == 4),
                capacity if divisor > 4 else 0,
            )
            groups[feature] += 1
        # State: (used, small2, small3, small4, large_capacity).
        states = {(0, 0, 0, 0, 0): 1}
        for feature, alphabet_size in groups.items():
            unit2, unit3, unit4, unit_large = feature
            additions = []
            for state, multiplicity in states.items():
                used, count2, count3, count4, large = state
                for copies in range(1, 5 - used):
                    coefficient = comb(
                        alphabet_size + copies - 1, copies
                    )
                    additions.append(
                        (
                            (
                                used + copies,
                                count2 + copies * unit2,
                                count3 + copies * unit3,
                                count4 + copies * unit4,
                                large + copies * unit_large,
                            ),
                            multiplicity * coefficient,
                        )
                    )
            for state, multiplicity in additions:
                states[state] = states.get(state, 0) + multiplicity
        for state, multiplicity in states.items():
            used, count2, count3, count4, large = state
            if used == 4:
                result[((count2, count3, count4), large)] += sign * multiplicity
    require(
        all(multiplicity >= 0 for multiplicity in result.values()),
        ("negative weighted Mobius coefficient", D),
    )
    result += Counter()  # discard exact zero coefficients after cancellation
    expected = lcm_multiset_shapes(D, 4)
    require(sum(result.values()) == expected, ("shape GF changed", D))
    return result


def brute_shape_controls():
    cases = 0
    for D in range(1, 61):
        alphabet = tuple(d for d in divisors_of(D) if d > 1)
        brute = Counter()
        for values in combinations_with_replacement(alphabet, 4):
            if lcm(*values) != D:
                continue
            pattern = (
                values.count(2),
                values.count(3),
                values.count(4),
            )
            large = sum(
                (D // d) * ((d + 6) // 7)
                for d in values
                if d > 4
            )
            brute[(pattern, large)] += 1
        require(brute == shape_distribution(D), ("weighted GF failed", D))
        cases += 1
    return cases


def small_vectors(pattern, D, row_loads=None):
    """Return labelled small-mask weights and exact activity marginals."""
    weights = []
    marginals = []
    for divisor, copies in zip((2, 3, 4), pattern):
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


def suffix_counter(counter):
    """Sorted keys and suffix sums for integer threshold queries."""
    keys = sorted(counter)
    suffix = [0] * (len(keys) + 1)
    for index in range(len(keys) - 1, -1, -1):
        suffix[index] = suffix[index + 1] + counter[keys[index]]
    return keys, suffix


def count_at_least(keys, suffix, threshold):
    return suffix[bisect_left(keys, threshold)]


def main(max_D=None, progress=False):
    activity_control_count = activity_controls()
    shape_control_cases = brute_shape_controls()

    by_divisor = defaultdict(list)
    body_count = 0
    body_divisor_rows = 0
    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(body)
        for D in divisors_of(L):
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
            for divisor in (2, 3, 4):
                if D % divisor == 0:
                    histogram = combined.residue_load_histogram(arcs, divisor)
                    small_loads[divisor] = combined.top_class_load(histogram, 1)
            by_divisor[D].append(
                (support_count, body, L, small_loads)
            )

    require(body_count == 3003, "body universe changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")
    if max_D is None:
        require(sum(map(len, by_divisor.values())) == EXPECTED_ROWS, "k3 row count changed")
        require(len(by_divisor) == EXPECTED_DIVISORS, "k3 divisor alphabet changed")

    stage_names = ("raw", "scalar", "ambient_status", "support_status")
    occurrences = Counter()
    shapes = Counter()
    rows = {stage: set() for stage in stage_names}
    bodies = {stage: set() for stage in stage_names}
    divisors = {stage: set() for stage in stage_names}
    pattern_occurrences = {stage: Counter() for stage in stage_names}
    pattern_shapes = {stage: Counter() for stage in stage_names}
    divisor_occurrences = {stage: Counter() for stage in stage_names}
    status_limit_cache_summary = Counter()
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
            total_pattern_shapes = suffix[0]
            ambient_weights, marginals = small_vectors(pattern, D)
            ambient_small_capacity = sum(ambient_weights)
            ambient_limit = status_need_limit(ambient_weights, marginals)
            status_limit_cache_summary[
                (pattern, ambient_small_capacity, ambient_limit)
            ] += 1

            # A shape is counted at a stage if it survives at least one row.
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
            for support_count, body, L, small_loads in records:
                local_weights, local_marginals = small_vectors(
                    pattern, D, small_loads
                )
                require(local_marginals == marginals, "activity marginals changed")
                local_limit = status_need_limit(local_weights, marginals)
                support_requirements.append(support_count - local_limit)
            minimum_requirements["support_status"] = min(support_requirements)

            for stage in stage_names:
                stage_shape_count = count_at_least(
                    keys, suffix, minimum_requirements[stage]
                )
                shapes[stage] += stage_shape_count
                if stage_shape_count:
                    pattern_shapes[stage][pattern] += stage_shape_count

            for row_index, (support_count, body, L, small_loads) in enumerate(records):
                local_weights, local_marginals = small_vectors(
                    pattern, D, small_loads
                )
                local_limit = status_need_limit(local_weights, local_marginals)
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

    print("LRC14 k=3 aligned / four-drift divisor-status GF scout")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"support_cutoff={SUPPORT_CUTOFF}")
    print(f"three_aligned_safe_floor={THREE_SAFE_FLOOR}")
    print(f"activity_control_count={activity_control_count}")
    print(f"weighted_shape_control_cases={shape_control_cases}")
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
            f"stage={stage},divisor_occurrences_top="
            f"{divisor_occurrences[stage].most_common(30)}"
        )
    print(f"ambient_status_limit_types={len(status_limit_cache_summary)}")
    print(
        "ambient_status_limit_top="
        f"{status_limit_cache_summary.most_common(30)}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--max-D", type=int, default=None)
    parser.add_argument("--progress", action="store_true")
    arguments = parser.parse_args()
    main(max_D=arguments.max_D, progress=arguments.progress)
