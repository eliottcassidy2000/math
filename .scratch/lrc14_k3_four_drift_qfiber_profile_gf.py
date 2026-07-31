#!/usr/bin/env python3
"""Bounded exact q-fibre profile scout for the k=3/four-drift sector.

This complements ``lrc14_k3_four_drift_divisor_status_gf.py``.  It restricts
to the structurally important residual in which every drift denominator is
at least five, then aggregates denominator four-multisets by their exact
q-fibre load profile rather than enumerating them row by row.

For d,q|D, write

    ell=ceil(d/7), g=gcd(d,q), ell=A*g+r,
    H=D/lcm(d,q).

An enlarged denominator-d needle has q-fibre load

    H * (A + high_i(b)),

where ``high_i`` has exactly ``R_i=q*r/g`` true fibres.  For four masks the
baseline loads add, and every target threshold defines an upward event in
the high statuses.  The exact real relaxation from the general status
theorem is

    max heavy fibres =
      min(q, min_w sum_i R_i*w_i),

with w a fractional cover of the minimal true status sets.  An integer
target-heavy count larger than this rational upper bound is impossible.

The multiset generating function is propagated with lcm state and truncated
at arity four.  Its coefficient key is only

    (ambient capacity, q-baseline, multiset of (increment, high count)).

Thus the program can audit a chosen resolving denominator and several q
quotients without materializing its denominator tuples.  Screens for
different q are not intersected: their rowwise minimum is reported honestly
as an upper bound for the all-q intersection.
"""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / ".scratch" / "lrc14_k3_four_drift_divisor_status_gf.py"
spec = spec_from_file_location("lrc14_k3_status_gf", BASE_PATH)
base = module_from_spec(spec)
spec.loader.exec_module(base)
support = base.support
combined = base.combined


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
    require(0 <= high_count <= q, "high marginal count changed")
    return low, increment, high_count


def no_small_shape_count(D):
    """Mobius count of lcm-D four-multisets with every denominator >4."""
    total = 0
    for E in base.divisors_of(D):
        alphabet = sum(d > 4 for d in base.divisors_of(E))
        total += base.mobius(D // E) * comb(alphabet + 3, 4)
    require(total >= 0, ("negative no-small shape count", D))
    return total


@lru_cache(maxsize=None)
def profile_distribution(D, q):
    """Aggregate all lcm-D four-multisets with denominators d>4."""
    require(D % q == 0, "q does not divide D")
    # State: used, lcm, ambient, baseline, sorted (increment,R) profile.
    states = {(0, 1, 0, 0, ()): 1}
    for divisor in (d for d in base.divisors_of(D) if d > 4):
        ambient = (D // divisor) * ((divisor + 6) // 7)
        low, increment, high_count = q_type(D, divisor, q)
        additions = []
        for state, multiplicity in states.items():
            used, current_lcm, old_ambient, baseline, profile = state
            for copies in range(1, 5 - used):
                new_profile = list(profile)
                if increment:
                    new_profile.extend([(increment, high_count)] * copies)
                    new_profile.sort()
                additions.append(
                    (
                        (
                            used + copies,
                            lcm(current_lcm, divisor),
                            old_ambient + copies * ambient,
                            baseline + copies * low,
                            tuple(new_profile),
                        ),
                        multiplicity,
                    )
                )
        for state, multiplicity in additions:
            states[state] = states.get(state, 0) + multiplicity
    result = Counter()
    for state, multiplicity in states.items():
        used, current_lcm, ambient, baseline, profile = state
        if used == 4 and current_lcm == D:
            result[(ambient, baseline, profile)] += multiplicity
    return result


def target_thresholds(histogram):
    """Pairs (load, number of fibres with target load at least load)."""
    result = []
    running = 0
    for load, count in reversed(histogram):
        running += count
        if load:
            result.append((load, running))
    return tuple(result)


def exact_needle(D, d, phase, step):
    require(D % d == 0, "control denominator does not divide D")
    require(gcd(step, d) == 1, "control step is not a unit")
    width = (d + 6) // 7
    classes = {
        (phase + index * step) % d for index in range(width)
    }
    return sum(1 << x for x in range(D) if x % d in classes)


def bitset_arcs(mask, D):
    arcs = []
    left = None
    for point in range(D + 1):
        occupied = point < D and bool((mask >> point) & 1)
        if occupied and left is None:
            left = point
        elif not occupied and left is not None:
            arcs.append((left, point))
            left = None
    return tuple(arcs)


def positive_profile_controls():
    """Every literal enlarged-needle union must pass its own relaxation."""
    D = 84
    denominators = (6, 7, 12, 14)
    phases = (1, 2, 3, 4)
    steps = (1, 2, 5, 3)
    require(lcm(*denominators) == D, "control lcm changed")
    masks = tuple(
        exact_needle(D, d, phase, step)
        for d, phase, step in zip(denominators, phases, steps)
    )
    union = 0
    for mask in masks:
        union |= mask
    arcs = bitset_arcs(union, D)
    checks = 0
    for q in base.divisors_of(D):
        baseline = 0
        profile = []
        for d in denominators:
            low, increment, high_count = q_type(D, d, q)
            baseline += low
            if increment:
                profile.append((increment, high_count))
        profile.sort()
        thresholds = target_thresholds(
            combined.residue_load_histogram(arcs, q)
        )
        require(
            profile_passes(q, baseline, tuple(profile), thresholds),
            ("positive q-profile control rejected", q),
        )
        checks += 1
    return checks, union.bit_count()


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
        maximum_count = q * maximum_mass
        if Q(target_count) > maximum_count:
            return False
    return True


def main(D, local_moduli):
    base.activity_controls()
    positive_checks, positive_union_size = positive_profile_controls()
    rows = []
    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        if L % D:
            continue
        support_count = support.support_size_bitset(D, ranges)
        if Q(support_count, D) > base.SUPPORT_CUTOFF:
            continue
        arcs = combined.projected_support_arcs(D, ranges)
        require(
            sum(right - left for left, right in arcs) == support_count,
            ("support projection changed", body, D),
        )
        rows.append((support_count, body, L, arcs))
    require(rows, ("chosen D has no k3 support rows", D))
    rows.sort(key=lambda row: (row[0], row[1], row[2]))

    q_values = sorted(
        {
            D // modulus
            for modulus in local_moduli
            if modulus > 0 and D % modulus == 0
        }
        | {1}
    )
    scalar_shape_count = None
    scalar_occurrences = 0
    row_q_counts = {index: {} for index in range(len(rows))}
    q_summaries = []
    raw_no_small_shapes = None

    for q in q_values:
        distribution = profile_distribution(D, q)
        shape_count = sum(distribution.values())
        require(
            shape_count == no_small_shape_count(D),
            ("q-profile shape coefficient changed", D, q),
        )
        if raw_no_small_shapes is None:
            raw_no_small_shapes = shape_count
        require(
            shape_count == raw_no_small_shapes,
            ("q changed no-small shape universe", D, q),
        )
        stage_features = set()
        occurrences = 0
        passing_rows = 0
        feature_count = len(distribution)
        for row_index, (support_count, _body, _L, arcs) in enumerate(rows):
            histogram = combined.residue_load_histogram(arcs, q)
            thresholds = target_thresholds(histogram)
            row_count = 0
            for feature, multiplicity in distribution.items():
                ambient, baseline, profile = feature
                if ambient < support_count:
                    continue
                if not profile_passes(q, baseline, profile, thresholds):
                    continue
                row_count += multiplicity
                stage_features.add(feature)
            row_q_counts[row_index][q] = row_count
            occurrences += row_count
            passing_rows += bool(row_count)
        if q == 1:
            scalar_shape_count = sum(
                distribution[feature] for feature in stage_features
            )
            scalar_occurrences = occurrences
        q_summaries.append(
            (
                q,
                D // q,
                feature_count,
                occurrences,
                sum(distribution[feature] for feature in stage_features),
                passing_rows,
            )
        )

    best_row_counts = {
        row_index: min(counts.values())
        for row_index, counts in row_q_counts.items()
    }
    best_single_q_upper = sum(best_row_counts.values())
    best_rows = sum(count > 0 for count in best_row_counts.values())
    best_q_histogram = Counter()
    for row_index, counts in row_q_counts.items():
        best = min(counts.values())
        best_q = min(q for q, count in counts.items() if count == best)
        best_q_histogram[(best_q, D // best_q)] += 1

    print("LRC14 k=3/four-drift bounded q-fibre profile GF scout")
    print(f"D={D}")
    print(f"rows={len(rows)}")
    print(f"positive_profile_checks={positive_checks}")
    print(f"positive_union_size={positive_union_size}")
    print(f"raw_no_small_shapes={raw_no_small_shapes}")
    print(f"q_values={q_values}")
    print(f"scalar_shapes={scalar_shape_count}")
    print(f"scalar_occurrences={scalar_occurrences}")
    print(
        "q_summaries=(q,local_modulus,features,occurrences,"
        f"shapes,rows)={q_summaries}"
    )
    print(f"best_single_q_occurrence_upper={best_single_q_upper}")
    print(f"best_single_q_rows={best_rows}")
    print(f"best_q_histogram={best_q_histogram}")
    scalar_survivor_rows = [
        (
            rows[row_index][1],
            rows[row_index][2],
            rows[row_index][0],
            counts,
        )
        for row_index, counts in row_q_counts.items()
        if counts[1]
    ]
    print(f"scalar_survivor_rows={scalar_survivor_rows}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--D", type=int, required=True)
    parser.add_argument(
        "--local-moduli",
        type=str,
        default="5,6,7,10,12,13,14,20,21,28,35,40,42,55,60,65,84",
    )
    arguments = parser.parse_args()
    main(
        arguments.D,
        tuple(
            int(value)
            for value in arguments.local_moduli.split(",")
            if value
        ),
    )
