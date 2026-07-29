#!/usr/bin/env python3
"""Exact mixed three-drift activity/fibre sieve for THM-2928.

This candidate referee reconstructs the frozen 544,571
Lambda_1+Lambda_2+Lambda_3 occurrences for four aligned combs and three
drifts, then applies two independent necessary mechanisms.

1. Divisor-fibre overload.  For d,q|D, ell=ceil(d/7),
   g=gcd(d,q), and H=D/lcm(d,q), an enlarged denominator-d needle meets
   each q-fibre in at most

       H * ceil(ell/g)

   points.  Every q|D is tested.  This kills 125,060 occurrences, including
   the entire diagonal-denominator sector.

2. Common-u activity.  For d<7, a fixed-u denominator-d mask is nonempty
   exactly when

       ||c u|| < d/14,

   because multiplication by the unit c permutes the d residue classes.
   This activity set has Haar measure d/7.  The four-aligned carrier has
   measure at least 558/1183, which is strictly greater than d/7 for
   d=2,3.  Therefore, if the other two masks' total ambient capacities are
   smaller than |S_D|, the d-mask would have to be active for every carrier
   phase, an impossibility.

   A residual double-denominator-2 occurrence has a sharper parity version.
   Each active d=2 mask occupies exactly one parity class.  If the other
   activity screen has not already killed the row, the q=2 fibre load of
   S_D still exceeds the exact per-parity capacity of the third mask in
   every case.  The two parity masks would consequently have to be active
   on opposite parities for every carrier phase, forcing the carrier into
   a single d=2 activity set of measure 2/7.  This is impossible.

The script is numerator-free only in the first mechanism.  The second
mechanism deliberately restores the shared carrier phase u that independent
address-mask relaxations forget.  Survivors are not counterexamples: they
only survive these necessary screens.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
COMBINED_PATH = (
    HERE / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
COMBINED_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_three_drift_body_projection_fiber_thm2928.out"
)
EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
EXPECTED_COMBINED_OUTPUT_SHA256 = (
    "2e211620ad7064ea06f7544b5fbac709d6d52d9a0e261b464ae26b595f09b669"
)
EXPECTED_ALL_TOP_ATLAS_SHA256 = (
    "88ac1b6c34c11eadd286d02f930f69a44a82bd429fe90248417fd35caa947960"
)
EXPECTED_Q_SURVIVOR_SHA256 = (
    "87fc8e09d7e69c9d91cf7ebc895703393496e9b1d4b5a1623ca6c112c6cbb4fc"
)
EXPECTED_FINAL_RESIDUAL_SHA256 = (
    "6a0b57e9da265d9f5aab55b8dd6166e7439d759a0f4ba5f1ef6e43b4d212cab5"
)

FOUR_SAFE_FLOOR = Q(558, 1183)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(COMBINED_PATH) == EXPECTED_COMBINED_SHA256,
    "frozen combined referee changed",
)
require(
    file_sha256(COMBINED_OUTPUT_PATH) == EXPECTED_COMBINED_OUTPUT_SHA256,
    "frozen combined output changed",
)
require(FOUR_SAFE_FLOOR > Q(3, 7), "activity measure gap changed")

spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


def torus_distance(value):
    residue = value % 1
    return min(residue, 1 - residue)


def small_denominator_mask(d, c, u):
    """Literal residue classes hit at a fixed rational phase u."""
    require(gcd(c, d) == 1, "small-denominator numerator is not a unit")
    return tuple(
        r
        for r in range(d)
        if torus_distance(Q(c) * (r + u) / d) < Q(1, 14)
    )


def bitset_arcs(mask, D):
    """Singleton arcs for small exact controls."""
    return tuple((x, x + 1) for x in range(D) if (mask >> x) & 1)


def exact_needle(D, d, phase, step):
    """An enlarged lifted unit needle, used only by finite controls."""
    require(D % d == 0, "needle denominator does not divide ambient modulus")
    require(gcd(step, d) == 1, "needle step is not a unit")
    width = (d + 6) // 7
    classes = {
        (phase + index * step) % d for index in range(width)
    }
    return sum(1 << x for x in range(D) if x % d in classes)


def fibre_cap(D, d, q):
    """Sharp phase-free maximum of one enlarged needle on one q-fibre."""
    ell = (d + 6) // 7
    g = gcd(d, q)
    H = D // lcm(d, q)
    return H * ((ell + g - 1) // g)


def controls():
    """Positive and hostile controls for both proof mechanisms."""
    # Exhaustive rational samples of the small-denominator activity law.
    activity_samples = 0
    active_samples = 0
    for d in (2, 3):
        for c in range(1, 13):
            if gcd(c, d) != 1:
                continue
            for numerator in range(97):
                u = Q(numerator, 97)
                mask = small_denominator_mask(d, c, u)
                predicted = torus_distance(Q(c) * u) < Q(d, 14)
                require(bool(mask) == predicted, "activity equivalence failed")
                require(len(mask) <= 1, "small mask hit two residue classes")
                activity_samples += 1
                active_samples += bool(mask)

    # A genuine union of three enlarged needles must pass every fibre cap.
    positive_D = 420
    positive_ds = (20, 28, 30)
    require(lcm(*positive_ds) == positive_D, "positive lcm changed")
    positive_masks = (
        exact_needle(positive_D, positive_ds[0], 3, 3),
        exact_needle(positive_D, positive_ds[1], 5, 5),
        exact_needle(positive_D, positive_ds[2], 7, 7),
    )
    positive_union = positive_masks[0] | positive_masks[1] | positive_masks[2]
    positive_arcs = bitset_arcs(positive_union, positive_D)
    for q in support.divisors(positive_D):
        histogram = combined.residue_load_histogram(positive_arcs, q)
        maximum = max(load for load, count in histogram if count)
        cap = sum(fibre_cap(positive_D, d, q) for d in positive_ds)
        require(maximum <= cap, ("positive fibre control rejected", q))

    # Four points in one q=4 fibre defeat three diagonal D=28 sections.
    hostile_D = 28
    hostile_q = 4
    hostile = sum(1 << x for x in (0, 4, 8, 12))
    hostile_histogram = combined.residue_load_histogram(
        bitset_arcs(hostile, hostile_D), hostile_q
    )
    hostile_maximum = max(
        load for load, count in hostile_histogram if count
    )
    hostile_cap = 3 * fibre_cap(hostile_D, hostile_D, hostile_q)
    require(
        (hostile_maximum, hostile_cap) == (4, 3),
        "hostile fibre control lost overload",
    )
    return (
        activity_samples,
        active_samples,
        positive_union.bit_count(),
        hostile_maximum,
        hostile_cap,
    )


def repeat_type(d1, d2, d3):
    if d1 == d3:
        return "diagonal"
    if d1 == d2:
        return "low_pair"
    if d2 == d3:
        return "high_pair"
    return "distinct"


def quantiles(values):
    values = sorted(values)
    return {
        q: values[min(len(values) - 1, len(values) * q // 100)]
        for q in (0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100)
    }


def main():
    (
        activity_samples,
        active_samples,
        positive_union_size,
        hostile_maximum,
        hostile_cap,
    ) = controls()

    by_divisor = {}
    body_count = 0
    divisor_rows = 0
    support_hard_rows = 0
    for F in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(F)
        for D in support.divisors(L):
            divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > combined.SUPPORT_CUTOFF:
                continue
            support_hard_rows += 1
            by_divisor.setdefault(D, []).append(
                (support_count, F, L, tuple(ranges))
            )

    require(body_count == 3003, "body universe changed")
    require(divisor_rows == 251536, "body/divisor universe changed")
    require(support_hard_rows == 13778, "support-hard universe changed")
    require(len(by_divisor) == 206, "support-hard divisor alphabet changed")

    all_top_occurrences = 0
    all_top_diagonal = 0
    q_fibre_kills = 0
    q_fibre_survivors = 0
    q_fibre_diagonal_survivors = 0
    q_kill_local_modulus = Counter()
    q_best_margin_histogram = Counter()
    q_repeat_survivors = Counter()
    general_activity_kills = 0
    general_activity_reasons = Counter()
    double_two_q_survivors = 0
    double_two_general_overlap = 0
    double_two_after_general = 0
    double_two_parity_kills = 0
    double_two_parity_margin = Counter()
    final_records = []
    final_rows = set()
    final_shapes = set()
    final_bodies = set()
    final_Ds = set()
    final_repeat = Counter()
    final_minimum_denominator = Counter()
    all_top_hash = sha256()
    q_survivor_hash = sha256()
    final_hash = sha256()

    for D in sorted(by_divisor):
        divisors = tuple(d for d in support.divisors(D) if d > 1)
        capacities = tuple(
            (D // d) * ((d + 6) // 7) for d in divisors
        )
        rows = sorted(
            by_divisor[D],
            key=lambda record: (record[0], record[1], record[2]),
        )
        supports = tuple(record[0] for record in rows)
        top_loads_by_row = []
        max_fibre_by_row = []
        parity_load_by_row = []
        for support_count, F, _L, ranges in rows:
            arcs = combined.projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support arc mismatch", F, D),
            )
            tops = []
            maxima = []
            parity_load = None
            for q in divisors:
                histogram = combined.residue_load_histogram(arcs, q)
                tops.append(
                    combined.top_class_load(histogram, (q + 6) // 7)
                )
                maxima.append(
                    max(load for load, count in histogram if count)
                )
                if q == 2:
                    require(
                        histogram == ((support_count // 2, 2),),
                        ("reflection parity balance changed", F, D, histogram),
                    )
                    parity_load = support_count // 2
            top_loads_by_row.append(tuple(tops))
            max_fibre_by_row.append(tuple(maxima))
            parity_load_by_row.append(parity_load)

        for first_index, d1 in enumerate(divisors):
            for second_index in range(first_index, len(divisors)):
                d2 = divisors[second_index]
                first_second_lcm = lcm(d1, d2)
                for third_index in range(second_index, len(divisors)):
                    d3 = divisors[third_index]
                    if lcm(first_second_lcm, d3) != D:
                        continue
                    full_capacity = (
                        capacities[first_index]
                        + capacities[second_index]
                        + capacities[third_index]
                    )
                    capacity_row_count = bisect_right(
                        supports, full_capacity
                    )
                    all_top_rows = []
                    for row_index in range(capacity_row_count):
                        top_loads = top_loads_by_row[row_index]
                        top_sum = (
                            top_loads[first_index]
                            + top_loads[second_index]
                            + top_loads[third_index]
                        )
                        if supports[row_index] <= top_sum:
                            all_top_rows.append(row_index)
                    if not all_top_rows:
                        continue

                    fibre_caps = tuple(
                        sum(fibre_cap(D, d, q) for d in (d1, d2, d3))
                        for q in divisors
                    )
                    for row_index in all_top_rows:
                        support_count, F, L, _ranges = rows[row_index]
                        all_top_occurrences += 1
                        all_top_diagonal += d1 == d3
                        all_top_hash.update(
                            (
                                f"{F}|{L}|{D}|{support_count}|"
                                f"{d1}|{d2}|{d3}\n"
                            ).encode()
                        )

                        margins = tuple(
                            maximum - cap
                            for maximum, cap in zip(
                                max_fibre_by_row[row_index],
                                fibre_caps,
                            )
                        )
                        best_margin = max(margins)
                        q_best_margin_histogram[best_margin] += 1
                        if best_margin > 0:
                            q_fibre_kills += 1
                            best_q = min(
                                q
                                for q, margin in zip(divisors, margins)
                                if margin == best_margin
                            )
                            q_kill_local_modulus[D // best_q] += 1
                            continue

                        q_fibre_survivors += 1
                        q_fibre_diagonal_survivors += d1 == d3
                        q_repeat_survivors[
                            repeat_type(d1, d2, d3)
                        ] += 1
                        q_record = (
                            F,
                            L,
                            D,
                            support_count,
                            d1,
                            d2,
                            d3,
                            best_margin,
                        )
                        q_survivor_hash.update(f"{q_record}\n".encode())

                        ds = (d1, d2, d3)
                        caps = (
                            capacities[first_index],
                            capacities[second_index],
                            capacities[third_index],
                        )
                        activity_reasons = tuple(
                            (index, d)
                            for index, d in enumerate(ds)
                            if d in (2, 3)
                            and support_count
                            > sum(
                                caps[other]
                                for other in range(3)
                                if other != index
                            )
                        )
                        is_double_two = d1 == d2 == 2
                        double_two_q_survivors += is_double_two
                        if activity_reasons:
                            general_activity_kills += 1
                            general_activity_reasons[activity_reasons] += 1
                            double_two_general_overlap += is_double_two
                            continue

                        if is_double_two:
                            double_two_after_general += 1
                            require(
                                parity_load_by_row[row_index] is not None,
                                ("double-two row lacks parity quotient", F, D),
                            )
                            parity_load = parity_load_by_row[row_index]
                            third_parity_cap = fibre_cap(D, d3, 2)
                            require(
                                parity_load > third_parity_cap,
                                (
                                    "double-two parity survivor",
                                    F,
                                    D,
                                    d3,
                                    parity_load,
                                    third_parity_cap,
                                ),
                            )
                            double_two_parity_kills += 1
                            double_two_parity_margin[
                                parity_load - third_parity_cap
                            ] += 1
                            continue

                        final_records.append(q_record)
                        final_rows.add((F, D))
                        final_shapes.add((D, d1, d2, d3))
                        final_bodies.add(F)
                        final_Ds.add(D)
                        final_repeat[repeat_type(d1, d2, d3)] += 1
                        final_minimum_denominator[d1] += 1
                        final_hash.update(f"{q_record}\n".encode())

    require(all_top_occurrences == 544571, "all-top occurrence ledger changed")
    require(all_top_diagonal == 2636, "diagonal baseline changed")
    require(
        all_top_hash.hexdigest() == EXPECTED_ALL_TOP_ATLAS_SHA256,
        "all-top semantic atlas changed",
    )
    require(q_fibre_kills == 125060, "q-fibre kill count changed")
    require(q_fibre_survivors == 419511, "q-fibre survivor count changed")
    require(
        q_fibre_diagonal_survivors == 0,
        "q-fibre screen left a diagonal occurrence",
    )
    require(
        q_survivor_hash.hexdigest() == EXPECTED_Q_SURVIVOR_SHA256,
        "q-fibre survivor atlas changed",
    )
    require(
        general_activity_kills == 383391,
        "general activity kill count changed",
    )
    require(
        double_two_q_survivors == 15379,
        "double-two q-survivor count changed",
    )
    require(
        double_two_general_overlap == 8623,
        "double-two activity overlap changed",
    )
    require(
        double_two_after_general == 6756,
        "double-two post-activity count changed",
    )
    require(
        double_two_parity_kills == 6756,
        "double-two parity completion changed",
    )
    require(len(final_records) == 29364, "final occurrence count changed")
    require(len(final_rows) == 2974, "final row count changed")
    require(len(final_shapes) == 4993, "final shape count changed")
    require(len(final_bodies) == 2878, "final body count changed")
    require(len(final_Ds) == 110, "final divisor count changed")
    require(
        final_repeat
        == Counter({"distinct": 20946, "high_pair": 5984, "low_pair": 2434}),
        "final repetition atlas changed",
    )
    require(
        final_hash.hexdigest() == EXPECTED_FINAL_RESIDUAL_SHA256,
        "final residual semantic atlas changed",
    )

    smallest = sorted(
        final_records,
        key=lambda record: (
            record[2],
            record[4],
            record[5],
            record[6],
            record[0],
        ),
    )[:12]
    D_profile = quantiles(record[2] for record in final_records)

    print("LRC14 mixed three-drift activity/fibre exact referee")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"combined_output_sha256={file_sha256(COMBINED_OUTPUT_PATH)}")
    print(f"four_comb_safe_floor={FOUR_SAFE_FLOOR}")
    print(f"activity_measure_d2={Q(2, 7)}")
    print(f"activity_measure_d3={Q(3, 7)}")
    print(f"activity_control_samples={activity_samples}")
    print(f"activity_control_active_samples={active_samples}")
    print(f"positive_fibre_union_size={positive_union_size}")
    print(
        "hostile_fibre_control="
        f"maximum={hostile_maximum},capacity={hostile_cap}"
    )
    print(f"all_top_occurrences={all_top_occurrences}")
    print(f"all_top_diagonal={all_top_diagonal}")
    print(f"all_top_semantic_sha256={all_top_hash.hexdigest()}")
    print(f"q_fibre_kills={q_fibre_kills}")
    print(f"q_fibre_survivors={q_fibre_survivors}")
    print(
        f"q_fibre_diagonal_survivors={q_fibre_diagonal_survivors}"
    )
    print(f"q_fibre_best_margin_histogram={q_best_margin_histogram}")
    print(
        "q_fibre_kill_local_modulus_top="
        f"{q_kill_local_modulus.most_common(60)}"
    )
    print(f"q_fibre_repeat_survivors={q_repeat_survivors}")
    print(f"q_survivor_semantic_sha256={q_survivor_hash.hexdigest()}")
    print(f"general_activity_kills={general_activity_kills}")
    print(f"general_activity_reasons={general_activity_reasons}")
    print(
        "post_general_activity_survivors="
        f"{q_fibre_survivors-general_activity_kills}"
    )
    print(f"double_two_q_survivors={double_two_q_survivors}")
    print(
        f"double_two_general_overlap={double_two_general_overlap}"
    )
    print(f"double_two_after_general={double_two_after_general}")
    print(f"double_two_parity_kills={double_two_parity_kills}")
    print(
        "double_two_parity_margin="
        f"distinct={len(double_two_parity_margin)},"
        f"minimum={min(double_two_parity_margin)},"
        f"maximum={max(double_two_parity_margin)},"
        f"top={double_two_parity_margin.most_common(60)}"
    )
    print(f"final_occurrences={len(final_records)}")
    print(f"final_rows={len(final_rows)}")
    print(f"final_shapes={len(final_shapes)}")
    print(f"final_bodies={len(final_bodies)}")
    print(f"final_divisors={len(final_Ds)}")
    print(f"final_repeat_types={final_repeat}")
    print(
        "final_minimum_denominator_histogram="
        f"{Counter(dict(sorted(final_minimum_denominator.items())))}"
    )
    print(f"final_D_quantiles={D_profile}")
    print(f"final_smallest_examples={smallest}")
    print(f"final_semantic_sha256={final_hash.hexdigest()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
