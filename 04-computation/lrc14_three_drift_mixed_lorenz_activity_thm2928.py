#!/usr/bin/env python3
"""Exact divisor-fiber Lorenz/activity sieve for THM-2928's three drifts.

This referee imports the frozen combined four-aligned/three-drift referee
and reconstructs exactly its final
Lambda_1+Lambda_2+Lambda_3 occurrence ledger.

For d,q|D and ell=ceil(d/7), the enlarged denominator-d needle has, on
the q quotient fibers, the exact phase-independent load multiset

    H(A+1), repeated r*q/g times,
    H*A,   repeated q-r*q/g times,

where g=gcd(d,q), H=D/lcm(d,q), and ell=A*g+r.  Consequently its largest
s q-fiber loads sum to

    H * (A*s + min(s, r*q/g)).

If S_D is covered by three such needles, then for every q|D and every s
the sum of the s largest q-fiber loads of S_D is at most the sum of the
three displayed needle Lorenz bounds.  It suffices to test breakpoints of
the four piecewise-linear functions.

The Lorenz screen is deliberately phase-free: each needle may choose its best
q-fibers independently.  Hence every rejection is a rigorous obstruction,
while every survivor merely defeats this relaxation.

The final cross-tab also applies the proved common-u activity obstruction:
if d_i is 2 or 3 and the support exceeds the other two masks' total
capacities, then the d_i mask would have to be active on the full aligned
carrier, impossible because d_i/7 < 558/1183.  The separate double-(2,2)
capacity inequality is also retained.
"""

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
COMBINED_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_three_drift_body_projection_fiber_thm2928.py"
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
EXPECTED_LORENZ_SEMANTIC_SHA256 = (
    "cca2cdc865edced5335de67486003f20dd38e299efc5d00eb9e52b1d134c6431"
)
EXPECTED_COMBINED_SEMANTIC_SHA256 = (
    "cc8dea01ecf116f2514ea1fbbe917be736c429c5e41ad932b17158ffc27e7a52"
)


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
spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


def support_lorenz_profile(arcs, q):
    """Return descending-load cumulative breakpoints for q fibers."""
    histogram = combined.residue_load_histogram(arcs, q)
    cumulative_count = 0
    cumulative_mass = 0
    counts = []
    masses = []
    slopes = []
    for load, count in reversed(histogram):
        require(count > 0, "zero-size support histogram class")
        cumulative_count += count
        cumulative_mass += load * count
        counts.append(cumulative_count)
        masses.append(cumulative_mass)
        slopes.append(load)
    require(counts[-1] == q, "support Lorenz profile lost fibers")
    return (tuple(counts), tuple(masses), tuple(slopes))


def support_lorenz_value(profile, s):
    """Sum of the s largest support-fiber loads."""
    counts, masses, slopes = profile
    require(0 <= s <= counts[-1], "Lorenz argument outside quotient")
    if s == 0:
        return 0
    index = bisect_left(counts, s)
    previous_count = counts[index - 1] if index else 0
    previous_mass = masses[index - 1] if index else 0
    return previous_mass + (s - previous_count) * slopes[index]


def needle_lorenz_parameters(D, d, q):
    """Return base slope, exceptional height, exceptional-fiber count."""
    ell = (d + 6) // 7
    g = gcd(d, q)
    H = D // lcm(d, q)
    A, r = divmod(ell, g)
    high_count = r * (q // g)
    require(
        H * (A * q + high_count) == (D // d) * ell,
        ("needle profile lost mass", D, d, q),
    )
    return (H * A, H, high_count)


def phase_free_best_margin(
    q,
    support_profile,
    first_profile,
    second_profile,
    third_profile,
):
    """Return max(L_support-L_needles), with an attaining breakpoint."""
    profiles = (first_profile, second_profile, third_profile)
    base = sum(profile[0] for profile in profiles)
    breakpoints = set(support_profile[0])
    breakpoints.update(
        profile[2] for profile in profiles if 0 < profile[2] < q
    )
    breakpoints.add(1)
    breakpoints.add(q)
    best = None
    for s in breakpoints:
        left = support_lorenz_value(support_profile, s)
        right = base * s + sum(
            height * min(s, high_count)
            for _base, height, high_count in profiles
        )
        margin = left - right
        if best is None or margin > best[0]:
            best = (margin, s, left, right)
    require(best is not None, "empty Lorenz breakpoint set")
    return best


def phase_free_violation(
    q,
    support_profile,
    first_profile,
    second_profile,
    third_profile,
):
    """Return the strongest positive Lorenz margin, or None."""
    best = phase_free_best_margin(
        q,
        support_profile,
        first_profile,
        second_profile,
        third_profile,
    )
    return best if best[0] > 0 else None


def exact_needle(D, d, phase, step):
    """The enlarged lifted needle as a D-bit word (small controls only)."""
    require(gcd(step, d) == 1, "control step is not a unit")
    ell = (d + 6) // 7
    classes = {(phase + index * step) % d for index in range(ell)}
    return sum(1 << x for x in range(D) if x % d in classes)


def bitset_arcs(mask, D):
    """Convert a small non-wrapping bitset to singleton half-open arcs."""
    return tuple((point, point + 1) for point in range(D) if mask >> point & 1)


def capacity(D, d):
    return (D // d) * ((d + 6) // 7)


def activity_status(D, support_count, denominators):
    """Return the common-u/double-2 status of a Lorenz survivor."""
    capacities = tuple(capacity(D, d) for d in denominators)
    for index, d in enumerate(denominators):
        if d in (2, 3) and support_count > (
            sum(capacities) - capacities[index]
        ):
            return "common-u"
    if denominators[0] == denominators[1] == 2:
        # Reflection x -> D-1-x swaps parity and preserves S_D, so each
        # parity load is exactly support_count/2.  Compare that load with
        # the third needle's sharp q=2 fibre cap, not with half of its
        # global capacity (the two quantities need not agree).
        require(support_count % 2 == 0, "double-2 support lost parity balance")
        third = denominators[2]
        ell = (third + 6) // 7
        common = gcd(third, 2)
        height = D // lcm(third, 2)
        third_parity_cap = height * ((ell + common - 1) // common)
        if support_count // 2 > third_parity_cap:
            return "double-2"
    return "survive"


def controls():
    """Positive exact-cover and hostile overloaded-fiber controls."""
    D = 420
    ds = (20, 28, 30)
    require(lcm(*ds) == D, "positive-control lcm changed")
    masks = (
        exact_needle(D, ds[0], 3, 3),
        exact_needle(D, ds[1], 5, 5),
        exact_needle(D, ds[2], 7, 7),
    )
    union = masks[0] | masks[1] | masks[2]
    arcs = bitset_arcs(union, D)
    for q in support.divisors(D):
        profile = support_lorenz_profile(arcs, q)
        needle_profiles = tuple(
            needle_lorenz_parameters(D, d, q) for d in ds
        )
        require(
            phase_free_violation(q, profile, *needle_profiles) is None,
            ("exact-cover positive control rejected", q),
        )

    # Four support points in one mod-4 fiber cannot be covered by three
    # diagonal 28-needles: each is a section over Z/4Z.
    D = 28
    q = 4
    ds = (28, 28, 28)
    hostile = sum(1 << point for point in (0, 4, 8, 12))
    profile = support_lorenz_profile(bitset_arcs(hostile, D), q)
    needle_profiles = tuple(
        needle_lorenz_parameters(D, d, q) for d in ds
    )
    witness = phase_free_violation(q, profile, *needle_profiles)
    require(witness == (1, 1, 4, 3), "hostile control lost unit overload")
    require(
        activity_status(42, 13, (2, 7, 42)) == "common-u",
        "common-u hostile control failed",
    )
    require(
        activity_status(42, 12, (2, 7, 42)) == "survive",
        "common-u boundary control failed",
    )
    require(
        activity_status(28, 10, (2, 2, 28)) == "double-2",
        "double-2 hostile control failed",
    )
    require(
        activity_status(28, 4, (2, 2, 28)) == "survive",
        "double-2 equality boundary control failed",
    )
    return union.bit_count(), witness


def valuation(number, prime):
    result = 0
    while number % prime == 0:
        result += 1
        number //= prime
    return result


def main():
    positive_size, hostile_witness = controls()

    # Reconstruct the support-hard universe exactly as in the frozen referee.
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

    baseline_occurrences = 0
    baseline_shapes = set()
    baseline_rows = set()
    baseline_diagonal = 0
    single_fiber_killed = 0
    single_fiber_diagonal_killed = 0
    higher_lorenz_killed = 0
    killed_occurrences = 0
    surviving_occurrences = 0
    killed_shapes = set()
    surviving_shapes = set()
    killed_rows = set()
    surviving_rows = set()
    first_kill_q = Counter()
    best_margin_histogram = Counter()
    survivor_v7_patterns = Counter()
    survivor_equality_patterns = Counter()
    survivor_full_D_counts = Counter()
    survivor_pair_lcm_counts = Counter()
    survivor_best_q = Counter()
    survivor_best_height = Counter()
    survivor_best_slack = Counter()
    survivor_gcd_patterns = Counter()
    survivor_ratio_patterns = Counter()
    survivor_Ds = Counter()
    survivor_records = []
    kill_witnesses = []
    semantic_hash = sha256()
    activity_status_counts = Counter()
    combined_shapes = set()
    combined_rows = set()
    combined_best_height = Counter()
    combined_best_slack = Counter()
    combined_denominator_prefixes = Counter()
    combined_equality_patterns = Counter()
    combined_full_D_counts = Counter()
    combined_best_q = Counter()
    combined_Ds = Counter()
    combined_records = []
    combined_semantic_hash = sha256()

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
        arcs_by_row = []
        top_loads_by_row = []
        for support_count, F, _L, ranges in rows:
            arcs = combined.projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support arc mismatch", F, D),
            )
            arcs_by_row.append(arcs)
            top_loads_by_row.append(
                tuple(
                    combined.top_class_load(
                        combined.residue_load_histogram(arcs, d),
                        (d + 6) // 7,
                    )
                    for d in divisors
                )
            )

        # q=D cannot improve cardinality for a set-valued support.  q=1 is
        # cardinality itself.  Test every proper nontrivial divisor.
        quotients = tuple(q for q in divisors if q < D)
        quotient_profiles = {
            q: tuple(
                support_lorenz_profile(arcs, q) for arcs in arcs_by_row
            )
            for q in quotients
        }
        needle_profiles = {
            q: tuple(
                needle_lorenz_parameters(D, d, q) for d in divisors
            )
            for q in quotients
        }

        for first_index, first in enumerate(divisors):
            for second_index in range(first_index, len(divisors)):
                second = divisors[second_index]
                first_second_lcm = lcm(first, second)
                for third_index in range(second_index, len(divisors)):
                    third = divisors[third_index]
                    if lcm(first_second_lcm, third) != D:
                        continue
                    # The frozen final stage tests these separate top-class
                    # relaxations.  Rows are sorted only for the coarse cap;
                    # the all-top inequality itself is row-specific.
                    full_capacity = (
                        capacities[first_index]
                        + capacities[second_index]
                        + capacities[third_index]
                    )
                    row_count = bisect_right(supports, full_capacity)
                    if row_count == 0:
                        continue
                    shape = (D, first, second, third)
                    any_baseline = False
                    for row_index in range(row_count):
                        support_count = supports[row_index]
                        top_loads = top_loads_by_row[row_index]
                        if support_count > (
                            top_loads[first_index]
                            + top_loads[second_index]
                            + top_loads[third_index]
                        ):
                            continue
                        any_baseline = True
                        baseline_occurrences += 1
                        baseline_shapes.add(shape)
                        row = (D, rows[row_index][1])
                        baseline_rows.add(row)
                        if first == third:
                            baseline_diagonal += 1

                        single_fiber_violation = None
                        closest_surviving_margin = None
                        # Divisors with larger q tend to expose thin lifted
                        # fibers first; ordering changes only the diagnostic
                        # first-kill ledger, never the exhaustive result.
                        for q in reversed(quotients):
                            support_profile = quotient_profiles[q][row_index]
                            profiles = (
                                needle_profiles[q][first_index],
                                needle_profiles[q][second_index],
                                needle_profiles[q][third_index],
                            )
                            left = support_lorenz_value(
                                support_profile,
                                1,
                            )
                            right = sum(
                                base
                                + (height if high_count else 0)
                                for base, height, high_count in profiles
                            )
                            if left > right:
                                single_fiber_violation = (
                                    q,
                                    (left - right, 1, left, right),
                                )
                                break

                        if single_fiber_violation is not None:
                            killed_occurrences += 1
                            single_fiber_killed += 1
                            if first == third:
                                single_fiber_diagonal_killed += 1
                            killed_shapes.add(shape)
                            killed_rows.add(row)
                            first_kill_q[single_fiber_violation[0]] += 1
                            best_margin_histogram[
                                single_fiber_violation[1][0]
                            ] += 1
                            if len(kill_witnesses) < 20:
                                kill_witnesses.append(
                                    (
                                        "s=1",
                                        rows[row_index][1],
                                        D,
                                        (first, second, third),
                                        single_fiber_violation,
                                    )
                                )
                            continue

                        first_violation = None
                        for q in reversed(quotients):
                            witness_any_sign = phase_free_best_margin(
                                q,
                                quotient_profiles[q][row_index],
                                needle_profiles[q][first_index],
                                needle_profiles[q][second_index],
                                needle_profiles[q][third_index],
                            )
                            if (
                                closest_surviving_margin is None
                                or witness_any_sign[0]
                                > closest_surviving_margin[1][0]
                                or (
                                    witness_any_sign[0]
                                    == closest_surviving_margin[1][0]
                                    and q > closest_surviving_margin[0]
                                )
                            ):
                                closest_surviving_margin = (
                                    q,
                                    witness_any_sign,
                                )
                            if witness_any_sign[0] > 0:
                                first_violation = (q, witness_any_sign)
                                break

                        if first_violation is not None:
                            killed_occurrences += 1
                            higher_lorenz_killed += 1
                            killed_shapes.add(shape)
                            killed_rows.add(row)
                            first_kill_q[first_violation[0]] += 1
                            best_margin_histogram[first_violation[1][0]] += 1
                            if len(kill_witnesses) < 20:
                                kill_witnesses.append(
                                    (
                                        "s>1",
                                        rows[row_index][1],
                                        D,
                                        (first, second, third),
                                        first_violation,
                                    )
                                )
                            continue

                        surviving_occurrences += 1
                        require(
                            closest_surviving_margin is not None,
                            ("survivor has no proper quotient", D),
                        )
                        surviving_shapes.add(shape)
                        surviving_rows.add(row)
                        v7_pattern = tuple(
                            valuation(d, 7) for d in (first, second, third)
                        )
                        if first == third:
                            equality_pattern = "all_equal"
                        elif first == second:
                            equality_pattern = "first_pair_equal"
                        elif second == third:
                            equality_pattern = "last_pair_equal"
                        else:
                            equality_pattern = "all_distinct"
                        full_D_count = sum(
                            d == D for d in (first, second, third)
                        )
                        pair_lcm_count = sum(
                            lcm(a, b) == D
                            for a, b in (
                                (first, second),
                                (first, third),
                                (second, third),
                            )
                        )
                        common_gcd = gcd(gcd(first, second), third)
                        gcd_pattern = (
                            common_gcd,
                            gcd(first, second),
                            gcd(first, third),
                            gcd(second, third),
                        )
                        ratio_pattern = (
                            D // first,
                            D // second,
                            D // third,
                        )
                        survivor_v7_patterns[v7_pattern] += 1
                        survivor_equality_patterns[equality_pattern] += 1
                        survivor_full_D_counts[full_D_count] += 1
                        survivor_pair_lcm_counts[pair_lcm_count] += 1
                        best_q, best_margin = closest_surviving_margin
                        best_height = D // best_q
                        best_slack = -best_margin[0]
                        survivor_best_q[best_q] += 1
                        survivor_best_height[best_height] += 1
                        survivor_best_slack[best_slack] += 1
                        survivor_gcd_patterns[gcd_pattern] += 1
                        survivor_ratio_patterns[ratio_pattern] += 1
                        survivor_Ds[D] += 1
                        record = (
                            rows[row_index][1],
                            rows[row_index][2],
                            D,
                            support_count,
                            (first, second, third),
                            best_q,
                            best_height,
                            best_slack,
                            best_margin[1],
                            best_margin[2],
                            best_margin[3],
                        )
                        survivor_records.append(record)
                        semantic_hash.update(f"{record}\n".encode())
                        status = activity_status(
                            D,
                            support_count,
                            (first, second, third),
                        )
                        activity_status_counts[status] += 1
                        if status == "survive":
                            combined_shapes.add(shape)
                            combined_rows.add(row)
                            combined_best_height[best_height] += 1
                            combined_best_slack[best_slack] += 1
                            combined_denominator_prefixes[
                                (first, second)
                            ] += 1
                            combined_equality_patterns[
                                equality_pattern
                            ] += 1
                            combined_full_D_counts[full_D_count] += 1
                            combined_best_q[best_q] += 1
                            combined_Ds[D] += 1
                            combined_records.append(record)
                            combined_semantic_hash.update(
                                f"{record}\n".encode()
                            )
                    require(
                        any_baseline == (shape in baseline_shapes),
                        "shape baseline flag drifted",
                    )

    require(
        baseline_occurrences == 544571,
        "frozen all-top occurrence ledger not reconstructed",
    )
    require(len(baseline_shapes) == 36614, "baseline shape ledger changed")
    require(len(baseline_rows) == 13577, "baseline row ledger changed")
    require(baseline_diagonal == 2636, "baseline diagonal ledger changed")
    require(
        single_fiber_killed == 125060,
        "independent single-fiber cap ledger changed",
    )
    require(
        single_fiber_diagonal_killed == baseline_diagonal,
        "single-fiber cap did not empty diagonal ledger",
    )
    require(
        killed_occurrences + surviving_occurrences == baseline_occurrences,
        "Lorenz screen lost occurrences",
    )
    require(higher_lorenz_killed == 147, "higher Lorenz ledger changed")
    require(killed_occurrences == 125207, "total Lorenz kill ledger changed")
    require(
        surviving_occurrences == 419364,
        "total Lorenz residual ledger changed",
    )
    require(
        semantic_hash.hexdigest() == EXPECTED_LORENZ_SEMANTIC_SHA256,
        "Lorenz residual semantic ledger changed",
    )
    require(
        activity_status_counts
        == Counter(
            {"common-u": 383389, "double-2": 6756, "survive": 29219}
        ),
        "Lorenz/activity cross-tab changed",
    )
    require(len(combined_shapes) == 4970, "combined shape ledger changed")
    require(len(combined_rows) == 2974, "combined row ledger changed")
    require(
        combined_semantic_hash.hexdigest()
        == EXPECTED_COMBINED_SEMANTIC_SHA256,
        "combined residual semantic ledger changed",
    )

    print("LRC14 mixed three-drift divisor-fiber Lorenz/activity referee")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"combined_output_sha256={file_sha256(COMBINED_OUTPUT_PATH)}")
    print(f"positive_exact_cover_union_size={positive_size}")
    print(f"hostile_overload_witness={hostile_witness}")
    print(f"baseline_occurrences={baseline_occurrences}")
    print(f"baseline_shapes={len(baseline_shapes)}")
    print(f"baseline_rows={len(baseline_rows)}")
    print(f"baseline_diagonal_occurrences={baseline_diagonal}")
    print(f"single_fiber_killed_occurrences={single_fiber_killed}")
    print(
        "single_fiber_diagonal_killed="
        f"{single_fiber_diagonal_killed}"
    )
    print(f"higher_lorenz_killed_occurrences={higher_lorenz_killed}")
    print(f"lorenz_killed_occurrences={killed_occurrences}")
    print(f"lorenz_surviving_occurrences={surviving_occurrences}")
    print(f"killed_shapes={len(killed_shapes)}")
    print(f"surviving_shapes={len(surviving_shapes)}")
    print(f"killed_rows={len(killed_rows)}")
    print(f"surviving_rows={len(surviving_rows)}")
    print(f"first_kill_q_top={first_kill_q.most_common(20)}")
    print(
        "first_positive_margin_histogram="
        f"{sorted(best_margin_histogram.items())}"
    )
    print(
        "survivor_v7_patterns="
        f"{survivor_v7_patterns.most_common(20)}"
    )
    print(
        "survivor_equality_patterns="
        f"{survivor_equality_patterns.most_common()}"
    )
    print(
        "survivor_full_D_counts="
        f"{sorted(survivor_full_D_counts.items())}"
    )
    print(
        "survivor_pair_lcm_counts="
        f"{sorted(survivor_pair_lcm_counts.items())}"
    )
    print(f"survivor_best_q_top={survivor_best_q.most_common(20)}")
    print(
        "survivor_best_fiber_height="
        f"{sorted(survivor_best_height.items())}"
    )
    print(
        "survivor_min_slack_histogram="
        f"{sorted(survivor_best_slack.items())}"
    )
    print(
        "survivor_gcd_patterns_top="
        f"{survivor_gcd_patterns.most_common(20)}"
    )
    print(
        "survivor_ratio_patterns_top="
        f"{survivor_ratio_patterns.most_common(20)}"
    )
    print(f"survivor_Ds_top={survivor_Ds.most_common(20)}")
    print(f"kill_witnesses={kill_witnesses}")
    print(f"survivor_first={survivor_records[:5]}")
    print(f"survivor_last={survivor_records[-5:]}")
    print(f"survivor_semantic_sha256={semantic_hash.hexdigest()}")
    print(f"activity_status_counts={activity_status_counts}")
    print(f"combined_residual_shapes={len(combined_shapes)}")
    print(f"combined_residual_rows={len(combined_rows)}")
    print(
        "combined_best_fiber_height="
        f"{sorted(combined_best_height.items())}"
    )
    print(
        "combined_best_slack="
        f"{sorted(combined_best_slack.items())}"
    )
    print(
        "combined_denominator_prefixes_top="
        f"{combined_denominator_prefixes.most_common(20)}"
    )
    print(
        "combined_equality_patterns="
        f"{combined_equality_patterns.most_common()}"
    )
    print(
        "combined_full_D_counts="
        f"{sorted(combined_full_D_counts.items())}"
    )
    print(f"combined_best_q_top={combined_best_q.most_common(20)}")
    print(f"combined_Ds_top={combined_Ds.most_common(20)}")
    print(f"combined_first={combined_records[:5]}")
    print(f"combined_last={combined_records[-5:]}")
    print(
        "combined_semantic_sha256="
        f"{combined_semantic_hash.hexdigest()}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
