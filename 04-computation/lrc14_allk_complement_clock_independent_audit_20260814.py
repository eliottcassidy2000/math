#!/usr/bin/env python3
"""Independent exact audit of THM-3366's pool-14 complement-clock census.

This audit deliberately avoids the S174 recursive bitset solver.  It imports
only the pinned THM-2928 safe-cell geometry, then rebuilds the pool-14 census
from scratch by:

* recomputing each projected support count from merged cyclic arcs;
* constructing exact unsupported-cell targets from strict endpoints plus open
  arrangement atoms;
* enumerating every subset of clocks {1,...,14} of size <= 5 (and <= 4 for
  the k=0 budget) to classify each unique target by its least exact cover; and
* recomputing denominator-multiset occurrence weights from an independently
  coded divisor-Mobius formula checked against brute-force multisets.

All guards use RuntimeError and remain active under ``python -O``.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "04-computation/lrc14_two_drift_body_projection_support_thm2928.py"
TARGET_SCRIPT_PATH = (
    ROOT / "04-computation/lrc14_allk_universal_complement_clock_scan_kps_s174.py"
)
TARGET_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_allk_universal_complement_clock_scan_kps_s174.out"
)

EXPECTED_SUPPORT_SHA256 = "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a"
EXPECTED_TARGET_SCRIPT_SHA256 = (
    "372f1b0d2bd8c1cc453080e9ab55880480352c5277a9a216739a7262d198efba"
)
EXPECTED_TARGET_OUTPUT_SHA256 = (
    "1104900cd805a2af05e3c7252b7d7bacd417b36b0f015a90d1976b6913fc91c5"
)

POOL = tuple(range(1, 15))
POINT_BUDGET = 5
K0_BUDGET = 4

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
CUTOFFS = {k: (Q(1) - SAFE_FLOORS[7 - k]) / SAFE_FLOORS[k] for k in range(8)}
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

EXPECTED_INPUT_ROWS = (27210, 27240, 27163, 26970, 13778, 10976, 6237)
EXPECTED_INPUT_OCCURRENCES = (
    1504842061942849,
    38954725590760,
    951545890235,
    21357714101,
    298255882,
    3066274,
    6237,
)
EXPECTED_CLOSED_ROWS = (17561, 19273, 19198, 19053, 5964, 3334, 570)
EXPECTED_CLOSED_OCCURRENCES = (
    113853177071279,
    5116993586195,
    182987150531,
    5887257171,
    47788532,
    528176,
    570,
)
EXPECTED_SIZE_HISTOGRAMS = (
    ((1, 12662), (2, 2766), (3, 1006), (4, 1127)),
    ((1, 12665), (2, 2766), (3, 1007), (4, 1142), (5, 1693)),
    ((1, 12662), (2, 2764), (3, 998), (4, 1106), (5, 1668)),
    ((1, 12659), (2, 2764), (3, 976), (4, 1052), (5, 1602)),
    ((2, 2404), (3, 959), (4, 1040), (5, 1561)),
    ((2, 398), (3, 692), (4, 780), (5, 1464)),
    ((4, 173), (5, 397)),
)
EXPECTED_OCCURRENCE_HISTOGRAMS = (
    ((1, 89654307406587), (2, 6978600287173), (3, 3540422163243), (4, 13679847214276)),
    ((1, 3427334397071), (2, 332132448500), (3, 179613246178), (4, 497575939015), (5, 680337555431)),
    ((1, 122664738423), (2, 14373925905), (3, 7408833370), (4, 16153083294), (5, 22386569539)),
    ((1, 3971732547), (2, 550517943), (3, 264376426), (4, 451833167), (5, 648797088)),
    ((2, 13640121), (3, 7236404), (4, 10653866), (5, 16258141)),
    ((2, 59665), (3, 97596), (4, 121920), (5, 248995)),
    ((4, 173), (5, 397)),
)
EXPECTED_DISTINCT_TARGETS = 11180
EXPECTED_ARRANGEMENT_POINTS = 206
EXPECTED_ATOMS = 205


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def normalized_sha256(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_support():
    require(
        normalized_sha256(SUPPORT_PATH) == EXPECTED_SUPPORT_SHA256,
        "THM-2928 support dependency changed",
    )
    spec = spec_from_file_location("lrc14_support_thm2928", SUPPORT_PATH)
    require(spec is not None and spec.loader is not None, "support import failed")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def divisors(number):
    low = []
    high = []
    candidate = 1
    while candidate * candidate <= number:
        if number % candidate == 0:
            low.append(candidate)
            quotient = number // candidate
            if quotient != candidate:
                high.append(quotient)
        candidate += 1
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


def exact_lcm_multiset_count(modulus, arity):
    total = 0
    for divisor in divisors(modulus):
        divisor_count = len(divisors(divisor)) - 1
        total += mobius(modulus // divisor) * comb(divisor_count + arity - 1, arity)
    require(total >= 0, ("negative exact-lcm count", modulus, arity, total))
    return total


def brute_mobius_controls():
    cases = 0
    for modulus in range(1, 41):
        allowed = tuple(d for d in divisors(modulus) if d > 1)
        for arity in range(1, 8):
            brute = sum(
                1
                for values in combinations_with_replacement(allowed, arity)
                if lcm_tuple(values) == modulus
            )
            exact = exact_lcm_multiset_count(modulus, arity)
            require(
                exact == brute,
                ("Mobius multiset mismatch", modulus, arity, exact, brute),
            )
            cases += 1
    return cases


def lcm_tuple(values):
    current = 1
    for value in values:
        current = (current * value) // gcd(current, value)
    return current


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def arrangement_points(clocks):
    points = {Q(0), Q(1)}
    for clock in clocks:
        for tooth in range(clock + 1):
            for sign in (-1, 1):
                point = Q(14 * tooth + sign, 14 * clock)
                if Q(0) <= point <= Q(1):
                    points.add(point)
    return tuple(sorted(points))


def circular_distance(value):
    residue = value - (value.numerator // value.denominator)
    return residue if residue <= 1 - residue else 1 - residue


def danger(clock, point):
    return circular_distance(clock * point) < Q(1, 14)


def projected_support_arcs(modulus, safe_ranges):
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
    pieces.sort()
    merged = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def unsupported_gaps(modulus, arcs):
    gaps = []
    cursor = 0
    for left, right in arcs:
        if cursor < left:
            gaps.append((Q(cursor, modulus), Q(left, modulus)))
        cursor = right
    if cursor < modulus:
        gaps.append((Q(cursor, modulus), Q(1)))
    return tuple(gaps)


def make_clock_masks(points, atoms):
    masks = {}
    atom_offset = len(points)
    for clock in POOL:
        mask = 0
        for index, point in enumerate(points):
            if danger(clock, point):
                mask |= 1 << index
        for index, (left, right) in enumerate(atoms):
            midpoint = (left + right) / 2
            if danger(clock, midpoint):
                mask |= 1 << (atom_offset + index)
        masks[clock] = mask
    return masks


def target_masks(modulus, gaps, points, atoms):
    strict_endpoints = 0
    strict_endpoints_with_grid = 0
    atom_mask = 0
    for index, point in enumerate(points):
        if any(left < point < right for left, right in gaps):
            strict_endpoints_with_grid |= 1 << index
            if (modulus * point).denominator != 1:
                strict_endpoints |= 1 << index
    atom_offset = len(points)
    for index, (left, right) in enumerate(atoms):
        if any(max(left, gap_left) < min(right, gap_right) for gap_left, gap_right in gaps):
            atom_mask |= 1 << (atom_offset + index)
    return (
        strict_endpoints | atom_mask,
        atom_mask,
        strict_endpoints_with_grid | atom_mask,
        strict_endpoints,
    )


def enumerate_subset_masks(clock_masks):
    by_size = {size: [] for size in range(POINT_BUDGET + 1)}
    by_size[0].append(((), 0))
    for size in range(1, POINT_BUDGET + 1):
        for subset in combinations(POOL, size):
            mask = 0
            for clock in subset:
                mask |= clock_masks[clock]
            by_size[size].append((subset, mask))
    return by_size


def least_cover(mask, subset_masks_by_size, max_size, cache):
    key = (mask, max_size)
    if key in cache:
        return cache[key]
    for size in range(max_size + 1):
        for subset, cover_mask in subset_masks_by_size[size]:
            if mask & ~cover_mask == 0:
                cache[key] = subset
                return subset
    cache[key] = None
    return None


def first_hostile_endpoint_witness(rows, exact_cache, atom_cache, subset_masks_by_size):
    for row in rows:
        atom_subset = least_cover(row["atom_target"], subset_masks_by_size, POINT_BUDGET, atom_cache)
        exact_subset = exact_cache[row["target"]]
        if atom_subset is None:
            continue
        if exact_subset is None or len(atom_subset) < len(exact_subset):
            return {
                "body": row["body"],
                "L": row["L"],
                "D": row["D"],
                "support": row["support"],
                "atom_only": atom_subset,
                "exact": exact_subset,
                "strict_endpoint_bits": bit_count(row["endpoint_mask"]),
            }
    raise RuntimeError("failed to find midpoint-only hostile witness")


def first_grid_ownership_witness(rows, exact_cache, wrong_cache, subset_masks_by_size):
    for row in rows:
        if row["wrong_target"] == row["target"]:
            continue
        exact_subset = exact_cache[row["target"]]
        wrong_subset = least_cover(row["wrong_target"], subset_masks_by_size, POINT_BUDGET, wrong_cache)
        if exact_subset is not None and (wrong_subset is None or len(wrong_subset) != len(exact_subset)):
            return {
                "body": row["body"],
                "L": row["L"],
                "D": row["D"],
                "support": row["support"],
                "extra_grid_bits": bit_count(row["wrong_target"] ^ row["target"]),
                "exact": exact_subset,
                "without_ownership": wrong_subset,
            }
    for row in rows:
        if row["wrong_target"] != row["target"]:
            return {
                "body": row["body"],
                "L": row["L"],
                "D": row["D"],
                "support": row["support"],
                "extra_grid_bits": bit_count(row["wrong_target"] ^ row["target"]),
                "exact": exact_cache[row["target"]],
                "without_ownership": least_cover(row["wrong_target"], subset_masks_by_size, POINT_BUDGET, wrong_cache),
            }
    raise RuntimeError("failed to find D-grid ownership witness")


def first_budget_witness(rows, exact_cache, subset_masks_by_size, budget):
    smaller_cache = {}
    for row in rows:
        exact_subset = exact_cache[row["target"]]
        if exact_subset is None or len(exact_subset) != budget:
            continue
        insufficient = least_cover(row["target"], subset_masks_by_size, budget - 1, smaller_cache)
        if insufficient is None:
            return {
                "body": row["body"],
                "L": row["L"],
                "D": row["D"],
                "support": row["support"],
                "exact": exact_subset,
                "failed_budget": budget - 1,
            }
    raise RuntimeError(("failed to find insufficient-budget witness", budget))


def bit_count(number):
    return number.bit_count()


def format_histogram(counter):
    return tuple(sorted(counter.items()))


def main():
    require(CUTOFFS == EXPECTED_CUTOFFS, ("cutoffs changed", CUTOFFS))
    require(
        normalized_sha256(TARGET_SCRIPT_PATH) == EXPECTED_TARGET_SCRIPT_SHA256,
        "THM-3366 script changed",
    )
    require(
        normalized_sha256(TARGET_OUTPUT_PATH) == EXPECTED_TARGET_OUTPUT_SHA256,
        "THM-3366 output changed",
    )

    support = load_support()
    mobius_control_cases = brute_mobius_controls()

    points = arrangement_points(POOL)
    atoms = tuple(zip(points, points[1:]))
    require(len(points) == EXPECTED_ARRANGEMENT_POINTS, "arrangement point count changed")
    require(len(atoms) == EXPECTED_ATOMS, "arrangement atom count changed")

    clock_masks = make_clock_masks(points, atoms)
    subset_masks_by_size = enumerate_subset_masks(clock_masks)

    rows = []
    unique_targets = {}
    total_body_divisor_rows = 0
    body_count = 0
    arc_projection_controls = 0

    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, safe_ranges = support.safe_cell_ranges(body)
        for D in divisors(L):
            total_body_divisor_rows += 1
            arcs = projected_support_arcs(D, safe_ranges)
            support_count = sum(right - left for left, right in arcs)
            control_count = support.support_size_bitset(D, safe_ranges)
            require(
                support_count == control_count,
                ("support projection mismatch", body, L, D, support_count, control_count),
            )
            arc_projection_controls += 1
            density = Q(support_count, D)
            eligible = tuple(k for k in range(7) if density <= CUTOFFS[k])
            if not eligible:
                continue
            gaps = unsupported_gaps(D, arcs)
            target, atom_target, wrong_target, endpoint_mask = target_masks(
                D, gaps, points, atoms
            )
            row = {
                "body": body,
                "L": L,
                "D": D,
                "support": support_count,
                "eligible": eligible,
                "target": target,
                "atom_target": atom_target,
                "wrong_target": wrong_target,
                "endpoint_mask": endpoint_mask,
            }
            rows.append(row)
            if target not in unique_targets:
                unique_targets[target] = row

    require(body_count == 3003, "body universe changed")
    require(total_body_divisor_rows == 251536, "body/divisor universe changed")
    require(
        arc_projection_controls == total_body_divisor_rows,
        "arc projection control count changed",
    )
    require(len(unique_targets) == EXPECTED_DISTINCT_TARGETS, "distinct target count changed")

    cover_cache = {}
    exact_subsets = {}
    semantic_targets = sha256()
    target_size_histogram = Counter()
    uncovered_targets = 0
    for target in sorted(unique_targets):
        subset = least_cover(target, subset_masks_by_size, POINT_BUDGET, cover_cache)
        exact_subsets[target] = subset
        semantic_targets.update(f"{target}|{subset}\n".encode())
        if subset is None:
            uncovered_targets += 1
            target_size_histogram[None] += 1
        else:
            target_size_histogram[len(subset)] += 1

    input_rows = [0] * 7
    input_occurrences = [0] * 7
    closed_rows = [0] * 7
    closed_occurrences = [0] * 7
    size_histograms = [Counter() for _ in range(7)]
    occurrence_histograms = [Counter() for _ in range(7)]
    semantics = [sha256() for _ in range(7)]
    occurrence_cache = {}

    for row in rows:
        exact_subset = exact_subsets[row["target"]]
        for k in row["eligible"]:
            p = 7 - k
            weight_key = (row["D"], p)
            if weight_key not in occurrence_cache:
                occurrence_cache[weight_key] = exact_lcm_multiset_count(row["D"], p)
            weight = occurrence_cache[weight_key]
            input_rows[k] += 1
            input_occurrences[k] += weight
            budget = K0_BUDGET if k == 0 else POINT_BUDGET
            if exact_subset is None or len(exact_subset) > budget:
                continue
            closed_rows[k] += 1
            closed_occurrences[k] += weight
            size_histograms[k][len(exact_subset)] += 1
            occurrence_histograms[k][len(exact_subset)] += weight
            semantics[k].update(
                (
                    f"{row['body']}|{row['L']}|{row['D']}|{row['support']}|"
                    f"{exact_subset}\n"
                ).encode()
            )

    require(tuple(input_rows) == EXPECTED_INPUT_ROWS, ("input row mismatch", tuple(input_rows)))
    require(
        tuple(input_occurrences) == EXPECTED_INPUT_OCCURRENCES,
        ("input occurrence mismatch", tuple(input_occurrences)),
    )
    require(tuple(closed_rows) == EXPECTED_CLOSED_ROWS, ("closed row mismatch", tuple(closed_rows)))
    require(
        tuple(closed_occurrences) == EXPECTED_CLOSED_OCCURRENCES,
        ("closed occurrence mismatch", tuple(closed_occurrences)),
    )
    require(
        tuple(format_histogram(counter) for counter in size_histograms)
        == EXPECTED_SIZE_HISTOGRAMS,
        "minimal-size histogram mismatch",
    )
    require(
        tuple(format_histogram(counter) for counter in occurrence_histograms)
        == EXPECTED_OCCURRENCE_HISTOGRAMS,
        "occurrence histogram mismatch",
    )

    atom_cache = {}
    wrong_cache = {}
    hostile_endpoint = first_hostile_endpoint_witness(rows, exact_subsets, atom_cache, subset_masks_by_size)
    hostile_grid = first_grid_ownership_witness(rows, exact_subsets, wrong_cache, subset_masks_by_size)
    hostile_budget5 = first_budget_witness(rows, exact_subsets, subset_masks_by_size, POINT_BUDGET)
    hostile_budget4 = first_budget_witness(
        [row for row in rows if 0 in row["eligible"]],
        exact_subsets,
        subset_masks_by_size,
        K0_BUDGET,
    )

    print("LRC14 ALL-K COMPLEMENT-CLOCK INDEPENDENT AUDIT")
    print("status=FINITE-EXACT independent subset-enumeration replay of THM-3366 pool-14 census")
    print(
        f"support_sha256={normalized_sha256(SUPPORT_PATH)};"
        f"target_script_sha256={normalized_sha256(TARGET_SCRIPT_PATH)};"
        f"target_output_sha256={normalized_sha256(TARGET_OUTPUT_PATH)}"
    )
    print(
        f"arrangement_points={len(points)};atoms={len(atoms)};"
        f"body_count={body_count};body_divisor_rows={total_body_divisor_rows};"
        f"distinct_targets={len(unique_targets)};mobius_control_cases={mobius_control_cases};"
        f"arc_projection_controls={arc_projection_controls}"
    )
    print(
        f"target_min_size_histogram={tuple(sorted(target_size_histogram.items(), key=lambda item: (99 if item[0] is None else item[0], item[1])))};"
        f"target_semantic_sha256={semantic_targets.hexdigest()};"
        f"uncovered_targets_at_budget5={uncovered_targets}"
    )
    for k in range(7):
        budget = K0_BUDGET if k == 0 else POINT_BUDGET
        print(
            f"k={k};p={7-k};budget={budget};cutoff={CUTOFFS[k]};"
            f"input_rows={input_rows[k]};closed_rows={closed_rows[k]};remaining_rows={input_rows[k]-closed_rows[k]};"
            f"input_occurrences={input_occurrences[k]};closed_occurrences={closed_occurrences[k]};"
            f"remaining_occurrences={input_occurrences[k]-closed_occurrences[k]};"
            f"size_histogram={format_histogram(size_histograms[k])};"
            f"occurrence_histogram={format_histogram(occurrence_histograms[k])};"
            f"audit_semantic_sha256={semantics[k].hexdigest()};"
            "thm3366_compare=PASS"
        )
    print(
        "hostile_midpoint_endpoint_control="
        f"{hostile_endpoint}"
    )
    print(
        "hostile_D_grid_ownership_control="
        f"{hostile_grid}"
    )
    print(
        "hostile_insufficient_budget_control="
        f"{hostile_budget5}"
    )
    print(
        "hostile_k0_budget_control="
        f"{hostile_budget4}"
    )
    print("discrepancy=NONE")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
