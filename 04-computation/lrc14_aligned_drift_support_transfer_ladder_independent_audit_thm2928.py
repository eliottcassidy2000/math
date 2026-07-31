#!/usr/bin/env python3
"""Independent audit of the all-k support-transfer ladder referee.

This deliberately avoids the referee's two central implementations:

* body support is projected by merged cyclic arcs, never by a bitset; and
* exact-lcm denominator multisets are counted by divisor-poset recurrence,
  never by an explicitly evaluated Moebius sum.

The imported referee is used only as the object under audit: its pinned source,
constants, and Moebius formula are compared against independently generated
data.  All guards remain active under ``python -O``.
"""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TARGET_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_aligned_drift_support_transfer_ladder_thm2928.py"
)
EXPECTED_TARGET_SHA256 = (
    "8db781fb3e7dc8fdc4df2bf3c6d83869a9ffe52f41c7d70c25bbd0a9b0122bea"
)

# These are copied from their proved theorem sources, rather than read from the
# target at run time.  Index j means the universal safe floor for j combs.
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
EXPECTED_ROWS = (27210, 27240, 27163, 26970, 13778, 10976, 6237, 0)
EXPECTED_BODIES = (3003, 3003, 3003, 3003, 3003, 3003, 3003, 0)
EXPECTED_DIVISORS = (219, 219, 219, 217, 206, 193, 171, 0)
EXPECTED_SHAPES = (
    161535777082757,
    3095010121875,
    50874159718,
    694921995,
    7483350,
    56419,
    171,
    0,
)
EXPECTED_OCCURRENCES = (
    1504842061942849,
    38954725590760,
    951545890235,
    21357714101,
    298255882,
    3066274,
    6237,
    0,
)
EXPECTED_SEMANTIC_SHA256 = (
    "e7ecc3d7aed7b692564a03734f0f240d4f9c6ca9d68d012a66f5a4a626d0a7c4",
    "931e04a95f70e4b16a3004496df981a0678f0839dbfa85bc6194d89290a4a76b",
    "2c2909d3d1e514b4e8f09166d23145ac6b577fe6e10f39d39a205094e7656dd5",
    "3738fce411911b8cfb64596a6ed595a01d54f9093b502b4a4161fcd63ac91b2a",
    "86109feddd36a3ece9068d4e50d091b8423930a226f3c50c5024578df41572b0",
    "38aa87113f5acf83273d756e1ad9c261f745d1f3ae106b6d1d1aa5d02ef03739",
    "ee5bd68fe275192c6eeb6b09081982c226251b0305268fdbb91df894da250941",
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(TARGET_PATH) == EXPECTED_TARGET_SHA256,
    "support-transfer referee changed since this audit was written",
)
spec = spec_from_file_location("support_transfer_ladder_under_audit", TARGET_PATH)
target = module_from_spec(spec)
spec.loader.exec_module(target)


def divisors(number):
    """Increasing positive divisors, generated without target code."""
    low = []
    high = []
    candidate = 1
    while candidate * candidate <= number:
        quotient, remainder = divmod(number, candidate)
        if remainder == 0:
            low.append(candidate)
            if quotient != candidate:
                high.append(quotient)
        candidate += 1
    return tuple(low + high[::-1])


def body_safe_ranges(body):
    """Exact full-safe cell ranges from a fresh integer endpoint sweep."""
    L = 14 * lcm(*body)
    forbidden = []
    for speed in body:
        period = L // speed
        radius = period // 14
        require(14 * radius == period, ("nonintegral tooth radius", body, speed))
        for tooth in range(speed + 1):
            center = tooth * period
            forbidden.append(
                (max(0, center - radius), min(L, center + radius))
            )
    forbidden.sort()

    # A cell j is forbidden precisely when its interior meets a strict danger
    # interval, namely left <= j < right.  Thus touching a danger endpoint is
    # retained as safe, which is the required open-danger convention.
    merged = []
    for left, right in forbidden:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))

    safe = []
    cursor = 0
    for left, right in merged:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < L:
        safe.append((cursor, L))
    return L, tuple(safe)


def projected_support_by_arcs(modulus, safe_ranges):
    """Cardinality of the image modulo ``modulus`` via cyclic arc merging."""
    arcs = []
    for left, right in safe_ranges:
        length = right - left
        if length >= modulus:
            return modulus
        start = left % modulus
        end = start + length
        if end <= modulus:
            arcs.append((start, end))
        else:
            arcs.append((start, modulus))
            arcs.append((0, end - modulus))
    arcs.sort()
    merged = []
    for left, right in arcs:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return sum(right - left for left, right in merged)


@lru_cache(maxsize=None)
def exact_lcm_multisets(modulus, arity):
    """Exact-lcm multisets by downward divisor-poset subtraction."""
    require(arity >= 1, ("recurrence is only for positive arity", arity))
    divisor_list = divisors(modulus)
    available = len(divisor_list) - 1  # divisors > 1
    all_dividing = comb(available + arity - 1, arity)
    exact = all_dividing
    for proper in divisor_list[:-1]:
        exact -= exact_lcm_multisets(proper, arity)
    require(exact >= 0, ("negative exact-lcm count", modulus, arity, exact))
    return exact


def brute_shape_controls():
    """Direct multiset enumeration on a wider hostile range than the target."""
    cases = 0
    for modulus in range(1, 85):
        allowed = divisors(modulus)[1:]
        for arity in range(1, 6):
            brute = sum(
                lcm(*values) == modulus
                for values in combinations_with_replacement(allowed, arity)
            )
            recurrent = exact_lcm_multisets(modulus, arity)
            require(
                recurrent == brute,
                ("shape recurrence versus brute force", modulus, arity),
            )
            cases += 1

    # Hostile identities catch accidental inclusion of the forbidden divisor 1.
    for prime in (2, 3, 5, 7, 11, 13, 17, 19):
        for arity in range(1, 8):
            require(
                exact_lcm_multisets(prime, arity) == 1,
                ("prime exact-lcm identity failed", prime, arity),
            )
    for modulus in range(2, 85):
        require(
            exact_lcm_multisets(modulus, 1) == 1,
            ("arity-one exact-lcm identity failed", modulus),
        )
    for arity in range(1, 8):
        require(
            exact_lcm_multisets(1, arity) == 0,
            ("D=1 should have no positive-denominator tuple", arity),
        )
    return cases


def endpoint_controls():
    """Check the safe-range endpoint convention on canonical hostile bodies."""
    bodies = (
        (1, 2, 3, 4, 5, 6),
        (1, 4, 8, 10, 12, 13),
        (2, 6, 8, 10, 11, 14),
        (9, 10, 11, 12, 13, 14),
    )
    checks = 0
    for body in bodies:
        L, ranges = body_safe_ranges(body)
        target_L, target_ranges = target.support.safe_cell_ranges(body)
        require(L == target_L, ("body ruler mismatch", body, L, target_L))
        require(
            tuple(ranges) == tuple(target_ranges),
            ("safe-range mismatch", body),
        )
        for left, right in ranges:
            require(0 <= left < right <= L, ("bad safe range", body, left, right))
            for cell in (left, right - 1):
                for speed in body:
                    period = L // speed
                    radius = period // 14
                    # Full closed cell [cell,cell+1] avoids every strict tooth.
                    for tooth in range(speed + 1):
                        center = tooth * period
                        require(
                            cell + 1 <= center - radius
                            or cell >= center + radius,
                            ("unsafe endpoint cell retained", body, cell, speed),
                        )
                checks += 1
    return checks


def main():
    shape_control_cases = brute_shape_controls()
    endpoint_control_cases = endpoint_controls()

    cutoffs = {
        k: (Q(1) - SAFE_FLOORS[7 - k]) / SAFE_FLOORS[k]
        for k in range(8)
    }
    require(cutoffs == EXPECTED_CUTOFFS, ("cutoff arithmetic changed", cutoffs))
    require(target.SAFE_FLOORS == SAFE_FLOORS, "target safe floors differ")
    require(target.EXPECTED_CUTOFFS == EXPECTED_CUTOFFS, "target cutoffs differ")

    rows_by_k_D = {k: Counter() for k in range(8)}
    bodies_by_k = {k: set() for k in range(8)}
    semantic = {k: sha256() for k in range(8)}
    equality_rows = {k: [] for k in range(8)}
    all_moduli = set()
    total_rows = 0
    body_count = 0

    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, safe_ranges = body_safe_ranges(body)
        # Comparing the independently swept ranges on all 3,003 bodies detects
        # endpoint and ruler drift before projection is attempted.
        target_L, target_ranges = target.support.safe_cell_ranges(body)
        require(L == target_L, ("all-body ruler mismatch", body))
        require(tuple(safe_ranges) == tuple(target_ranges), ("all-body ranges", body))

        for modulus in divisors(L):
            total_rows += 1
            all_moduli.add(modulus)
            support_count = projected_support_by_arcs(modulus, safe_ranges)
            require(
                0 < support_count <= modulus,
                ("projected support out of range", body, modulus, support_count),
            )
            for k, cutoff in cutoffs.items():
                density = Q(support_count, modulus)
                if density <= cutoff:
                    rows_by_k_D[k][modulus] += 1
                    bodies_by_k[k].add(body)
                    semantic[k].update(
                        f"{body}|{L}|{modulus}|{support_count}\n".encode()
                    )
                    if density == cutoff and len(equality_rows[k]) < 5:
                        equality_rows[k].append((body, L, modulus, support_count))

    require(body_count == 3003, ("body universe", body_count))
    require(total_rows == 251536, ("body/divisor universe", total_rows))

    # Compare the target's Moebius formula to the independent recurrence on
    # every modulus that actually occurs, for every positive drift arity.
    formula_cases = 0
    for modulus in sorted(all_moduli):
        for arity in range(1, 8):
            recurrent = exact_lcm_multisets(modulus, arity)
            moebius = target.lcm_multiset_shapes(modulus, arity)
            require(
                recurrent == moebius,
                ("Moebius formula mismatch", modulus, arity, recurrent, moebius),
            )
            # The zeta/Moebius defining identity is checked explicitly.
            all_dividing = comb(len(divisors(modulus)) + arity - 2, arity)
            zeta_sum = sum(
                exact_lcm_multisets(d, arity) for d in divisors(modulus)
            )
            require(
                zeta_sum == all_dividing,
                ("divisor-poset zeta identity", modulus, arity),
            )
            formula_cases += 1

    actual_rows = []
    actual_bodies = []
    actual_divisors = []
    actual_shapes = []
    actual_occurrences = []
    actual_semantic = []

    print("Independent audit: LRC14 aligned/drift support-transfer ladder")
    print(f"target_script_sha256={file_sha256(TARGET_PATH)}")
    print(f"shape_bruteforce_cases={shape_control_cases}")
    print(f"safe_endpoint_checks={endpoint_control_cases}")
    print(f"actual_moduli={len(all_moduli)}")
    print(f"recurrence_moebius_cases={formula_cases}")
    print(f"body_count={body_count}")
    print(f"body_divisor_rows={total_rows}")
    print(f"support_cutoffs={cutoffs}")

    for k in range(8):
        arity = 7 - k
        row_count = sum(rows_by_k_D[k].values())
        divisor_count = len(rows_by_k_D[k])
        if arity:
            shape_count = sum(
                exact_lcm_multisets(modulus, arity)
                for modulus in rows_by_k_D[k]
            )
            occurrence_count = sum(
                multiplicity * exact_lcm_multisets(modulus, arity)
                for modulus, multiplicity in rows_by_k_D[k].items()
            )
        else:
            # The all-aligned branch is outside the positive-drift multiset
            # formula and is already empty under its zero support cutoff.
            shape_count = 0
            occurrence_count = 0
        digest = semantic[k].hexdigest()
        actual_rows.append(row_count)
        actual_bodies.append(len(bodies_by_k[k]))
        actual_divisors.append(divisor_count)
        actual_shapes.append(shape_count)
        actual_occurrences.append(occurrence_count)
        actual_semantic.append(digest)
        print(
            f"k={k},p={arity},cutoff={cutoffs[k]},rows={row_count},"
            f"bodies={len(bodies_by_k[k])},divisors={divisor_count},"
            f"denominator_shapes={shape_count},raw_occurrences={occurrence_count},"
            f"equality_examples={equality_rows[k]},semantic_sha256={digest}"
        )

    require(tuple(actual_rows) == EXPECTED_ROWS, "row ladder mismatch")
    require(tuple(actual_bodies) == EXPECTED_BODIES, "body ladder mismatch")
    require(tuple(actual_divisors) == EXPECTED_DIVISORS, "divisor ladder mismatch")
    require(tuple(actual_shapes) == EXPECTED_SHAPES, "shape ladder mismatch")
    require(
        tuple(actual_occurrences) == EXPECTED_OCCURRENCES,
        "occurrence ladder mismatch",
    )
    require(
        tuple(actual_semantic) == EXPECTED_SEMANTIC_SHA256,
        "semantic ladder mismatch",
    )
    require(tuple(actual_rows) == target.EXPECTED_ROWS, "target row constants differ")
    require(
        tuple(actual_shapes) == target.EXPECTED_SHAPES,
        "target shape constants differ",
    )
    require(
        tuple(actual_occurrences) == target.EXPECTED_OCCURRENCES,
        "target occurrence constants differ",
    )
    require(
        tuple(actual_semantic) == target.EXPECTED_SEMANTIC_SHA256,
        "target semantic constants differ",
    )
    print("scope=necessary_body_divisor_rows_and_unordered_denominator_multisets")
    print("scope_not=realized_covers_or_physical_slope_phase_packets")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
