#!/usr/bin/env python3
"""Exact coarse address sieve for four aligned combs and three drifts.

Let F be a literal six-speed body in {1,...,14}, let L=L_F, let J be
the exact safe 1/L-cell word, and put S_D=J mod D for D|L.  If four
aligned multiplier combs leave normalized safe mass u_A and three drift
combs cover the residual, then

    (|S_D|/D) u_A <= measure(D_a1 union D_a2 union D_a3).

THM-2928 (25a), from THM-1166, gives u_A>=558/1183.  THM-1166 S77 gives
the sharp union cap 36/91 for three distinct danger combs.  Hence

    |S_D|/D <= (36/91)/(558/1183) = 26/31.

For d|D put k_d=ceil(d/7), let lambda_d(s) count S_D in the residue
class s mod d, and let Lambda_d be the sum of the k_d largest class
loads.  An actual fixed-phase d-mask is contained in the lift of k_d
classes, so its S_D load is at most Lambda_d.  For each unordered
denominator triple d_1<=d_2<=d_3, d_i>1, lcm(d_1,d_2,d_3)=D, this script
checks the nested necessary relaxations

    |S_D| <= C_1+C_2+C_3,
    |S_D| <= Lambda_1+C_2+C_3,
    |S_D| <= Lambda_1+Lambda_2+C_3,
    |S_D| <= Lambda_1+Lambda_2+Lambda_3,

where C_i=(D/d_i)k_(d_i) is the full ambient capacity.

The diagonal sector has an additional finite-ring Kakeya obstruction.
Every capacity-admissible diagonal row has D=7k.  A k-term unit-step
block in Z/DZ is a section of Z/DZ -> Z/kZ, so three such blocks meet
each fiber at most three times.  A support fiber of multiplicity at
least four therefore kills the diagonal denominator triple (D,D,D).

The saturated three-height rows have a second quotient.  They all satisfy
D=49m.  On a fixed residue modulo m, a full diagonal block restricts to a
seven-point unit-step arithmetic progression in Z/49Z.  We again relax by
allowing the three local needles to vary independently from slice to slice.
There are 1,029 distinct local needles.  An exact depth-three set-cover
test finds a noncoverable slice in every saturated row, emptying the full
diagonal denominator sector.

The body-support referee is hash-pinned.  Support projections are checked
by exact merged cyclic arcs.  Residue-load histograms are computed by an
event sweep, not an O(d)-sized array.  All checks use explicit RuntimeError
guards and remain active under optimized Python.
"""

from bisect import bisect_right
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SUPPORT_PATH = HERE / "lrc14_two_drift_body_projection_support_thm2928.py"
SUPPORT_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_two_drift_body_projection_support_thm2928.out"
)
EXPECTED_SUPPORT_SHA256 = (
    "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a"
)
EXPECTED_SUPPORT_OUTPUT_SHA256 = (
    "648327d3b9b5b9a50c7760f0afd89a7a33161f57fa98c1b9e181d6b5b791a25f"
)
EXPECTED_SUPPORT_SEMANTIC_SHA256 = (
    "e5d991029cb066546755681d3fe312c1e1b6521d9f1f06d5a3ab8214443761f2"
)
EXPECTED_ALL_TOP_SEMANTIC_SHA256 = (
    "8dfccbb9277486621ec4dde48616be63aa32faee3e33ffee4a72582bf8e02a95"
)
EXPECTED_DIAGONAL_SEMANTIC_SHA256 = (
    "58d9b576035d0be020510318654a517e9f1e0400a0c8c2fef7415324c7ef4315"
)
EXPECTED_LOCAL_WITNESS_SEMANTIC_SHA256 = (
    "b91a9e93724d94f2d9eaa8cda4ef1e4a5a21441000b4e167f44a2e06303a03ab"
)

FOUR_SAFE_FLOOR = Q(558, 1183)
THREE_UNION_CAP = Q(36, 91)
SUPPORT_CUTOFF = THREE_UNION_CAP / FOUR_SAFE_FLOOR


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(SUPPORT_PATH) == EXPECTED_SUPPORT_SHA256,
    "support-script dependency hash changed",
)
require(
    file_sha256(SUPPORT_OUTPUT_PATH) == EXPECTED_SUPPORT_OUTPUT_SHA256,
    "support-output dependency hash changed",
)
require(SUPPORT_CUTOFF == Q(26, 31), "three-drift support cutoff changed")
spec = spec_from_file_location("lrc14_body_projection_support", SUPPORT_PATH)
support_module = module_from_spec(spec)
spec.loader.exec_module(support_module)


def projected_support_arcs(D, ranges):
    """Merge the image modulo D of half-open integer ranges."""
    pieces = []
    for left, right in ranges:
        length = right - left
        if length >= D:
            return ((0, D),)
        residue = left % D
        endpoint = residue + length
        if endpoint <= D:
            pieces.append((residue, endpoint))
        else:
            pieces.append((residue, D))
            pieces.append((0, endpoint - D))
    pieces.sort()
    merged = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return tuple(merged)


def residue_load_histogram(arcs, d):
    """Histogram of lambda_d via a cyclic interval-event sweep."""
    common = 0
    events = defaultdict(int)
    events[0] = 0
    events[d] = 0

    def add_interval(left, right):
        if left < right:
            events[left] += 1
            events[right] -= 1

    for left, right in arcs:
        quotient, remainder = divmod(right - left, d)
        common += quotient
        start = left % d
        endpoint = start + remainder
        if endpoint <= d:
            add_interval(start, endpoint)
        else:
            add_interval(start, d)
            add_interval(0, endpoint - d)

    histogram = defaultdict(int)
    points = sorted(events)
    running = common
    for left, right in zip(points, points[1:]):
        running += events[left]
        histogram[running] += right - left
    result = tuple(sorted(histogram.items()))
    require(sum(count for _load, count in result) == d, "class count changed")
    require(
        sum(load * count for load, count in result)
        == sum(right - left for left, right in arcs),
        "class load lost support mass",
    )
    return result


def top_class_load(histogram, width):
    """Sum the width largest loads represented by a histogram."""
    remaining = width
    result = 0
    for load, count in reversed(histogram):
        take = min(remaining, count)
        result += take * load
        remaining -= take
        if remaining == 0:
            break
    require(remaining == 0, "top-class width exceeded denominator")
    return result


def slice_support_mask(arcs, modulus, residue):
    """Support heights in residue+modulus*(Z/49Z), as a 49-bit word."""
    result = 0
    arc_index = 0
    for height in range(49):
        point = residue + height * modulus
        while arc_index < len(arcs) and arcs[arc_index][1] <= point:
            arc_index += 1
        if (
            arc_index < len(arcs)
            and arcs[arc_index][0] <= point < arcs[arc_index][1]
        ):
            result |= 1 << height
    return result


def density_text(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    support_hash = sha256()
    by_divisor = {}
    body_count = 0
    divisor_rows = 0
    support_killed_rows = 0
    support_hard_rows = 0
    maximum_hard_divisors_per_body = 0
    minimum_hard_D = None
    minimum_hard_rows = []

    for F in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support_module.safe_cell_ranges(F)
        hard_this_body = 0
        for D in support_module.divisors(L):
            divisor_rows += 1
            support_count = support_module.support_size_bitset(D, ranges)
            support_hash.update(
                (
                    f"{','.join(map(str, F))}|{L}|{D}|{support_count}\n"
                ).encode()
            )
            density = Q(support_count, D)
            if density > SUPPORT_CUTOFF:
                support_killed_rows += 1
                continue
            support_hard_rows += 1
            hard_this_body += 1
            record = (support_count, F, L, tuple(ranges))
            by_divisor.setdefault(D, []).append(record)
            if minimum_hard_D is None or D < minimum_hard_D:
                minimum_hard_D = D
                minimum_hard_rows = [
                    (F, L, D, support_count, density)
                ]
            elif D == minimum_hard_D:
                minimum_hard_rows.append(
                    (F, L, D, support_count, density)
                )
        maximum_hard_divisors_per_body = max(
            maximum_hard_divisors_per_body,
            hard_this_body,
        )

    require(body_count == 3003, "body universe changed")
    require(divisor_rows == 251536, "body/divisor universe changed")
    require(support_killed_rows == 237758, "support kill count changed")
    require(support_hard_rows == 13778, "support frontier changed")
    require(len(by_divisor) == 206, "support divisor alphabet changed")
    require(
        maximum_hard_divisors_per_body == 8,
        "per-body support frontier changed",
    )
    require(
        minimum_hard_rows
        == [
            (
                (1, 2, 3, 4, 6, 12),
                168,
                28,
                22,
                Q(11, 14),
            )
        ],
        "minimum support-hard row changed",
    )
    require(
        support_hash.hexdigest() == EXPECTED_SUPPORT_SEMANTIC_SHA256,
        "body-support semantic ledger changed",
    )

    denominator_triple_shapes = 0
    raw_triple_occurrences = 0
    stage_occurrences = [0, 0, 0, 0]
    stage_rows = [0, 0, 0, 0]
    stage_shapes = [0, 0, 0, 0]
    all_top_hash = sha256()

    diagonal_candidates = 0
    diagonal_kills = 0
    diagonal_survivors = []
    diagonal_survivor_arcs = []
    diagonal_maximum_histogram = Counter()
    diagonal_bodies = set()
    diagonal_divisors = set()

    for D in sorted(by_divisor):
        divisors = [
            divisor
            for divisor in support_module.divisors(D)
            if divisor > 1
        ]
        capacities = [
            (D // divisor) * ((divisor + 6) // 7)
            for divisor in divisors
        ]
        rows = sorted(
            by_divisor[D],
            key=lambda record: (record[0], record[1], record[2]),
        )
        supports = [record[0] for record in rows]
        arcs_by_row = []
        top_loads_by_row = []
        for support_count, F, _L, ranges in rows:
            arcs = projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support arc mismatch", F, D),
            )
            arcs_by_row.append(arcs)
            top_loads = []
            for divisor in divisors:
                width = (divisor + 6) // 7
                histogram = residue_load_histogram(arcs, divisor)
                top_loads.append(top_class_load(histogram, width))
            top_loads_by_row.append(tuple(top_loads))

        row_flags = [[False] * len(rows) for _stage in range(4)]
        divisor_shape_count = 0
        for first_index, first in enumerate(divisors):
            for second_index in range(first_index, len(divisors)):
                second = divisors[second_index]
                first_second_lcm = lcm(first, second)
                for third_index in range(second_index, len(divisors)):
                    third = divisors[third_index]
                    if lcm(first_second_lcm, third) != D:
                        continue
                    divisor_shape_count += 1
                    denominator_triple_shapes += 1
                    full_capacity = (
                        capacities[first_index]
                        + capacities[second_index]
                        + capacities[third_index]
                    )
                    capacity_row_count = bisect_right(
                        supports,
                        full_capacity,
                    )
                    raw_triple_occurrences += len(rows)
                    if capacity_row_count == 0:
                        continue
                    stage_shapes[0] += 1
                    stage_occurrences[0] += capacity_row_count
                    any_one_top = False
                    any_two_top = False
                    any_all_top = False
                    all_top_row_mask = 0
                    for row_index in range(capacity_row_count):
                        row_flags[0][row_index] = True
                        support_count = supports[row_index]
                        top_loads = top_loads_by_row[row_index]
                        if support_count > (
                            top_loads[first_index]
                            + capacities[second_index]
                            + capacities[third_index]
                        ):
                            continue
                        stage_occurrences[1] += 1
                        row_flags[1][row_index] = True
                        any_one_top = True
                        if support_count > (
                            top_loads[first_index]
                            + top_loads[second_index]
                            + capacities[third_index]
                        ):
                            continue
                        stage_occurrences[2] += 1
                        row_flags[2][row_index] = True
                        any_two_top = True
                        if support_count > (
                            top_loads[first_index]
                            + top_loads[second_index]
                            + top_loads[third_index]
                        ):
                            continue
                        stage_occurrences[3] += 1
                        row_flags[3][row_index] = True
                        any_all_top = True
                        all_top_row_mask |= 1 << row_index
                    stage_shapes[1] += any_one_top
                    stage_shapes[2] += any_two_top
                    stage_shapes[3] += any_all_top
                    if all_top_row_mask:
                        mask_bytes = all_top_row_mask.to_bytes(
                            (len(rows) + 7) // 8,
                            "little",
                        )
                        all_top_hash.update(
                            (
                                f"{D}|{first}|{second}|{third}|"
                            ).encode()
                        )
                        all_top_hash.update(mask_bytes)

        for stage in range(4):
            stage_rows[stage] += sum(row_flags[stage])

        diagonal_width = (D + 6) // 7
        diagonal_capacity = 3 * diagonal_width
        diagonal_row_count = bisect_right(supports, diagonal_capacity)
        for row_index in range(diagonal_row_count):
            support_count, F, L, _ranges = rows[row_index]
            diagonal_candidates += 1
            diagonal_bodies.add(F)
            diagonal_divisors.add(D)
            require(
                D == 7 * diagonal_width,
                ("capacity-admissible diagonal is not septimal", F, D),
            )
            sparse_histogram = residue_load_histogram(
                arcs_by_row[row_index],
                diagonal_width,
            )
            maximum = max(
                load for load, count in sparse_histogram if count
            )
            sparse_dictionary = dict(sparse_histogram)
            histogram = tuple(
                (load, sparse_dictionary.get(load, 0))
                for load in range(maximum + 1)
            )
            diagonal_maximum_histogram[maximum] += 1
            diagonal_record = (
                F,
                L,
                D,
                support_count,
                Q(support_count, D),
                diagonal_width,
                histogram,
            )
            if maximum > 3:
                diagonal_kills += 1
                continue
            require(maximum == 3, "unexpected low diagonal fiber maximum")
            require(L == D, "surviving diagonal is not canonical-ruler")
            require(7 in F and 14 not in F, "surviving body lost speed seven")
            require(
                all(
                    load in (0, 2, 3) or count == 0
                    for load, count in histogram
                ),
                "surviving fiber alphabet changed",
            )
            diagonal_survivors.append(diagonal_record)
            diagonal_survivor_arcs.append(arcs_by_row[row_index])

    require(
        denominator_triple_shapes == 7483350,
        "denominator triple-shape universe changed",
    )
    require(
        raw_triple_occurrences == 298255882,
        "raw triple occurrence universe changed",
    )
    require(
        stage_occurrences
        == [143852683, 44573157, 1385991, 544571],
        "nested top-class occurrence ledger changed",
    )
    require(
        stage_rows == [13697, 13590, 13577, 13577],
        "nested top-class row ledger changed",
    )
    require(
        stage_shapes == [7313788, 6629229, 79451, 36614],
        "nested top-class shape ledger changed",
    )
    require(
        all_top_hash.hexdigest() == EXPECTED_ALL_TOP_SEMANTIC_SHA256,
        "all-top semantic ledger changed",
    )
    require(diagonal_candidates == 2636, "diagonal universe changed")
    require(diagonal_kills == 2601, "diagonal fiber kill count changed")
    require(len(diagonal_survivors) == 35, "diagonal residual count changed")
    require(len(diagonal_bodies) == 2617, "diagonal body count changed")
    require(len(diagonal_divisors) == 109, "diagonal divisor count changed")
    require(
        diagonal_maximum_histogram
        == Counter({6: 1496, 4: 702, 5: 403, 3: 35}),
        "diagonal maximum-fiber histogram changed",
    )

    diagonal_hash = sha256()
    for record in diagonal_survivors:
        diagonal_hash.update(f"{record}\n".encode())
    require(
        diagonal_hash.hexdigest() == EXPECTED_DIAGONAL_SEMANTIC_SHA256,
        "diagonal residual semantic ledger changed",
    )

    local_needles = tuple(
        sorted(
            {
                sum(
                    1 << ((start + index * step) % 49)
                    for index in range(7)
                )
                for start in range(49)
                for step in range(1, 49)
                if gcd(step, 49) == 1
            }
        )
    )
    require(len(local_needles) == 1029, "local needle alphabet changed")
    needles_through = [[] for _residue in range(49)]
    for needle in local_needles:
        for residue in range(49):
            if (needle >> residue) & 1:
                needles_through[residue].append(needle)
    require(
        {len(bank) for bank in needles_through} == {147},
        "local needle point degree changed",
    )

    @lru_cache(maxsize=None)
    def one_needle_coverable(target):
        return any(target & ~needle == 0 for needle in local_needles)

    @lru_cache(maxsize=None)
    def two_needle_coverable(target):
        if target == 0:
            return True
        if target.bit_count() > 14:
            return False
        first_point = (target & -target).bit_length() - 1
        return any(
            one_needle_coverable(target & ~needle)
            for needle in needles_through[first_point]
        )

    @lru_cache(maxsize=None)
    def three_needle_coverable(target):
        if target == 0:
            return True
        if target.bit_count() > 21:
            return False
        first_point = (target & -target).bit_length() - 1
        return any(
            two_needle_coverable(target & ~needle)
            for needle in needles_through[first_point]
        )

    local_witnesses = []
    local_witness_size_histogram = Counter()
    for record, arcs in zip(diagonal_survivors, diagonal_survivor_arcs):
        F, _L, D, _support, _density, _width, _histogram = record
        require(D % 49 == 0, ("saturated diagonal is not 49-adic", F, D))
        slice_modulus = D // 49
        witness = None
        for residue in range(slice_modulus):
            target = slice_support_mask(arcs, slice_modulus, residue)
            if not three_needle_coverable(target):
                heights = tuple(
                    height
                    for height in range(49)
                    if (target >> height) & 1
                )
                witness = (
                    F,
                    D,
                    slice_modulus,
                    residue,
                    len(heights),
                    heights,
                    residue + 1,
                )
                break
        require(witness is not None, ("locally coverable diagonal row", F, D))
        local_witnesses.append(witness)
        local_witness_size_histogram[witness[4]] += 1
    require(len(local_witnesses) == 35, "local needle kill count changed")
    require(
        all(14 <= witness[4] <= 18 for witness in local_witnesses),
        "local hostile target-size range changed",
    )
    require(
        max(witness[3] for witness in local_witnesses) == 792,
        "latest first hostile slice changed",
    )
    require(
        local_witness_size_histogram
        == Counter({16: 18, 15: 8, 18: 5, 14: 4}),
        "local hostile target-size histogram changed",
    )
    local_witness_hash = sha256()
    for witness in local_witnesses:
        local_witness_hash.update(f"{witness}\n".encode())
    require(
        local_witness_hash.hexdigest()
        == EXPECTED_LOCAL_WITNESS_SEMANTIC_SHA256,
        "local witness semantic ledger changed",
    )
    require(
        (
            one_needle_coverable.cache_info().hits,
            one_needle_coverable.cache_info().misses,
            two_needle_coverable.cache_info().hits,
            two_needle_coverable.cache_info().misses,
            three_needle_coverable.cache_info().hits,
            three_needle_coverable.cache_info().misses,
        )
        == (371639, 76417, 4089, 3302, 6313, 71),
        "local set-cover traversal changed",
    )

    print("LRC14 four-aligned/three-drift body-projection exact referee")
    print(f"support_script_sha256={file_sha256(SUPPORT_PATH)}")
    print(f"support_output_sha256={file_sha256(SUPPORT_OUTPUT_PATH)}")
    print(f"four_comb_safe_floor={density_text(FOUR_SAFE_FLOOR)}")
    print(f"three_comb_union_cap={density_text(THREE_UNION_CAP)}")
    print(f"support_cutoff={density_text(SUPPORT_CUTOFF)}")
    print(f"body_count={body_count}")
    print(f"divisor_rows={divisor_rows}")
    print(f"support_killed_rows={support_killed_rows}")
    print(f"support_hard_rows={support_hard_rows}")
    print(f"support_hard_divisors={len(by_divisor)}")
    print(
        "maximum_hard_divisors_per_body="
        f"{maximum_hard_divisors_per_body}"
    )
    print(f"minimum_hard_D={minimum_hard_D}")
    print(f"minimum_hard_rows={minimum_hard_rows}")
    print(f"support_semantic_sha256={support_hash.hexdigest()}")
    print(f"denominator_triple_shapes={denominator_triple_shapes}")
    print(f"raw_triple_occurrences={raw_triple_occurrences}")
    print(
        "nested_relaxations="
        "[C1+C2+C3,Lambda1+C2+C3,"
        "Lambda1+Lambda2+C3,Lambda1+Lambda2+Lambda3]"
    )
    print(f"stage_occurrences={stage_occurrences}")
    print(f"stage_rows={stage_rows}")
    print(f"stage_shapes={stage_shapes}")
    print(f"all_top_semantic_sha256={all_top_hash.hexdigest()}")
    print(f"diagonal_candidates={diagonal_candidates}")
    print(f"diagonal_candidate_bodies={len(diagonal_bodies)}")
    print(f"diagonal_candidate_divisors={len(diagonal_divisors)}")
    print(f"diagonal_maximum_histogram={diagonal_maximum_histogram}")
    print(f"diagonal_fiber_kills={diagonal_kills}")
    print(f"diagonal_fiber_survivors={len(diagonal_survivors)}")
    print(f"diagonal_survivor_first={diagonal_survivors[0]}")
    print(f"diagonal_survivor_last={diagonal_survivors[-1]}")
    print(
        "diagonal_structure=all survivors have D=L_F, 7 in F, "
        "14 not in F, and nonempty fiber sizes in {2,3}"
    )
    print(f"diagonal_semantic_sha256={diagonal_hash.hexdigest()}")
    print(f"local_needle_count={len(local_needles)}")
    print(f"local_needle_point_degree={len(needles_through[0])}")
    print(f"local_witness_size_histogram={local_witness_size_histogram}")
    print(f"local_witness_count={len(local_witnesses)}")
    print(f"local_witness_first={local_witnesses[0]}")
    print(
        "local_witness_latest="
        f"{max(local_witnesses, key=lambda witness: witness[3])}"
    )
    print(f"local_witness_last={local_witnesses[-1]}")
    print(f"local_one_cover_cache={one_needle_coverable.cache_info()}")
    print(f"local_two_cover_cache={two_needle_coverable.cache_info()}")
    print(f"local_three_cover_cache={three_needle_coverable.cache_info()}")
    print(f"local_witness_semantic_sha256={local_witness_hash.hexdigest()}")
    print("local_needle_kills=35")
    print("final_diagonal_survivors=0")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
