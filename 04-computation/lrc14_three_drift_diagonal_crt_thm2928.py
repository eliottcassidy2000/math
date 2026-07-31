#!/usr/bin/env python3
"""Exact diagonal referee for THM-2928's three-drift address frontier.

Let F be a six-element subset of {1,...,14}, let L=14*lcm(F), and let J
be the exact safe 1/L-cell word of the literal body.  For D|L put
S_D=J mod D.  The four-aligned safe floor 558/1183 and the sharp union cap
36/91 force |S_D|/D<=26/31.  A diagonal denominator triple (D,D,D) then
has fixed-phase cardinality capacity 3*ceil(D/7).

This referee exhausts that exact universe and applies two nested fiber
obstructions.

First, every cardinality-admissible diagonal row has D=7k.  A full
k-term unit-step block in Z/DZ is a section modulo k, so three blocks
meet each k-fiber at most three times.

The 35 rows surviving that screen all have D=49m with gcd(m,49)=1.  In
one m-fiber, a full D/7=7m term unit-step block becomes a seven-point
unit-step arithmetic progression in Z/49Z.  There are 1,029 distinct
such "needles".  We grant the three global blocks independently chosen
needles in every m-fiber, which is a relaxation.  An exact depth-three
set-cover recursion finds a noncoverable support slice on every row.

The recursion precomputes all 54,587 subsets coverable by one needle.
For depth>1 it pivots on the least target bit and exhausts every needle
through that bit.  This is complete because some member of any cover
must contain the pivot.  All checks use RuntimeError and remain active
under optimized Python.
"""

from bisect import bisect_right
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
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
require(SUPPORT_CUTOFF == Q(26, 31), "support cutoff changed")
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
    """Return the exact histogram of S_D loads modulo d."""
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

    histogram = Counter()
    points = sorted(events)
    running = common
    for left, right in zip(points, points[1:]):
        running += events[left]
        histogram[running] += right - left
    result = tuple(sorted(histogram.items()))
    require(sum(count for _load, count in result) == d, "fiber count changed")
    require(
        sum(load * count for load, count in result)
        == sum(right - left for left, right in arcs),
        "fiber load lost support mass",
    )
    return result


def support_slice(arcs, starts, m, residue):
    """Return {y in Z/49Z: residue+m*y belongs to S_D} as a bit mask."""
    result = 0
    for height in range(49):
        point = residue + m * height
        index = bisect_right(starts, point) - 1
        if index >= 0 and point < arcs[index][1]:
            result |= 1 << height
    return result


def needle_universe():
    raw = []
    for phase in range(49):
        for step in range(49):
            if gcd(step, 49) != 1:
                continue
            raw.append(
                sum(1 << ((phase + index * step) % 49) for index in range(7))
            )
    needles = tuple(sorted(set(raw)))
    by_bit = []
    for bit in range(49):
        by_bit.append(tuple(mask for mask in needles if mask & (1 << bit)))
    one_coverable = set()
    for needle in needles:
        subset = needle
        while True:
            one_coverable.add(subset)
            if subset == 0:
                break
            subset = (subset - 1) & needle
    return tuple(raw), needles, tuple(by_bit), frozenset(one_coverable)


RAW_NEEDLES, NEEDLES, NEEDLES_BY_BIT, ONE_COVERABLE = needle_universe()


@lru_cache(maxsize=None)
def coverable(mask, depth):
    """Whether mask is contained in a union of at most depth needles."""
    if mask == 0:
        return True
    if mask in ONE_COVERABLE:
        return True
    if depth <= 1:
        return False
    pivot = (mask & -mask).bit_length() - 1
    for needle in NEEDLES_BY_BIT[pivot]:
        if coverable(mask & ~needle, depth - 1):
            return True
    return False


def mask_elements(mask):
    return tuple(bit for bit in range(49) if mask & (1 << bit))


def main():
    require(len(RAW_NEEDLES) == 2058, "raw needle count changed")
    require(len(NEEDLES) == 1029, "unique needle count changed")
    require(len(ONE_COVERABLE) == 54587, "one-needle submask count changed")
    require(
        all(mask.bit_count() == 7 for mask in NEEDLES),
        "needle cardinality changed",
    )
    require(
        all(len(NEEDLES_BY_BIT[bit]) == 147 for bit in range(49)),
        "needle point degree changed",
    )

    # Positive controls for recursion depth and a cardinality-negative control.
    control_one = NEEDLES[0]
    control_two = NEEDLES[0] | NEEDLES[len(NEEDLES) // 2]
    control_three = control_two | NEEDLES[-1]
    require(coverable(control_one, 1), "one-needle positive control failed")
    require(coverable(control_two, 2), "two-needle positive control failed")
    require(coverable(control_three, 3), "three-needle positive control failed")
    require(
        not coverable((1 << 49) - 1, 3),
        "full-circle negative control unexpectedly covered",
    )

    body_count = 0
    divisor_rows = 0
    support_hard_rows = 0
    support_killed_rows = 0
    diagonal_cardinality = 0
    first_fiber_kills = 0
    first_fiber_survivors = []

    for F in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support_module.safe_cell_ranges(F)
        for D in support_module.divisors(L):
            divisor_rows += 1
            support_count = support_module.support_size_bitset(D, ranges)
            density = Q(support_count, D)
            if density > SUPPORT_CUTOFF:
                support_killed_rows += 1
                continue
            support_hard_rows += 1
            width = (D + 6) // 7
            if support_count > 3 * width:
                continue
            diagonal_cardinality += 1
            require(D % 7 == 0, ("diagonal row not divisible by 7", F, D))
            require(width == D // 7, ("diagonal width changed", F, D))
            arcs = projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support arc mismatch", F, D),
            )
            first_histogram = residue_load_histogram(arcs, width)
            first_maximum = first_histogram[-1][0]
            record = (
                F,
                L,
                D,
                support_count,
                density,
                first_maximum,
                first_histogram,
                arcs,
            )
            if first_maximum > 3:
                first_fiber_kills += 1
            else:
                first_fiber_survivors.append(record)

    require(body_count == 3003, "body universe changed")
    require(divisor_rows == 251536, "body/divisor universe changed")
    require(support_killed_rows == 237758, "support kill count changed")
    require(support_hard_rows == 13778, "support frontier changed")
    require(diagonal_cardinality == 2636, "diagonal cardinality changed")
    require(first_fiber_kills == 2601, "first fiber kill count changed")
    require(len(first_fiber_survivors) == 35, "first fiber frontier changed")

    witness_hash = sha256()
    witness_records = []
    for (
        F,
        L,
        D,
        support_count,
        density,
        first_maximum,
        first_histogram,
        arcs,
    ) in first_fiber_survivors:
        require(D == L, ("first-fiber survivor is not full ruler", F, L, D))
        require(7 in F, ("first-fiber survivor body omits seven", F, D))
        require(D % 49 == 0, ("first-fiber survivor lacks 49", F, D))
        m = D // 49
        require(gcd(m, 49) == 1, ("CRT factors not coprime", F, D, m))
        require(
            set(load for load, count in first_histogram if count)
            <= {0, 2, 3},
            ("unexpected first-fiber alphabet", F, D, first_histogram),
        )
        starts = tuple(left for left, _right in arcs)
        smallest_witness = None
        for residue in range(m):
            target = support_slice(arcs, starts, m, residue)
            if not coverable(target, 3):
                smallest_witness = (residue, target)
                break
        require(
            smallest_witness is not None,
            ("three needles cover every relaxed CRT slice", F, D),
        )
        residue, target = smallest_witness
        elements = mask_elements(target)
        record = (
            F,
            D,
            m,
            support_count,
            density,
            residue,
            elements,
        )
        witness_records.append(record)
        witness_hash.update(
            (
                f"{','.join(map(str, F))}|{D}|{m}|{support_count}|"
                f"{density.numerator}/{density.denominator}|{residue}|"
                f"{','.join(map(str, elements))}\n"
            ).encode()
        )

    require(len(witness_records) == 35, "CRT witness count changed")
    require(
        len({(F, D) for F, D, *_rest in witness_records}) == 35,
        "duplicate CRT row",
    )

    print("THM-2928 literal three-drift diagonal CRT referee")
    print(f"support_cutoff={SUPPORT_CUTOFF.numerator}/{SUPPORT_CUTOFF.denominator}")
    print(f"body_count={body_count}")
    print(f"divisor_rows={divisor_rows}")
    print(f"support_killed_rows={support_killed_rows}")
    print(f"support_hard_rows={support_hard_rows}")
    print(f"diagonal_cardinality_survivors={diagonal_cardinality}")
    print(f"first_fiber_kills={first_fiber_kills}")
    print(f"first_fiber_survivors={len(first_fiber_survivors)}")
    print("first_fiber_survivors_all_D_eq_L=True")
    print("first_fiber_survivors_all_body_contains_7=True")
    print("first_fiber_survivors_all_D_eq_49m_coprime=True")
    print(f"raw_Z49_needles={len(RAW_NEEDLES)}")
    print(f"unique_Z49_needles={len(NEEDLES)}")
    print(f"one_needle_coverable_submasks={len(ONE_COVERABLE)}")
    print(f"needles_through_each_point={len(NEEDLES_BY_BIT[0])}")
    print(
        "recursion_controls="
        f"{int(coverable(control_one,1))},"
        f"{int(coverable(control_two,2))},"
        f"{int(coverable(control_three,3))},"
        f"{int(coverable((1<<49)-1,3))}"
    )
    print(f"second_CRT_slice_kills={len(witness_records)}")
    print("diagonal_survivors=0")
    print(f"witness_sha256={witness_hash.hexdigest()}")
    print("smallest_noncoverable_slice_ledger:")
    for F, D, m, support_count, density, residue, elements in witness_records:
        print(
            "F="
            + ",".join(map(str, F))
            + f";D={D};m={m};support={support_count};"
            + f"density={density.numerator}/{density.denominator};"
            + f"x={residue};size={len(elements)};"
            + "heights="
            + ",".join(map(str, elements))
        )


if __name__ == "__main__":
    main()
