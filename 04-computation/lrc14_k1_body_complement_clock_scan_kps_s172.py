#!/usr/bin/env python3
"""Exact scout for small-LRC complement-clock terminals in the k=1 ledger.

For a body E with ruler L and divisor D, let S_D be its safe-address support.
The projected carrier's complement is contained in the aligned comb D_(D*a)
plus the open D-cells outside S_D.  A body speed e divisible by M=L/D descends
to the integer quotient clock c=e/M.  This scout asks whether at most five of
those descended body clocks cover every unsupported open D-cell.  If so, six
drift clocks, the aligned clock, and the complement clocks would globally
cover the circle with at most twelve distinct integer combs, contradicting the
cited LRC theorem through twelve speeds.
"""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, floor
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "04-computation/lrc14_two_drift_body_projection_support_thm2928.py"
EXPECTED_SUPPORT_SHA256 = "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a"
CUTOFF = Q(106, 117)
INPUT_OCCURRENCES = 38_954_725_590_760


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def file_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_support():
    require(file_hash(SUPPORT_PATH) == EXPECTED_SUPPORT_SHA256, "support dependency changed")
    spec = spec_from_file_location("kps_s172_support", SUPPORT_PATH)
    require(spec is not None and spec.loader is not None, "support import")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def all_divisors(number):
    low = []
    high = []
    d = 1
    while d * d <= number:
        if number % d == 0:
            low.append(d)
            if d * d != number:
                high.append(number // d)
        d += 1
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


@lru_cache(maxsize=None)
def denominator_shapes(D):
    return sum(
        mobius(D // divisor) * comb(len(all_divisors(divisor)) + 4, 6)
        for divisor in all_divisors(D)
    )


def residue_arcs(D, ranges):
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
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def unsupported_gaps(D, arcs):
    gaps = []
    cursor = 0
    for left, right in arcs:
        if cursor < left:
            gaps.append((Q(cursor, D), Q(left, D)))
        cursor = right
    if cursor < D:
        gaps.append((Q(cursor, D), Q(1)))
    return tuple(gaps)


def circular_distance(value):
    residue = value - floor(value)
    return min(residue, 1 - residue)


def danger(clock, x):
    return circular_distance(clock * x) < Q(1, 14)


def arrangement_points(clocks):
    points = {Q(0), Q(1)}
    for clock in clocks:
        for k in range(clock + 1):
            for sign in (-1, 1):
                x = Q(14 * k + sign, 14 * clock)
                if Q(0) <= x <= Q(1):
                    points.add(x)
    return tuple(sorted(points))


def safe_for(clocks, x):
    return all(not danger(clock, x) for clock in clocks)


@lru_cache(maxsize=None)
def safe_structure(clocks):
    points = arrangement_points(clocks)
    safe_points = tuple(x for x in points if safe_for(clocks, x))
    safe_atoms = tuple(
        (left, right)
        for left, right in zip(points, points[1:])
        if safe_for(clocks, (left + right) / 2)
    )
    return safe_points, safe_atoms


def covers_unsupported_cells(clocks, D, gaps):
    safe_points, safe_atoms = safe_structure(clocks)
    # Equality points of the open danger combs are safe.  A D-grid point is
    # nevertheless owned by the aligned clock D_(D*alpha); every other safe
    # arrangement endpoint must lie outside the unsupported open cells.
    for x in safe_points:
        if (D * x).denominator == 1:
            continue
        if any(left < x < right for left, right in gaps):
            return False
    # On each open arrangement atom the danger mask is constant.  Check the
    # whole atom against unsupported cell arcs, not only its midpoint's cell.
    for left, right in safe_atoms:
        if any(max(left, gap_left) < min(right, gap_right) for gap_left, gap_right in gaps):
            return False
    return True


@lru_cache(maxsize=None)
def least_completion(candidates, D, gaps):
    for size in range(min(5, len(candidates)) + 1):
        for clocks in combinations(candidates, size):
            if covers_unsupported_cells(clocks, D, gaps):
                return clocks
    return None


def main():
    support = load_support()
    rows = 0
    closed = []
    size_histogram = Counter()
    divisor_histogram = Counter()
    occurrence_size_histogram = Counter()
    semantic = sha256()
    input_occurrences = 0
    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            count = support.support_size_bitset(D, ranges)
            if Q(count, D) > CUTOFF:
                continue
            rows += 1
            input_occurrences += denominator_shapes(D)
            arcs = residue_arcs(D, ranges)
            require(sum(right - left for left, right in arcs) == count, (body, D, count, arcs))
            gaps = unsupported_gaps(D, arcs)
            M = L // D
            candidates = tuple(e // M for e in body if e % M == 0)
            require(len(candidates) == len(set(candidates)), (body, D, candidates))
            completion = least_completion(candidates, D, gaps)
            if completion is not None:
                row = (body, L, D, count, candidates, completion)
                closed.append(row)
                size_histogram[len(completion)] += 1
                divisor_histogram[D] += 1
                occurrence_size_histogram[len(completion)] += denominator_shapes(D)
                semantic.update((repr(row) + "\n").encode())

    require(rows == 27_240, rows)
    require(input_occurrences == INPUT_OCCURRENCES, input_occurrences)
    require(denominator_shapes(14) == 26, denominator_shapes(14))
    require(all(len(row[-1]) <= 5 for row in closed), "completion size")
    closed_occurrences = sum(denominator_shapes(row[2]) for row in closed)
    print("LRC14 K=1 BODY-COMPLEMENT CLOCK SCOUT")
    print("status=FINITE-EXACT support-level sufficient terminals; no physical cover claimed beyond listed rows")
    print(f"input_rows={rows};closed_rows={len(closed)};remaining_rows={rows-len(closed)}")
    print(f"completion_size_histogram={tuple(sorted(size_histogram.items()))}")
    print(f"closed_occurrences={closed_occurrences};remaining_occurrences={INPUT_OCCURRENCES-closed_occurrences}")
    print(f"occurrence_size_histogram={tuple(sorted(occurrence_size_histogram.items()))}")
    print(f"closed_divisor_histogram={tuple(sorted(divisor_histogram.items()))}")
    print(f"first_rows={tuple(closed[:20])}")
    print(f"last_rows={tuple(closed[-5:])}")
    print(f"support_dependency_sha256={file_hash(SUPPORT_PATH)}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
