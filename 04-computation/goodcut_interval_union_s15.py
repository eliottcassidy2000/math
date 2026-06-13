#!/usr/bin/env python3
"""
goodcut_interval_union_s15.py

opus-2026-05-29-S15

Good-cut buckets are an interval-cover gas on the cut path.

For an n-vertex base-path tiling, there are N = n-1 legal cuts.  A tile
(hi, lo) with hi >= lo+2 contributes the interval {lo+1, ..., hi}; these
are exactly the intervals of the N-cut path with length at least 2.

This script derives the bucket-count polynomial

    B_N(x) = sum_g #{tilings with g good cuts} x^g

from one-dimensional interval covers, verifies it against direct tiling
enumeration for small n, and prints larger exact counts using the recurrence.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache


def run_lengths(mask: int, length: int) -> list[int]:
    """Lengths of 1-runs in a bit mask of the given length."""
    runs: list[int] = []
    cur = 0
    for i in range(length):
        if mask & (1 << i):
            cur += 1
        elif cur:
            runs.append(cur)
            cur = 0
    if cur:
        runs.append(cur)
    return runs


@lru_cache(maxsize=None)
def run_cover_count(length: int) -> int:
    """Number of subsets of intervals of length >=2 covering a full run.

    The interval family on a run of length L has C(L, 2) members.  Inclusion-
    exclusion over uncovered cut positions gives

        c_L = sum_{A subset [L]} (-1)^|A| 2^{sum_R C(|R|, 2)},

    where R ranges over the 1-runs of [L] \\ A.
    """
    if length < 2:
        return 0
    total = 0
    full = (1 << length) - 1
    for missing in range(1 << length):
        available = full ^ missing
        interval_slots = sum(r * (r - 1) // 2 for r in run_lengths(available, length))
        term = 1 << interval_slots
        if missing.bit_count() & 1:
            total -= term
        else:
            total += term
    return total


def bucket_polynomial(num_cuts: int) -> Counter[int]:
    """Return B_N as a Counter degree -> coefficient."""
    polys: list[Counter[int]] = [Counter() for _ in range(num_cuts + 1)]
    polys[0][0] = 1

    for ncuts in range(1, num_cuts + 1):
        poly = Counter(polys[ncuts - 1])

        for run_len in range(2, ncuts + 1):
            cover = run_cover_count(run_len)
            if run_len == ncuts:
                tail = Counter({0: 1})
            else:
                # The position after the initial run must be absent, hence
                # the remaining suffix has length ncuts-run_len-1.
                tail = polys[ncuts - run_len - 1]
            for degree, coeff in tail.items():
                poly[degree + run_len] += cover * coeff

        polys[ncuts] = poly

    return polys[num_cuts]


def polynomial_from_good_sets(num_cuts: int) -> Counter[int]:
    """Independent check: sum component-cover weights over all good-cut sets."""
    counts: Counter[int] = Counter()
    for mask in range(1 << num_cuts):
        weight = 1
        for run_len in run_lengths(mask, num_cuts):
            c = run_cover_count(run_len)
            if c == 0:
                weight = 0
                break
            weight *= c
        if weight:
            counts[mask.bit_count()] += weight
    return counts


def tile_intervals(num_vertices: int) -> list[tuple[int, int]]:
    """All tile intervals as cut endpoints [lo+1, hi], 1-indexed."""
    intervals: list[tuple[int, int]] = []
    for lo in range(num_vertices):
        for hi in range(lo + 2, num_vertices):
            intervals.append((lo + 1, hi))
    return intervals


def brute_tiling_counts(num_vertices: int) -> Counter[int]:
    """Directly enumerate tilings and count good cuts."""
    intervals = tile_intervals(num_vertices)
    counts: Counter[int] = Counter()
    for bits in range(1 << len(intervals)):
        good: set[int] = set()
        for k, (start, end) in enumerate(intervals):
            if bits & (1 << k):
                good.update(range(start, end + 1))
        counts[len(good)] += 1
    return counts


def fmt_poly(poly: Counter[int]) -> str:
    return "{" + ", ".join(f"{k}: {poly[k]}" for k in sorted(poly)) + "}"


def verify_small(max_vertices: int = 8) -> None:
    print("SMALL-N VERIFICATION AGAINST DIRECT TILING ENUMERATION")
    for n in range(3, max_vertices + 1):
        formula = bucket_polynomial(n - 1)
        by_good_sets = polynomial_from_good_sets(n - 1)
        brute = brute_tiling_counts(n)
        ok = formula == by_good_sets == brute
        print(f"  n={n}: ok={ok}, counts={fmt_poly(formula)}")
        if not ok:
            print(f"    recurrence: {fmt_poly(formula)}")
            print(f"    good sets : {fmt_poly(by_good_sets)}")
            print(f"    brute     : {fmt_poly(brute)}")
            raise SystemExit(1)


def print_run_counts(max_run: int = 12) -> None:
    print("\nCONNECTED RUN COVER COUNTS c_L")
    for length in range(2, max_run + 1):
        slots = length * (length - 1) // 2
        print(f"  L={length:2d}: c_L={run_cover_count(length)}  (all interval subsets=2^{slots})")


def print_bucket_table(max_vertices: int = 13) -> None:
    print("\nBUCKET POLYNOMIALS FROM THE INTERVAL RECURRENCE")
    for n in range(3, max_vertices + 1):
        num_cuts = n - 1
        interval_count = num_cuts * (num_cuts - 1) // 2
        poly = bucket_polynomial(num_cuts)
        total = sum(poly.values())
        expected = 1 << interval_count
        top = poly[num_cuts]
        print(
            f"  n={n:2d}, cuts={num_cuts:2d}, tiles={interval_count:2d}, "
            f"total_ok={total == expected}, top={top}, counts={fmt_poly(poly)}"
        )


def main() -> None:
    print("GOOD-CUT INTERVAL-UNION ENUMERATION")
    print("opus-2026-05-29-S15")
    print("B_N = B_{N-1} + sum_{L=2..N} c_L x^L B_{N-L-1}, with B_{-1}=B_0=1")
    verify_small()
    print_run_counts()
    print_bucket_table()


if __name__ == "__main__":
    main()
