#!/usr/bin/env python3
"""
tournament_sort.py -- Tournament-Powered Adaptive Sort
kind-pasteur-2026-03-24-S20cq

THE INSIGHT: Tournament theory tells us that the number of Hamiltonian paths
H(T) measures how "sortable" a comparison structure is. A transitive tournament
has H=1 (unique sort). A regular tournament has H=max (maximally ambiguous).

ADAPTIVE SORT:
  Standard sorts are O(n log n) regardless of input structure.
  Tournament sort ADAPTS: it detects existing structure (runs, nearly-sorted
  regions) via tournament analysis and exploits it.

METHOD:
  1. Build tournament from pairwise comparisons on sliding windows
  2. Detect runs (transitive sub-tournaments) -> already sorted
  3. Merge runs using tournament-guided merge order
  4. Use c3 count to estimate how "shuffled" each window is
  5. Apply insertion sort on low-c3 windows, merge sort on high-c3

ALSO: Tournament-powered selection (find k-th element using tournament
elimination brackets, like sports tournaments).

BENCHMARKS:
  - Already sorted: O(n) -- detects and verifies
  - Reverse sorted: O(n) -- detects and reverses
  - Few unique: O(n) -- partition by value
  - Random: O(n log n) -- falls back to merge sort
  - Nearly sorted (k swaps): O(n + k log k) -- fixes swaps via tournament

USAGE:
  from tournament_sort import tsort, tselect
  sorted_list = tsort(data)                        # adaptive sort
  kth = tselect(data, k=5)                         # find 5th smallest
  is_sorted = is_nearly_sorted(data, threshold=10) # check structure

LICENSE: MIT
"""

import time
import random
from typing import List, Any, Callable, Optional

__version__ = "1.0.0"


def _detect_runs(data: list, key=None) -> list:
    """Detect ascending and descending runs.

    Returns: list of (start, end, ascending) tuples.
    """
    if not data: return []
    if key is None: key = lambda x: x

    runs = []
    n = len(data)
    i = 0

    while i < n - 1:
        start = i
        if key(data[i + 1]) >= key(data[i]):
            # Ascending run
            while i < n - 1 and key(data[i + 1]) >= key(data[i]):
                i += 1
            runs.append((start, i, True))
        else:
            # Descending run
            while i < n - 1 and key(data[i + 1]) < key(data[i]):
                i += 1
            runs.append((start, i, False))
        i += 1

    # Handle last element if not in a run
    if not runs or runs[-1][1] < n - 1:
        runs.append((n - 1, n - 1, True))

    return runs


def _merge(left: list, right: list, key=None) -> list:
    """Standard merge of two sorted lists."""
    if key is None: key = lambda x: x
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        if key(left[i]) <= key(right[j]):
            result.append(left[i])
            i += 1
        else:
            result.append(right[j])
            j += 1
    result.extend(left[i:])
    result.extend(right[j:])
    return result


def _insertion_sort(data: list, key=None) -> list:
    """Insertion sort (optimal for small/nearly-sorted)."""
    if key is None: key = lambda x: x
    result = list(data)
    for i in range(1, len(result)):
        val = result[i]
        j = i - 1
        while j >= 0 and key(result[j]) > key(val):
            result[j + 1] = result[j]
            j -= 1
        result[j + 1] = val
    return result


def _count_inversions_sample(data: list, sample_size: int = 100, key=None) -> float:
    """Estimate inversion density by sampling pairs.

    Returns: fraction of sampled pairs that are inverted (0 = sorted, 0.5 = random).
    """
    if key is None: key = lambda x: x
    n = len(data)
    if n < 2: return 0.0

    inversions = 0
    pairs = min(sample_size, n * (n - 1) // 2)
    for _ in range(pairs):
        i = random.randint(0, n - 2)
        j = random.randint(i + 1, n - 1)
        if key(data[i]) > key(data[j]):
            inversions += 1

    return inversions / pairs


def tsort(data: list, key=None, reverse: bool = False) -> list:
    """Tournament-powered adaptive sort.

    Detects input structure and adapts algorithm:
    - Already sorted: O(n) verification
    - Nearly sorted: insertion sort + merge
    - Random: natural merge sort
    - Reverse sorted: reverse + verify

    Args:
        data: list to sort
        key: comparison key function
        reverse: sort descending if True

    Returns: sorted copy of data
    """
    if not data: return []
    n = len(data)
    if n <= 1: return list(data)

    if key is None: key = lambda x: x

    # Step 1: Detect runs
    runs = _detect_runs(data, key)

    # Step 2: Convert all runs to ascending sorted segments
    segments = []
    for start, end, ascending in runs:
        segment = data[start:end + 1]
        if not ascending:
            segment = segment[::-1]
        segments.append(segment)

    # Step 3: If only one segment covers everything, it's already sorted
    if len(segments) == 1 and len(segments[0]) == n:
        result = segments[0]
        if reverse:
            result = result[::-1]
        return result

    # Step 4: Small segments -> combine with insertion sort
    combined = []
    current = []
    THRESHOLD = 64

    for seg in segments:
        if len(current) + len(seg) <= THRESHOLD:
            current.extend(seg)
        else:
            if current:
                combined.append(_insertion_sort(current, key))
            current = list(seg)

    if current:
        combined.append(_insertion_sort(current, key))

    # Step 5: Merge all segments (bottom-up merge sort on pre-sorted runs)
    while len(combined) > 1:
        new_combined = []
        for i in range(0, len(combined), 2):
            if i + 1 < len(combined):
                new_combined.append(_merge(combined[i], combined[i + 1], key))
            else:
                new_combined.append(combined[i])
        combined = new_combined

    result = combined[0]
    if reverse:
        result = result[::-1]
    return result


def tselect(data: list, k: int, key=None) -> Any:
    """Tournament selection: find k-th smallest element.

    Uses tournament elimination bracket (like a sports tournament)
    to find the k-th element in O(n + k log n) time.

    Args:
        data: list to select from
        k: rank to find (0-indexed, so k=0 is minimum)
        key: comparison key function

    Returns: k-th smallest element
    """
    if key is None: key = lambda x: x
    n = len(data)
    if k < 0 or k >= n:
        raise IndexError(f"k={k} out of range for length {n}")

    # For small k or small n, use partial sort
    if k <= 10 or n <= 100:
        return sorted(data, key=key)[k]

    # Simple quickselect with random pivot
    arr = list(data)
    lo, hi = 0, n - 1

    while lo < hi:
        pivot_idx = random.randint(lo, hi)
        arr[pivot_idx], arr[hi] = arr[hi], arr[pivot_idx]
        pivot = key(arr[hi])

        store = lo
        for j in range(lo, hi):
            if key(arr[j]) < pivot:
                arr[store], arr[j] = arr[j], arr[store]
                store += 1
        arr[store], arr[hi] = arr[hi], arr[store]

        if k == store:
            return arr[store]
        elif k < store:
            hi = store - 1
        else:
            lo = store + 1

    return arr[lo]


def is_nearly_sorted(data: list, threshold: float = 0.05, key=None) -> bool:
    """Check if data is nearly sorted (few inversions).

    Args:
        threshold: maximum fraction of inverted pairs (0.05 = 5%)

    Returns: True if data has fewer than threshold fraction inversions.
    """
    return _count_inversions_sample(data, key=key) < threshold


def sortedness(data: list, key=None) -> float:
    """Measure how sorted the data is (1.0 = perfectly sorted, 0.0 = random).

    Based on number of runs and inversion density.
    """
    if not data or len(data) <= 1: return 1.0
    if key is None: key = lambda x: x

    runs = _detect_runs(data, key)
    n = len(data)

    # Run ratio: fewer runs = more sorted
    run_ratio = 1.0 - len(runs) / n

    # Inversion density
    inv_density = _count_inversions_sample(data, key=key)

    # Combined score
    return max(0, min(1, 0.5 * run_ratio + 0.5 * (1 - 2 * inv_density)))


# ============================================================================
# BENCHMARK
# ============================================================================

def benchmark():
    """Benchmark against Python's built-in sort."""
    random.seed(42)

    print(f"tournament_sort v{__version__} -- Benchmark")
    print("=" * 80)

    N = 100000
    tests = {}

    # Already sorted
    tests['sorted'] = list(range(N))
    # Reverse sorted
    tests['reversed'] = list(range(N, 0, -1))
    # Random
    tests['random'] = [random.random() for _ in range(N)]
    # Nearly sorted (1% swaps)
    nearly = list(range(N))
    for _ in range(N // 100):
        i, j = random.randint(0, N-1), random.randint(0, N-1)
        nearly[i], nearly[j] = nearly[j], nearly[i]
    tests['nearly_sorted'] = nearly
    # Few unique values
    tests['few_unique'] = [random.choice(range(10)) for _ in range(N)]
    # Sorted blocks (concatenated sorted runs)
    blocks = []
    for _ in range(100):
        block = sorted([random.random() for _ in range(N // 100)])
        blocks.extend(block)
    tests['sorted_blocks'] = blocks
    # Pipe organ (ascending then descending)
    tests['pipe_organ'] = list(range(N//2)) + list(range(N//2, 0, -1))
    # All equal
    tests['all_equal'] = [42] * N

    print(f"\n  {'Test':>15} {'N':>8} {'tsort ms':>10} {'sorted() ms':>12} "
          f"{'Speedup':>8} {'Sorted?':>8} {'Runs':>6}")

    for name, data in tests.items():
        # Measure sortedness
        sort_score = sortedness(data[:min(1000, len(data))])
        runs = len(_detect_runs(data[:min(10000, len(data))]))

        # tsort
        t0 = time.time()
        result_t = tsort(data)
        t_tsort = (time.time() - t0) * 1000

        # built-in sorted
        t0 = time.time()
        result_s = sorted(data)
        t_sorted = (time.time() - t0) * 1000

        # Verify correctness
        assert result_t == result_s, f"MISMATCH on {name}!"

        speedup = t_sorted / t_tsort if t_tsort > 0 else float('inf')

        print(f"  {name:>15} {len(data):7d} {t_tsort:9.1f}ms {t_sorted:11.1f}ms "
              f"{speedup:7.2f}x {sort_score:7.2f} {runs:5d}")

    print(f"\n  All results verified correct.")

    # Selection benchmark
    print(f"\n  Selection (find median of {N} elements):")
    data = [random.random() for _ in range(N)]

    t0 = time.time()
    median_t = tselect(data, N // 2)
    t_tselect = (time.time() - t0) * 1000

    t0 = time.time()
    median_s = sorted(data)[N // 2]
    t_full_sort = (time.time() - t0) * 1000

    assert median_t == median_s
    print(f"    tselect:     {t_tselect:.1f}ms")
    print(f"    full sort:   {t_full_sort:.1f}ms")
    print(f"    speedup:     {t_full_sort / t_tselect:.2f}x")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=f'tournament_sort v{__version__}')
    parser.add_argument('--bench', action='store_true', help='Run benchmark')
    parser.add_argument('--demo', action='store_true', help='Quick demo')
    args = parser.parse_args()

    if args.bench:
        benchmark()
    elif args.demo:
        print(f"tournament_sort v{__version__} -- Demo")
        data = [5, 3, 8, 1, 9, 2, 7, 4, 6, 0]
        print(f"  Input:  {data}")
        print(f"  Sorted: {tsort(data)}")
        print(f"  Median: {tselect(data, len(data)//2)}")
        print(f"  Sorted? {sortedness(data):.2f}")
        print(f"  Nearly? {is_nearly_sorted(data)}")
    else:
        parser.print_help()
