#!/usr/bin/env python3
"""
Benchmark: Bitmask DP vs OCF for computing H(T)
=================================================
Compares hamiltonian_path_count(T) [bitmask DP, O(2^n * n^2)]
against hamiltonian_path_count_ocf(T) [odd-cycle / independence polynomial].

For each n in {5,6,7,8,9,10,11}, generates 10 random tournaments,
times both methods, and reports average time.

Either method is skipped for subsequent n values once a single call
exceeds 30 seconds.
"""

import sys
import time
import random

sys.path.insert(0, 'C:/Users/Eliott/Documents/GitHub/math/03-artifacts/code')
from tournament_lib import random_tournament, hamiltonian_path_count, hamiltonian_path_count_ocf

NUM_TRIALS = 10
TIMEOUT = 30.0  # seconds — skip method for larger n if a single call exceeds this

def benchmark_single(T, method):
    """Time a single call. Returns (result, elapsed_seconds)."""
    t0 = time.perf_counter()
    result = method(T)
    t1 = time.perf_counter()
    return result, t1 - t0

def run_benchmarks():
    sizes = [5, 6, 7, 8, 9, 10, 11]
    rng = random.Random(42)  # reproducible

    rows = []
    dp_skipped = False
    ocf_skipped = False

    for n in sizes:
        print(f"\n--- n = {n} ---")
        tournaments = [random_tournament(n, rng) for _ in range(NUM_TRIALS)]

        # --- DP timing ---
        dp_avg = None
        dp_results = [None] * NUM_TRIALS
        if dp_skipped:
            dp_note = "skipped"
            print(f"  DP   skipped (previous n exceeded {TIMEOUT}s)")
        else:
            dp_times = []
            timed_out = False
            for i, T in enumerate(tournaments):
                res, elapsed = benchmark_single(T, hamiltonian_path_count)
                dp_times.append(elapsed)
                dp_results[i] = res
                print(f"  DP   trial {i+1:2d}: H={res:>12d}  time={elapsed:.6f}s")
                if elapsed > TIMEOUT:
                    timed_out = True
                    print(f"  DP exceeded {TIMEOUT}s — will skip DP for larger n")
                    break
            dp_avg = sum(dp_times) / len(dp_times)
            dp_note = f"timed out (>{TIMEOUT}s)" if timed_out else ""
            if timed_out:
                dp_skipped = True

        # --- OCF timing ---
        ocf_avg = None
        ocf_results = [None] * NUM_TRIALS
        if ocf_skipped:
            ocf_note = "skipped"
            print(f"  OCF  skipped (previous n exceeded {TIMEOUT}s)")
        else:
            ocf_times = []
            timed_out = False
            for i, T in enumerate(tournaments):
                res, elapsed = benchmark_single(T, hamiltonian_path_count_ocf)
                ocf_times.append(elapsed)
                ocf_results[i] = res
                print(f"  OCF  trial {i+1:2d}: H={res:>12d}  time={elapsed:.6f}s")
                if elapsed > TIMEOUT:
                    timed_out = True
                    print(f"  OCF exceeded {TIMEOUT}s — will skip OCF for larger n")
                    break
            ocf_avg = sum(ocf_times) / len(ocf_times)
            ocf_note = f"timed out (>{TIMEOUT}s)" if timed_out else ""
            if timed_out:
                ocf_skipped = True

        # --- Verify agreement where both ran ---
        for i in range(NUM_TRIALS):
            if dp_results[i] is not None and ocf_results[i] is not None:
                if dp_results[i] != ocf_results[i]:
                    print(f"  *** MISMATCH trial {i+1}: DP={dp_results[i]}, OCF={ocf_results[i]} ***")

        # Build note
        notes = []
        if dp_note: notes.append(f"DP {dp_note}")
        if ocf_note: notes.append(f"OCF {ocf_note}")
        note = "; ".join(notes)

        rows.append((n, dp_avg, ocf_avg, note))

    # --- Print summary table ---
    print("\n" + "=" * 80)
    print("BENCHMARK SUMMARY: Bitmask DP vs OCF  (H(T) = I(Omega(T), 2))")
    print(f"Trials per size: {NUM_TRIALS}, RNG seed: 42")
    print("=" * 80)
    print(f"{'n':>3s}  {'DP avg (s)':>12s}  {'OCF avg (s)':>12s}  {'DP/OCF':>10s}  {'Notes'}")
    print("-" * 80)

    for n, dp_avg, ocf_avg, note in rows:
        dp_str = f"{dp_avg:.6f}" if dp_avg is not None else "---"
        ocf_str = f"{ocf_avg:.6f}" if ocf_avg is not None else "---"
        if dp_avg is not None and ocf_avg is not None and ocf_avg > 0:
            ratio = dp_avg / ocf_avg
            ratio_str = f"{ratio:.2f}x"
        else:
            ratio_str = "---"
        print(f"{n:>3d}  {dp_str:>12s}  {ocf_str:>12s}  {ratio_str:>10s}  {note}")

    print("-" * 80)
    print("DP/OCF > 1 means OCF is faster; < 1 means DP is faster.")
    print()


if __name__ == "__main__":
    run_benchmarks()
