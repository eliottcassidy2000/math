#!/usr/bin/env python3
"""
Benchmark: Bitmask DP vs OCF for computing H(T)
=================================================
Compares hamiltonian_path_count(T) [bitmask DP, O(2^n * n^2)]
against hamiltonian_path_count_ocf(T) [odd-cycle / independence polynomial].

For each n in {5,6,7,8,9,10,11}, generates 10 random tournaments,
times both methods, and reports average times.

DP is skipped for n >= 12 or if a single call exceeds 30 seconds.
"""

import sys
import time
import random

sys.path.insert(0, 'C:/Users/Eliott/Documents/GitHub/math/03-artifacts/code')
from tournament_lib import random_tournament, hamiltonian_path_count, hamiltonian_path_count_ocf

NUM_TRIALS = 10
DP_TIMEOUT = 30.0  # seconds — skip DP for larger n if exceeded

def benchmark_single(T, method):
    """Time a single call. Returns (result, elapsed_seconds)."""
    t0 = time.perf_counter()
    result = method(T)
    t1 = time.perf_counter()
    return result, t1 - t0

def run_benchmarks():
    sizes = [5, 6, 7, 8, 9, 10, 11]
    rng = random.Random(42)  # reproducible

    # Column widths for table
    hdr_n    = "n"
    hdr_dp   = "DP avg (s)"
    hdr_ocf  = "OCF avg (s)"
    hdr_ratio= "Speedup"
    hdr_note = "Notes"

    rows = []
    dp_skipped = False  # once True, skip DP for all larger n

    for n in sizes:
        print(f"\n--- n = {n} ---")
        tournaments = [random_tournament(n, rng) for _ in range(NUM_TRIALS)]

        # --- OCF timing ---
        ocf_times = []
        ocf_results = []
        for i, T in enumerate(tournaments):
            res, elapsed = benchmark_single(T, hamiltonian_path_count_ocf)
            ocf_times.append(elapsed)
            ocf_results.append(res)
            print(f"  OCF  trial {i+1:2d}: H={res:>12d}  time={elapsed:.6f}s")
        ocf_avg = sum(ocf_times) / len(ocf_times)

        # --- DP timing (skip if previously timed out) ---
        if dp_skipped:
            dp_avg = None
            note = "DP skipped (previous n exceeded timeout)"
            print(f"  DP   skipped (previous n exceeded {DP_TIMEOUT}s timeout)")
        else:
            dp_times = []
            dp_results = []
            timed_out = False
            for i, T in enumerate(tournaments):
                res, elapsed = benchmark_single(T, hamiltonian_path_count)
                dp_times.append(elapsed)
                dp_results.append(res)
                print(f"  DP   trial {i+1:2d}: H={res:>12d}  time={elapsed:.6f}s")

                # Verify agreement
                if res != ocf_results[i]:
                    print(f"  *** MISMATCH at trial {i+1}: DP={res}, OCF={ocf_results[i]} ***")

                # Check timeout on this single call
                if elapsed > DP_TIMEOUT:
                    timed_out = True
                    print(f"  DP exceeded {DP_TIMEOUT}s — will skip DP for larger n")
                    break

            if timed_out:
                dp_avg = sum(dp_times) / len(dp_times)
                note = f"DP timed out (>{DP_TIMEOUT}s)"
                dp_skipped = True
            else:
                dp_avg = sum(dp_times) / len(dp_times)
                note = ""

        rows.append((n, dp_avg, ocf_avg, note))

    # --- Print summary table ---
    print("\n" + "=" * 80)
    print("BENCHMARK SUMMARY: Bitmask DP vs OCF  (H(T) = I(Omega(T), 2))")
    print(f"Trials per size: {NUM_TRIALS}, RNG seed: 42")
    print("=" * 80)
    print(f"{'n':>3s}  {'DP avg (s)':>12s}  {'OCF avg (s)':>12s}  {'Speedup':>10s}  {'Notes'}")
    print("-" * 80)

    for n, dp_avg, ocf_avg, note in rows:
        dp_str = f"{dp_avg:.6f}" if dp_avg is not None else "---"
        ocf_str = f"{ocf_avg:.6f}"
        if dp_avg is not None and ocf_avg > 0:
            ratio = dp_avg / ocf_avg
            ratio_str = f"{ratio:.2f}x"
        else:
            ratio_str = "---"
        print(f"{n:>3d}  {dp_str:>12s}  {ocf_str:>12s}  {ratio_str:>10s}  {note}")

    print("-" * 80)
    print("Speedup = DP_time / OCF_time  (>1 means OCF is faster)")
    print()


if __name__ == "__main__":
    run_benchmarks()
