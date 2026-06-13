#!/usr/bin/env python3
"""
Comprehensive OCF verification: H(T) = I(Omega(T), 2) for tournaments.

Exhaustive for n <= 6 (33,864 tournaments), random sampling for n >= 7.

Usage:
    python3 verify_ocf_sweep.py [--max-n N] [--seed S]
"""

from tournament_lib import *
import random
import time
import sys


def main():
    max_n = 10
    seed = 42
    for arg in sys.argv[1:]:
        if arg.startswith('--max-n='):
            max_n = int(arg.split('=')[1])
        elif arg.startswith('--seed='):
            seed = int(arg.split('=')[1])

    rng = random.Random(seed)
    print(f"OCF Verification Sweep (seed={seed}, max_n={max_n})")
    print("=" * 70)

    sample_sizes = {7: 5000, 8: 500, 9: 100, 10: 30, 11: 10, 12: 5}

    for n in range(3, max_n + 1):
        t0 = time.time()
        tested = failures = 0

        if n <= 6:
            for T in all_tournaments(n):
                ok, h, i = verify_ocf(T)
                tested += 1
                if not ok:
                    failures += 1
                    print(f"  FAILURE at n={n}: H={h}, I={i}")
            mode = "exhaustive"
        else:
            count = sample_sizes.get(n, 5)
            for _ in range(count):
                T = random_tournament(n, rng)
                ok, h, i = verify_ocf(T)
                tested += 1
                if not ok:
                    failures += 1
                    print(f"  FAILURE at n={n}: H={h}, I={i}")
            mode = f"random({count})"

        elapsed = time.time() - t0
        status = "PASS" if failures == 0 else f"FAIL({failures})"
        print(f"  n={n:2d}: {status} | {tested:6d} tested | {mode:15s} | {elapsed:.1f}s")

    print("=" * 70)
    print("Done.")


if __name__ == "__main__":
    main()
