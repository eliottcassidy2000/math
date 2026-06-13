#!/usr/bin/env python3
"""
Verification Runner
====================
Exhaustively or randomly verify Claims A, B, and Redei's theorem.

Usage:
    python3 03-artifacts/code/verify.py --self-test
    python3 03-artifacts/code/verify.py --redei --n 5
    python3 03-artifacts/code/verify.py --claim-a --n 4
    python3 03-artifacts/code/verify.py --claim-a --n 5
    python3 03-artifacts/code/verify.py --claim-a --n 6          # ~15-30 min
    python3 03-artifacts/code/verify.py --claim-a --n 7 --sample 500
    python3 03-artifacts/code/verify.py --claim-b --n 5
    python3 03-artifacts/code/verify.py --all --n 5              # Claims A, B, and Redei
"""

import sys
import os
import time
import argparse
import random

# Allow running from repo root or from code directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tournament_lib import (
    all_tournaments, random_tournament, hamiltonian_path_count,
    verify_claim_a, verify_claim_b, verify_redei, self_test,
    find_odd_cycles, delete_vertex,
)


def fmt_time(seconds):
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"


def verify_exhaustive(n, claim, verbose=False):
    """Exhaustively verify a claim for all tournaments on n vertices."""
    m = n * (n - 1) // 2
    total_tournaments = 1 << m
    total_pairs = total_tournaments * n

    print(f"\nVerifying {claim} for ALL tournaments on n={n} vertices")
    print(f"Tournaments: {total_tournaments:,} | Pairs (T,v): {total_pairs:,}")
    print()

    verify_fn = verify_claim_a if claim == "Claim A" else verify_claim_b
    failures = 0
    pairs = 0
    t0 = time.time()
    last_report = t0

    for T in all_tournaments(n):
        for v in range(n):
            ok, lhs, rhs = verify_fn(T, v)
            pairs += 1
            if not ok:
                failures += 1
                if verbose or failures <= 10:
                    print(f"  FAILURE #{failures}: v={v}, LHS={lhs}, RHS={rhs}")

        # Progress report every 5 seconds
        now = time.time()
        if now - last_report > 5:
            elapsed = now - t0
            rate = pairs / elapsed if elapsed > 0 else 0
            eta = (total_pairs - pairs) / rate if rate > 0 else 0
            pct = 100 * pairs / total_pairs
            print(f"  [{pct:5.1f}%] {pairs:,}/{total_pairs:,} pairs | "
                  f"{rate:.0f} pairs/s | elapsed {fmt_time(elapsed)} | "
                  f"ETA {fmt_time(eta)} | failures: {failures}")
            last_report = now

    elapsed = time.time() - t0
    print(f"\nResult: {failures} failures out of {pairs:,} pairs. "
          f"({fmt_time(elapsed)})")
    if failures == 0:
        print(f"{claim} HOLDS for n={n}.")
    else:
        print(f"{claim} FAILS for n={n}: {failures} counterexamples found.")
    return failures


def verify_sampled(n, claim, sample_size, verbose=False):
    """Randomly verify a claim for sampled tournaments on n vertices."""
    print(f"\nVerifying {claim} for {sample_size} random tournaments on n={n}")
    print()

    verify_fn = verify_claim_a if claim == "Claim A" else verify_claim_b
    rng = random.Random(42)  # reproducible
    failures = 0
    pairs = 0
    t0 = time.time()

    for i in range(sample_size):
        T = random_tournament(n, rng)
        for v in range(n):
            ok, lhs, rhs = verify_fn(T, v)
            pairs += 1
            if not ok:
                failures += 1
                if verbose or failures <= 10:
                    print(f"  FAILURE #{failures}: tournament {i}, v={v}, "
                          f"LHS={lhs}, RHS={rhs}")

        if (i + 1) % max(1, sample_size // 20) == 0:
            elapsed = time.time() - t0
            print(f"  [{100*(i+1)/sample_size:5.1f}%] {i+1}/{sample_size} "
                  f"tournaments | failures: {failures} | {fmt_time(elapsed)}")

    elapsed = time.time() - t0
    print(f"\nResult: {failures} failures out of {pairs:,} pairs "
          f"({sample_size} tournaments). ({fmt_time(elapsed)})")
    if failures == 0:
        print(f"{claim} HOLDS for all sampled tournaments at n={n}.")
    else:
        print(f"{claim} FAILS: {failures} counterexamples found.")
    return failures


def verify_redei_all(n):
    """Verify Redei's theorem for all tournaments on n vertices."""
    m = n * (n - 1) // 2
    total = 1 << m
    print(f"\nVerifying Redei's theorem for all {total:,} tournaments on n={n}")

    failures = 0
    for T in all_tournaments(n):
        ok, h = verify_redei(T)
        if not ok:
            failures += 1
            if failures <= 5:
                print(f"  FAILURE: H(T) = {h} (even)")

    if failures == 0:
        print(f"Redei's theorem HOLDS for n={n}: all H(T) are odd.")
    else:
        print(f"Redei's theorem FAILS for n={n}: {failures} counterexamples.")
    return failures


def print_tournament_stats(n, sample=None):
    """Print some basic statistics about tournaments of size n."""
    print(f"\n--- Tournament statistics for n={n} ---")
    m = n * (n - 1) // 2
    total = 1 << m

    if sample is not None and total > sample:
        rng = random.Random(42)
        tournaments = [random_tournament(n, rng) for _ in range(sample)]
        label = f"(sampled {sample})"
    else:
        tournaments = list(all_tournaments(n))
        label = "(exhaustive)"

    h_values = [hamiltonian_path_count(T) for T in tournaments]
    cycle_counts = [len(find_odd_cycles(T)) for T in tournaments]

    print(f"  Count {label}: {len(tournaments)}")
    print(f"  H(T) range: [{min(h_values)}, {max(h_values)}]")
    print(f"  H(T) mean:  {sum(h_values)/len(h_values):.1f}")
    print(f"  Odd cycles range: [{min(cycle_counts)}, {max(cycle_counts)}]")
    print(f"  Odd cycles mean:  {sum(cycle_counts)/len(cycle_counts):.1f}")
    print(f"  All H(T) odd: {all(h % 2 == 1 for h in h_values)}")


def main():
    parser = argparse.ArgumentParser(description="Verify tournament claims")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--self-test", action="store_true",
                       help="Run quick self-tests")
    group.add_argument("--claim-a", action="store_true",
                       help="Verify Claim A")
    group.add_argument("--claim-b", action="store_true",
                       help="Verify Claim B")
    group.add_argument("--redei", action="store_true",
                       help="Verify Redei's theorem")
    group.add_argument("--all", action="store_true",
                       help="Verify Claims A, B, and Redei")
    group.add_argument("--stats", action="store_true",
                       help="Print tournament statistics")

    parser.add_argument("--n", type=int, help="Number of vertices")
    parser.add_argument("--sample", type=int, default=None,
                        help="Number of random tournaments to sample "
                             "(default: exhaustive)")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    if not args.n:
        parser.error("--n is required (except for --self-test)")

    if args.stats:
        print_tournament_stats(args.n, sample=args.sample)
        return

    if args.claim_a or args.all:
        claim = "Claim A"
        if args.sample:
            verify_sampled(args.n, claim, args.sample, args.verbose)
        else:
            verify_exhaustive(args.n, claim, args.verbose)

    if args.claim_b or args.all:
        claim = "Claim B"
        if args.sample:
            verify_sampled(args.n, claim, args.sample, args.verbose)
        else:
            verify_exhaustive(args.n, claim, args.verbose)

    if args.redei or args.all:
        verify_redei_all(args.n)


if __name__ == "__main__":
    main()
