#!/usr/bin/env python3
"""
n=8 OCF proof v3: direct H(T) = I(Omega(T), 2) check.

Optimized cycle finding + independence polynomial.
Key insight: at n=8, max independent set size in Omega is 2
(three mutually VD odd cycles need 9+ vertices).

So I(Omega, 2) = 1 + 2*c + 4*p where:
  c = # odd cycles
  p = # vertex-disjoint odd cycle PAIRS

This avoids full independent set enumeration.

Instance: opus-2026-03-05-S4
"""

import time
import random
from itertools import combinations

N = 8
FULL = (1 << N) - 1


def ham_count(T):
    """H(T) via bitmask DP. T is flat array T[a*N+b]."""
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << N):
        for last in range(N):
            c = dp[mask][last]
            if c == 0:
                continue
            for nxt in range(N):
                if mask & (1 << nxt):
                    continue
                if T[last * N + nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[FULL])


def find_all_odd_cycles_bitmask(T):
    """Find all directed odd cycles. Returns list of vertex bitmasks."""
    cycles = []
    for length in range(3, N + 1, 2):
        for combo in combinations(range(N), length):
            vmask = 0
            for v in combo:
                vmask |= (1 << v)

            # Check if there's a directed Hamiltonian cycle in T[combo]
            start = combo[0]
            others = list(combo[1:])
            m = len(others)
            if m == 0:
                continue

            # DP: paths from start through subsets of others
            dp = {}
            for i, o in enumerate(others):
                if T[start * N + o]:
                    dp[(1 << i, i)] = True

            ofull = (1 << m) - 1
            for omask in range(1, ofull + 1):
                for li in range(m):
                    if not (omask & (1 << li)):
                        continue
                    if (omask, li) not in dp:
                        continue
                    for ni in range(m):
                        if omask & (1 << ni):
                            continue
                        if T[others[li] * N + others[ni]]:
                            dp[(omask | (1 << ni), ni)] = True

            # Check if any full path returns to start
            has_cycle = False
            for li in range(m):
                if (ofull, li) in dp and T[others[li] * N + start]:
                    has_cycle = True
                    break

            if has_cycle:
                cycles.append(vmask)

    return cycles


def indpoly2_from_cycles(cycle_masks):
    """Compute I(Omega, 2) from cycle bitmasks.
    At n=8: max independent set size = 2.
    I = 1 + 2*c + 4*p where c=#cycles, p=#VD pairs."""
    c = len(cycle_masks)
    p = 0
    for i in range(c):
        for j in range(i + 1, c):
            if cycle_masks[i] & cycle_masks[j] == 0:  # vertex-disjoint
                p += 1
    return 1 + 2 * c + 4 * p


def test_sample(sample_size=100):
    arc_pairs = [(a, b) for a in range(N) for b in range(a+1, N) if (a,b) != (0,1)]
    random.seed(42)
    print(f"Testing {sample_size} random n=8 tournaments (direct OCF check)...")

    t0 = time.time()
    passes = 0
    for trial in range(sample_size):
        mask = random.randint(0, (1 << 27) - 1)
        T = [0] * (N * N)
        T[0 * N + 1] = 1
        for bit, (a, b) in enumerate(arc_pairs):
            if mask & (1 << bit):
                T[a * N + b] = 1
            else:
                T[b * N + a] = 1

        ht = ham_count(T)
        cycles = find_all_odd_cycles_bitmask(T)
        it = indpoly2_from_cycles(cycles)

        if ht == it:
            passes += 1
        else:
            print(f"  FAIL trial {trial}: H(T)={ht}, I(Omega,2)={it}, diff={ht-it}")
            if trial - passes > 5:
                print("Too many failures, stopping.")
                break

    elapsed = time.time() - t0
    rate = sample_size / elapsed if elapsed > 0 else 0
    print(f"Result: {passes}/{sample_size} pass ({elapsed:.1f}s, {rate:.1f}/s)")
    est = (1 << 27) / rate if rate > 0 else float('inf')
    print(f"Estimated full: {est:.0f}s = {est/3600:.1f}hr")
    return passes == sample_size


def full_verification():
    arc_pairs = [(a, b) for a in range(N) for b in range(a+1, N) if (a,b) != (0,1)]
    total = 1 << 27
    chunk = 1 << 15  # 32768

    print(f"=== n=8 FULL OCF PROOF: H(T) = I(Omega(T), 2) ===")
    print(f"Total: {total:,} configs")

    start = time.time()
    checked = 0
    fails = 0

    for cs in range(0, total, chunk):
        ce = min(cs + chunk, total)
        t0 = time.time()
        cf = 0

        for mask in range(cs, ce):
            T = [0] * (N * N)
            T[0 * N + 1] = 1
            for bit, (a, b) in enumerate(arc_pairs):
                if mask & (1 << bit):
                    T[a * N + b] = 1
                else:
                    T[b * N + a] = 1

            ht = ham_count(T)
            cycles = find_all_odd_cycles_bitmask(T)
            it = indpoly2_from_cycles(cycles)

            if ht != it:
                cf += 1
                if fails + cf <= 3:
                    print(f"  FAIL mask {mask}: H={ht}, I={it}")

        el = time.time() - t0
        checked += (ce - cs)
        fails += cf
        rate = (ce - cs) / el if el > 0 else 0
        eta = (total - checked) / rate if rate > 0 else 0
        print(f"  [{100*checked/total:5.1f}%] {checked:>12,}/{total:,} | "
              f"{cf} fails | {rate:.0f}/s | ETA {eta/3600:.2f}hr")

        if fails > 10:
            print("Aborting.")
            return

    print(f"\nTOTAL: {checked:,} checked, {fails} fails, {time.time()-start:.0f}s")
    if fails == 0:
        print("*** PROVED: H(T) = I(Omega(T), 2) for ALL n=8 tournaments ***")


if __name__ == "__main__":
    import sys
    if "--full" in sys.argv:
        full_verification()
    else:
        test_sample(100)
