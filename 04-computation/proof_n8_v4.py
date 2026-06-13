#!/usr/bin/env python3
"""
n=8 OCF proof v4: correctly count individual directed odd cycles.

A vertex set of size L can have up to (L-1)!/2 distinct directed cycles.
The conflict graph Omega has one vertex per directed cycle, not per vertex set.
Two cycles are adjacent iff they share a vertex.

At n=8, max indep set size = 2 (need 9+ verts for 3 VD odd cycles).

Instance: opus-2026-03-05-S4
"""

import time
import random
from itertools import combinations

N = 8
FULL = (1 << N) - 1


def ham_count(T):
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


def count_directed_cycles(T, combo):
    """Count # directed Hamiltonian cycles in T[combo].
    A directed cycle on L vertices: (L-1)! orderings / L rotations = (L-1)!/1
    But direction matters, so we fix start vertex and count directed paths back.
    Each cycle counted once (fix start, count paths returning)."""
    m = len(combo)
    if m < 3:
        return 0
    start = combo[0]
    others = list(combo[1:])
    mo = len(others)
    ofull = (1 << mo) - 1

    # DP: count paths from start through subsets of others
    dp = [[0] * mo for _ in range(1 << mo)]
    for i, o in enumerate(others):
        if T[start * N + o]:
            dp[1 << i][i] = 1

    for omask in range(1, ofull + 1):
        for li in range(mo):
            c = dp[omask][li]
            if c == 0:
                continue
            for ni in range(mo):
                if omask & (1 << ni):
                    continue
                if T[others[li] * N + others[ni]]:
                    dp[omask | (1 << ni)][ni] += c

    # Count paths that return to start
    total = 0
    for li in range(mo):
        if T[others[li] * N + start]:
            total += dp[ofull][li]

    return total


def compute_ocf(T):
    """Compute I(Omega(T), 2) correctly.

    At n=8, max independent set size in Omega is 2.
    I = 1 + 2*total_cycles + 4*VD_pairs

    For VD pairs: two cycles are VD if their vertex sets don't intersect.
    Cycles on the same vertex set are always adjacent (share all vertices).
    So VD pairs only come from different vertex sets.
    """
    # For each odd-size vertex subset, count directed cycles
    cycle_info = []  # (vmask, count)

    for length in range(3, N + 1, 2):
        for combo in combinations(range(N), length):
            vmask = 0
            for v in combo:
                vmask |= (1 << v)
            cnt = count_directed_cycles(T, combo)
            if cnt > 0:
                cycle_info.append((vmask, cnt))

    total_cycles = sum(cnt for _, cnt in cycle_info)

    # VD pairs: for cycles on vertex sets A and B with A∩B=∅
    # Number of VD cycle pairs = sum over VD set pairs (A,B) of count_A * count_B
    vd_pairs = 0
    nc = len(cycle_info)
    for i in range(nc):
        for j in range(i + 1, nc):
            if cycle_info[i][0] & cycle_info[j][0] == 0:
                vd_pairs += cycle_info[i][1] * cycle_info[j][1]

    return 1 + 2 * total_cycles + 4 * vd_pairs


def test_sample(sample_size=200):
    arc_pairs = [(a, b) for a in range(N) for b in range(a+1, N) if (a,b) != (0,1)]
    random.seed(42)
    print(f"Testing {sample_size} random n=8 tournaments...")

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
        it = compute_ocf(T)

        if ht == it:
            passes += 1
        else:
            if trial - passes <= 3:
                print(f"  FAIL trial {trial}: H={ht}, I={it}, diff={ht-it}")
            if trial - passes > 10:
                print(f"Stopping early: {passes}/{trial+1}")
                break

    elapsed = time.time() - t0
    rate = (trial + 1) / elapsed if elapsed > 0 else 0
    print(f"Result: {passes}/{trial+1} pass ({elapsed:.1f}s, {rate:.1f}/s)")
    est = (1 << 27) / rate if rate > 0 else float('inf')
    print(f"Estimated full: {est:.0f}s = {est/3600:.1f}hr")
    return passes == sample_size


def full_verification():
    arc_pairs = [(a, b) for a in range(N) for b in range(a+1, N) if (a,b) != (0,1)]
    total = 1 << 27
    chunk = 1 << 15

    print(f"=== n=8 FULL OCF PROOF ===")
    print(f"Total: {total:,}")

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
            it = compute_ocf(T)

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

    print(f"\nTOTAL: {checked:,}, {fails} fails, {time.time()-start:.0f}s")
    if fails == 0:
        print("*** PROVED: H(T) = I(Omega(T), 2) for ALL n=8 tournaments ***")


if __name__ == "__main__":
    import sys
    if "--full" in sys.argv:
        full_verification()
    else:
        test_sample(200)
