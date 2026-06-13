#!/usr/bin/env python3
"""
n=8 OCF proof using multiprocessing.

Direct check: H(T) = I(Omega(T), 2) for all 2^27 tournaments.
Uses 4 worker processes.

Instance: opus-2026-03-05-S4
"""

import time
import multiprocessing as mp
from itertools import combinations

N = 8
FULL = (1 << N) - 1

ARC_PAIRS = [(a, b) for a in range(N) for b in range(a+1, N) if (a,b) != (0,1)]
assert len(ARC_PAIRS) == 27


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
    m = len(combo)
    if m < 3:
        return 0
    start = combo[0]
    others = list(combo[1:])
    mo = len(others)
    ofull = (1 << mo) - 1
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
    total = 0
    for li in range(mo):
        if T[others[li] * N + start]:
            total += dp[ofull][li]
    return total


def compute_ocf(T):
    cycle_info = []
    for length in range(3, N + 1, 2):
        for combo in combinations(range(N), length):
            vmask = 0
            for v in combo:
                vmask |= (1 << v)
            cnt = count_directed_cycles(T, combo)
            if cnt > 0:
                cycle_info.append((vmask, cnt))

    total_cycles = sum(cnt for _, cnt in cycle_info)
    vd_pairs = 0
    nc = len(cycle_info)
    for i in range(nc):
        for j in range(i + 1, nc):
            if cycle_info[i][0] & cycle_info[j][0] == 0:
                vd_pairs += cycle_info[i][1] * cycle_info[j][1]

    return 1 + 2 * total_cycles + 4 * vd_pairs


def verify_chunk(args):
    start_mask, end_mask = args
    fails = 0
    fail_examples = []
    for mask in range(start_mask, end_mask):
        T = [0] * (N * N)
        T[0 * N + 1] = 1
        for bit, (a, b) in enumerate(ARC_PAIRS):
            if mask & (1 << bit):
                T[a * N + b] = 1
            else:
                T[b * N + a] = 1

        ht = ham_count(T)
        it = compute_ocf(T)
        if ht != it:
            fails += 1
            if len(fail_examples) < 3:
                fail_examples.append((mask, ht, it))

    return end_mask - start_mask, fails, fail_examples


def main():
    total = 1 << 27
    n_workers = 4
    chunk_size = 1 << 15  # 32768 per chunk

    chunks = [(s, min(s + chunk_size, total)) for s in range(0, total, chunk_size)]
    n_chunks = len(chunks)

    print(f"=== n=8 FULL OCF PROOF (multiprocessing) ===")
    print(f"Total: {total:,}, Workers: {n_workers}, Chunks: {n_chunks}")

    start = time.time()
    checked = 0
    total_fails = 0

    with mp.Pool(n_workers) as pool:
        for i, (count, fails, examples) in enumerate(pool.imap_unordered(verify_chunk, chunks)):
            checked += count
            total_fails += fails
            for mask, ht, it in examples:
                print(f"  FAIL mask {mask}: H={ht}, I={it}")

            if (i + 1) % 64 == 0 or checked == total:
                elapsed = time.time() - start
                rate = checked / elapsed
                eta = (total - checked) / rate if rate > 0 else 0
                print(f"  [{100*checked/total:5.1f}%] {checked:>12,}/{total:,} | "
                      f"{total_fails} fails | {rate:.0f}/s | ETA {eta/3600:.2f}hr | "
                      f"{elapsed:.0f}s elapsed")

            if total_fails > 10:
                print("Too many failures, aborting.")
                pool.terminate()
                return

    elapsed = time.time() - start
    print(f"\n{'='*60}")
    print(f"TOTAL: {checked:,} checked, {total_fails} failures, {elapsed:.0f}s")
    if total_fails == 0:
        print(f"\n*** PROVED: H(T) = I(Omega(T), 2) for ALL n=8 tournaments ***")


if __name__ == "__main__":
    main()
