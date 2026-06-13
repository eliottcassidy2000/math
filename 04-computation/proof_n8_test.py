#!/usr/bin/env python3
"""Quick test of n=8 OCF verifier on a small sample."""

import time, random
from itertools import combinations, permutations

N = 8
FULL = (1 << N) - 1


def ham_path_count(T):
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << N):
        for last in range(N):
            if dp[mask][last] == 0:
                continue
            c = dp[mask][last]
            for nxt in range(N):
                if mask & (1 << nxt):
                    continue
                if T[last * N + nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[FULL])


def find_odd_cycles(T):
    cycles = set()
    for length in range(3, N + 1, 2):
        for combo in combinations(range(N), length):
            start = combo[0]
            others = list(combo[1:])
            _dfs(T, start, others, [start], combo, length, cycles)
    return cycles


def _dfs(T, start, remaining, path, combo, length, cycles):
    if len(path) == length:
        if T[path[-1] * N + start]:
            cycles.add(frozenset(combo))
        return
    for i, v in enumerate(remaining):
        if T[path[-1] * N + v]:
            _dfs(T, start, remaining[:i] + remaining[i+1:], path + [v], combo, length, cycles)


def indpoly2(T):
    cycles = list(find_odd_cycles(T))
    nc = len(cycles)
    if nc == 0:
        return 1
    adj = [0] * nc  # bitmask adjacency
    for a in range(nc):
        for b in range(a + 1, nc):
            if cycles[a] & cycles[b]:
                adj[a] |= (1 << b)
                adj[b] |= (1 << a)
    total = 1
    for size in range(1, nc + 1):
        for combo in combinations(range(nc), size):
            ok = True
            for idx in range(len(combo)):
                mask = 0
                for idx2 in range(idx + 1, len(combo)):
                    mask |= (1 << combo[idx2])
                if adj[combo[idx]] & mask:
                    ok = False
                    break
            if ok:
                total += (1 << size)
    return total


def test():
    arc_pairs = [(a, b) for a in range(N) for b in range(a+1, N) if (a,b) != (0,1)]
    assert len(arc_pairs) == 27

    random.seed(42)
    sample_size = 200
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

        Tp = T[:]
        Tp[0 * N + 1] = 0
        Tp[1 * N + 0] = 1

        ht = ham_path_count(T)
        htp = ham_path_count(Tp)
        it = indpoly2(T)
        itp = indpoly2(Tp)

        if ht - htp == it - itp:
            passes += 1
        else:
            print(f"  FAIL trial {trial}: dH={ht-htp}, dI={it-itp}")

    elapsed = time.time() - t0
    rate = sample_size / elapsed
    print(f"Result: {passes}/{sample_size} pass ({elapsed:.1f}s, {rate:.1f}/s)")
    print(f"Estimated full time: {(1<<27)/rate:.0f}s = {(1<<27)/rate/3600:.1f}hr")


if __name__ == "__main__":
    test()
