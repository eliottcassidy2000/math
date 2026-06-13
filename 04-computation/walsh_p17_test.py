#!/usr/bin/env python3
"""
walsh_p17_test.py -- Extend Walsh sign analysis to p=17

At p=17 (p=1 mod 4), m=8, 256 circulant tournaments.
Tests whether the phase transition pattern continues.

Key questions:
1. Is deg-4 still positive (helping Interval) at p=17?
2. Does deg-6 appear with significant energy?
3. Is Interval still the maximizer at p=17 (p=1 mod 4)?

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from math import comb
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def walsh_decomposition(H_dict, m):
    n = 1 << m
    h_hat = {}
    for bits in range(n):
        S = [i for i in range(m) if bits & (1 << i)]
        coeff = 0
        for sigma_bits in range(n):
            sigma = [(1 if sigma_bits & (1 << i) else -1) for i in range(m)]
            H = H_dict[tuple(sigma)]
            prod = 1
            for i in S:
                prod *= sigma[i]
            coeff += H * prod
        h_hat[tuple(S)] = coeff / n
    return h_hat


def main():
    p = 17
    m = (p - 1) // 2
    print(f"p={p}, m={m}, p mod 4 = {p % 4}")
    print(f"Computing H for {1 << m} circulant tournaments...")

    pairs = [(s, p - s) for s in range(1, m + 1)]

    t0 = time.time()
    H_dict = {}
    for bits in range(1 << m):
        sigma = []
        S = []
        for i, (a, b) in enumerate(pairs):
            if bits & (1 << i):
                S.append(a)
                sigma.append(+1)
            else:
                S.append(b)
                sigma.append(-1)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_dict[tuple(sigma)] = H

        if (bits + 1) % 32 == 0:
            elapsed = time.time() - t0
            done = bits + 1
            total = 1 << m
            eta = elapsed / done * (total - done)
            print(f"  {done}/{total} computed ({elapsed:.1f}s elapsed, ETA {eta:.0f}s)")

    t1 = time.time()
    print(f"\nAll {1 << m} H values computed in {t1-t0:.1f}s")

    # Distinct H values
    distinct = sorted(set(H_dict.values()), reverse=True)
    print(f"\n{len(distinct)} distinct H values:")
    for h in distinct[:10]:
        count = sum(1 for v in H_dict.values() if v == h)
        print(f"  H={h}: {count} orientations")
    if len(distinct) > 10:
        print(f"  ... ({len(distinct)} total)")
    for h in distinct[-3:]:
        count = sum(1 for v in H_dict.values() if v == h)
        print(f"  H={h}: {count} orientations")

    # Mean H
    mean_H = sum(H_dict.values()) / len(H_dict)

    # Interval
    sigma_int = tuple([1] * m)
    H_int = H_dict[sigma_int]
    print(f"\nMean H = {mean_H:.2f}")
    print(f"H(Interval) = {H_int}")
    print(f"Excess = {H_int - mean_H:.2f}")

    max_H = max(H_dict.values())
    print(f"Max H = {max_H}")
    print(f"Interval is maximizer: {H_int == max_H}")

    # Walsh decomposition
    print(f"\nComputing Walsh decomposition ({1 << m} x {1 << m} transform)...")
    t2 = time.time()
    h_hat = walsh_decomposition(H_dict, m)
    t3 = time.time()
    print(f"Walsh transform computed in {t3-t2:.1f}s")

    # Energy by degree
    energy = defaultdict(float)
    contrib = defaultdict(float)
    for S, coeff in h_hat.items():
        deg = len(S)
        energy[deg] += coeff**2
        contrib[deg] += coeff  # at all-+1

    total_E = sum(energy.values())
    print(f"\nWalsh energy by degree:")
    print(f"  {'deg':>4} {'energy':>22} {'pct':>10} {'contrib_int':>16}")
    for deg in sorted(energy):
        E = energy[deg]
        pct = 100 * E / total_E if total_E > 0 else 0
        c = contrib[deg]
        if E > 1e-10 or abs(c) > 1e-10:
            print(f"  {deg:>4} {E:>22.1f} {pct:>10.6f}% {c:>+16.2f}")

    # Sign pattern of degree-2 and degree-4
    print(f"\nSign pattern:")
    for target_deg in [2, 4, 6]:
        items = [(S, c) for S, c in h_hat.items() if len(S) == target_deg and abs(c) > 1e-10]
        if not items:
            continue
        n_pos = sum(1 for _, c in items if c > 0)
        n_neg = sum(1 for _, c in items if c < 0)
        mags = sorted(set(round(abs(c), 2) for _, c in items))
        total = sum(c for _, c in items)
        print(f"  deg {target_deg}: {n_pos}+ {n_neg}- coefficients, "
              f"magnitudes={mags}, sum={total:.2f}")

    # Phase transition summary
    print(f"\n{'='*60}")
    print(f"PHASE TRANSITION SUMMARY for p={p} (mod 4 = {p%4}):")
    print(f"  deg-2 contribution to Interval: {contrib.get(2, 0):>+16.2f}")
    print(f"  deg-4 contribution to Interval: {contrib.get(4, 0):>+16.2f}")
    print(f"  deg-6 contribution to Interval: {contrib.get(6, 0):>+16.2f}")
    print(f"  deg-8 contribution to Interval: {contrib.get(8, 0):>+16.2f}")
    print(f"  Total non-constant contribution: {sum(c for d,c in contrib.items() if d > 0):>+16.2f}")
    print(f"  Interval excess over mean: {H_int - mean_H:>+16.2f}")


if __name__ == '__main__':
    main()
