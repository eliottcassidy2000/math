#!/usr/bin/env python3
"""
CORRECTION: Paley tournaments exist ONLY at p = 3 mod 4.

At p = 1 mod 4, -1 is a QR, so QR is closed under negation,
and the "Paley digraph" has bidirectional edges — NOT a tournament.

Previous computations of "T_13" and "T_17" used QR mod 13 = {1,3,4,9,10,12}
and QR mod 17 respectively. These are NOT tournaments.

Valid Paley tournaments: p = 3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, ...

For n=13 (verification of THM-052), we should use a proper circulant tournament.
A valid generator set S for n=13 must have |S| = 6 and for each d in {1,...,12},
exactly one of d, 13-d is in S.

kind-pasteur-2026-03-06-S25d
"""

import time

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def is_valid_tournament(T, n):
    """Check that T is actually a tournament (exactly one direction per pair)."""
    for i in range(n):
        for j in range(i+1, n):
            if T[(i,j)] == T[(j,i)]:
                return False
    return True

def compute_f_single_dp(T, n, d):
    """Compute f(d,j) for a single d using DP."""
    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not ((mask >> last) & 1): continue
            cnt = dp.get((mask, last), 0)
            if cnt == 0: continue
            for nxt in range(n):
                if (mask >> nxt) & 1: continue
                if T.get((last, nxt), 0) == 0: continue
                nkey = (mask | (1 << nxt), nxt)
                dp[nkey] = dp.get(nkey, 0) + cnt
    H = sum(dp.get((full, v), 0) for v in range(n))

    dp_suf = {}
    for v in range(n):
        dp_suf[(1 << v, v)] = 1
    for popcount in range(2, n + 1):
        for mask in range(1, 1 << n):
            if bin(mask).count('1') != popcount: continue
            for first in range(n):
                if not ((mask >> first) & 1): continue
                total = 0
                for nxt in range(n):
                    if nxt == first or not ((mask >> nxt) & 1): continue
                    if T.get((first, nxt), 0) == 0: continue
                    total += dp_suf.get((mask & ~(1 << first), nxt), 0)
                if total > 0:
                    dp_suf[(mask, first)] = total

    f = [0] * (n - 1)
    for x, y in [(0, d), (d, 0)]:
        if T.get((x, y), 0) == 0: continue
        for j in range(n - 1):
            pc = j + 1
            for prefix_mask in range(1, 1 << n):
                if bin(prefix_mask).count('1') != pc: continue
                if not ((prefix_mask >> x) & 1): continue
                if (prefix_mask >> y) & 1: continue
                pcnt = dp.get((prefix_mask, x), 0)
                if pcnt == 0: continue
                scnt = dp_suf.get((full ^ prefix_mask, y), 0)
                if scnt == 0: continue
                f[j] += pcnt * scnt
    return f, H


# ============================================================
# Show which previous computations were on actual tournaments
# ============================================================
print("=" * 70)
print("VERIFICATION: Which previous computations were valid tournaments?")
print("=" * 70)

test_cases = [
    (5, {1, 4}, "Paley T_5 (p=5, 1 mod 4)"),
    (7, {1, 2, 4}, "Paley T_7 (p=7, 3 mod 4)"),
    (7, {1, 2, 3}, "Circulant n=7 {1,2,3}"),
    (9, {1, 2, 4, 6}, "Circulant n=9"),
    (11, {1, 3, 4, 5, 9}, "Paley T_11 (p=11, 3 mod 4)"),
    (13, {1, 3, 4, 9, 10, 12}, "Paley T_13 (p=13, 1 mod 4)"),
    (15, {1, 2, 3, 4, 5, 6, 7}, "Circulant n=15"),
    (17, {1, 2, 4, 8, 9, 13, 15, 16}, "Paley T_17 (p=17, 1 mod 4)"),
]

for n, S, label in test_cases:
    T = circulant_tournament(n, S)
    valid = is_valid_tournament(T, n)
    print(f"  {label}: {'TOURNAMENT' if valid else 'NOT A TOURNAMENT!'}")

# ============================================================
# Re-verify n=13 with a VALID tournament generator set
# ============================================================
print("\n" + "=" * 70)
print("n=13: Re-verification with valid circulant tournament")
print("=" * 70)

n = 13
# Valid generator: first half {1,2,3,4,5,6}
S = {1, 2, 3, 4, 5, 6}
T = circulant_tournament(n, S)
valid = is_valid_tournament(T, n)
print(f"  S = {sorted(S)}, valid tournament: {valid}")

t0 = time.time()
for d in [1, 2, 6, 7]:
    f, H = compute_f_single_dp(T, n, d)
    is_pal = all(f[j] == f[n-2-j] for j in range(n-1))
    alt_sum = sum((-1)**j * f[j] for j in range(n-1))
    print(f"  d={d}: H={H}, pal={is_pal}, alt={alt_sum}")
    if d == 1:
        print(f"    f(1) = {f}")

elapsed = time.time() - t0
print(f"  Total time: {elapsed:.1f}s")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
