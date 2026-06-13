#!/usr/bin/env python3
"""
Test palindromic N for a non-circulant VT tournament at n=21.

The first non-circulant VT tournaments appear at n=21 (22 of them).
If palindromic N holds for these, our theorem extends beyond circulant.

We can't easily get McKay's data, so let's construct a known non-circulant
VT tournament. The Paley tournament T_7 on Z/7Z has Aut(T_7) = AGL(1,7)
restricted to QR multipliers. But T_7 IS circulant.

For n=21 = 3*7: we can try a tournament on Z/3Z x Z/7Z with a non-cyclic
automorphism group. But constructing such a tournament requires care.

Instead, let's verify that our DP approach works at n=15 and n=17 (all circulant).

kind-pasteur-2026-03-06-S25d
"""

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def compute_f_dp(T, n, target_d=None):
    """Compute f(d,j) = N(0,d,j) using DP.
    If target_d specified, only compute for that d (faster)."""
    full = (1 << n) - 1

    # dp[mask][last] = # Ham paths visiting mask ending at last
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for mask in range(1, 1 << n):
        for last in range(n):
            if not ((mask >> last) & 1):
                continue
            cnt = dp.get((mask, last), 0)
            if cnt == 0:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if T.get((last, nxt), 0) == 0:
                    continue
                nkey = (mask | (1 << nxt), nxt)
                dp[nkey] = dp.get(nkey, 0) + cnt

    H = sum(dp.get((full, v), 0) for v in range(n))

    # dp_suf[mask][first]
    dp_suf = {}
    for v in range(n):
        dp_suf[(1 << v, v)] = 1
    for popcount in range(2, n + 1):
        for mask in range(1, 1 << n):
            if bin(mask).count('1') != popcount:
                continue
            for first in range(n):
                if not ((mask >> first) & 1):
                    continue
                total = 0
                for nxt in range(n):
                    if nxt == first or not ((mask >> nxt) & 1):
                        continue
                    if T.get((first, nxt), 0) == 0:
                        continue
                    total += dp_suf.get((mask & ~(1 << first), nxt), 0)
                if total > 0:
                    dp_suf[(mask, first)] = total

    ds = [target_d] if target_d else range(1, n)
    f = {}
    for d in ds:
        f[d] = [0] * (n - 1)
        for x, y in [(0, d), (d, 0)]:
            if T.get((x, y), 0) == 0:
                continue
            for j in range(n - 1):
                pc = j + 1
                for prefix_mask in range(1, 1 << n):
                    if bin(prefix_mask).count('1') != pc:
                        continue
                    if not ((prefix_mask >> x) & 1):
                        continue
                    if (prefix_mask >> y) & 1:
                        continue
                    pcnt = dp.get((prefix_mask, x), 0)
                    if pcnt == 0:
                        continue
                    scnt = dp_suf.get((full ^ prefix_mask, y), 0)
                    if scnt == 0:
                        continue
                    f[d][j] += pcnt * scnt

    return f, H


# Test at n=15
import time

print("=" * 70)
print("n=15: Palindromic N for a circulant tournament")
print("=" * 70)

n = 15
# Use a simple generator set
S = set()
for d in range(1, (n+1)//2):
    S.add(d)  # Add first half
# S = {1,2,3,4,5,6,7}
print(f"  S = {sorted(S)}")

T = circulant_tournament(n, S)
t0 = time.time()
f, H = compute_f_dp(T, n, target_d=1)
t1 = time.time()
print(f"  H = {H}, computed in {t1-t0:.1f}s")

d = 1
fvals = f[d]
is_pal = all(fvals[j] == fvals[n-2-j] for j in range(n-1))
alt_sum = sum((-1)**j * fvals[j] for j in range(n-1))
print(f"  f(1) = {fvals}")
print(f"  palindromic: {is_pal}, alt_sum: {alt_sum}")

# Check a few more d values
for d in [2, 3, 7]:
    f2, _ = compute_f_dp(T, n, target_d=d)
    fvals = f2[d]
    is_pal = all(fvals[j] == fvals[n-2-j] for j in range(n-1))
    alt_sum = sum((-1)**j * fvals[j] for j in range(n-1))
    print(f"  f({d}): pal={is_pal}, alt_sum={alt_sum}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
