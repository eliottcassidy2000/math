#!/usr/bin/env python3
"""
Test palindromic N(a,b,j) for Paley T_11.

Uses DP for efficiency.

kind-pasteur-2026-03-06-S25d
"""

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def compute_N_dp(T, n):
    """Compute N(a,b,j) using DP."""
    full = (1 << n) - 1

    # Build dp_prefix[mask][last] = # paths visiting mask ending at last
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

    # Build dp_suffix[mask][first] = # paths visiting mask starting at first
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
                    rest_mask = mask & ~(1 << first)
                    total += dp_suf.get((rest_mask, nxt), 0)
                if total > 0:
                    dp_suf[(mask, first)] = total

    H = sum(dp.get((full, v), 0) for v in range(n))

    # Compute N(0, b, j) for all b (by circulant symmetry, enough to check vertex 0)
    # N_directed(a->b at position j): a at pos j, b at pos j+1
    N_0 = [[0]*(n-1) for _ in range(n)]  # N_0[b][j]

    for b in range(1, n):
        for direction in [(0, b), (b, 0)]:
            x, y = direction
            if T.get((x, y), 0) == 0:
                continue
            for j in range(n-1):
                pc = j + 1  # popcount of prefix
                for prefix_mask in range(1, 1 << n):
                    if bin(prefix_mask).count('1') != pc:
                        continue
                    if not ((prefix_mask >> x) & 1):
                        continue
                    if (prefix_mask >> y) & 1:
                        continue
                    prefix_count = dp.get((prefix_mask, x), 0)
                    if prefix_count == 0:
                        continue
                    suffix_mask = full ^ prefix_mask
                    suffix_count = dp_suf.get((suffix_mask, y), 0)
                    if suffix_count == 0:
                        continue
                    N_0[b][j] += prefix_count * suffix_count

    return N_0, H


print("=" * 70)
print("n=11 Paley: Palindromic N test")
print("=" * 70)

n = 11
# QR mod 11 = {1, 3, 4, 5, 9}
T = circulant_tournament(n, {1, 3, 4, 5, 9})

print(f"  Computing N(0,b,j) for Paley T_{n}...")
N_0, H = compute_N_dp(T, n)
print(f"  H = {H}")

is_pal = True
for b in range(1, n):
    nab = N_0[b]
    for j in range(n-1):
        if nab[j] != nab[n-2-j]:
            is_pal = False
            alt = sum((-1)**j2 * nab[j2] for j2 in range(n-1))
            print(f"  N(0,{b}) NOT palindromic: {nab}, alt_sum={alt}")
            break

if is_pal:
    print(f"  ALL N(0,b,j) are palindromic!")
    print(f"  => M[0,b] = 0 for all b != 0")
    print(f"  => By circulant symmetry, M = (H/n)*I = ({H}/{n})*I = {H//n}*I")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
