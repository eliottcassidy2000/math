# ⚠️ WARNING: This script uses QR mod p for p ≡ 1 (mod 4), which does NOT
# produce a tournament (S ∩ (-S) ≠ ∅ gives bidirectional edges).
# Results for those primes are INVALID. See MISTAKE-011b.
# Valid Paley tournaments require p ≡ 3 (mod 4).

#!/usr/bin/env python3
"""
Verify palindromic N at n=17 (Paley tournament T_17).
QR mod 17 = {1, 2, 4, 8, 9, 13, 15, 16}.

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

def compute_f_single_dp(T, n, d):
    """Compute f(d,j) for a single d using DP."""
    full = (1 << n) - 1

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

    f = [0] * (n - 1)
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
                f[j] += pcnt * scnt

    return f, H


n = 17
# QR mod 17
QR = set()
for a in range(1, n):
    QR.add((a*a) % n)

print(f"n={n} Paley: QR = {sorted(QR)}")

T = circulant_tournament(n, QR)

t0 = time.time()
# Check 3 representative d values
for d in [1, 2, 3]:
    f, H = compute_f_single_dp(T, n, d)
    is_pal = all(f[j] == f[n-2-j] for j in range(n-1))
    alt_sum = sum((-1)**j * f[j] for j in range(n-1))
    elapsed = time.time() - t0
    print(f"  d={d}: H={H}, pal={is_pal}, alt={alt_sum}, f={f}")
    print(f"    elapsed: {elapsed:.1f}s")

print("DONE")
