#!/usr/bin/env python3
"""Auditor prior-art (S92).  Exhaustive H-spectrum of all 2^21 labelled tournaments on 7 vertices (numba), and a
random search at order 8 and 9 for the prime labels of the cycle of -17 that are not reached at order 7.
Purpose: an own check that every label |4n + 1| of the 29 arcs of the five known cycles of 3n + 1 on Z other than 7
is a Hamiltonian path count (S92 note section 7 uses THM-1370 and THM-4875 for this)."""
import sys

import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


@njit(cache=True)
def hp(n, adj):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp[mask, v]
            if c != 0:
                out = adj[v] & (~mask) & full
                while out != 0:
                    low = out & (-out)
                    w = 0
                    t = low
                    while t > 1:
                        t >>= 1
                        w += 1
                    dp[mask | low, w] += c
                    out ^= low
    s = 0
    for v in range(n):
        s += dp[full, v]
    return s


@njit(parallel=True, cache=True)
def spectrum7(maxH):
    n = 7
    npairs = 21
    seen = np.zeros((12, maxH + 1), dtype=np.uint8)
    total = 1 << npairs
    chunk = total // 12 + 1
    for t in prange(12):
        lo = t * chunk
        hi = min(total, lo + chunk)
        adj = np.zeros(n, dtype=np.int64)
        for bits in range(lo, hi):
            for v in range(n):
                adj[v] = 0
            k = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (bits >> k) & 1:
                        adj[i] |= 1 << j
                    else:
                        adj[j] |= 1 << i
                    k += 1
            h = hp(n, adj)
            seen[t, h] = 1
    out = np.zeros(maxH + 1, dtype=np.uint8)
    for t in range(12):
        for h in range(maxH + 1):
            if seen[t, h]:
                out[h] = 1
    return out


@njit(cache=True)
def random_search(n, targets, tries, seed):
    np.random.seed(seed)
    found = np.zeros(len(targets), dtype=np.int64)
    wit = np.zeros((len(targets), n), dtype=np.int64)
    adj = np.zeros(n, dtype=np.int64)
    for _ in range(tries):
        for v in range(n):
            adj[v] = 0
        for i in range(n):
            for j in range(i + 1, n):
                if np.random.random() < 0.5:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
        h = hp(n, adj)
        for t in range(len(targets)):
            if found[t] == 0 and targets[t] == h:
                found[t] = 1
                for v in range(n):
                    wit[t, v] = adj[v]
    return found, wit


spec = spectrum7(200)
vals = [int(h) for h in range(201) if spec[h]]
print("order 7: %d distinct H values, max %d" % (len(vals), max(vals)))
print("order 7 odd values missing in [1, 189]:", [h for h in range(1, 190, 2) if not spec[h]])
labels = [1, 5, 17, 9, 3, 7, 19, 55, 27, 79, 39, 67, 199, 99, 295, 147, 439, 219, 655, 327, 163, 487, 243, 727, 363,
          1087, 543, 271, 135]
small = [l for l in labels if l <= 189]
print("labels <= 189 not attained at order 7:", [l for l in small if not spec[l]])
big = np.array(sorted(l for l in labels if l > 189 or not spec[l]), dtype=np.int64)
big = big[big != 7]
print("labels to search at orders 8, 9, 10:", big.tolist())
remaining = big.copy()
for n, tries in ((8, 3_000_000), (9, 3_000_000), (10, 2_000_000)):
    if len(remaining) == 0:
        break
    found, wit = random_search(n, remaining, tries, 12345 + n)
    for t in range(len(remaining)):
        if found[t]:
            w = wit[t]
            # recount the witness
            assert hp(n, w) == remaining[t]
            print("   H = %4d attained at order %d, out-neighbour masks %s" % (remaining[t], n, w.tolist()))
    remaining = remaining[found == 0]
print("labels without a witness after the searches:", remaining.tolist())
# products for the rest (H is multiplicative over the ordered join)
att = set(vals)
print("ALL LABELS EXCEPT 7 ARE PATH COUNTS" if len(remaining) == 0 else "SEE REMAINING")
