#!/usr/bin/env python3
"""Scratch probe (not canon): the block identity of the one-apex count and the collisions by common prefix.
x(h) = 2h + 1 - 2q(h), q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2).
Identity (h = A 2^m + B, 0 <= B < 2^m):  x(h) = A c(B) + f(B) - 2 q(A) - B_(m-1) (A mod 2),
   c(B) = 2^(m+1) - 2 R_m(B) - B_(m-1),  R_m(B) = sum_{i <= m-2} B_i 2^(m-2-i),  f(B) = 2B + 1 - 2 q_m(B),
   q_m(B) = q of B read as an m-bit block (prefixes inside the block only).
Collisions x(h) = x(h') among n-bit h (top bit 1): histogram by the length of the common prefix and suffix.
"""
import sys
from collections import Counter

import numpy as np


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    size = 1 << n
    q = np.zeros(size, np.int64)
    for h in range(1, size):
        hh = h >> 1
        q[h] = q[hh] + (((hh + 1) >> 1) if (h & 1) else 0)
    x = 2 * np.arange(size, dtype=np.int64) + 1 - 2 * q
    # identity
    ok = True
    rng = np.random.default_rng(1)
    for _ in range(200000):
        m = int(rng.integers(1, 12))
        A = int(rng.integers(1, 1 << 9))
        B = int(rng.integers(0, 1 << m))
        h = (A << m) + B
        if h >= size:
            continue
        # q_m(B): prefixes inside the block
        qm = 0
        for i in range(m):
            if (B >> i) & 1:
                pre = B >> (i + 1)
                qm += (pre + 1) >> 1
        R = 0
        for i in range(0, m - 1):
            if (B >> i) & 1:
                R += 1 << (m - 2 - i)
        top = (B >> (m - 1)) & 1
        c = (1 << (m + 1)) - 2 * R - top
        f = 2 * B + 1 - 2 * qm
        if int(x[h]) != A * c + f - 2 * int(q[A]) - top * (A & 1):
            ok = False
            print("IDENTITY FAILS", A, m, B)
            break
    print("block identity x(A 2^m + B) = A c(B) + f(B) - 2 q(A) - B_(m-1) (A mod 2):", ok)
    # c(B) is a bijection from m-bit blocks onto (2^m, 2^(m+1)]
    m = 10
    cs = set()
    for B in range(1 << m):
        R = 0
        for i in range(0, m - 1):
            if (B >> i) & 1:
                R += 1 << (m - 2 - i)
        cs.add((1 << (m + 1)) - 2 * R - ((B >> (m - 1)) & 1))
    print("c is a bijection onto (2^m, 2^(m+1)] for m = 10:", cs == set(range((1 << m) + 1, (1 << (m + 1)) + 1)))
    # collisions among n-bit words
    lo, hi = 1 << (n - 1), 1 << n
    vals = x[lo:hi]
    order = np.argsort(vals, kind="stable")
    sv = vals[order]
    hist = Counter()
    pairs = 0
    i = 0
    N = len(sv)
    maxmult = 0
    while i < N:
        j = i + 1
        while j < N and sv[j] == sv[i]:
            j += 1
        if j - i > 1:
            grp = order[i:j] + lo
            maxmult = max(maxmult, j - i)
            for a in range(len(grp)):
                for b in range(a + 1, len(grp)):
                    h1, h2 = int(grp[a]), int(grp[b])
                    d = h1 ^ h2
                    pre = n - d.bit_length()                 # common prefix length
                    suf = (d & -d).bit_length() - 1          # common suffix length
                    hist[(pre, suf)] += 1
                    pairs += 1
        i = j
    distinct = len(np.unique(sv))
    print(f"n = {n}: {N} words, {distinct} distinct values ({distinct / N:.4f}), {pairs} colliding pairs ({pairs / N:.4f} per word), max multiplicity {maxmult}")
    by_pre = Counter()
    diff = Counter()
    for (pre, suf), c in hist.items():
        by_pre[pre] += c
        diff[pre - suf] += c
    print("colliding pairs by common prefix length:", sorted(by_pre.items()))
    print("prefix minus suffix:", sorted(diff.items()))
    q4 = sum(c for (pre, suf), c in hist.items() if pre >= n / 4)
    print(f"pairs with common prefix >= n/4: {q4} ({q4 / max(pairs, 1):.4f} of all)")


if __name__ == "__main__":
    main()
