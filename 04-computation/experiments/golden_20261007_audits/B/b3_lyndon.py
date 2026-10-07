#!/usr/bin/env python3
"""Independent Lyndon-word census of integer Terras cycles for L <= LMAX, all k (no j >= L/2 pruning),
plus rational-cycle denominators for the Ellison shapes and the E(L,k) ranking for 12 <= L <= 1500.
Own code: Duval's algorithm for Lyndon words with a fixed number of ones (via simple recursion with
canonical-rotation test), parity-vector cycle formula x = d/(2^L - 3^k)."""
import sys
from math import comb, gcd, log
from itertools import combinations

def mobius(n):
    r, q = 1, 2
    while q * q <= n:
        if n % q == 0:
            n //= q
            if n % q == 0:
                return 0
            r = -r
        q += 1
    return -r if n > 1 else r

def lyn_count(L, k):
    g = gcd(L, k) if k else L
    s = sum(mobius(e) * comb(L // e, k // e) for e in range(1, g + 1) if g % e == 0)
    assert s % L == 0
    return s // L

def lyndon_words(L, k):
    """All binary Lyndon words of length L with k ones (generate by FKM restricted to k ones)."""
    a = [0] * (L + 1)
    out = []
    def gen(t, p, ones):
        if ones > k or ones + (L - t + 1) < k:
            return
        if t > L:
            if p == L and ones == k:
                out.append(tuple(a[1:]))
            return
        a[t] = a[t - p]
        gen(t + 1, p, ones + a[t])
        if a[t - p] == 0:
            a[t] = 1
            gen(t + 1, t, ones + 1)
            a[t] = 0
    if L == 1:
        return [(k,)] if k in (0, 1) else []
    gen(1, 1, 0)
    return out

def cycle_value(w):
    """x with Terras parity vector w (periodic), as (numerator d, denominator 2^L - 3^k)."""
    d, k = 0, 0
    for t, b in enumerate(w):
        if b:
            d = 3 * d + 2 ** t; k += 1
    return d, 2 ** len(w) - 3 ** k

def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2

if __name__ == "__main__":
    LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    print(f"--- integer cycles via Lyndon words, all k, L <= {LMAX} ---")
    tot = 0
    for L in range(1, LMAX + 1):
        for k in range(0, L + 1):
            W = lyndon_words(L, k)
            assert len(W) == lyn_count(L, k), (L, k, len(W), lyn_count(L, k))
            tot += len(W)
            for w in W:
                d, D = cycle_value(w)
                if d % D == 0:
                    x = d // D
                    # verify by iteration
                    y = x
                    for _ in range(L):
                        y = T(y)
                    assert y == x
                    print(f"  L={L} k={k} x={x} word={''.join(map(str, w))}")
    print("  words enumerated:", tot)
