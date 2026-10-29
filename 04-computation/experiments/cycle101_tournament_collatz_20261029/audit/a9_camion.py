#!/usr/bin/env python3
"""Spot check of the results note's Camion claims (section 1a, 1d, 6): s(m) = strong labelled tournaments,
c(m) = 2^((m-1)/2) circulant count for odd m; the strong Camion test is s(n) = c(n) (mod n).
Direct recurrence modulo n for the small cases (no Korselt reduction), plus s(19) - 2^9 mod 101."""
from math import comb

def s_mod(n_max, mod):
    s = [0] * (n_max + 1)
    t = [pow(2, j * (j - 1) // 2, mod) for j in range(n_max + 1)]
    row = [1]  # Pascal row m-1 mod mod
    for m in range(1, n_max + 1):
        # row for m
        new = [1] * (m + 1)
        for k in range(1, m):
            new[k] = (row[k - 1] + row[k]) % mod
        row = new
        acc = t[m]
        for k in range(1, m):
            acc -= row[k] * s[k] % mod * t[m - k]
        s[m] = acc % mod
    return s

def s_exact(n):
    s = [0] * (n + 1)
    for m in range(1, n + 1):
        s[m] = 2 ** comb(m, 2) - sum(comb(m, k) * s[k] * 2 ** comb(m - k, 2) for k in range(1, m))
    return s

S = s_exact(25)
print("s(1..8) =", S[1:9], "(A054946: 1, 0, 2, 24, 544, 22320, 1677488, 236522496)")
print("s(19) - 2^9 mod 101 =", (S[19] - 2 ** 9) % 101)
for n in (169, 1313, 1501, 1919, 3211):
    sm = s_mod(n, n)[n]
    c = pow(2, (n - 1) // 2, n)
    print(f"n = {n}: s(n) mod n = {sm}, c(n) mod n = {c}, strong Camion test passes: {sm == c}")
