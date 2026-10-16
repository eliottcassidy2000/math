#!/usr/bin/env python3
"""S37 audit, part 7b: brute-force circulant orbit counts on Z_(p^2) for the non-generous prime p = 40487.
g = 5 (least primitive root mod p, 5^(p-1) = 1 mod p^2), g' = 5^((p-1)/d') = 25 (d' = odd part of p - 1 = 20243):
the companion note predicts 2E orbits with E = (p+1)(p-1)/(2d') = 40488.  Control: the least primitive root h of
p^2, h' = h^2, predicts E = (p-1)/d' = 2.  Bitmap traversal of all p^2 - 1 = 1,639,197,168 nonzero residues.
"""
import time
import numpy as np
import numba as nb
import sympy as sp

t00 = time.time()


@nb.njit(cache=False)
def count_orbits(m, g):
    vis = np.zeros((m >> 3) + 1, np.uint8)
    norb = 0
    selfneg = 0
    for x in range(1, m):
        if (vis[x >> 3] >> (x & 7)) & 1:
            continue
        norb += 1
        y = x
        neg = m - x
        found = False
        while True:
            vis[y >> 3] |= np.uint8(1 << (y & 7))
            if y == neg:
                found = True
            y = (y * g) % m
            if y == x:
                break
        if found:
            selfneg += 1
    return norb, selfneg


p = 40487
m = p * p
d = p - 1
while d % 2 == 0:
    d //= 2
which = __import__('sys').argv[1] if len(__import__('sys').argv) > 1 else 'g'
if which == 'g':
    gp = pow(5, (p - 1) // d, m)
    norb, selfneg = count_orbits(m, gp)
    print(f"p = 40487, g = 5, g' = {gp}: {norb} orbits on the nonzero residues mod p^2 (self-negative {selfneg}) -> "
          f"E = {norb // 2}; formula (p+1)(p-1)/(2d') = {(p + 1) * (p - 1) // (2 * d)}  [{time.time() - t00:.1f} s]", flush=True)
else:
    h = 2
    while not (sp.n_order(h % p, p) == p - 1 and pow(h, p - 1, m) != 1):
        h += 1
    hp = pow(h, (p - 1) // d, m)
    norb, selfneg = count_orbits(m, hp)
    print(f"control: least primitive root of p^2 is {h}; h' = {hp}: {norb} orbits (self-negative {selfneg}) -> "
          f"E = {norb // 2}; formula (p-1)/d' = {(p - 1) // d}  [{time.time() - t00:.1f} s]", flush=True)
