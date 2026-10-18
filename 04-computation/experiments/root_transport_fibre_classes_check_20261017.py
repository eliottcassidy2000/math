#!/usr/bin/env python3
"""Brute-force check of THM-4667 statement 8 (sibling factorisation of the odd-transport fibres), S41.

Seeds s = 1 (mod 8) (heads, a = 2) and 2s+1 (tails, a = 1) have the same image; every other seed (x = 7 mod 8 or
x = 11 mod 16) is a singleton x = (2u - 1)/3 with u = 2 (mod 3) and a = 1. Checked for all seeds below M = 2 x 10^6
(pairs for heads below M/2). The independent S41 audit (s41_audit_20261017/a2_fibre_classes.py) extends this to 4 x 10^6.
Usage: python root_transport_fibre_classes_check_20261017.py
Session opus-2026-10-17-S41.
"""
import numpy as np
M = 2_000_000
s = np.arange(3, M + 1, 2, dtype=np.int64); s = s[s % 8 != 5]
u = 3 * s + 1; a = np.where(u % 4 == 0, 2, 1); u = u >> a
m = u % 8 == 5
while m.any():
    u[m] = (u[m] - 1) // 4; m = u % 8 == 5
# pair property: s = 1 mod 8 -> 2s+1 has the same image, a(s)=2, a(2s+1)=1
img = dict()
pos = {int(x): i for i, x in enumerate(s[:200000])}
ok_pair = True
for i in range(len(s)):
    x = int(s[i])
    if x > M // 2: break
    if x % 8 == 1:
        y = 2 * x + 1
        # compute image of y directly
        v = 3 * y + 1; aa = 0
        while v % 2 == 0: v //= 2; aa += 1
        while v % 8 == 5: v = (v - 1) // 4
        if v != u[i] or aa != 1 or a[i] != 2:
            ok_pair = False; print('pair fails', x); break
print('pair property for s = 1 mod 8, s <= M/2:', ok_pair)
# singleton property: elements with s = 3 or 7 mod 8 that are NOT 2t+1 for t = 1 mod 8 must be (2u-1)/3 with u = 2 mod 3
x = s; t = (x - 1) // 2
is_partner = (x % 16 == 3)   # 2t+1 with t = 1 mod 8  <=> x = 3 mod 16
single = (x % 8 != 1) & (~is_partner)
uu = u[single]; xx = x[single]
print('singletons:', single.sum(), ' all equal (2u-1)/3 with u = 2 mod 3:', bool(np.all((3 * xx == 2 * uu - 1) & (uu % 3 == 2))))
print('a = 1 for all singletons:', bool(np.all(a[single] == 1)))
