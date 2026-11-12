#!/usr/bin/env python3
"""Verifier C: cross-check that {x(h)} is the one-apex spectrum recorded in the cited canon.

THM-4878 (5): 521,915 of the 1,048,576 odd numbers below 2^21 are one-apex counts.
THM-4878 (6): G(1,0) attains 560 of the 1024 odd numbers below 2048.
Session 191 one_apex/PROOF.md section 7: 13 is not a one-apex count; strong counts at chain lengths <= 4: {3,5,9,11}.
THM-4890 (E): g(Y_j) is contained in [2^(j-2) + 1, 2^j - 1]  (Y_j = even j-bit words), i.e. for an n-bit word h
(j = n + 1): 2^(n-1) + 1 <= x(h) <= 2^(n+1) - 1.
"""
import numpy as np

L = 1 << 21
q = np.zeros(L, dtype=np.int64)
for k in range(1, 22):                      # fill by bit length: q[h] = q[h >> 1] + (h & 1) * ceil((h >> 1)/2)
    lo, hi = 1 << (k - 1), 1 << k
    h = np.arange(lo, hi, dtype=np.int64)
    q[lo:hi] = q[h >> 1] + (h & 1) * (((h >> 1) + 1) >> 1)
h = np.arange(L, dtype=np.int64)
x = 2 * h + 1 - 2 * q
print("all x(h) odd:", bool(np.all(x & 1 == 1)))
vals = np.unique(x[x < (1 << 21)])
# x(h) >= h + 2 for h >= 1, so every value below 2^21 comes from some h < 2^21
print("distinct x(h) below 2^21:", vals.size, "(THM-4878 (5): 521,915)")
print("distinct x(h) below 2048:", int(np.sum(vals < 2048)), "(THM-4878 (6): 560)")
print("13 in spectrum:", bool(np.any(vals == 13)), "| 7 in spectrum:", bool(np.any(vals == 7)),
      "| 21 in spectrum:", bool(np.any(vals == 21)))
print("values below 40:", [int(v) for v in vals[vals < 40]])
ok = True
for n in range(1, 22):
    lo, hi = 1 << (n - 1), 1 << n
    xs = x[lo:hi]
    mn, mx = int(xs.min()), int(xs.max())
    if not ((1 << (n - 1)) + 1 <= mn and mx <= (1 << (n + 1)) - 1 and mn >= lo + 2):
        ok = False
    if n <= 8 or n == 21:
        print("n = %2d: x on n-bit words in [%d, %d]; 2^(n-1) + 2 = %d, 2^(n+1) - 1 = %d" % (n, mn, mx, lo + 2, 2 * hi - 1))
print("range checks for n <= 21:", ok)
