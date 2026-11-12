#!/usr/bin/env python3
"""Lane twoapex - added by the FINISHING agent.  Why variance / mean = 2.00 for the primitive multiplicities:
the converse map w -> w* (reverse the word, letter (a, b) -> (1 - b, 1 - a); H(w*) = H(w), checked in
a1_verify_transfer.py) is an involution of the primitive words of each length, so prim(N) is even unless a
self-converse word has value N.  Share of odd prim(N) per dyadic window, and the dispersion of prim(N) / 2.

usage: python z3_pairing_check.py STATE_FILE
"""
import sys

import numpy as np

st = np.load(sys.argv[1])
prim = st["last"].astype(np.int64) - st["prev"].astype(np.int64)
L = int(st["upto"])
print(f"primitive words of length <= {L}")
print("  window        share of odd prim(N)   mean of prim/2   var/mean of prim/2 (blocks of 256)")
B = 256
for k in range(16, 26):
    lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
    seg = prim[lo:hi]
    half = seg.astype(np.float64) / 2
    nb = len(half) // B
    blk = half[: nb * B].reshape(nb, B)
    var = float(((blk - blk.mean(axis=1, keepdims=True)) ** 2).sum() / (nb * (B - 1)))
    print(f"  [2^{k}, 2^{k + 1})   {float((seg & 1).mean()):10.6f}        {half.mean():10.2f}      {var / half.mean():8.4f}")
