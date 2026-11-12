#!/usr/bin/env python3
"""Lane twoapex, deliverable (c): residue classes BY VALUE.  For a dyadic window [2^k, 2^(k+1)) of odd N and a
modulus m, the share of the class c among the representations (words of length <= L with H = N in the window),
times the number of admissible classes; with the multiplicity cum(N) (all words of length <= L) and with the
primitive multiplicity prim(N) = c_L(N) - c_(L-1)(N) (words with first letter not 0 and last letter not 3).
Also the same shares for the SET of attained values (which is uniform as soon as every odd N of the window is
attained).

usage: python c5_byvalue_windows.py BOUND_EXP STATE_FILE
"""
import sys

import numpy as np


def main():
    bexp, sfile = int(sys.argv[1]), sys.argv[2]
    st = np.load(sfile)
    cum, last, prev, L = st["cum"].astype(np.float64), st["last"].astype(np.float64), st["prev"].astype(np.float64), int(st["upto"])
    prim = last - prev
    N = 2 * np.arange(len(cum), dtype=np.int64) + 1
    print(f"words of length <= {L}; by-value class factors in dyadic windows (factor 1 = the uniform law on the"
          f" admissible classes)")
    for m in (3, 5, 7, 9, 11, 13, 4, 8, 16):
        odd_only = m % 2 == 0
        active = [c for c in range(m) if (c % 2 == 1 or not odd_only)]
        print(f"===== modulus {m} =====")
        for label, wt in (("all words", cum), ("primitive words", prim)):
            print(f"  weights: {label}")
            for k in range(14, bexp):
                lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
                w = wt[lo:hi]
                r = N[lo:hi] % m
                tot = w.sum()
                if tot == 0:
                    continue
                sh = np.bincount(r, weights=w, minlength=m)[active] / tot * len(active)
                # statistical scale if the representations were independent: sqrt(#classes / total)
                scale = np.sqrt(len(active) / tot)
                print(f"    [2^{k}, 2^{k + 1}): " + " ".join(f"{c}:{v:.5f}" for c, v in zip(active, sh))
                      + f" | max |f-1| = {np.abs(sh - 1).max():.5f} (independent-sample scale {scale:.5f})")


if __name__ == "__main__":
    main()
