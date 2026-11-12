#!/usr/bin/env python3
"""Lane twoapex, deliverable (b): window statistics of the census state written by b1_census.py.

usage: python b2_windows.py BOUND_EXP STATE_DIR
For the words of length <= L (L = the state's last length): per dyadic window [2^k, 2^(k+1)) of odd N,
  multiplicity cum(N) = number of words of length <= L with H = N: min, median, mean, number of zeros;
  primitive multiplicity prim(N) = c_L(N) - c_(L-1)(N) = number of words of length <= L with H = N, first letter
  not 0 and last letter not 3 (no padding source or sink): min, median, mean;
  the largest minimal length over the attained N of the window;
and the list of the odd N not attained.
"""
import sys

import numpy as np


def main():
    bexp, sdir = int(sys.argv[1]), sys.argv[2]
    st = np.load(f"{sdir}/census_state_b{bexp}.npz")
    cum, last, prev, minlen, L = st["cum"], st["last"], st["prev"], st["minlen"], int(st["upto"])
    prim = last.astype(np.int64) - prev.astype(np.int64)
    assert prim.min() >= 0
    print(f"words of length <= {L} (4^0 + ... + 4^{L} = {sum(4 ** i for i in range(L + 1))} words); odd N < 2^{bexp}")
    print(f"total multiplicity recorded below the bound: {int(cum.sum(dtype=np.int64))}; primitive words recorded: "
          f"{int(prim.sum())}")
    print()
    print("  window            odd N    zeros      min   median       mean | prim: min  median      mean | max minimal length")
    first_bad = None
    for k in range(0, bexp):
        lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
        seg = cum[lo:hi]
        ps = prim[lo:hi]
        ml = minlen[lo:hi]
        zeros = int((seg == 0).sum())
        att = seg != 0
        mxlen = int(ml[att].max()) if att.any() else -1
        print(f"  [2^{k:2d}, 2^{k + 1:2d})  {hi - lo:9d} {zeros:8d} {int(seg.min()):8d} {int(np.median(seg)):8d} "
              f"{seg.mean():10.1f} | {int(ps.min()):9d} {int(np.median(ps)):7d} {ps.mean():9.1f} | {mxlen}")
        if zeros and first_bad is None and k >= 5:
            first_bad = k
    miss = np.nonzero(cum == 0)[0]
    missN = 2 * miss + 1
    print()
    print(f"odd N < 2^{bexp} not attained by words of length <= {L}: {len(missN)} in all")
    print("the first 60:", missN[:60].tolist())
    for e in range(10, bexp + 1):
        c = int((missN < (1 << e)).sum())
        print(f"  not attained below 2^{e}: {c}" + (f"  -> {missN[missN < (1 << e)].tolist()}" if c <= 12 else ""))
    # minimal length profile: number of odd N with minimal length exactly n
    print()
    print("number of odd N < bound whose minimal word length is n (attained N only):")
    att = cum != 0
    for n in range(0, L + 1):
        c = int(((minlen == n) & att).sum())
        if c:
            idx = np.nonzero((minlen == n) & att)[0]
            print(f"  n = {n:2d}: {c:9d} values; smallest {2 * int(idx[0]) + 1}, largest {2 * int(idx[-1]) + 1}")
    # largest N such that every odd value in [23, N] is attained
    bad = missN[missN > 21]
    top = int(bad[0]) - 2 if len(bad) else (1 << bexp) - 1
    print()
    print(f"every odd N in [23, {top}] is attained by a word of length <= {L}; log2 = {np.log2(top):.4f}")


if __name__ == "__main__":
    main()
