#!/usr/bin/env python3
"""Lane twoapex - added by the FINISHING agent: the collision energy of the two-apex words of exact length n,

    E_n = sum_N c_n(N)^2 = #{(w, w') : |w| = |w'| = n, H(w) = H(w')},

from the census states of b1_census.py (arrays last = c_L, prev = c_(L-1); the states after L = 14, 16, 17 were
kept).  "Linear energy" would be E_n <= K (16/3)^n (the least possible order, since 4^n words take at most 3^n
odd values); the table gives E_n (3/16)^n and the diagonal share 4^n / E_n.
For n = 16 the 667 words with H >= 2^26 are missing from the histogram and for n = 17 the 301,885,533 words with
H >= 2^26 are missing (the numbers for these two lengths are lower bounds; for n = 17 the bound is not sharp).

usage: python z2_energy_from_census.py STATE_DIR
"""
import sys

import numpy as np


def main():
    sdir = sys.argv[1]
    done = set()
    print("  n   words counted            E_n        E_n (3/16)^n    4^n / E_n   max_N c_n(N)")
    for tag in ("_L14", "_L16", "_L17"):
        st = np.load(f"{sdir}/census_state_b26{tag}.npz")
        L = int(st["upto"])
        for name, n in (("prev", L - 1), ("last", L)):
            if n in done:
                continue
            done.add(n)
            c = st[name].astype(np.float64)
            tot = int(c.sum())
            e = float((c * c).sum())
            note = "" if tot == 4 ** n else f"   (incomplete: {4 ** n - tot} words above 2^26 not counted)"
            print(f"  {n:2d}  {tot:14d}  {e:18.6e}  {e * (3 / 16) ** n:14.5f}  {4.0 ** n / e:10.6f}   {int(c.max())}{note}")
            del c


if __name__ == "__main__":
    main()
