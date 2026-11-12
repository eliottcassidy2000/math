#!/usr/bin/env python3
"""Lane twoapex, deliverable (e): two-block sub-languages L1^p L2^q with two-letter alphabets L1, L2
(one free binary digit per position; two free integer parameters, one per block).

For each ordered pair (L1, L2) of two-letter alphabets among 01, 02, 03, 12, 13, 23:
  * restricted seam rank: the exact rank of [H(w1 w2)], w1 in L1^(<=4), w2 in L2^(<=4)  (the number of bilinear
    terms f_k(w1) g_k(w2) needed for a closed form on the family);
  * the largest value at p + q = 16 (growth) and the number of distinct values there;
  * coverage: the share of the odd N of the window [2^15, 2^16) attained by the words with p + q <= PMAX.
Letters: 0 = (a, b) = (0, 0), 1 = (1, 0), 2 = (0, 1), 3 = (1, 1); a = [u -> s], b = [v -> s], u -> v.
The family F1 of THM-4875 (one-apex x one-apex joins with the reversed arc between the apexes) is (02, 23).
"""
import itertools
import sys

import numpy as np
from numba import njit

from ta_fast import FIN24, INIT24, M24, NS
from ta_lattice import rank_frac


def alpha_lang(letters, n):
    cur = INIT24[None, :].copy()
    out = [cur]
    for _ in range(n):
        cur = np.concatenate([cur @ M24[l] for l in letters])
        out.append(cur)
    return out            # out[k]: all words of length k over the alphabet


def beta_lang(letters, n):
    cur = FIN24[None, :].copy()
    out = [cur]
    for _ in range(n):
        cur = np.concatenate([cur @ M24[l].T for l in letters])
        out.append(cur)
    return out


@njit(cache=True)
def mark(al, be, lo, hi, seen):
    mx = 0
    for i in range(al.shape[0]):
        for j in range(be.shape[0]):
            h = 0
            for s in range(al.shape[1]):
                h += al[i, s] * be[j, s]
            if h > mx:
                mx = h
            if lo <= h < hi:
                seen[(h - lo) >> 1] = 1
    return mx


def main():
    pmax = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    alph = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    lo, hi = 1 << 15, 1 << 16
    half = pmax // 2 + 1
    print(f"window [2^15, 2^16): {(hi - lo) // 2} odd N; words L1^p L2^q with p + q <= {pmax}, p, q <= 18")
    print(" L1  L2 | seam rank | max H at p+q=16 (log2) | distinct values at p+q=16 of 2^16*17 words | coverage of the window")
    rows = []
    for L1 in alph:
        A = alpha_lang(L1, min(pmax, 18))
        for L2 in alph:
            B = beta_lang(L2, min(pmax, 18))
            pre = [a for k in range(5) for a in A[k]]
            suf = [b for k in range(5) for b in B[k]]
            G = [[int(a @ b) for b in suf] for a in pre]
            rk = rank_frac(G)
            seen = np.zeros((hi - lo) >> 1, np.uint8)
            mx16 = 0
            vals16 = set()
            for p in range(0, min(pmax, 18) + 1):
                for q in range(0, min(pmax - p, 18) + 1):
                    al = np.ascontiguousarray(A[p])
                    be = np.ascontiguousarray(B[q])
                    m = mark(al, be, lo, hi, seen)
                    if p + q == 16:
                        mx16 = max(mx16, int(m))
            cov = float(seen.mean())
            rows.append((L1, L2, rk, mx16, cov))
            print(f" {L1[0]}{L1[1]}  {L2[0]}{L2[1]} | {rk:9d} | {mx16:10d} ({np.log2(max(mx16, 1)):5.2f}) | {cov:8.4f}",
                  flush=True)
    print()
    print("sorted by coverage:")
    for L1, L2, rk, mx16, cov in sorted(rows, key=lambda r: -r[4]):
        print(f"  ({L1[0]}{L1[1]}, {L2[0]}{L2[1]}): rank {rk}, growth 2^{np.log2(max(mx16, 1)) / 16:.3f} per letter, "
              f"coverage {cov:.4f}")


if __name__ == "__main__":
    main()
