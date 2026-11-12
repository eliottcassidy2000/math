#!/usr/bin/env python3
"""Independent verifier, parts (A)2, (A)3 and (A)4 (direct form) of THM-5141 (draft).

For each w the full table  r_w(a, b) = #{(u, u') in [0,2^w)^2 : u - u' = a, rev_w u - rev_w u' = b}
is built DIRECTLY from the definition (all 4^w ordered pairs, no signed strings involved).  Then:
  (i)    number of nonzero cells = 3^w                      [image of an injective map on {-1,0,1}^w]
  (ii)   sum r = 4^w                                         [sanity]
  (iii)  sum r^2 = 6^w                                       [(A)3]
  (iv)   every nonzero r is a power of 2, and the number of cells with r = 2^z is C(w, z) 2^(w-z)
  (v)    for every eps in {-1,0,1}^w (independent enumeration): r(D(eps), D*(eps)) = 2^z(eps), each cell is hit by
         exactly one eps, and no nonzero cell is left over   [(A)2 in full: injective, r = 2^z, r = 0 off the image]
  (vi)   (A)4 directly: sum_{a != 0} sum_b r(a,b)^2 |a|^(-1/2)  compared with 42 * 6^w * 2^(-w/2)
         (row sums of r^2 exact; the final weighted sum in 50-digit mpmath).
Separately, injectivity of eps -> (D, D*) for larger w by sorting the 3^w keys.
Usage: python a2_r_identity.py [WMAX table (default 12)] [WMAX sort (default 15)]
"""
import sys
import time
from math import comb
import numpy as np
from numba import njit
import mpmath as mp

mp.mp.dps = 50


@njit(cache=False)
def fill_table(w, T):
    n = 1 << w
    rev = np.zeros(n, np.int64)
    for x in range(n):
        r = 0
        for i in range(w):
            if (x >> i) & 1:
                r |= 1 << (w - 1 - i)
        rev[x] = r
    for u in range(n):
        ru = rev[u]
        for v in range(n):
            T[u - v + n - 1, ru - rev[v] + n - 1] += 1


@njit(cache=False)
def table_stats(w, T):
    n = 1 << w
    size = 2 * n - 1
    nnz = 0
    s1 = 0
    s2 = 0
    rowsq = np.zeros(size, np.int64)
    hist = np.zeros(w + 1, np.int64)      # hist[z] = number of cells with r = 2^z
    notpow2 = 0
    for i in range(size):
        acc = 0
        for j in range(size):
            v = np.int64(T[i, j])
            if v != 0:
                nnz += 1
                s1 += v
                acc += v * v
                if v & (v - 1) != 0:
                    notpow2 += 1
                else:
                    z = 0
                    while (1 << z) < v:
                        z += 1
                    hist[z] += 1
        rowsq[i] = acc
        s2 += acc
    return nnz, s1, s2, rowsq, hist, notpow2


@njit(cache=False)
def check_eps(w, T):
    """for every eps: T[cell] must equal 2^z; then zero the cell (a second hit is detected as a mismatch).
    returns (mismatches, leftover nonzero cells)."""
    n = 1 << w
    size = 2 * n - 1
    total = 1
    for _ in range(w):
        total *= 3
    digits = np.empty(w, np.int64)
    D = 0
    Ds = 0
    z = 0
    for i in range(w):
        digits[i] = -1
        D -= (1 << i)
        Ds -= (1 << (w - 1 - i))
    bad = 0
    for _ in range(total):
        v = np.int64(T[D + n - 1, Ds + n - 1])
        if v != (1 << z):
            bad += 1
        T[D + n - 1, Ds + n - 1] = 0
        i = 0
        while i < w:
            if digits[i] < 1:
                if digits[i] == -1:
                    z += 1          # -1 -> 0
                else:
                    z -= 1          # 0 -> 1
                digits[i] += 1
                D += (1 << i)
                Ds += (1 << (w - 1 - i))
                break
            else:
                digits[i] = -1      # 1 -> -1 : z unchanged
                D -= 2 * (1 << i)
                Ds -= 2 * (1 << (w - 1 - i))
                i += 1
    left = 0
    for i in range(size):
        for j in range(size):
            if T[i, j] != 0:
                left += 1
    return bad, left


@njit(cache=False)
def eps_keys(w):
    total = 1
    for _ in range(w):
        total *= 3
    keys = np.empty(total, np.int64)
    digits = np.empty(w, np.int64)
    D = 0
    Ds = 0
    for i in range(w):
        digits[i] = -1
        D -= (1 << i)
        Ds -= (1 << (w - 1 - i))
    off = 1 << w
    for t in range(total):
        keys[t] = (D + off) * (off << 2) + (Ds + off)
        i = 0
        while i < w:
            if digits[i] < 1:
                digits[i] += 1
                D += (1 << i)
                Ds += (1 << (w - 1 - i))
                break
            else:
                digits[i] = -1
                D -= 2 * (1 << i)
                Ds -= 2 * (1 << (w - 1 - i))
                i += 1
    return keys


def main():
    WT = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    WS = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    ok = True
    print("=== direct table r_w(a,b) from all 4^w pairs ===")
    worst_ratio = mp.mpf(0)
    for w in range(1, WT + 1):
        t0 = time.time()
        n = 1 << w
        size = 2 * n - 1
        T = np.zeros((size, size), np.uint16 if w <= 15 else np.uint32)
        fill_table(w, T)
        nnz, s1, s2, rowsq, hist, notpow2 = table_stats(w, T)
        c_i = (nnz == 3 ** w)
        c_ii = (s1 == 4 ** w)
        c_iii = (s2 == 6 ** w)
        c_iv = (notpow2 == 0) and all(int(hist[z]) == comb(w, z) * 2 ** (w - z) for z in range(w + 1))
        # (A)4 direct
        tot = mp.mpf(0)
        for i in range(size):
            a = i - (n - 1)
            if a != 0 and rowsq[i]:
                tot += mp.mpf(int(rowsq[i])) / mp.sqrt(abs(a))
        ratio = tot / (mp.mpf(6) ** w * mp.mpf(2) ** (-mp.mpf(w) / 2))
        worst_ratio = max(worst_ratio, ratio)
        c_vi = ratio <= 42
        row0 = int(rowsq[n - 1])            # a = 0 row: only (0,0) should contribute: 4^w
        c_row0 = (row0 == 4 ** w)
        bad, left = check_eps(w, T)
        c_v = (bad == 0 and left == 0)
        allok = c_i and c_ii and c_iii and c_iv and c_v and c_vi and c_row0
        ok &= allok
        print(f"  w = {w:2d}: nonzero cells = {nnz} (3^w: {c_i}); sum r = 4^w: {c_ii}; sum r^2 = {s2} (= 6^w: {c_iii}); "
              f"powers of 2 with binomial histogram: {c_iv}; r(D,D*) = 2^z for all eps, each cell hit once, nothing left: {c_v} "
              f"(mismatches {bad}, leftover {left}); row a = 0 carries exactly 4^w: {c_row0}; "
              f"(A)4 sum / (6^w 2^(-w/2)) = {mp.nstr(ratio, 8)} (<= 42: {c_vi})   [{time.time()-t0:.1f}s]", flush=True)
        del T
    print(f"  largest (A)4 ratio over w <= {WT}: {mp.nstr(worst_ratio, 8)}")
    print("=== injectivity of eps -> (D, D*) by sorting all 3^w keys ===")
    for w in range(1, WS + 1):
        t0 = time.time()
        k = eps_keys(w)
        k.sort()
        inj = bool(np.all(k[1:] != k[:-1])) if k.size > 1 else True
        ok &= inj
        print(f"  w = {w:2d}: {k.size} strings, all images distinct: {inj}   [{time.time()-t0:.1f}s]", flush=True)
        del k
    print("ALL CHECKS PASSED" if ok else "*** SOME CHECK FAILED ***")


if __name__ == "__main__":
    main()
