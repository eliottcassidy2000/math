#!/usr/bin/env python3
"""Independent verifier, part (A)4 of THM-5141 (draft) and every step of its proof.

Notation: eps in {-1,0,1}^w, D = sum eps_i 2^i, z = number of zero digits, weight 4^z.
  W_w(a) = sum of 4^z over eps with D(eps) = a   (so sum_a W_w(a) = 6^w, W_w(0) = 4^w).
  F_j = sum of 4^z over eps with |D| < 2^(w-j),   0 <= j <= w.
Steps checked:
  S1  j-digit signed strings of value 0 / +1 / -1: exact list, by enumeration (pure python), against the claimed
      families; weights 4^j and (4^j - 1)/3.
  S2  for every eps and every j with |D(eps)| < 2^(w-j): the value T of the top j digits has |T| <= 1, and
      |L| <= 2^(w-j) - 1 for the low part (brute force over all 3^w strings, w <= WB).
  S3  F_j <= ((5*4^j - 2)/3) * 6^(w-j) <= (5/3)(2/3)^j 6^w  for all 0 <= j <= w  (exact integers):
      brute force for w <= WB, and an independent exact DP (top-digit recursion on a half array) for w <= WD;
      brute force and DP must agree on the whole distribution W_w for w <= WB.
  S4  dyadic decomposition: every eps != 0 has 1 <= |D| <= 2^w - 1, hence a unique j in [0, w-1] with
      2^(w-j-1) <= |D| < 2^(w-j)  (brute force); D(eps) = 0 only for eps = 0.
  S5  the chain  Q_w <= U_w <= V_w <= G_w < 20 + 15 sqrt 2 < 42  (all divided by 6^w 2^(-w/2)), where
        Q_w = sum_{eps != 0} 4^z |D|^(-1/2),  U_w = sum_j 2^(-(w-j-1)/2) (F_j - F_(j+1)),
        V_w = sum_{j<w} 2^(-(w-j-1)/2) F_j,   G_w = (5 sqrt 2/3) 6^w 2^(-w/2) sum_{j<w} (2 sqrt 2/3)^j.
      20 + 15 sqrt 2 < 42 is checked in exact integers: (42 - 20)^2 = 484 > 450 = (15 sqrt 2)^2.
  S6  (extra, not in the theorem) the exact recursion Q_(w+1) = rho Q_w + (sqrt 2/3) J_w, rho = 2 sqrt 2/3,
      J_w = E (1 + D_w/2^w)^(-1/2), which gives the true size of the constant: lim Q_w = (sqrt 2/3) J/(1 - rho).
Usage: python a4_item4.py [WB (default 14)] [WD (default 24)]
"""
import sys
import time
import itertools
from fractions import Fraction
import numpy as np
from numba import njit
import mpmath as mp

mp.mp.dps = 40


# ---------------------------------------------------------------- S1
def s1_check(jmax):
    ok = True
    for j in range(0, jmax + 1):
        found = {0: [], 1: [], -1: []}
        for e in itertools.product((-1, 0, 1), repeat=j):
            v = sum(e[i] << i for i in range(j))
            if v in found:
                found[v].append(e)
        claim1 = sorted(tuple([-1] * i + [1] + [0] * (j - 1 - i)) for i in range(j))
        claimm1 = sorted(tuple(-x for x in t) for t in claim1)
        wt = lambda lst: sum(4 ** sum(1 for x in e if x == 0) for e in lst)
        c0 = (found[0] == [tuple([0] * j)])
        c1 = (sorted(found[1]) == claim1)
        cm = (sorted(found[-1]) == claimm1)
        w0 = (wt(found[0]) == 4 ** j)
        w1 = (3 * wt(found[1]) == 4 ** j - 1) and (3 * wt(found[-1]) == 4 ** j - 1)
        tot = wt(found[0]) + wt(found[1]) + wt(found[-1])
        wtot = (3 * tot == 5 * 4 ** j - 2) and (3 * tot < 5 * 4 ** j)
        good = c0 and c1 and cm and w0 and w1 and wtot
        ok &= good
        print(f"  j = {j:2d}: value 0: {len(found[0])} string(s) (only the zero string: {c0}, weight 4^j: {w0}); "
              f"value +1: {len(found[1])} (= the family (-1)^i 1 0^(j-1-i), i < j: {c1}); value -1: {len(found[-1])} (negatives: {cm}); "
              f"weights (4^j-1)/3: {w1}; total (5*4^j-2)/3 < (5/3)4^j: {wtot}")
    return ok


# ---------------------------------------------------------------- S2-S4 brute force
@njit(cache=False)
def brute(w):
    """returns Wabs (weight by |D|, index 0..2^w-1), violations of |T| <= 1, violations of |L| bound,
    number of nonzero eps with D = 0, number of eps != 0 with no dyadic j, max |D|."""
    n = 1 << w
    total = 1
    for _ in range(w):
        total *= 3
    Wabs = np.zeros(n, np.int64)
    digits = np.empty(w, np.int64)
    for i in range(w):
        digits[i] = -1
    badT = 0
    badL = 0
    zeroD = 0
    nodyad = 0
    maxD = 0
    for t in range(total):
        D = 0
        z = 0
        nzero = True
        for i in range(w):
            D += digits[i] * (1 << i)
            if digits[i] == 0:
                z += 1
            else:
                nzero = False
        aD = D if D >= 0 else -D
        if aD > maxD:
            maxD = aD
        Wabs[aD] += (1 << (2 * z))
        if D == 0 and not nzero:
            zeroD += 1
        if not nzero:
            # dyadic index
            found = 0
            for j in range(0, w):
                if (1 << (w - j - 1)) <= aD and aD < (1 << (w - j)):
                    found += 1
            if found != 1:
                nodyad += 1
        for j in range(0, w + 1):
            # low part L = digits 0..w-j-1, top part T = digits w-j..w-1 (as a j-digit string)
            L = 0
            for i in range(w - j):
                L += digits[i] * (1 << i)
            T = 0
            for i in range(w - j, w):
                T += digits[i] * (1 << (i - (w - j)))
            if D != (1 << (w - j)) * T + L:
                badL += 1
            aL = L if L >= 0 else -L
            if aL > (1 << (w - j)) - 1:
                badL += 1
            if aD < (1 << (w - j)):
                if T > 1 or T < -1:
                    badT += 1
        i = 0
        while i < w:
            if digits[i] < 1:
                digits[i] += 1
                break
            else:
                digits[i] = -1
                i += 1
    return Wabs, badT, badL, zeroD, nodyad, maxD


# ---------------------------------------------------------------- DP (half array, exact int64 up to 6^24)
def dp_levels(wmax):
    Wh = np.array([1], dtype=np.int64)         # w = 0
    yield 0, Wh
    for w in range(0, wmax):
        n = 1 << w
        new = np.zeros(2 * n, np.int64)
        new[:n] = 4 * Wh
        if n > 1:
            new[1:n] += Wh[1:][::-1]           # W_w(a - 2^w) = W_w(2^w - a), 0 < a < 2^w
        new[n:] = Wh                           # W_w(a - 2^w), 2^w <= a < 2^(w+1)
        Wh = new
        yield w + 1, Wh


def F_values(Wh, w):
    """exact F_j, j = 0..w, as python ints, from the half array."""
    cs = np.cumsum(Wh, dtype=np.int64)         # cs[a] = sum_{0<=b<=a} Wh[b]
    out = []
    for j in range(0, w + 1):
        m = 1 << (w - j)                       # |a| < m
        half = int(cs[m - 1])                  # a = 0..m-1
        out.append(2 * half - int(Wh[0]))
    return out


def main():
    WB = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    WD = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    ok = True
    print("=== S1: signed strings of value 0, +1, -1 ===")
    ok &= s1_check(12)

    print("=== exact constant: 20 + 15 sqrt 2 < 42 ===")
    c = (42 - 20) ** 2 > 15 * 15 * 2
    ok &= c
    rho = 2 * mp.sqrt(2) / 3
    closed = (5 * mp.sqrt(2) / 3) / (1 - rho)
    print(f"  (5 sqrt2/3) sum_j (2 sqrt2/3)^j = (5 sqrt2/3)/(1 - 2 sqrt2/3) = {mp.nstr(closed, 12)}; "
          f"20 + 15 sqrt2 = {mp.nstr(20 + 15 * mp.sqrt(2), 12)}; identity holds (to 1e-30): {abs(closed - (20 + 15 * mp.sqrt(2))) < mp.mpf(10) ** -30}; "
          f"(42-20)^2 = 484 > 450 = (15 sqrt2)^2: {c}")

    print("=== S2-S4 brute force over all 3^w strings ===")
    brute_W = {}
    for w in range(1, WB + 1):
        t0 = time.time()
        Wabs, badT, badL, zeroD, nodyad, maxD = brute(w)
        brute_W[w] = Wabs
        tot = int(Wabs.sum())
        good = (badT == 0 and badL == 0 and zeroD == 0 and nodyad == 0 and maxD == (1 << w) - 1 and tot == 6 ** w and int(Wabs[0]) == 4 ** w)
        ok &= good
        print(f"  w = {w:2d}: total weight 6^w: {tot == 6 ** w}; W(0) = 4^w (D = 0 only for eps = 0: {zeroD == 0}); max |D| = 2^w - 1: {maxD == (1 << w) - 1}; "
              f"|T| <= 1 whenever |D| < 2^(w-j) (violations {badT}); D = 2^(w-j) T + L with |L| <= 2^(w-j) - 1 (violations {badL}); "
              f"eps != 0 without a unique dyadic j in [0, w-1]: {nodyad}   [{time.time()-t0:.1f}s]", flush=True)

    print("=== S3, S5, S6: exact DP ===")
    sq2 = mp.sqrt(2)
    prevQ = None
    prevJ = None
    worst_Fratio = Fraction(0)
    worst_Fj = None
    maxQ = mp.mpf(0)
    for w, Wh in dp_levels(WD):
        if w == 0:
            prevQ = mp.mpf(0)
            prevJ = mp.mpf(1)
            continue
        t0 = time.time()
        n = 1 << w
        tot = int(Wh[0]) + 2 * int(Wh[1:].sum(dtype=np.int64))
        c_tot = (tot == 6 ** w)
        c_brute = True
        if w in brute_W:
            # brute force weight by |D|: Wabs[0] = W(0), Wabs[a] = 2 W(a)
            bw = brute_W[w]
            c_brute = (int(bw[0]) == int(Wh[0])) and bool(np.all(bw[1:] == 2 * Wh[1:]))
        F = F_values(Wh, w)
        c_F0 = (F[0] == 6 ** w) and (F[w] == 4 ** w)
        c_F = True
        c_Fmid = True
        for j in range(0, w + 1):
            # F_j <= ((5*4^j - 2)/3) 6^(w-j)  and  <= (5/3)(2/3)^j 6^w   (exact)
            if 3 * F[j] > (5 * 4 ** j - 2) * 6 ** (w - j):
                c_Fmid = False
            if F[j] * 3 ** (j + 1) > 5 * 2 ** j * 6 ** w:
                c_F = False
            r = Fraction(F[j] * 3 ** j, 2 ** j * 6 ** w)
            if r > worst_Fratio:
                worst_Fratio = r
                worst_Fj = (w, j)
        # Q_w (float64 pairwise sum, then mp), U_w, V_w, G_w
        Qf = 0.0
        Jf = float(Wh[0])
        CH = 1 << 20
        for lo in range(1, n, CH):
            hi = min(n, lo + CH)
            a = np.arange(lo, hi, dtype=np.float64)
            wv = Wh[lo:hi].astype(np.float64)
            Qf += 2.0 * float(np.sum(wv / np.sqrt(a)))
            x = a / float(n)
            Jf += float(np.sum(wv * (1.0 / np.sqrt(1.0 + x) + 1.0 / np.sqrt(1.0 - x))))
        Jf /= float(6 ** w)
        norm = mp.mpf(6) ** w * mp.mpf(2) ** (-mp.mpf(w) / 2)
        Q = mp.mpf(Qf) / norm
        U = sum(mp.mpf(2) ** (-mp.mpf(w - j - 1) / 2) * (F[j] - F[j + 1]) for j in range(w)) / norm
        V = sum(mp.mpf(2) ** (-mp.mpf(w - j - 1) / 2) * F[j] for j in range(w)) / norm
        G = (5 * sq2 / 3) * sum(rho ** j for j in range(w))
        c_chain = (Q <= U * (1 + mp.mpf(10) ** -12)) and (U <= V) and (V <= G) and (G < 20 + 15 * sq2) and (G < 42)
        # J_w and the recursion
        J = mp.mpf(Jf)
        rec = rho * prevQ + (sq2 / 3) * prevJ
        c_rec = abs(rec - Q) < mp.mpf(10) ** -9 * max(1, Q)
        good = c_tot and c_brute and c_F0 and c_F and c_Fmid and c_chain and c_rec
        ok &= good
        maxQ = max(maxQ, Q)
        print(f"  w = {w:2d}: sum W = 6^w: {c_tot}; DP = brute force: {c_brute if w in brute_W else 'n/a'}; F_0 = 6^w, F_w = 4^w: {c_F0}; "
              f"F_j <= ((5*4^j-2)/3) 6^(w-j) all j: {c_Fmid}; F_j <= (5/3)(2/3)^j 6^w all j: {c_F}; "
              f"Q = {mp.nstr(Q, 9)} <= U = {mp.nstr(U, 7)} <= V = {mp.nstr(V, 7)} <= G = {mp.nstr(G, 7)} < 41.2132 < 42: {c_chain}; "
              f"J_w = {mp.nstr(J, 9)}; recursion Q_w = rho Q_(w-1) + (sqrt2/3) J_(w-1): {c_rec}   [{time.time()-t0:.1f}s]", flush=True)
        prevQ, prevJ = Q, J
    lim = (sq2 / 3) * prevJ / (1 - rho)
    print(f"  largest F_j (3/2)^j / 6^w over all j <= w <= {WD}: {float(worst_Fratio):.6f} at (w, j) = {worst_Fj}  (bound 5/3 = 1.666667)")
    print(f"  largest Q_w for w <= {WD}: {mp.nstr(maxQ, 9)}  (theorem's constant: 42; proof's constant: 41.2132)")
    print(f"  limit of Q_w (from the recursion, J = J_{WD}): (sqrt2/3) J / (1 - 2 sqrt2/3) = {mp.nstr(lim, 9)}")
    print("ALL CHECKS PASSED" if ok else "*** SOME CHECK FAILED ***")


if __name__ == "__main__":
    main()
