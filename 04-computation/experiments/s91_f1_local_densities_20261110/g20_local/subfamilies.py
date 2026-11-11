#!/usr/bin/env python3
"""S91 helper (g20_local), part B(ii): the eight sub-languages A1* A2* of the two-apex grammar in which the first part
of the word is seen by one apex only and the second part by the other one, and the rank-2 bias of F1 inside G(2,0).

Letters: bit 0 = (a -> vertex), bit 1 = (b -> vertex); a -> b.
   {0,1}*{1,3}*   a sees C1 (C1 -> b), b sees C2 (a -> C2):   ordinal sum T1 => T2, H = x1 x2          (seam rank 1)
   {0,2}*{2,3}*   b sees C1 (C1 -> a), a sees C2 (b -> C2):   F1 (THM-4875), one reversed arc a -> b   (seam rank 3:
                  bilinear rank 2 in (x, y) plus the flag [s = 0])
   and six further families with |C1| or |C2| reversed arcs.
For each family:
   1. the seam rank over Q of the pairing (w1 in A1*, w2 in A2*) -> H(w1 w2)  (exact, from the minimal representation,
      words of length <= 5 on each side);
   2. the exact residue distribution of H over all words of the language of length n (each word once), n = 8, 12, 16,
      20, by a depth-first walk (36-state transfer): class factors for p = 3, 5, 7, 11, 13.
For F1 and the product family, the same distributions are predicted from the exact digit recursion for the law of
(x, y) mod p over h < 2^k (the recursion of box_identity.py) and the bilinear class law; they must agree EXACTLY with
the enumeration for n <= 20 (check), and are then evaluated at n = 40, 80, 160, 320 (floating point) against the limits
     F1:       L(0) = L(-1) = 1 - 1/p^2 + 1/(2p),  other classes 1 - 1/p^2      (sheets s = 0 and s >= 1 have limiting
               word-share 1/2 each)
     product:  L(0) = 2 - 1/p,  other classes 1 - 1/p.

Usage: python3 subfamilies.py [NMAX_BRUTE=20]
"""
import itertools
import json
import os
import sys
from fractions import Fraction

import numpy as np
from numba import njit

import g20_core as g
import minrep
from lattice_tools import rank_q

MOD = 720720
PRIMES = [3, 5, 7, 11, 13]
FAMILIES = [
    ("{0,1}*{1,3}*", (0, 1), (1, 3), "a sees C1, C1 -> b; b sees C2, a -> C2: ordinal sum (product)"),
    ("{0,2}*{2,3}*", (0, 2), (2, 3), "b sees C1, C1 -> a; a sees C2, b -> C2: F1"),
    ("{0,1}*{0,2}*", (0, 1), (0, 2), "a sees C1, C1 -> b; b sees C2, C2 -> a"),
    ("{2,3}*{1,3}*", (2, 3), (1, 3), "a sees C1, b -> C1; b sees C2, a -> C2"),
    ("{2,3}*{0,2}*", (2, 3), (0, 2), "a sees C1, b -> C1; b sees C2, C2 -> a"),
    ("{1,3}*{2,3}*", (1, 3), (2, 3), "b sees C1, a -> C1; a sees C2, b -> C2"),
    ("{0,2}*{0,1}*", (0, 2), (0, 1), "b sees C1, C1 -> a; a sees C2, C2 -> b"),
    ("{1,3}*{0,1}*", (1, 3), (0, 1), "b sees C1, a -> C1; a sees C2, C2 -> b"),
]


@njit(cache=True)
def nb_lang_hist(n, in1, in2, mod, t_to, t_cnt, acc, init):
    """histogram of H mod `mod` over the words of length n of the language A1* A2* (each word once);
    in1[c], in2[c] = 1 iff the letter c belongs to A1, A2."""
    hist = np.zeros(mod, dtype=np.int64)
    st = np.zeros((n + 1, 36), dtype=np.int64)
    for j in range(36):
        st[0, j] = init[j]
    if n == 0:
        hist[1 % mod] += 1
        return hist, 1
    letter = np.zeros(n + 1, dtype=np.int64)
    phase = np.zeros(n + 2, dtype=np.int64)
    d = 0
    words = 0
    while True:
        if letter[d] == 4:
            d -= 1
            if d < 0:
                break
            letter[d] += 1
            continue
        c = letter[d]
        ph = phase[d]
        if ph == 0:
            if in1[c] == 1:
                nph = 0
            elif in2[c] == 1:
                nph = 1
            else:
                letter[d] += 1
                continue
        else:
            if in2[c] == 1:
                nph = 1
            else:
                letter[d] += 1
                continue
        g.nb_step(st[d], c, st[d + 1], t_to, t_cnt)
        if d + 1 == n:
            tot = 0
            for j in range(36):
                if acc[j]:
                    tot += st[n, j]
            hist[tot % mod] += 1
            words += 1
            letter[d] += 1
        else:
            d += 1
            letter[d] = 0
            phase[d] = nph
    return hist, words


def exact_counts(p, kmax):
    """cnt[k][x][y] = #{h < 2^k : x(h) = x, y(h) = y mod p} (Python integers), k = 0..kmax; digit recursion on the
    state (h mod p, h mod 2, q mod p) with x = y - 2 q, y = 2 h + 1."""
    inv2 = pow(2, -1, p)
    cur = {(0, 0, 0): 1}
    out = []
    for k in range(kmax + 1):
        tab = [[0] * p for _ in range(p)]
        for (eta, beta, r), c in cur.items():
            y = (2 * eta + 1) % p
            x = (y - 2 * r) % p
            tab[x][y] += c
        out.append(tab)
        nxt = {}
        for (eta, beta, r), c in cur.items():
            k0 = ((2 * eta) % p, 0, r)
            nxt[k0] = nxt.get(k0, 0) + c
            k1 = ((2 * eta + 1) % p, 1, (r + (eta + beta) * inv2) % p)
            nxt[k1] = nxt.get(k1, 0) + c
        cur = nxt
    return out


def f1_exact_hist(p, n, cnts):
    """exact class counts of N = x1 x2 + 2^s y1 y2 - [s = 0] mod p over the words alpha 2^s beta of length n:
    |alpha| = a (empty or bits(h1) 0, h1 < 2^(a-1)), |beta| = b likewise, a + s + b = n."""
    out = [0] * p
    xs = np.arange(p)
    for a in range(0, n + 1):
        ka = max(a - 1, 0)
        D1 = cnts[ka]
        for b in range(0, n - a + 1):
            s = n - a - b
            kb = max(b - 1, 0)
            D2 = cnts[kb]
            u = pow(2, s, p)
            shift = 1 if s == 0 else 0
            for x1 in range(p):
                for y1 in range(p):
                    w1 = D1[x1][y1]
                    if not w1:
                        continue
                    for x2 in range(p):
                        for y2 in range(p):
                            w2 = D2[x2][y2]
                            if w2:
                                out[(x1 * x2 + u * y1 * y2 - shift) % p] += w1 * w2
    return out


def f1_float_factors(p, n, cnts):
    """class factors of the F1 word family of length n in floating point (vectorised)"""
    idx = {}
    xs = np.arange(p)
    X1, Y1, X2, Y2 = np.meshgrid(xs, xs, xs, xs, indexing="ij")
    base1 = (X1 * X2) % p
    base2 = (Y1 * Y2) % p
    Dn = [np.array(c, dtype=np.float64) / float(1 << k) for k, c in enumerate(cnts)]
    tot = np.zeros(p)
    wsum = 0.0
    for a in range(0, n + 1):
        ka = max(a - 1, 0)
        for b in range(0, n - a + 1):
            s = n - a - b
            kb = max(b - 1, 0)
            u = pow(2, s, p)
            key = (u, s == 0)
            if key not in idx:
                idx[key] = ((base1 + u * base2 - (1 if s == 0 else 0)) % p).ravel()
            # weight of the box relative to 2^n: 2^(ka + kb - n)
            w = 2.0 ** (ka + kb - n)
            outer = np.multiply.outer(Dn[ka], Dn[kb]).ravel()
            tot += w * np.bincount(idx[key], weights=outer, minlength=p)
            wsum += w
    return p * tot / wsum


def product_exact_hist(p, n, cnts):
    """exact class counts of x1 x2 over the words of {0,1}*{1,3}* of length n.
    A word is alpha beta with alpha in {0,1}^a the maximal prefix in {0,1}*, beta in {1,3}^b empty or starting with 3.
    H = H1(alpha) H1(b-bits of beta), H1(word) = x(h) with the word, stripped of trailing 1s, equal to bits(h) 0."""
    # law of H1 over all binary words of length a: strip trailing ones: t trailing ones, then a 0, then free bits
    def h1_counts(length):
        out = [0] * p
        # all ones: H = 1
        out[1 % p] += 1
        for t in range(0, length):       # t trailing ones, a 0 before them, length - t - 1 free leading bits
            k = length - t - 1
            row = cnts[k]
            for x in range(p):
                out[x] += sum(row[x])
        return out
    # beta: b-bits; beta empty, or beta = 3 followed by b-1 letters of {1,3}: b-bit word = 1 followed by b-1 free bits
    def h1_counts_first_one(length):
        if length == 0:
            o = [0] * p
            o[1 % p] = 1
            return o
        # words 1 w' with w' free of length length-1, value x(h) for stripped word
        o = [0] * p
        m = length - 1
        o[1 % p] += 1                      # all ones
        for t in range(0, m):              # t trailing ones, a 0, then the leading part 1 + (m - t - 1 free bits)
            k = m - t - 1
            # h = 2^k + h', h' < 2^k : counts over [2^k, 2^(k+1)) = cnts[k+1] - cnts[k]
            for x in range(p):
                o[x] += sum(cnts[k + 1][x]) - sum(cnts[k][x])
        return o
    out = [0] * p
    for a in range(0, n + 1):
        b = n - a
        ca = h1_counts(a)
        cb = h1_counts_first_one(b)
        for x1 in range(p):
            if ca[x1]:
                for x2 in range(p):
                    if cb[x2]:
                        out[(x1 * x2) % p] += ca[x1] * cb[x2]
    return out


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    ok = True

    def check(name, cond):
        nonlocal ok
        print(("PASS  " if cond else "FAIL  ") + name)
        ok &= bool(cond)

    rep = minrep.build()
    R = rep["R"]
    # 1. seam ranks
    print("1. seam ranks over Q (pairing of first-part words with second-part words, lengths <= 5)")
    seam = {}
    for name, A1, A2, desc in FAMILIES:
        lw = [rep["l_of"](w) for k in range(0, 6) for w in itertools.product(A1, repeat=k)]
        rw = [rep["r_of"](w) for k in range(0, 6) for w in itertools.product(A2, repeat=k)]
        # rank of the matrix (l . r)
        # reduce: basis of span of lw and rw first
        def basis(vs):
            b = []
            for v in vs:
                if rank_q(b + [v]) > len(b):
                    b.append(v)
            return b
        bl, br = basis(lw), basis(rw)
        Gm = [[sum(x * y for x, y in zip(l, r)) for r in br] for l in bl]
        seam[name] = (len(bl), len(br), rank_q(Gm))
        print(f"   {name}: left span {len(bl)}, right span {len(br)}, seam rank {seam[name][2]}   ({desc})")
    check("seam rank of the product family is 1 and of F1 is 3", seam["{0,1}*{1,3}*"][2] == 1 and seam["{0,2}*{2,3}*"][2] == 3)

    # 2. enumeration
    hists = {}
    lens = [n for n in (8, 12, 16, 20) if n <= nmax]
    for name, A1, A2, desc in FAMILIES:
        in1 = np.array([1 if c in A1 else 0 for c in range(4)], dtype=np.int64)
        in2 = np.array([1 if c in A2 else 0 for c in range(4)], dtype=np.int64)
        for n in lens:
            h, words = nb_lang_hist(n, in1, in2, MOD, g.T_TO, g.T_CNT, g.ACC, g.INIT)
            hists[(name, n)] = h
            ok &= words == (n + 2) * 2 ** (n - 1)
    check("every family has (n + 2) 2^(n-1) words of length n", ok)

    def factors(h, p):
        res = np.arange(MOD) % p
        c = np.bincount(res, weights=h.astype(np.float64), minlength=p)
        return p * c / c.sum()

    print("2. class factors p * share over the words of length n (enumeration, exact counts)")
    for name, A1, A2, desc in FAMILIES:
        print(f"   {name}  ({desc})")
        for n in lens:
            line = f"      n = {n:2d}: "
            for p in PRIMES:
                f = factors(hists[(name, n)], p)
                # report class 0, class -1 and the mean of the others
                others = [f[c] for c in range(p) if c not in (0, p - 1)]
                line += f"p={p}: L(0)={f[0]:.4f} L(-1)={f[p - 1]:.4f} rest={np.mean(others):.4f} "
                line += f"[{min(others):.4f},{max(others):.4f}] | "
            print(line)

    # 3. F1 and product: exact predictions
    print("3. F1 and the product family: prediction from the exact law of (x, y) mod p and the bilinear class law")
    allgood = True
    for p in PRIMES:
        cnts = exact_counts(p, 330)
        for n in lens:
            pred = f1_exact_hist(p, n, cnts)
            h = hists[("{0,2}*{2,3}*", n)]
            meas = np.bincount(np.arange(MOD) % p, weights=h.astype(np.float64), minlength=p)
            good = all(int(meas[c]) == pred[c] for c in range(p))
            allgood &= good
            predp = product_exact_hist(p, n, cnts)
            h = hists[("{0,1}*{1,3}*", n)]
            meas = np.bincount(np.arange(MOD) % p, weights=h.astype(np.float64), minlength=p)
            goodp = all(int(meas[c]) == predp[c] for c in range(p))
            allgood &= goodp
        line = f"   p = {p:2d}: F1 class factors (L(0), L(-1), mean of the rest): "
        for n in (20, 40, 80, 160, 320):
            f = f1_float_factors(p, n, cnts)
            others = [f[c] for c in range(p) if c not in (0, p - 1)]
            line += f"n={n}: {f[0]:.4f} {f[p - 1]:.4f} {np.mean(others):.4f} | "
        line += f"limit: {1 - 1 / p ** 2 + 1 / (2 * p):.4f} {1 - 1 / p ** 2 + 1 / (2 * p):.4f} {1 - 1 / p ** 2:.4f}"
        print(line)
    check("F1 and product: predicted class counts = enumerated class counts for every p and n <= 20 (exact integers)",
          allgood)
    for p in PRIMES:
        cnts = exact_counts(p, 330)
        line = f"   p = {p:2d}: product family (L(0), mean of the rest): "
        for n in (20, 80, 320):
            pr = product_exact_hist(p, n, cnts)
            tot = sum(pr)
            f0 = float(Fraction(p * pr[0], tot))
            rest = float(Fraction(p * (tot - pr[0]), tot * (p - 1)))
            line += f"n={n}: {f0:.4f} {rest:.4f} | "
        line += f"limit: {2 - 1 / p:.4f} {1 - 1 / p:.4f}"
        print(line)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
