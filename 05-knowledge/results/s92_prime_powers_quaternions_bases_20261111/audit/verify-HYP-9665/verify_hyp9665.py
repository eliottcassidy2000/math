#!/usr/bin/env python3
"""Verifier of HYP-9665 after the audit corrections (S92, 2026-11-11). Read-only checks.

V1  the bound 2 R_m <= R_prod <= 2 R_m + 2 F R_flat (exact rational arithmetic on the audit histograms
    hypotheses/even_moduli_hist_E.json, E = 20..28), and every number the file quotes from those histograms;
V2  Lemma M: T T^* = 2 I for an even modulus and an odd character; T computes the character sum of q;
V3  the 2-adic algebra of (A): N = const - 2 x1 q2 (mod 2^j) on a cell with j - 1 bottom digits prescribed;
    joint counts of (q mod 2^(j-1) m, h mod m) on such cells;
V4  the algebra of (C): N_b = 1 + alpha_1 b + alpha_2 b^2; alpha_1 = const - q2 (mod 2^a) on a cell with a bottom
    digits prescribed; the affine form modulo an odd b; the mixed case;
V5  (D): the states of type p as the language of an automaton on (q mod p, h mod 2p); first states of type p; the
    primes below 400 whose fan is the first state of type p; divergence of the uncentred sum; partial sums of A_p;
V6  (B): the mean of S_sigma over a class modulo m is L_m(c - c_sigma); arithmetic of the quoted ratios;
V7  the file itself: front matter, paths, labels, notation, ASCII, backticks.
Exit status 1 on failure.
"""
import json
import math
import os
import re
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(AUDIT, "..", "..", "..", ".."))
HYP = os.path.join(ROOT, "05-knowledge", "hypotheses",
                   "HYP-9665-f1-densities-beyond-odd-moduli-even-moduli-single-values-primes-of-the-base.md")
OK = True


def check(cond, text):
    global OK
    OK &= bool(cond)
    print(("  ok   " if cond else "  FAIL ") + text)


def build(n):
    """psi, q for h < n (q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2))."""
    q = np.zeros(n, dtype=np.int64)
    k = 1
    while (1 << (k - 1)) < n:
        lo, hi = 1 << (k - 1), min(n, 1 << k)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
        k += 1
    h = np.arange(n, dtype=np.int64)
    return h - q, q


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def L_m(m, e):
    g = math.gcd(m, e % m) if e % m else m
    tot = Fraction(0)
    for d in divisors(g):
        md = m // d
        pr = Fraction(1)
        for p in (3, 5, 7, 11, 13):
            if md % p == 0:
                pr *= 1 - Fraction(1, p * p)
        tot += pr / d
    return tot


# ---------------------------------------------------------------------------------------------------------------
def v1():
    print("V1 even moduli: the audit histograms (exact rational arithmetic)")
    DIV315 = divisors(315)
    csig = {0: -1, 1: 0}
    moduli = ((2, 3), (3, 3), (2, 9), (3, 9), (2, 5), (3, 5), (2, 7), (3, 7), (2, 15), (3, 15), (2, 45), (3, 45), (3, 315))
    table = {}
    worst_low, worst_up, loose_fail = None, None, []
    for E in (20, 22, 24, 26, 28):
        hist = np.array(json.load(open(os.path.join(AUDIT, "hypotheses", "even_moduli_hist_%d.json" % E)))["hist"],
                        dtype=object)
        for (j, m) in moduli:
            Mod = (1 << j) * m
            for sh in (0, 1):
                T = int(sum(int(v) for row in hist[sh] for v in row))
                cnt = [0] * Mod
                cntm = [0] * m
                D = {d: 0 for d in divisors(m)}
                for di, dd in enumerate(DIV315):
                    row = hist[sh][di]
                    D[math.gcd(dd, m)] += int(sum(int(v) for v in row))
                    for c in range(2520):
                        v = int(row[c])
                        if v:
                            cnt[c % Mod] += v
                            cntm[c % m] += v
                pred_m = {}
                for c in range(m):
                    e = (c - csig[sh]) % m
                    gg = math.gcd(m, e) if e else m
                    pred_m[c] = sum(Fraction(d, m) * D[d] for d in divisors(gg))
                Rm = max(Fraction(m) * abs(cntm[c] - pred_m[c]) / T for c in range(m))
                F = max(Fraction(m * cntm[c], T) for c in range(m))
                Rprod = Fraction(0)
                Rflat = Fraction(0)
                for c in range(1, Mod, 2):
                    Rprod = max(Rprod, Fraction(Mod) * abs(cnt[c] - pred_m[c % m] / (1 << (j - 1))) / T)
                    Rflat = max(Rflat, abs(Fraction((1 << (j - 1)) * cnt[c], cntm[c % m]) - 1))
                assert all(cnt[c] == 0 for c in range(0, Mod, 2))
                low = Rprod - 2 * Rm
                up = 2 * Rm + 2 * F * Rflat - Rprod
                if worst_low is None or low < worst_low[0]:
                    worst_low = (low, Mod, E, sh)
                if worst_up is None or up < worst_up[0]:
                    worst_up = (up, Mod, E, sh)
                if Rprod > 2 * Rm + Rflat:
                    loose_fail.append((Mod, E, sh, float(Rprod), float(2 * Rm + Rflat)))
                table[(Mod, E, sh)] = (float(Rprod), float(Rm), float(Rflat), float(F), T, j, m, cnt, cntm)
    check(worst_low[0] >= 0, "2 R_m <= R_prod at the 13 moduli, X = 2^20..2^28, both sheets (130 cases); least slack "
          "%.6f at %d, 2^%d, sheet %d" % (float(worst_low[0]), worst_low[1], worst_low[2], worst_low[3]))
    check(worst_up[0] >= 0, "R_prod <= 2 R_m + 2 F R_flat in the same 130 cases; least slack %.2e at %d, 2^%d, sheet %d"
          % (float(worst_up[0]), worst_up[1], worst_up[2], worst_up[3]))
    print("    the looser bound R_prod <= 2 R_m + R_flat fails in %d of the 130 cases, e.g. " % len(loose_fail)
          + "; ".join("%d at 2^%d sheet %d: %.4f > %.4f" % t for t in loose_fail[:4]))
    check(any(t[0] == 20 and t[1] == 22 and t[2] == 0 for t in loose_fail),
          "it fails at 20, 2^22, sheet 0 (the case quoted by the editor)")

    def fmt4(Mod, E):
        return "%.4f / %.4f" % (table[(Mod, E, 0)][2], table[(Mod, E, 1)][2])

    quoted = {12: ["0.0021 / 0.0014", "0.0009 / 0.0005", "0.0005 / 0.0001", "0.0002 / 0.0002"],
              20: ["0.0028 / 0.0018", "0.0009 / 0.0007", "0.0006 / 0.0003", "0.0001 / 0.0002"],
              36: ["0.0040 / 0.0024", "0.0020 / 0.0011", "0.0005 / 0.0004", "0.0005 / 0.0005"],
              72: ["0.0053 / 0.0048", "0.0036 / 0.0019", "0.0013 / 0.0009", "0.0006 / 0.0006"]}
    for Mod in (12, 20, 36, 72):
        got = [fmt4(Mod, E) for E in (22, 24, 26, 28)]
        check(got == quoted[Mod], "R_flat at %d, X = 2^22..2^28: %s" % (Mod, ", ".join(got)))
    scales = []
    for Mod in (12, 20, 36, 72):
        for E in (22, 24, 26, 28):
            for sh in (0, 1):
                Rprod, Rm, Rflat, F, T, j, m, _, _ = table[(Mod, E, sh)]
                scales.append((Rflat / math.sqrt((1 << (j - 1)) * m / T), Mod, E, sh))
    lo, hi = min(scales), max(scales)
    print("    R_flat in sampling scales sqrt(2^(j-1) m/#T), moduli 12, 20, 36, 72, X = 2^22..2^28: min %.2f (%d, 2^%d, "
          "sheet %d), max %.2f (%d, 2^%d, sheet %d)" % (lo + hi))
    check(round(lo[0], 1) == 0.6 and round(hi[0], 1) == 4.0, "between 0.6 and 4.0 sampling scales at these sizes")
    R, _, Rf, _, T, j, m, _, _ = table[(72, 20, 0)]
    check(round(Rf / math.sqrt((1 << (j - 1)) * m / T), 1) == 4.3, "4.3 sampling scales at 72 at 2^20 (sheet 0): %.2f"
          % (Rf / math.sqrt((1 << (j - 1)) * m / T)))
    falls = []
    for Mod in (12, 20, 36, 72):
        for sh in (0, 1):
            falls.append((table[(Mod, 24, sh)][2] / table[(Mod, 28, sh)][2], Mod, sh,
                          table[(Mod, 24, sh)][2], table[(Mod, 28, sh)][2]))
    print("    fall of R_flat from 2^24 to 2^28: " + "; ".join("%d sheet %d: %.2f (%.5f -> %.5f)" % (t[1], t[2], t[0], t[3], t[4])
                                                              for t in falls))
    lo, hi = min(falls), max(falls)
    check(lo[1] == 36 and lo[2] == 1 and round(lo[0], 1) == 2.1, "least fall 2.1 at 36, sheet 1 (%.5f -> %.5f)" % (lo[3], lo[4]))
    check(hi[1] == 20 and hi[2] == 0 and 11.5 < hi[0] < 13.5, "largest fall about 12 at 20, sheet 0 (%.5f -> %.5f; exact ratio %.2f)"
          % (hi[3], hi[4], hi[0]))
    check(all(t[0] >= 2 for t in falls), "the factor 2 is reached at 12, 20, 36, 72 on both sheets")
    v36 = {}
    for E in (26, 28):
        for sh in (0, 1):
            Rprod, Rm, Rflat, F, T, j, m, _, _ = table[(36, E, sh)]
            v36[(E, sh)] = ("%.5f" % Rflat, "%.1f" % (Rflat / math.sqrt((1 << (j - 1)) * m / T)))
    check([v36[k][0] for k in ((26, 0), (26, 1), (28, 0), (28, 1))] == ["0.00052", "0.00045", "0.00050", "0.00051"]
          and [v36[k][1] for k in ((26, 0), (26, 1), (28, 0), (28, 1))] == ["1.7", "1.7", "3.5", "4.0"],
          "at 36 R_flat is stationary from 2^26 to 2^28: %s / %s -> %s / %s, that is %s / %s -> %s / %s sampling scales"
          % (v36[(26, 0)][0], v36[(26, 1)][0], v36[(28, 0)][0], v36[(28, 1)][0],
             v36[(26, 0)][1], v36[(26, 1)][1], v36[(28, 0)][1], v36[(28, 1)][1]))
    samp = [math.sqrt(table[(12, 28, sh)][4] / table[(12, 24, sh)][4]) for sh in (0, 1)]
    check(all(round(v, 1) == 4.3 for v in samp), "a pure sampling term falls by sqrt(#T(2^28)/#T(2^24)) = %.2f / %.2f" % tuple(samp))
    rp = ["%.3f" % table[(12, E, 0)][0] for E in (22, 24, 26, 28)]
    rm = ["%.4f" % table[(12, E, 0)][1] for E in (22, 24, 26, 28)]
    check(rp == ["0.060", "0.028", "0.036", "0.034"], "R_prod at 12, sheet 0: " + ", ".join(rp))
    rm3 = ["%.3f" % table[(12, E, 0)][1] for E in (22, 24, 26, 28)]
    print("    R_m at 3, sheet 0, four decimals: " + ", ".join(rm) + "; three decimals: " + ", ".join(rm3))
    check(rm3 == ["0.028", "0.014", "0.018", "0.017"], "remainder at 3, sheet 0, to three decimals: 0.028, 0.014, 0.018, 0.017")
    # class factors at 12 at 2^28
    out = []
    for sh in (0, 1):
        _, _, _, _, T, j, m, cnt, cntm = table[(12, 28, sh)]
        a = [12 * cnt[c] / T for c in range(1, 12, 2) if (c - csig[sh]) % 3 != 0]
        b = [12 * cnt[c] / T for c in range(1, 12, 2) if (c - csig[sh]) % 3 == 0]
        out.append((sum(a) / len(a), sum(b) / len(b), len(a), len(b)))
    check(["%.4f" % out[0][0], "%.4f" % out[1][0]] == ["1.6996", "1.6936"] and out[0][2] == 4,
          "class factors at 12, 2^28, four classes with 3 not dividing N - c_sigma: %.4f / %.4f" % (out[0][0], out[1][0]))
    check(["%.4f" % out[0][1], "%.4f" % out[1][1]] == ["2.6008", "2.6128"] and out[0][3] == 2,
          "the two other classes: %.4f / %.4f" % (out[0][1], out[1][1]))
    check("%.4f" % float(2 * L_m(3, 1)) == "1.7778" and "%.4f" % float(2 * L_m(3, 0)) == "2.4444",
          "limits 2 L_3 = %.4f and %.4f" % (float(2 * L_m(3, 1)), float(2 * L_m(3, 0))))
    # the withdrawn falsifier, on the proved odd law (census of the session, cum_24.out and cum_28.out are not parsed
    # here; the audit histograms give 9 directly)
    for (m, sh) in ((9, 0),):
        a = table[(36, 24, sh)][1]
        b = table[(36, 28, sh)][1]
        check(b > a, "remainder of THM-5134 (B) at 9, sheet 0: %.4f at 2^24 -> %.4f at 2^28 (does not decrease)" % (a, b))
    a, b = table[(12, 24, 0)][0], table[(12, 28, 0)][0]
    check(b > a, "remainder of the product law at 12, sheet 0: %.4f at 2^24 -> %.4f at 2^28 (the old falsifier was met)" % (a, b))


# ---------------------------------------------------------------------------------------------------------------
def v2(q):
    print("V2 Lemma M: an even modulus M2 and an odd character k1")
    good = True
    radii = []
    for M2 in (2, 4, 6, 8, 10, 12, 14, 18, 20, 30, 36):
        for k1 in range(1, M2, 2):
            n = 2 * M2
            T = np.zeros((n, n), dtype=complex)
            for w in range(n):
                c = (w + 1) // 2
                T[w, (2 * w) % n] += 1
                T[w, (2 * w + 1) % n] += np.exp(2j * np.pi * k1 * c / M2)
            good &= np.allclose(T @ T.conj().T, 2 * np.eye(n), atol=1e-10)
            radii.append(max(abs(np.linalg.eigvals(T))) / 2)
            # T computes the character sum of q over h < 2^nn
            nn = 10
            v = np.zeros(n, dtype=complex)
            v[0] = 1
            for _ in range(nn):
                v = v @ T
            direct = np.exp(2j * np.pi * k1 * q[: 1 << nn] / M2).sum()
            good &= abs(v.sum() - direct) < 1e-6
            good &= abs(direct) <= 2 ** ((nn + 1) / 2) * math.sqrt(n) + 1e-9
    check(good, "T T^* = 2 I for M2 in {2, ..., 36} even and every odd k1; e_0^T T^10 1 = sum_{h < 2^10} e(k1 q(h)/M2)")
    check(max(abs(r - 2 ** -0.5) for r in radii) < 1e-9, "rho(T)/2 = 2^(-1/2) = 0.7071 in every case (%d operators)" % len(radii))
    # an even character whose k is even but not divisible by the full power of 2: reduction k/M1 = k1/M2
    red = True
    for (jm1, m) in ((1, 3), (2, 3), (3, 3), (2, 5), (2, 9), (1, 15)):
        M1 = (1 << jm1) * m
        for k in range(1, M1):
            a = (k & -k).bit_length() - 1
            even_order = (M1 // math.gcd(k, M1)) % 2 == 0
            red &= even_order == (k % (1 << jm1) != 0)
            if even_order:
                k1, M2 = k >> a, M1 >> a
                red &= k1 % 2 == 1 and M2 % 2 == 0 and Fraction(k, M1) == Fraction(k1, M2)
            else:
                red &= (k >> jm1) << jm1 == k
    check(red, "a character k of Z/(2^(j-1) m) has even order iff 2^(j-1) does not divide k; then k/(2^(j-1) m) = k1/M2 "
          "with k1 odd, M2 even; otherwise it is a character of Z/m")

    # rates of the moduli M of q (operators on the states w mod 2M, as in the file), by parity of the order
    def op(k, M):
        n = 2 * M
        A = np.zeros((n, n), dtype=complex)
        for w in range(n):
            A[w, (2 * w) % n] += 1
            A[w, (2 * w + 1) % n] += np.exp(2j * np.pi * k * ((w + 1) // 2) / M)
        return A

    def rate(M):
        ev, od = 0.0, 0.0
        for k in range(1, M):
            r = max(abs(np.linalg.eigvals(op(k, M)))) / 2
            if (M // math.gcd(k, M)) % 2 == 0:
                ev = max(ev, r)
            else:
                od = max(od, r)
        return ev, od

    odd_rate = {3: rate(3)[1], 5: rate(5)[1], 9: rate(9)[1], 7: rate(7)[1], 15: rate(15)[1], 45: rate(45)[1]}
    want = {6: "0.8660", 12: "0.8660", 24: "0.8660", 18: "0.8660", 36: "0.8660", 10: "0.7418", 20: "0.7418", 40: "0.7418",
            14: "0.9320", 28: "0.9320", 30: "0.8922", 60: "0.8922", 90: "0.8922"}
    good = True
    rows = []
    for M in (6, 12, 24, 10, 20, 40, 18, 36, 14, 28, 30, 60, 90):
        ev, od = rate(M)
        mo = M
        while mo % 2 == 0:
            mo //= 2
        good &= abs(ev - 2 ** -0.5) < 1e-6 and "%.4f" % max(ev, od) == want[M] and abs(od - odd_rate[mo]) < 1e-6
        rows.append("%d: %.4f" % (M, max(ev, od)))
    check(good, "M = 6, ..., 90 (13 even moduli): every character of even order has rho(T_k)/2 = 0.7071; the rate of M is the "
          "rate of its odd part: " + ", ".join(rows))
    check(all(v > 2 ** -0.5 for v in odd_rate.values()), "the rates of the odd parts 3, 5, 9, 7, 15, 45 (%s) are above 2^(-1/2)"
          % ", ".join("%.4f" % odd_rate[mm] for mm in (3, 5, 9, 7, 15, 45)))


# ---------------------------------------------------------------------------------------------------------------
def v3(psi, q):
    print("V3 the 2-adic algebra of (A)")
    rng = np.random.default_rng(9665)
    good = True
    for _ in range(4000):
        j = int(rng.integers(2, 6))
        h1 = int(rng.integers(0, 2000))
        s = int(rng.integers(0, 6))
        low = int(rng.integers(0, 1 << (j - 1)))            # the j - 1 bottom digits of h2
        top = rng.integers(1, 1 << 12, size=8)
        x1, y1 = 2 * int(psi[h1]) + 1, 2 * h1 + 1
        const = None
        for t in top:
            h2 = (int(t) << (j - 1)) | low
            x2, y2, q2 = 2 * int(psi[h2]) + 1, 2 * h2 + 1, int(q[h2])
            N = x1 * x2 + (1 << s) * y1 * y2 - (1 if s == 0 else 0)
            good &= x2 == y2 - 2 * q2 and N % 2 == 1
            c = (N + 2 * x1 * q2) % (1 << j)
            const = c if const is None else const
            good &= c == const and y2 % (1 << j) == (2 * low + 1) % (1 << j)
    check(good, "x2 = y2 - 2 q2; on a cell with j - 1 bottom digits prescribed y2 mod 2^j is constant and "
          "N = const - 2 x1 q2 (mod 2^j), N odd (4000 random cells, j = 2..5)")
    bij = all(sorted((-2 * x1 * t) % (1 << j) for t in range(1 << (j - 1))) == list(range(0, 1 << j, 2))
              for j in range(2, 7) for x1 in range(1, 64, 2))
    check(bij, "q2 -> -2 x1 q2 is a bijection of Z/2^(j-1) onto 2Z/2^j for odd x1")
    # joint counts of (q mod 2^(j-1) m, h mod m) on the cells of B_n with j - 1 bottom digits prescribed
    for (j, m) in ((2, 3), (3, 3), (2, 5)):
        M1 = (1 << (j - 1)) * m
        row = []
        for n in (14, 18, 20):
            hs = np.arange(1 << (n - 1), 1 << n, dtype=np.int64)
            worst = 0.0
            for low in range(1 << (j - 1)):
                sel = hs[(hs & ((1 << (j - 1)) - 1)) == low]
                key = (q[sel] % M1) * m + (sel % m)
                cnt = np.bincount(key, minlength=M1 * m)
                worst = max(worst, abs(cnt * (M1 * m) / len(sel) - 1).max())
            row.append((n, worst))
        dec = (row[-1][1] / row[0][1]) ** (1 / (row[-1][0] - row[0][0]))
        check(row[-1][1] < row[0][1], "(q mod %d, h mod %d) on the cells of B_n with %d bottom digit(s) prescribed: relative "
              "deviation %s; decay %.3f per digit" % (M1, m, j - 1, ", ".join("%.4f (n = %d)" % (w, n) for n, w in row), dec))
    # the statistic of mixed_cells.out: cells with 3 top and 3 bottom digits prescribed, r = n - 6 free digits
    decs = []
    for (j, m) in ((2, 3), (3, 3), (2, 5), (2, 9), (2, 7), (2, 15)):
        M1 = (1 << (j - 1)) * m
        ncl = M1 * m
        rel = []
        for n in (14, 18, 20):
            i = 3
            hs = np.arange(1 << (n - 1), 1 << n, dtype=np.int64)
            cell = ((hs >> (n - i)) - (1 << (i - 1))) * (1 << i) + (hs & ((1 << i) - 1))
            ncell = (1 << (i - 1)) * (1 << i)
            cnt = np.bincount(cell * ncl + (q[hs] % M1) * m + (hs % m), minlength=ncell * ncl).reshape(ncell, ncl)
            rel.append(abs(cnt - (1 << (n - 2 * i)) / ncl).max() * ncl / (1 << (n - 2 * i)))
        decs.append(((rel[2] / rel[0]) ** (1 / 6), (1 << j) * m))
    check(all(0.70 < d[0] < 0.87 for d in decs), "cells with 3 top and 3 bottom digits prescribed, r = 8 to 14 free digits: decay "
          "of the relative deviation per free digit " + ", ".join("%.3f (modulus %d)" % d for d in decs)
          + " (mixed_cells.out, r = 8 to 16: 0.74 to 0.84)")


# ---------------------------------------------------------------------------------------------------------------
def v4(psi, q):
    print("V4 the algebra of (C)")
    rng = np.random.default_rng(96650)
    good = True
    goodc = True
    for _ in range(4000):
        b = int(rng.choice([2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 20]))
        h1, h2 = int(rng.integers(0, 4000)), int(rng.integers(0, 4000))
        s = int(rng.integers(0, 6))
        p1, p2, y1, y2 = int(psi[h1]), int(psi[h2]), 2 * h1 + 1, 2 * h2 + 1
        c3 = ((1 << s) * y1 * y2 - (1 if s == 0 else 0))
        good &= c3 % 2 == 0
        c3 //= 2
        Nb = (1 + b * p1) * (1 + b * p2) + b * c3
        a1, a2 = p1 + p2 + c3, p1 * p2
        good &= Nb == 1 + a1 * b + a2 * b * b and Nb % b == 1 and ((Nb - 1) // b) % b == a1 % b
        good &= p2 == h2 - int(q[h2])
        if b % 2 == 1:
            inv2 = pow(2, -1, b)
            x2 = 2 * p2 + 1
            if s >= 1:
                form = inv2 * x2 + ((1 << (s - 1)) * y1) * y2 + (p1 - inv2)
            else:
                form = inv2 * x2 + (inv2 * y1) * y2 + (p1 - 1)
            goodc &= (form - a1) % b == 0
    check(good, "N_b = x_b(h1) x_b(h2) + b (2^s y1 y2 - [s = 0])/2 = 1 + alpha_1 b + alpha_2 b^2, second digit alpha_1 mod b, "
          "psi2 = h2 - q2 (4000 random triples, eleven bases)")
    check(goodc, "odd b: alpha_1 = 2^(-1) x2 + (2^(s-1) y1) y2 + psi1 - 2^(-1) (s >= 1), 2^(-1) x2 + (2^(-1) y1) y2 + psi1 - 1 "
          "(s = 0) modulo b")
    good = True
    for _ in range(3000):
        a = int(rng.integers(1, 5))
        h1 = int(rng.integers(0, 2000))
        s = int(rng.integers(0, 6))
        low = int(rng.integers(0, 1 << a))                   # the a bottom digits of h2
        const = None
        for t in rng.integers(1, 1 << 12, size=8):
            h2 = (int(t) << a) | low
            y1, y2 = 2 * h1 + 1, 2 * h2 + 1
            a1 = int(psi[h1]) + int(psi[h2]) + ((1 << s) * y1 * y2 - (1 if s == 0 else 0)) // 2
            c = (a1 + int(q[h2])) % (1 << a)
            const = c if const is None else const
            good &= c == const
    check(good, "b = 2^a: on a cell with a bottom digits prescribed alpha_1 = const - q2 (mod 2^a) (3000 random cells, "
          "a = 1..4, s = 0..5)")
    good = True
    for _ in range(3000):
        a = int(rng.integers(1, 4))
        m = int(rng.choice([3, 5, 7, 9]))
        b = (1 << a) * m
        h1 = int(rng.integers(0, 2000))
        s = int(rng.integers(0, 6))
        h2 = int(rng.integers(1, 1 << 14))
        y1, y2 = 2 * h1 + 1, 2 * h2 + 1
        a1 = int(psi[h1]) + int(psi[h2]) + ((1 << s) * y1 * y2 - (1 if s == 0 else 0)) // 2
        # alpha_1 is a function of (q2 mod b, h2 mod b): affine, with coefficient -1 at q2
        coef = (1 + (1 << s) * y1) if s >= 1 else (1 + y1)
        const = int(psi[h1]) + ((1 << (s - 1)) * y1 if s >= 1 else h1)
        good &= (a1 - (-int(q[h2]) + coef * h2 + const)) % b == 0
    check(good, "every b: alpha_1 = -q2 + (1 + 2^s y1) h2 + const (mod b) (s = 0: coefficient 1 + y1), so for b = 2^a m it "
          "is uniform as soon as (q2 mod 2^a m, h2 mod m) is equidistributed on cells with a bottom digits prescribed")


# ---------------------------------------------------------------------------------------------------------------
def v5(psi, q):
    print("V5 (D): states of type p, fans, constants")
    n = len(q)
    h = np.arange(n, dtype=np.int64)
    x, y = 2 * psi + 1, 2 * h + 1
    g = np.gcd(x, y)
    check(bool(np.all(g == np.gcd(y, q))), "g(h) = gcd(x, y) = gcd(2h + 1, q(h)) for h < 2^%d" % (n.bit_length() - 1))
    good = True
    for p in (3, 5, 7, 11, 13):
        # automaton on (q mod p, h mod 2p), digits read from the top
        st_q = np.zeros(1 << 16, dtype=np.int64)
        st_h = np.zeros(1 << 16, dtype=np.int64)
        for hh in range(1, 1 << 16):
            par, bit = hh >> 1, hh & 1
            if bit == 0:
                st_q[hh] = st_q[par]
                st_h[hh] = (2 * st_h[par]) % (2 * p)
            else:
                st_q[hh] = (st_q[par] + (st_h[par] + 1) // 2) % p
                st_h[hh] = (2 * st_h[par] + 1) % (2 * p)
        acc = (st_q == 0) & ((2 * st_h + 1) % p == 0)
        good &= bool(np.all(acc == (g[: 1 << 16] % p == 0)))
        good &= len(set(zip(st_q.tolist(), st_h.tolist()))) <= 2 * p * p
    check(good, "p | g(h) iff 2h + 1 = 0 and q(h) = 0 (mod p); the automaton on (q mod p, h mod 2p) (at most 2 p^2 "
          "states, digits from the top) accepts exactly these h, h < 2^16, p = 3, 5, 7, 11, 13")
    bal = h[x == y]
    check(bal.tolist() == [0] + [1 << i for i in range(n.bit_length() - 1)], "x = y exactly at h = 0 and at the powers of 2 "
          "below 2^%d" % (n.bit_length() - 1))
    primes = [p for p in range(3, 400, 2) if all(p % r for r in range(3, int(p ** 0.5) + 1, 2))]
    first, fan_first, not_first, odd_order = {}, [], [], []
    for p in primes:
        idx = np.flatnonzero(g % p == 0)
        first[p] = int(idx[0]) if len(idx) else None
        o, v = 1, 2 % p
        while v != 1:
            v = (2 * v) % p
            o += 1
        if o % 2:
            odd_order.append(p)
            continue
        t = o // 2
        assert pow(2, t, p) == p - 1 and all(pow(2, u, p) != p - 1 for u in range(1, t))
        if first[p] is not None and first[p] == (1 << (t - 1)):
            fan_first.append(p)
        else:
            not_first.append(p)
    check(all(first[p] is not None for p in primes), "every prime below 400 has a state of type p below 2^%d (largest first "
          "state %d, at p = %d)" % (n.bit_length() - 1, max(first.values()), max(first, key=first.get)))
    check([p for p in fan_first if p < 200] == [3, 5, 11, 13, 17, 43, 113], "primes below 200 whose fan is the first state "
          "of type p: %s" % [p for p in fan_first if p < 200])
    check([p for p in fan_first if p > 200] == [241, 257, 331], "then, below 400: %s" % [p for p in fan_first if p > 200])
    check([first[p] for p in (19, 29, 37, 41)] == [218, 1435, 1350, 266] and all(p in not_first for p in (19, 29, 37, 41)),
          "not the first for 19, 29, 37, 41: first states 218, 1435, 1350, 266")
    check([first[p] for p in (7, 23, 31)] == [38, 103, 8664] and all(p in odd_order for p in (7, 23, 31)),
          "odd order (no balanced state of type p): 7, 23, 31 with first states 38, 103, 8664")
    f65537 = int(np.flatnonzero(g % 65537 == 0)[0])
    check(f65537 == 1 << 15 and first[257] == 128, "the fan is the first state of type p also for 257 (h = 128) and 65537 "
          "(h = 2^15), the two cases beyond 200 recorded by the S91 audit")
    big = [(p, first[p]) for p in primes if p < 200 and first[p] >= (1 << 16)]
    print("    primes below 200 whose first state of type p is beyond 2^16: %s" % big)
    check(big == [(127, 109283), (163, 68867), (193, 75173)], "these are 127, 163, 193 (so a search below 2^16 does not see "
          "them; none of them is a fan)")
    # divergence of the uncentred sums and partial sums of the centred ones
    xf, yf = x.astype(float), y.astype(float)
    w0 = np.log(2 * (xf + yf) / (xf + 2 * yf)) / xf
    w1 = np.log(1 + xf / (2 * yf)) / xf
    li2 = lambda z: sum(z ** k / k ** 2 for k in range(1, 400))
    l = [math.pi ** 2 / 12] + [-li2(-2.0 ** -jj) for jj in range(1, 4)]          # l_0 = -Li_2(-1) = pi^2/12
    kap = (l[0] - 2 * l[1] + l[2], l[1] - l[2])
    check(abs(kap[0] - 0.1615389173) < 1e-9 and abs(kap[1] - 0.2125139092) < 1e-9, "kappa_0, kappa_1 = %.10f, %.10f" % kap)
    for p, ref in ((3, (0.0556, 0.0789)), (5, (0.02943, 0.04250))):
        typ = (g % p == 0)
        cd = 1 / p ** 2
        res = []
        for sgm, w in ((0, w0), (1, w1)):
            cs = np.cumsum(np.where(typ, w, 0.0))
            H1, H2 = 1 << 12, n
            slope = (cs[H2 - 1] - cs[H1 - 1]) / (math.log(H2) - math.log(H1))
            cen = float(np.sum((typ - cd) * w))
            res.append((slope, cd * kap[sgm], cen))
        check(all(abs(r[0] / r[1] - 1) < 0.12 for r in res), "p = %d: sum over the states of type p below H of w^sigma grows "
              "like c_p kappa_sigma ln H: slope between 2^12 and 2^%d %.5f / %.5f against %.5f / %.5f"
              % (p, n.bit_length() - 1, res[0][0], res[1][0], res[0][1], res[1][1]))
        check(all(abs(r[2] - a) < 0.004 for r, a in zip(res, ref)), "p = %d: centred partial sums over h < 2^%d: %.4f / %.4f "
              "(THM-5128 (D4): %s / %s at 2^26): A_(p,p) is A_p" % (p, n.bit_length() - 1, res[0][2], res[1][2], ref[0], ref[1]))


# ---------------------------------------------------------------------------------------------------------------
def v6():
    print("V6 (B): class means of the series, and the quoted arithmetic")
    NMAX = 1 << 22
    sig = np.zeros(NMAX + 2, dtype=np.float64)                 # sigma_1(n)/n = sum_{d | n} 1/d, all n
    for d in range(1, NMAX + 2):
        sig[d::d] += 1.0 / d
    N = np.arange(1, NMAX, 2, dtype=np.int64)
    S = {}
    S[1] = (8 / math.pi ** 2) * sig[N]
    n0 = N + 1
    n0 = n0 // (n0 & -n0)
    S[0] = (8 / math.pi ** 2) * sig[n0]
    check(abs(S[0].mean() - 1) < 2e-3 and abs(S[1].mean() - 1) < 2e-3, "mean of S_sigma over the odd N < 2^22: %.5f / %.5f"
          % (S[0].mean(), S[1].mean()))
    worst = 0.0
    for m in (3, 9, 15, 45, 25):
        for sgm in (0, 1):
            csg = -1 if sgm == 0 else 0
            for c in range(m):
                sel = (N % m) == c
                worst = max(worst, abs(S[sgm][sel].mean() - float(L_m(m, c - csg))))
    check(worst < 4e-3, "mean of S_sigma(N) over the odd N < 2^22 of a class c modulo m equals L_m(c - c_sigma), "
          "m = 3, 9, 15, 45, 25, both sheets: largest difference %.5f" % worst)
    # arithmetic of the quoted figures (inputs: the outputs named in the file)
    check("%.2f" % (8.550 / 4.301) == "1.99" and "%.2f" % (10.766 / 5.438) == "1.98", "quintile ratios 8.550/4.301 = %.3f, "
          "10.766/5.438 = %.3f" % (8.550 / 4.301, 10.766 / 5.438))
    check("%.2f" % (1.609 / 1.001) == "1.61" and "%.2f" % (1.6094 / 1.0006) == "1.61", "abundancy ratio 1.609/1.001 = %.3f "
          "(1.6094/1.0006 = %.3f)" % (1.609 / 1.001, 1.6094 / 1.0006))
    check(round(100 * 0.4016 ** 2) == 16 and round(100 * 0.4275 ** 2) == 18, "squared correlations %.4f, %.4f"
          % (0.4016 ** 2, 0.4275 ** 2))
    check(round(100 * (0.4967 / 0.4542 - 1)) == 9 and round(100 * (0.4556 / 0.4142 - 1)) == 10,
          "L1 remainder over the Poisson-pair floor on [2^24, 2^25), from poisson_pair_floor.out: %.4f / %.4f, that is 9 %% / "
          "10 %% above the floor (verify_census.out, full precision: 1.0936 / 1.1001, so not within 10 %% on sheet 1)"
          % (0.4967 / 0.4542, 0.4556 / 0.4142))
    print("    the same ratio on [2^16, 2^17): %.4f / %.4f; on [2^20, 2^21): %.4f / %.4f"
          % (0.6051 / 0.5408, 0.5700 / 0.5063, 0.5396 / 0.5119, 0.5079 / 0.4464))
    d6 = ((0.02815 / 0.18410) ** (1 / 12), (0.02852 / 0.16307) ** (1 / 12))
    d3 = ((0.02374 / 0.13093) ** (1 / 12), (0.02279 / 0.12685) ** (1 / 12))
    check("%.3f / %.3f" % d6 == "0.855 / 0.865", "second digit, b = 6, decay per binary digit 2^16 -> 2^28: %.4f / %.4f" % d6)
    check("%.3f / %.3f" % d3 == "0.867 / 0.867", "b = 3: %.4f / %.4f; sqrt(3)/2 = %.4f" % (d3 + (math.sqrt(3) / 2,)))
    miss = (1 - 55892593 / 57818155, 1 - 35516685 / 35564915, 1 - 24945520 / 24947862, 1 - 18836356 / 18838778,
            1 - 14916908 / 14919260, 1 - 7632149 / 7634213)
    check("%.1f" % (100 * miss[0]) == "3.3" and max(miss[1:]) < 0.0014 and "%.2f" % (100 * miss[1]) == "0.14",
          "triples missed by the pre-audit census of every_base.py: " + ", ".join("%.3f %%" % (100 * v) for v in miss))


# ---------------------------------------------------------------------------------------------------------------
def v7():
    print("V7 the file")
    text = open(HYP, encoding="utf-8").read()
    check(all(ord(ch) < 128 for ch in text), "pure ASCII")
    check(text.count("`") % 2 == 0 and all(line.count("`") % 2 == 0 for line in text.split("\n")),
          "backticks balanced on every line")
    m = re.match(r"---\n(.*?)\n---\n", text, re.S)
    check(m is not None, "front matter delimited")
    try:
        import yaml
        fm = yaml.safe_load(m.group(1))
    except Exception as exc:                                  # noqa
        fm = None
        print("    YAML error: %s" % exc)
    check(isinstance(fm, dict) and set(fm) >= {"id", "title", "status", "session", "date", "source", "scripts", "related"},
          "front matter parses as YAML with id, title, status, session, date, source, scripts, related")
    paths = []
    if isinstance(fm, dict):
        check(fm["id"] == "HYP-9665" and str(fm["date"]) == "2026-11-11" and fm["session"] == "opus-2026-11-11-S92",
              "id, date, session")
        paths.append(fm["source"])
        for entry in list(fm["scripts"]) + list(fm["related"]):
            head = entry.split(" (")[0]
            paths += [p.strip() for p in head.split(", ")]
        missing = [p for p in paths if not os.path.exists(os.path.join(ROOT, p))]
        check(not missing, "%d paths of source / scripts / related exist%s" % (len(paths), "" if not missing else "; MISSING: %s" % missing))
        print("    scripts entries: %d, related entries: %d" % (len(fm["scripts"]), len(fm["related"])))
    body = text[m.end():] if m else text
    # every output named in the body is listed in scripts / related
    names = set(re.findall(r"[A-Za-z0-9_\-]+\.(?:out|txt|py|md|json)", body))
    listed = " ".join(paths)
    unl = sorted(nm for nm in names if nm not in listed)
    check(not unl, "every file named in the body is listed in the front matter%s" % ("" if not unl else "; NOT LISTED: %s" % unl))
    check("(H)" not in text and "PRE-AUDIT" not in text, "no (H), no PRE-AUDIT")
    check(not re.search(r"sigma\(n\)", text) and "sigma_1(n)" in text, "the divisor sum is written sigma_1(n)")
    check(not re.search(r"\bX1\b|\bX2\b|X_b", text), "no X1, X2, X_b (X is the bound only)")
    lab34 = set(re.findall(r"THM-5134 \(([A-Z])\)", text)) | set(re.findall(r"THM-5134 \(A\)-\(([A-Z])\)", text))
    lab28 = set(re.findall(r"THM-5128 \(([A-Z]\d?)\)", text))
    check(lab34 <= set("ABCDEFG") and lab28 <= {"A", "B", "B1", "B2", "B3", "B4", "C", "D", "D1", "D2", "D3", "D4", "E", "F", "G"},
          "labels cited: THM-5134 %s (statements (A)-(G) exist), THM-5128 %s" % (sorted(lab34), sorted(lab28)))
    check(not re.search(r"THM-5135 \((\d|[A-Za-z]\))", text), "THM-5135 is not cited by part in this file")
    check("9 % / 10 % above the counting-noise floor" in text and "within 10 %" not in text,
          "the L1 remainder is quoted as 9 % / 10 % above the floor (not as within 10 %)")
    check("admissible residues" in text and "65537" in text, "the repaired sentences are in place")
    check(set(re.findall(r"THM-5128 Remark (\d)", text)) <= {"1", "2", "3"}, "THM-5128 Remark 3 exists")
    lines = text.split("\n")
    prose = [ln for ln in lines[45:] if not ln.startswith("    ")]
    check(max(len(ln) for ln in prose) <= 120, "prose lines of the body at most 120 characters (longest %d; displayed "
          "formulas are longer, as in THM-5134)" % max(len(ln) for ln in prose))
    dup = [ln for ln in set(lines) if len(ln) > 40 and lines.count(ln) > 1]
    check(not dup, "no duplicated line")
    low = text.lower()
    check("chan-pak" not in low and "collatz" not in low, "no sentence about Chan-Pak or Collatz in the file")


def main():
    psi, q = build(1 << 20)
    v1()
    v2(q)
    v3(psi, q)
    v4(psi, q)
    v5(psi, q)
    v6()
    v7()
    print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
    sys.exit(0 if OK else 1)


if __name__ == "__main__":
    main()
