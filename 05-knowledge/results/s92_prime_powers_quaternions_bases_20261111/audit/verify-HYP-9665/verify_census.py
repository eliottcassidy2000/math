#!/usr/bin/env python3
"""Verifier of HYP-9665 (S92, 2026-11-11): independent censuses behind the numbers the file quotes from the audit.

Own enumeration (pairs h1 <= h2 with weight 2 off the diagonal, lower bound x(h) >= h + 1; no code of the session
or of the auditors is imported).
C1  F1 by value, joint histogram (sheet, type of the smaller state with respect to 315, N mod 2520) at X = 2^22 and
    2^24: compared entry by entry with the audit histograms hypotheses/even_moduli_hist_22.json, _24.json.
C2  the second base-b digit, both sheets, complete census: b = 6 at 2^16, 2^20, 2^24, 2^28; b = 4, 10 at 2^24 and
    2^28; b = 12 at 2^24; b = 2, 3, 5 at 2^24.
C3  single values on the odd N of [2^24, 2^25): the statistics quoted in (B).
Exit status 1 on failure.
"""
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.dirname(HERE)
OK = True


def check(cond, text):
    global OK
    OK &= bool(cond)
    print(("  ok   " if cond else "  FAIL ") + text, flush=True)


@njit
def build_psi(n):
    psi = np.zeros(n, dtype=np.int64)
    for h in range(1, n):
        k = h >> 1
        if h & 1:
            psi[h] = psi[k] + (k >> 1) + 1
        else:
            psi[h] = psi[k] + k
    return psi


@njit
def census_f1(X, psi, ty, hist):
    nH = psi.shape[0]
    for h1 in range(nH):
        x1 = 2 * psi[h1] + 1
        y1 = 2 * h1 + 1
        s = 0
        while True:
            p2 = np.int64(1) << s
            sub = 1 if s == 0 else 0
            sheet = 0 if s == 0 else 1
            if x1 * (h1 + 1) + p2 * y1 * y1 - sub > X:          # x(h2) >= h2 + 1 and h2 >= h1
                break
            for h2 in range(h1, nH):
                y2 = 2 * h2 + 1
                if x1 * (h2 + 1) + p2 * y1 * y2 - sub > X:
                    break
                N = x1 * (2 * psi[h2] + 1) + p2 * y1 * y2 - sub
                if N <= X:
                    hist[sheet, ty[h1], N % 2520] += 1 if h2 == h1 else 2
            s += 1


@njit
def census_digit(X, b, psi, dig):
    """dig[sheet, alpha_1 mod b] over the ordered triples with N_b <= X."""
    nH = psi.shape[0]
    short = False
    for h1 in range(nH):
        X1 = 1 + b * psi[h1]
        y1 = 2 * h1 + 1
        s = 0
        while True:
            sub = 1 if s == 0 else 0
            sheet = 0 if s == 0 else 1
            # c3 = (2^s y1 y2 - [s = 0]) / 2;  psi(h2) >= h2 // 2
            c3 = ((np.int64(1) << s) * y1 * y1 - sub) // 2
            if X1 * (1 + b * (h1 // 2)) + b * c3 > X:
                break
            closed = False
            for h2 in range(h1, nH):
                y2 = 2 * h2 + 1
                c3 = ((np.int64(1) << s) * y1 * y2 - sub) // 2
                if X1 * (1 + b * (h2 // 2)) + b * c3 > X:
                    closed = True
                    break
                Nb = X1 * (1 + b * psi[h2]) + b * c3
                if Nb <= X:
                    dig[sheet, ((Nb - 1) // b) % b] += 1 if h2 == h1 else 2
            if not closed:
                short = True
            s += 1
    return short


@njit
def census_r(lo, hi, psi, r):
    """r[sheet, (N - lo) // 2] for the odd N in [lo, hi)."""
    nH = psi.shape[0]
    X = hi - 1
    for h1 in range(nH):
        x1 = 2 * psi[h1] + 1
        y1 = 2 * h1 + 1
        s = 0
        while True:
            p2 = np.int64(1) << s
            sub = 1 if s == 0 else 0
            sheet = 0 if s == 0 else 1
            if x1 * (h1 + 1) + p2 * y1 * y1 - sub > X:
                break
            for h2 in range(h1, nH):
                y2 = 2 * h2 + 1
                if x1 * (h2 + 1) + p2 * y1 * y2 - sub > X:
                    break
                N = x1 * (2 * psi[h2] + 1) + p2 * y1 * y2 - sub
                if N <= X and N >= lo:
                    r[sheet, (N - lo) >> 1] += 1 if h2 == h1 else 2
            s += 1


@njit
def abundancy_odd(n):
    """a[i] = sigma_1(k)/k for the odd k = 2i + 1 < n."""
    a = np.zeros(n // 2, dtype=np.float64)
    for d in range(1, n, 2):
        inv = 1.0 / d
        k = d
        while k < n:
            a[k >> 1] += inv
            k += 2 * d
    return a


@njit
def types315(psi, lut):
    n = psi.shape[0]
    ty = np.zeros(n, dtype=np.int8)
    for h in range(n):
        a = 2 * psi[h] + 1
        b = 2 * h + 1
        while b:
            a, b = b, a % b
        c = 315
        while c:
            a, c = c, a % c
        ty[h] = lut[a]
    return ty


def c1():
    print("C1 F1 by value against the audit histograms")
    divs = [d for d in range(1, 316) if 315 % d == 0]
    for E in (20, 22, 24, 26, 28):
        t0 = time.time()
        X = 1 << E
        nH = (X - 1) // 3 + 2
        psi = build_psi(nH)
        lut = np.zeros(316, dtype=np.int64)
        for i, d in enumerate(divs):
            lut[d] = i
        hist = np.zeros((2, 12, 2520), dtype=np.int64)
        census_f1(X, psi, types315(psi, lut), hist)
        ref = np.array(json.load(open(os.path.join(AUDIT, "hypotheses", "even_moduli_hist_%d.json" % E)))["hist"],
                       dtype=np.int64)
        check(hist.shape == ref.shape and bool(np.all(hist == ref)),
              "X = 2^%d: %d / %d triples; the histogram (2 x 12 x 2520 entries) equals the audit histogram (%.1f s)"
              % (E, hist[0].sum(), hist[1].sum(), time.time() - t0))
        check(int(hist[:, :, 0::2].sum()) == 0, "X = 2^%d: every N is odd" % E)


def digits(b, E):
    X = 1 << E
    nH = (2 * X) // b + 8                                   # b psi(h) < N_b and psi(h) >= h/2
    psi = build_psi(nH)
    dig = np.zeros((2, b), dtype=np.int64)
    short = census_digit(X, b, psi, dig)
    assert not short
    return dig


def dev(dig, mod=None):
    b = dig.shape[1]
    out = []
    for sh in (0, 1):
        v = dig[sh].astype(float)
        if mod is not None:
            v = np.array([v[np.arange(b) % mod == c].sum() for c in range(mod)])
        out.append(float(np.abs(len(v) * v / v.sum() - 1).max()))
    return out


def c2():
    print("C2 the second base-b digit, complete census, max_c |b share(c) - 1|, sheet 0 / sheet 1")
    want6 = {16: ("0.184", "0.163"), 20: ("0.090", "0.098"), 24: ("0.048", "0.052"), 28: ("0.028", "0.029")}
    d6 = {}
    for E in (16, 20, 24, 28):
        t0 = time.time()
        dig = digits(6, E)
        d = dev(dig)
        d6[E] = (d, dev(dig, 2), dev(dig, 3), dig.sum(axis=1))
        check(("%.3f" % d[0], "%.3f" % d[1]) == want6[E], "b = 6, X = 2^%d: %d / %d triples, %.5f / %.5f; digit mod 2: "
              "%.5f / %.5f; digit mod 3: %.5f / %.5f (%.1f s)" % (E, d6[E][3][0], d6[E][3][1], d[0], d[1], d6[E][1][0],
                                                               d6[E][1][1], d6[E][2][0], d6[E][2][1], time.time() - t0))
    check(tuple(int(v) for v in d6[24][3]) == (9609816, 14919260) and tuple(int(v) for v in d6[28][3]) == (177517734, 280437250),
          "b = 6: triple counts at 2^24 and 2^28 as in second_digit_b6.out (sheet 1 at 2^24 as in every_base.out)")
    dec = tuple((d6[28][0][i] / d6[16][0][i]) ** (1 / 12) for i in (0, 1))
    check("%.3f / %.3f" % dec == "0.855 / 0.865", "b = 6: decay per binary digit from 2^16 to 2^28: %.4f / %.4f" % dec)
    check("%.5f / %.5f" % tuple(d6[28][1]) == "0.00007 / 0.00004", "b = 6 at 2^28: digit modulo 2 flat to %.5f / %.5f"
          % tuple(d6[28][1]))
    check(all(d6[28][2][i] > 0.95 * d6[28][0][i] for i in (0, 1)), "b = 6 at 2^28: the digit modulo 3 carries the deviation "
          "(%.5f / %.5f of %.5f / %.5f)" % (d6[28][2][0], d6[28][2][1], d6[28][0][0], d6[28][0][1]))
    want = {(4, 24): ("0.0007", "0.0009"), (4, 28): ("0.0001", "0.0003"), (10, 24): ("0.0051", "0.0034"),
            (10, 28): ("0.0015", "0.0007")}
    for (b, E) in ((4, 24), (4, 28), (10, 24), (10, 28)):
        t0 = time.time()
        dig = digits(b, E)
        d = dev(dig)
        check(("%.4f" % d[0], "%.4f" % d[1]) == want[(b, E)], "b = %d, X = 2^%d: %.5f / %.5f (%d / %d triples, %.1f s)"
              % (b, E, d[0], d[1], dig[0].sum(), dig[1].sum(), time.time() - t0))
    dig = digits(12, 24)
    d, d4 = dev(dig), dev(dig, 4)
    check(("%.3f" % d[0], "%.3f" % d[1]) == ("0.064", "0.062") and "%.4f" % max(d4) == "0.0015",
          "b = 12, X = 2^24: %.5f / %.5f; digit modulo 4 flat to %.5f / %.5f" % (d[0], d[1], d4[0], d4[1]))
    s1 = {}
    for b in (2, 3, 5):
        dig = digits(b, 24)
        s1[b] = (dev(dig), int(dig[1].sum()))
    check(["%.3f" % s1[b][0][1] for b in (2, 3, 5)] == ["0.000", "0.040", "0.002"]
          and [s1[b][1] for b in (2, 3, 5)] == [57818155, 35564915, 18838778],
          "sheet 1 at 2^24, b = 2, 3, 5: %s with %s triples (every_base.out)"
          % (", ".join("%.5f" % s1[b][0][1] for b in (2, 3, 5)), ", ".join(str(s1[b][1]) for b in (2, 3, 5))))
    worst = max(s1[2][0][1], s1[5][0][1], 0.00094, 0.00340)
    print("    sheet 1 at 2^24: largest deviation among b = 2, 4, 5, 10 is %.4f (b = 10)" % worst)


def c3():
    print("C3 single values on the odd N of [2^24, 2^25)")
    t0 = time.time()
    lo, hi = 1 << 24, 1 << 25
    psi = build_psi((hi - 1) // 3 + 2)
    r = np.zeros((2, (hi - lo) // 2), dtype=np.int64)
    census_r(lo, hi, psi, r)
    ab = abundancy_odd(hi + 2)
    N = np.arange(lo + 1, hi, 2, dtype=np.int64)
    kap = (0.1615389173, 0.2125139092)
    res = []
    for sh in (0, 1):
        if sh == 1:
            n = N
        else:
            n = N + 1
            n = n // (n & -n)
        a = ab[n >> 1]
        S = (8 / math.pi ** 2) * a
        rr = r[sh].astype(float)
        mean = rr.mean()
        model = 2 * kap[sh] * S * np.log(N)
        corr = float(np.corrcoef(rr, a)[0, 1])
        order = np.argsort(a, kind="stable")
        qs = np.array_split(order, 5)
        qr = [float(rr[q].mean()) for q in qs]
        qa = [float(a[q].mean()) for q in qs]
        res.append(dict(mean=mean, model=float(model.mean()), corr=corr, l1=float(np.abs(rr - model).mean() / mean),
                        l1c=float(np.abs(rr - mean).mean() / mean), miss=float((rr == 0).mean()),
                        disp=float(rr.var() / mean), qr=qr, qa=qa, trip=int(r[sh].sum()),
                        v3=[float(rr[(n % 3 ** (v + 1) != 0) & (n % 3 ** v == 0)].mean()) for v in range(3)],
                        v5=[float(rr[(n % 5 ** (v + 1) != 0) & (n % 5 ** v == 0)].mean()) for v in range(3)]))
    print("    %d / %d triples in the window (%.1f s)" % (res[0]["trip"], res[1]["trip"], time.time() - t0))
    check((res[0]["trip"], res[1]["trip"]) == (49644073, 62761167), "triples in the window as in the three audit outputs")
    check(["%.3f" % x["mean"] for x in res] == ["5.918", "7.482"], "mean multiplicities %.4f / %.4f" % (res[0]["mean"], res[1]["mean"]))
    check(["%.2f" % x["model"] for x in res] == ["5.50", "7.23"], "means of 2 kappa_sigma S_sigma(N) ln N: %.4f / %.4f"
          % (res[0]["model"], res[1]["model"]))
    check(["%.4f" % x["corr"] for x in res] == ["0.4016", "0.4275"], "correlation of r with sigma_1(n)/n: %.4f / %.4f"
          % (res[0]["corr"], res[1]["corr"]))
    check(["%.0f" % (100 * x["corr"] ** 2) for x in res] == ["16", "18"], "squared correlation: %.1f %% / %.1f %%"
          % (100 * res[0]["corr"] ** 2, 100 * res[1]["corr"] ** 2))
    check(["%.2f" % x["qr"][0] for x in res] == ["4.30", "5.44"] and ["%.2f" % x["qr"][4] for x in res] == ["8.55", "10.77"],
          "mean of r by quintile of sigma_1(n)/n: %s / %s" % (", ".join("%.3f" % v for v in res[0]["qr"]),
                                                              ", ".join("%.3f" % v for v in res[1]["qr"])))
    check(["%.2f" % (x["qr"][4] / x["qr"][0]) for x in res] == ["1.99", "1.98"], "ratios highest / lowest quintile: %.3f / %.3f"
          % (res[0]["qr"][4] / res[0]["qr"][0], res[1]["qr"][4] / res[1]["qr"][0]))
    check(all("%.3f" % x["qa"][0] == "1.001" and "%.3f" % x["qa"][4] == "1.609" and "%.2f" % (x["qa"][4] / x["qa"][0]) == "1.61"
              for x in res), "quintile means of sigma_1(n)/n: %.4f to %.4f (ratio %.3f) / %.4f to %.4f (ratio %.3f)"
          % (res[0]["qa"][0], res[0]["qa"][4], res[0]["qa"][4] / res[0]["qa"][0],
             res[1]["qa"][0], res[1]["qa"][4], res[1]["qa"][4] / res[1]["qa"][0]))
    v3 = [res[0]["v3"][i] / res[0]["v3"][0] for i in (1, 2)]
    v5 = [res[0]["v5"][i] / res[0]["v5"][0] for i in (1, 2)]
    check(["%.3f" % v for v in v3] == ["1.455", "1.722"] and ["%.3f" % v for v in v5] == ["1.318", "1.323"],
          "sheet 0: ratios by v_3(n) = 1, 2: %.3f, %.3f (limits 1.333, 1.444); by v_5(n): %.3f, %.3f (limits 1.200, 1.240)"
          % (v3[0], v3[1], v5[0], v5[1]))
    check(["%.3f" % x["l1"] for x in res] == ["0.497", "0.456"], "L1 remainder mean |r - 2 kappa S ln N| / mean r: %.4f / %.4f"
          % (res[0]["l1"], res[1]["l1"]))
    check(["%.3f" % x["l1c"] for x in res] == ["0.532", "0.505"], "constant predictor: %.4f / %.4f" % (res[0]["l1c"], res[1]["l1c"]))
    floor = []
    for x in res:
        lam = x["mean"] / 2
        k = int(math.floor(lam))
        floor.append(2 * math.exp(-lam) * lam ** (k + 1) / math.factorial(k) / lam)      # E|P - lam| / lam, P Poisson(lam)
    check(["%.3f" % v for v in floor] == ["0.454", "0.414"], "Poisson floor of r = 2 Poisson(mean/2): %.4f / %.4f" % tuple(floor))
    print("    L1 remainder over that floor: %.4f / %.4f" % (res[0]["l1"] / floor[0], res[1]["l1"] / floor[1]))
    check(["%.2f" % (100 * x["miss"]) for x in res] == ["8.07", "4.86"], "misses r = 0: %.3f %% / %.3f %%"
          % (100 * res[0]["miss"], 100 * res[1]["miss"]))
    check(["%.2f" % x["disp"] for x in res] == ["2.79", "3.06"], "Var(r)/mean(r): %.3f / %.3f" % (res[0]["disp"], res[1]["disp"]))


def main():
    c1()
    c2()
    c3()
    print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
    sys.exit(0 if OK else 1)


if __name__ == "__main__":
    main()
