#!/usr/bin/env python3
"""Verifier of THM-5134 (after the audit corrections): the analytic ingredients of the proof of (G) and the
constants quoted in the file. Own code (mpmath), written from the text of the theorem.

  A1  int_0^1 Z(v) dv = (X - k) w_s^b(h);  w_0^b = w^0_b;  sum_{s >= 1} w_s^b = w^1_b
  A2  the additive constant: |k| w_s^b <= ln 2, with the three sub-bounds of the proof
  A3  F(r) <= 2^(-s)/b, the range of r for h >= 1, the derivative of F
  A4  int_{1/2}^{1} F((b/2) rho) d rho = (4/b^2)(l_(s+1) - 2 l_(s+2) + l_(s+3)); kappa_0(b), kappa_1(b)
  A5  block sums of w_s^b over B_n against ln 2 * the integral
  A6  constants of (E), (F), (G): 8/pi^2, 9/pi^2, 25/(3 pi^2), e^gamma/2, 4 e^gamma/pi^2, S at 7 and 21
  A7  numbers derived from outputs: theta_m/2, noise, kappa slopes, binomial tail, rates per digit, r^2
Exit status 1 on failure.
"""
import math
import sys

import mpmath as mp

mp.mp.dps = 30
FAIL = []


def check(name, ok, extra=""):
    print(("  ok   " if ok else "  FAIL ") + name + ((" : " + extra) if extra else ""))
    if not ok:
        FAIL.append(name)


NB = 18
PSI = [0] * (1 << NB)
for h in range(1, 1 << NB):
    k = h >> 1
    PSI[h] = PSI[k] + (k if h % 2 == 0 else (k >> 1) + 1)


def xb(b, h):
    return 1 + b * PSI[h]


def y(h):
    return 2 * h + 1


def w_s(b, h, s):
    X_, Y_ = xb(b, h), y(h)
    return (2.0 / (b * X_)) * math.log((X_ + 2 ** s * Y_) / (X_ / 2.0 + 2 ** s * Y_))


def w0(b, h):
    X_, Y_ = xb(b, h), y(h)
    return (2.0 / (b * X_)) * math.log(2.0 * (X_ + Y_) / (X_ + 2.0 * Y_))


def w1(b, h):
    X_, Y_ = xb(b, h), y(h)
    return (2.0 / (b * X_)) * math.log(1.0 + X_ / (2.0 * Y_))


print("A1 the weights")
worst = 0.0
for b in (2, 3, 5, 6, 10):
    for h in (0, 1, 2, 3, 6, 45, 117, 1000, 54321):
        for s in (0, 1, 2, 5):
            a = b * xb(b, h) / 2.0
            bp = b * 2.0 ** (s - 1) * y(h)
            k = xb(b, h) * (1 - b / 2.0) - (b / 2.0 if s == 0 else 0.0)
            X = 1.0e7
            val = mp.quad(lambda v: (X - k) / (2 * (a * (1 - v / 2) + bp)), [0, 1])
            worst = max(worst, abs(float(val) / ((X - k) * w_s(b, h, s)) - 1))
check("int_0^1 Z(v) dv = (X - k) w_s^b(h), Z = (X - k)/(2 (a (1 - v/2) + b'))  (180 cases)", worst < 1e-12,
      "largest relative difference %.1e" % worst)
worst = 0.0
worst1 = 0.0
for b in (2, 3, 4, 5, 6, 10):
    for h in range(0, 3000, 7):
        worst = max(worst, abs(w_s(b, h, 0) - w0(b, h)) / w0(b, h))
        tot = sum(w_s(b, h, s) for s in range(1, 60))
        worst1 = max(worst1, abs(tot - w1(b, h)) / w1(b, h))
check("w_0^b = (2/(b x_b)) ln(2 (x_b + y)/(x_b + 2y))", worst < 1e-12, "%.1e" % worst)
check("sum_{s >= 1} w_s^b = (2/(b x_b)) ln(1 + x_b/(2y))", worst1 < 1e-12, "%.1e" % worst1)
ok = True
for b in (2, 3, 10):
    for h in (0, 5, 77):
        X_, Y_ = xb(b, h), y(h)
        for S in (1, 3, 10, 30):
            prod = 1.0
            for s in range(1, S + 1):
                prod *= 2.0 * (X_ + 2 ** s * Y_) / (X_ + 2 ** (s + 1) * Y_)
            closed = 2.0 ** S * (X_ + 2.0 * Y_) / (X_ + 2.0 ** (S + 1) * Y_)
            if abs(prod / closed - 1) > 1e-12:
                ok = False
check("product of the ratios over 1 <= s <= S is 2^S (x_b + 2y)/(x_b + 2^(S+1) y) -> 1 + x_b/(2y)", ok)
check("b = 2: these are the weights of THM-5128 (D2): w^0 = (1/x) ln(2 (x + y)/(x + 2y)), w^1 = (1/x) ln(1 + x/(2y))",
      all(abs(w0(2, h) - (1.0 / (2 * PSI[h] + 1)) * math.log(2.0 * (2 * PSI[h] + 1 + y(h)) / (2 * PSI[h] + 1 + 2 * y(h))))
          < 1e-15 for h in range(200)))

print("A2 the additive constant")
m1 = 0.0
m0 = 0.0
m00 = 0.0
okf = True
for b in range(2, 13):
    for h in range(0, 1 << 13):
        X_, Y_ = xb(b, h), y(h)
        for s in range(0, 7):
            k = X_ * (1 - b / 2.0) - (b / 2.0 if s == 0 else 0.0)
            val = abs(k) * w_s(b, h, s)
            if s >= 1:
                m1 = max(m1, val / ((1 - 2.0 / b) * math.log(2.0)) if b > 2 else 0.0)
                if b == 2 and val != 0.0:
                    okf = False
            else:
                form = (1 - 2.0 / b + 1.0 / X_) * math.log(2.0 * (X_ + Y_) / (X_ + 2.0 * Y_))
                if abs(val - form) > 1e-12:
                    okf = False
                if h == 0:
                    m00 = max(m00, val)
                else:
                    m0 = max(m0, val)
check("s >= 1: |k| w_s^b < (1 - 2/b) ln 2  (b <= 12, h < 2^13, s <= 6)", m1 < 1, "largest ratio %.6f" % m1)
check("s = 0: |k| w_0^b = (1 - 2/b + 1/x_b) ln(2 (x_b + y)/(x_b + 2y))", okf)
check("s = 0, h = 0: below 2 ln(4/3) = %.4f" % (2 * math.log(4 / 3.0)), m00 < 2 * math.log(4 / 3.0), "max %.4f" % m00)
check("s = 0, h >= 1: below ln 2 = %.4f" % math.log(2.0), m0 < math.log(2.0), "max %.4f" % m0)
check("2 ln(4/3) < ln 2, so |k| w_s^b <= ln 2 in every case", 2 * math.log(4 / 3.0) < math.log(2.0))

print("A3 regularity")


def F(b, r, s):
    return (2.0 / (b * r)) * math.log((r + 2.0 ** s) / (r / 2.0 + 2.0 ** s))


okF = True
okr = True
okw = True
rmin = {}
for b in range(2, 13):
    lo, hi = 1e9, 0.0
    for h in range(1, 1 << 14):
        r = xb(b, h) / float(y(h))
        rho = (2 * PSI[h] + 1) / float(y(h))
        if abs(r - (b / 2.0 * rho + (1 - b / 2.0) / y(h))) > 1e-12:
            okr = False
        lo, hi = min(lo, r), max(hi, r)
        for s in (0, 1, 4):
            if F(b, r, s) > 2.0 ** (-s) / b + 1e-15:
                okF = False
            if abs(F(b, r, s) / y(h) - w_s(b, h, s)) > 1e-15:
                okw = False
    rmin[b] = (lo, hi)
    if not (lo > (b + 4) / 12.0 and hi <= b / 2.0 + 1e-12):
        okr = False
check("w_s^b = F(r)/y, r = x_b/y = (b/2) rho + (1 - b/2)/y", okw and okr)
check("F(r) <= 2^(-s)/b", okF)
check("for h >= 1: (b + 4)/12 < r <= b/2  (b <= 12, h < 2^14)", okr,
      "; ".join("b=%d: [%.4f, %.4f] in (%.4f, %.1f]" % (b, rmin[b][0], rmin[b][1], (b + 4) / 12.0, b / 2.0)
                for b in (2, 3, 10)))
okd = True
for b in (2, 3, 5, 10):
    for s in (0, 1, 3, 6):
        for i in range(200):
            r = (b + 4) / 12.0 + (b / 2.0 - (b + 4) / 12.0) * i / 199.0
            d = (F(b, r + 1e-6, s) - F(b, r - 1e-6, s)) / 2e-6
            if abs(d) > (2.0 / b) * 2.0 ** (-s) / ((b + 4) / 12.0):
                okd = False
check("|F'(r)| <= (2/b) 2^(-s)/r_min on the range (derivative O_b(2^-s))", okd)

print("A4 the integral and the mass constants")


def l(b, j):
    return -mp.polylog(2, -mp.mpf(b) / 2 ** j)


worst = mp.mpf(0)
for b in (2, 3, 4, 5, 6, 10):
    for s in range(0, 6):
        num = mp.quad(lambda rho: (2 / (b * (mp.mpf(b) / 2) * rho))
                      * mp.log(((mp.mpf(b) / 2) * rho + 2 ** s) / ((mp.mpf(b) / 4) * rho + 2 ** s)), [mp.mpf(1) / 2, 1])
        closed = (mp.mpf(4) / b ** 2) * (l(b, s + 1) - 2 * l(b, s + 2) + l(b, s + 3))
        prim = (mp.mpf(4) / b ** 2) * ((mp.polylog(2, -mp.mpf(b) / 2 ** (s + 2)) - mp.polylog(2, -mp.mpf(b) / 2 ** (s + 1)))
                                      - (mp.polylog(2, -mp.mpf(b) / 2 ** (s + 3)) - mp.polylog(2, -mp.mpf(b) / 2 ** (s + 2))))
        worst = max(worst, abs(num - closed), abs(prim - closed))
check("int_{1/2}^{1} F((b/2) rho) d rho = (4/b^2)[Li_2(-b rho/2^(s+2)) - Li_2(-b rho/2^(s+1))] = "
      "(4/b^2)(l_(s+1) - 2 l_(s+2) + l_(s+3))", worst < mp.mpf(10) ** (-20), "largest difference %s" % mp.nstr(worst, 3))


def kappa0(b):
    return (mp.mpf(4) / b ** 2) * (l(b, 1) - 2 * l(b, 2) + l(b, 3))


def kappa1(b):
    return (mp.mpf(4) / b ** 2) * (l(b, 2) - l(b, 3))


tele = max(abs(sum((mp.mpf(4) / b ** 2) * (l(b, s + 1) - 2 * l(b, s + 2) + l(b, s + 3)) for s in range(1, 80))
               - kappa1(b)) for b in (2, 3, 10))
check("the sum over s >= 1 telescopes to (4/b^2)(l_2 - l_3)", tele < mp.mpf(10) ** (-20))
vals = {b: (float(kappa0(b)), float(kappa1(b))) for b in (2, 3, 4, 5, 10)}
print("       kappa_0(b), kappa_1(b): " + "; ".join("b = %d: %.10f, %.10f" % (b, vals[b][0], vals[b][1]) for b in vals))
check("b = 2: 0.1615389173, 0.2125139092 (THM-5128 (D1))",
      abs(vals[2][0] - 0.1615389173) < 5e-11 and abs(vals[2][1] - 0.2125139092) < 5e-11)
check("b = 3: 0.091805, 0.132470;  b = 10: 0.013661, 0.028354 (the values in (G))",
      "%.6f %.6f" % vals[3] == "0.091805 0.132470" and "%.6f %.6f" % vals[10] == "0.013661 0.028354",
      "%.6f %.6f | %.6f %.6f" % (vals[3] + vals[10]))
check("b = 2: kappa_0 = pi^2/12 + 2 Li_2(-1/2) - Li_2(-1/4), kappa_1 = Li_2(-1/4) - Li_2(-1/2)",
      abs(kappa0(2) - (mp.pi ** 2 / 12 + 2 * mp.polylog(2, -0.5) - mp.polylog(2, -0.25))) < 1e-25
      and abs(kappa1(2) - (mp.polylog(2, -0.25) - mp.polylog(2, -0.5))) < 1e-25)

print("A5 block sums")
for b in (3, 10):
    for s in (0, 1):
        integ = float((mp.mpf(4) / b ** 2) * (l(b, s + 1) - 2 * l(b, s + 2) + l(b, s + 3)))
        row = []
        okb = True
        for n in (10, 12, 14, 16, 18):
            tot = sum(w_s(b, h, s) for h in range(1 << (n - 1), 1 << n))
            row.append("n=%d: %.6f" % (n, tot))
            if abs(tot - math.log(2.0) * integ) > 3.0 * 2.0 ** (-s) * 2.0 ** (-n / 2.0):
                okb = False
        check("b = %d, s = %d: sum over B_n of w_s^b -> ln 2 * integral = %.6f" % (b, s, math.log(2.0) * integ), okb,
              "; ".join(row))

print("A6 constants")
pi2 = math.pi ** 2
check("8/pi^2 = 0.8106", "%.4f" % (8 / pi2) == "0.8106")
prod = 1.0
for p in range(3, 200000, 2):
    if all(p % d for d in range(3, int(p ** 0.5) + 1, 2)):
        prod *= 1 - 1.0 / (p * p)
check("prod_{p odd} (1 - p^-2) = 8/pi^2 (partial product over p < 200000: %.7f)" % prod, abs(prod - 8 / pi2) < 1e-6)
check("(8/pi^2)/(1 - 1/9) = 9/pi^2 and (8/pi^2)/(1 - 1/25) = 25/(3 pi^2)",
      abs((8 / pi2) / (1 - 1 / 9.0) - 9 / pi2) < 1e-15 and abs((8 / pi2) / (1 - 1 / 25.0) - 25 / (3 * pi2)) < 1e-15)
check("(6/pi^2)/(1 - 1/4) = 8/pi^2", abs((6 / pi2) / 0.75 - 8 / pi2) < 1e-15)
eg = math.exp(0.5772156649015329)
check("e^gamma/2 = 0.8905 and 4 e^gamma/pi^2 = 0.7218", "%.4f" % (eg / 2) == "0.8905" and "%.4f" % (4 * eg / pi2) == "0.7218")
k0, k1 = vals[2]
w_0, w_1 = k0 / (k0 + k1), k1 / (k0 + k1)


def sig1(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def S(N, sheet):
    n = odd_part(N + (1 if sheet == 0 else 0))
    return (8 / pi2) * sig1(n) / n


d7 = w_0 * S(7, 0) + w_1 * S(7, 1)
d21 = w_0 * S(21, 0) + w_1 * S(21, 1)
check("w_0 S_0 + w_1 S_1 = 0.876 at N = 7 and 1.084 at N = 21 (w_0 = %.5f, w_1 = %.5f)" % (w_0, w_1),
      "%.3f %.3f" % (d7, d21) == "0.876 1.084", "%.5f %.5f" % (d7, d21))
tot0 = sum(S(N, 0) for N in range(1, 200001, 2)) / 100000.0
tot1 = sum(S(N, 1) for N in range(1, 200001, 2)) / 100000.0
check("mean of S_0 and of S_1 over the odd N < 200000 is 1", abs(tot0 - 1) < 1e-3 and abs(tot1 - 1) < 1e-4,
      "%.5f %.5f" % (tot0, tot1))
check("S_1(N) > 8/pi^2 for odd 3 <= N < 20000; S_0(N) = 8/pi^2 exactly when N + 1 is a power of 2",
      all(S(N, 1) > 8 / pi2 + 1e-12 for N in range(3, 20000, 2))
      and all((abs(S(N, 0) - 8 / pi2) < 1e-12) == ((N + 1) & N == 0) for N in range(1, 20000, 2)))
check("limits at 9: 8/9, 32/27, 35/27 = 0.889 / 1.185 / 1.296",
      "%.3f %.3f %.3f" % (8 / 9.0, 32 / 27.0, 35 / 27.0) == "0.889 1.185 1.296")

print("A7 numbers derived from outputs")
rates = {"9, 27, 81": 0.866025, "25, 125": 0.815263, "multiples of 7": 0.932039, "45, 75": 0.892191}
th = {k: math.log(1 / r) / math.log(2 / r ** 2) / 2 for k, r in rates.items()}
print("       sup of theta_m/2: " + "; ".join("%s: %.4f" % (k, v) for k, v in th.items()))
check("suprema 0.0733, 0.0927, 0.0422, 0.0619; hence below 0.074, 0.093, 0.043, 0.062 (and NOT below 0.073, 0.042)",
      ["%.4f" % v for v in th.values()] == ["0.0733", "0.0927", "0.0422", "0.0619"]
      and th["9, 27, 81"] > 0.073 and th["multiples of 7"] > 0.042)
T0, T1 = 854415478, 1083261684
ratio = math.sqrt(T1 / float(T0))
print("       noise sqrt(m/#T) at 2^28, sheet 0 / sheet 1: " + "; ".join(
    "%d: %.5f / %.5f" % (m, math.sqrt(m / float(T0)), math.sqrt(m / float(T1))) for m in (9, 27, 81, 25, 125, 49, 45, 75)))
print("       ratio of the sheet-0 noise to the sheet-1 noise: %.4f (%.1f %% larger)" % (ratio, 100 * (ratio - 1)))
same = [m for m in (9, 27, 81, 25, 125, 49, 45, 75)
        if "%.4f" % math.sqrt(m / float(T0)) == "%.4f" % math.sqrt(m / float(T1))]
check("the two noise values round to the same 4-decimal entry except at 125 (0.0004 against 0.0003)",
      same == [9, 27, 81, 25, 49, 45, 75] and "%.4f" % math.sqrt(125 / float(T0)) == "0.0004"
      and "%.4f" % math.sqrt(125 / float(T1)) == "0.0003")
check("the sheet-0 noise is 12.6 %% larger (the file must not say 12 %%): %.1f" % (100 * (ratio - 1)),
      abs(100 * (ratio - 1) - 12.6) < 0.05)
# kappa slopes from the triple counts printed in proof-5134-FH/h5_kappa_base_b.out
counts = {
    (2, 0): (45886662, 198575546), (2, 1): (57818155, 251044267),
    (3, 0): (26257636, 113570879), (3, 1): (35564915, 154584109),
    (4, 0): (17415766, 75250116), (4, 1): (24947862, 108491355),
    (5, 0): (12578286, 54279580), (5, 1): (18838778, 81940677),
    (10, 0): (4489156, 19227277), (10, 1): (7634213, 33174398),
}
worst = 0.0
for (b, sh), (c24, c26) in counts.items():
    slope = (c26 / 2.0 ** 26 - c24 / 2.0 ** 24) / math.log(4.0)
    kap = float(kappa0(b) if sh == 0 else kappa1(b))
    worst = max(worst, abs(slope - kap))
check("slopes of #T^b/X against ln X between 2^24 and 2^26 agree with kappa_sigma(b) to within 0.00001 "
      "(b = 2, 3, 4, 5, 10, both sheets)", worst < 1e-5, "largest difference %.7f" % worst)
tail = sum(math.comb(12, j) for j in range(8, 13)) / 2.0 ** 12
check("P(at least 8 of 12) under fair coin flips = 0.19", "%.2f" % tail == "0.19", "%.4f" % tail)
r49 = (0.1689 / 0.2980) ** (1 / 8.0)
r3 = ((0.02374 / 0.13093) ** (1 / 12.0), (0.02279 / 0.12685) ** (1 / 12.0))
check("fixed state h1 = 0, s = 0 modulo 49: factor per binary digit 0.93 (0.2980 -> 0.1689 over 8 digits)",
      "%.2f" % r49 == "0.93", "%.4f" % r49)
check("second digit, b = 3: factor 0.867 per binary digit on both sheets (2^16 -> 2^28); sqrt(3)/2 = 0.8660",
      "%.3f %.3f" % r3 == "0.867 0.867" and "%.4f" % (math.sqrt(3) / 2) == "0.8660", "%.4f %.4f" % r3)
check("correlations 0.4016, 0.4275: r^2 = 16 %% and 18 %%", "%.0f %.0f" % (100 * 0.4016 ** 2, 100 * 0.4275 ** 2) == "16 18")
check("quintile ratios: 8.550/4.301 = 1.99, 10.766/5.438 = 1.98, 1.6094/1.0006 = 1.61",
      "%.2f %.2f %.2f" % (8.550 / 4.301, 10.766 / 5.438, 1.6094 / 1.0006) == "1.99 1.98 1.61")
ratios = [0.1520 / 0.040, 0.1512 / 0.040, 0.4483 / 0.040, 0.1682 / 0.040, 0.1582 / 0.040, 0.4645 / 0.040,
          0.1571 / 0.040, 0.1656 / 0.040, 0.4627 / 0.040, 0.1774 / 0.040, 0.1645 / 0.040, 0.4708 / 0.040,
          0.1571 / 0.045, 0.1728 / 0.045, 0.4698 / 0.045, 0.1774 / 0.045, 0.1688 / 0.045, 0.4752 / 0.045,
          0.1878 / 0.045, 0.1595 / 0.045, 0.6915 / 0.045, 0.2044 / 0.045, 0.1668 / 0.045, 0.7081 / 0.045]
print("       P1 at 9, 27, 81, 45: (deviation of the alternative law)/(tolerance) from numerics-prereg/eval_prereg.out: "
      "min %.2f, max %.2f" % (min(ratios), max(ratios)))
check("these ratios run from 3.5 to 15.7", "%.1f" % min(ratios) == "3.5" and "%.1f" % max(ratios) == "15.7")

print()
if FAIL:
    print("FAILURES: %d" % len(FAIL))
    for f in FAIL:
        print("   " + f)
    sys.exit(1)
print("ALL ANALYTIC CHECKS PASSED")
