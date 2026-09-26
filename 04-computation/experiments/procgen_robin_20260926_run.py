#!/usr/bin/env python3
"""procgen_robin_20260926_run -- runner for the robin lane (note:
05-knowledge/results/procgen_robin_20260926_robin_inequality.md).  Every printed claim is a check(...) that aborts on
failure.  Deterministic apart from the [time] lines.

Sections
  0  constants
  1  Lemma E (free exponential-trigonometric eigenfunctions); kappa < 1
  2  Lemma Z (the Robin zone condition at width m+2): rigorous (analytic bound + interval arithmetic)
  3  Theorem 1: N_m(L) <= 2^(HL) e^(|beta|(m-c)) kappa(pi/(m+2))^L   (exact DP to L = 180, float to L = 1024)
  4  Theorem 2: A_M(L) >= 2^(HL) kappa(pi/M)^L / P_2(M)  -- its lemmas (Robbins, cycle lemma, Hoeffding, bridges,
     end lemma) and the final inequality against exact counts
  5  Corollary 3: N_m(L) <= P(m) A_(m+2)(L); the actual ratios
  6  Corollaries 4-5: private price pi_L <= poly(L) rho^peak_L; sharp second order ln rho^peak, ln pi_L
  7  Theorem 6 (RM reduction of the sharp Conjecture R to a one-walk inequality): numerical RM constants
"""
import itertools
import math
import os
import random
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp
from mpmath import iv
import numpy as np

import procgen_robin_20260926_lib as R
import procgen_pairpeak_20260926_lib as PP
import procgen_peak_20260926_lib as PK

check = R.check
T0 = time.time()
c, s2, lam, HB, eta = R.C, R.SIG2, R.LAM, R.HB, R.ETA
LN2, LN3 = math.log(2), math.log(3)


def rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 * 1024) if sys.platform == "darwin" else r / 1024


def stamp(label):
    print("[time] %s: %.1fs elapsed, peak RSS %.0f MB" % (label, time.time() - T0, rss_mb()), flush=True)


def header(s):
    print("\n" + "=" * 100 + "\n" + s + "\n" + "=" * 100, flush=True)


# ------------------------------------------------------------------------------------------------------------------
header("0. Constants")
M1 = R.M1_threshold(3000)
print("   c = %.10f  sigma^2 = %.10f  lambda = %.10f  H = %.10f bits  1-H = %.6f" % (c, s2, lam, HB, eta))
print("   kappa_3 = (3/2)(pi^2 sigma^2)^(1/3)(ln 3)^(2/3) = %.6f ; theta* = (pi^2 sigma^2/ln 3)^(1/3) = %.6f" % (R.KAPPA3, R.THETA_STAR))
print("   v_up = (2c-1/2)(3/2-2c) = %.6f ; end-lemma threshold M1 = %d" % (R.V_UP, M1))
check(abs(eta - 0.050044) < 1e-6 and abs(R.KAPPA3 - 2.107580) < 1e-6, "1 - H = 0.050044 and kappa_3 = 2.107580")
check(abs(R.KAPPA3 - (R.THETA_STAR * LN3 + math.pi ** 2 * s2 / (2 * R.THETA_STAR ** 2))) < 1e-12,
      "kappa_3 = min_theta [theta ln 3 + pi^2 sigma^2/(2 theta^2)], attained at theta*")
check(M1 == 110 and all(R.n_of_M(M) >= math.floor((M - 2) ** 2 / (8 * math.log(M - 2))) >= 4 * M >= R.nmin_of_M(M)
                        for M in range(3000, 3101)),
      "M1 = 110: n(M) >= n_min(M) for all real M in [110, 3000] (interval-safe 1/4-grid); for M >= 3000, "
      "n(M) >= (M-2)^2/(8 ln(M-2)) >= 4M >= n_min(M) (checked at M = 3000..3100; both sides monotone)")
stamp("section 0")

# ------------------------------------------------------------------------------------------------------------------
header("1. Lemma E: f(s) = e^(beta s) sin(theta s + phi) satisfies c f(s+1-c) + (1-c) f(s-c) = kappa(theta) f(s)")
mp.mp.dps = 50
cM = mp.log(2) / mp.log(3)
rng = random.Random(9142)
worst = mp.mpf(0)
for _ in range(300):
    th = mp.mpf(rng.uniform(0.001, math.pi / c * 0.999))
    ph = mp.mpf(rng.uniform(-3, 3))
    s = mp.mpf(rng.uniform(-50, 50))
    b = mp.log((1 - cM) * mp.sin(th * cM) / (cM * mp.sin(th * (1 - cM))))
    k = cM * mp.e ** (b * (1 - cM)) * mp.cos(th * (1 - cM)) + (1 - cM) * mp.e ** (-b * cM) * mp.cos(th * cM)
    f = lambda x: mp.e ** (b * x) * mp.sin(th * x + ph)
    lhs = cM * f(s + 1 - cM) + (1 - cM) * f(s - cM)
    scale = abs(f(s + 1 - cM)) + abs(f(s - cM)) + abs(f(s))
    worst = max(worst, abs(lhs - k * f(s)) / scale)
check(worst < mp.mpf(10) ** -40, "300 random (theta in (0, pi/c), phi, s): identity holds to relative error %.1e (50-digit arithmetic)" % float(worst))
betas = [R.beta(math.pi / 2 * i / 4000) for i in range(1, 4001)]
check(max(betas) < 0, "beta(theta) < 0 on (0, pi/2] (grid of 4000 points; the note's proof: ln(sin x/x) is decreasing)")
iv.dps = 30
cI = iv.log(2) / iv.log(3)


def kappa_iv(t):
    bI = iv.log((1 - cI) * iv.sin(t * cI) / (cI * iv.sin(t * (1 - cI))))
    return cI * iv.exp(bI * (1 - cI)) * iv.cos(t * (1 - cI)) + (1 - cI) * iv.exp(-bI * cI) * iv.cos(t * cI)


def zoneF_iv(t, d=2):
    bI = iv.log((1 - cI) * iv.sin(t * cI) / (cI * iv.sin(t * (1 - cI))))
    kI = cI * iv.exp(bI * (1 - cI)) * iv.cos(t * (1 - cI)) + (1 - cI) * iv.exp(-bI * cI) * iv.cos(t * cI)
    return kI * iv.sin(d * t) - 2 * (1 - cI) * iv.exp(-bI * cI) * iv.sin((d + cI) * t)


def prove_positive(fun, a, b, maxdepth=40):
    """adaptive bisection with interval arithmetic: True iff the lower endpoint of fun([x,y]) > 0 on a cover of [a,b]"""
    stack = [(mp.mpf(a), mp.mpf(b), 0)]
    nboxes = 0
    while stack:
        x, y, d = stack.pop()
        val = fun(iv.mpf([x, y]))
        nboxes += 1
        if val.a > 0:
            continue
        if d >= maxdepth:
            return False, nboxes
        mid = (x + y) / 2
        stack.append((x, mid, d + 1))
        stack.append((mid, y, d + 1))
    return True, nboxes


# kappa < 1 and a quadratic bound ln kappa <= -gamma0 theta^2 on [0.15, pi/2] by intervals; analytic on (0, 0.15]
GAMMA0 = 0.1
TH1 = 0.15
PI2_UP = mp.mpf((iv.pi / 2).b)                    # rigorous upper end: an upper bound of pi/2 (a float would fall ~6e-17 short)
ok, nb = prove_positive(lambda t: -iv.log(kappa_iv(t)) - GAMMA0 * t * t, TH1, PI2_UP)
check(ok and PI2_UP >= mp.pi / 2, "interval arithmetic (%d boxes): ln kappa(theta) < -%.2f theta^2 on [%.2f, u] with u = an outward "
                                  "enclosure of pi/2 (so the closed interval up to pi/2 is covered)" % (nb, GAMMA0, TH1))
# analytic bound on (0, TH1]: |beta| <= b1 theta^2 with b1 = c^2 (1/6 + (c TH1)^2/150), and
# kappa <= 1 - sigma^2 theta^2/2 + C4 theta^4 (second-order Taylor bounds of exp and cos; the linear |beta| terms cancel)
x1 = c * TH1
b1 = c * c * (1 / 6 + x1 * x1 / 150)
check(all(-math.log(math.sin(x) / x) <= x * x / 6 + x ** 4 / 150 for x in [i / 1000 for i in range(1, 1001)]),
      "-ln(sin x/x) <= x^2/6 + x^4/150 on (0,1] (grid; series with negative coefficients, tail <= x^4/150)")
eb = math.exp(b1 * c * TH1 ** 2)
# upper Taylor bound: e^-y <= 1 - y + y^2/2, e^y <= 1 + y + y^2 e^y/2, cos u <= 1 - u^2/2 + u^4/24 (y = |beta| times c or 1-c);
# the terms linear in |beta| combine to sigma^2 |beta| [-(2c-1) theta^2/2 + (c^4-(1-c)^4) theta^4/24]
C4 = (c * (1 - c) ** 4 + (1 - c) * c ** 4) / 24 + s2 * b1 * (c ** 4 - (1 - c) ** 4) * TH1 ** 2 / 24 \
    + s2 * b1 * b1 * ((1 - c) + c * eb) / 2
check(s2 / 2 - C4 * TH1 ** 2 > GAMMA0, "on (0, 0.15]: kappa <= 1 - (sigma^2/2 - C4 theta^2) theta^2 with C4 = %.4f, "
      "so ln kappa <= -%.4f theta^2 <= -%.2f theta^2" % (C4, s2 / 2 - C4 * TH1 ** 2, GAMMA0))
worst_ub = max(R.kappa(t) - (1 - s2 * t * t / 2 + C4 * t ** 4) for t in [TH1 * i / 500 for i in range(1, 501)])
worst_lb = max((1 - s2 * t * t / 2 - s2 * (2 * c - 1) * b1 * t ** 4 / 2) - R.kappa(t) for t in [TH1 * i / 500 for i in range(1, 501)])
check(worst_ub <= 1e-15 and worst_lb <= 1e-15,
      "numerical confirmation of the two Taylor bounds 1 - s2 t^2/2 - s2(2c-1) b1 t^4/2 <= kappa(t) <= 1 - s2 t^2/2 + C4 t^4 on (0, 0.15]")
q4 = [(math.log(R.kappa(t)) + s2 * t * t / 2) / t ** 4 for t in (0.1, 0.03, 0.01)]
print("   (ln kappa(theta) + sigma^2 theta^2/2)/theta^4 at theta = 0.1, 0.03, 0.01: %.6f %.6f %.6f" % tuple(q4))
check(all(-0.0051 < v < -0.0049 for v in q4), "ln kappa(theta) = -sigma^2 theta^2/2 - 0.00496 theta^4 + O(theta^6) (EMPIRICAL coefficient)")
stamp("section 1")

# ------------------------------------------------------------------------------------------------------------------
header("2. Lemma Z: the Robin zone condition F(theta) = kappa sin(2 theta) - 2(1-c) e^(-beta c) sin((2+c) theta) > 0 on (0, pi/4]")
G0 = 2 - 2 * (1 - c) * (2 + c)
G1 = s2 + 4 / 3 + 2 * (1 - c) * (2 + c) * b1 * c * eb
eps4 = s2 * (2 * c - 1) * b1 / 2
print("   analytic part: F(theta)/theta >= G0 - G1 theta^2 on (0, 0.15], G0 = 2(c + c^2 - 1) = %.6f, G1 = %.6f" % (G0, G1))
check(s2 / 2 * 4 / 3 - 2 * eps4 >= 0, "the dropped theta^4 coefficient (sigma^2/2)(4/3) - 2 eps4 = %.5f is >= 0" % (s2 / 2 * 4 / 3 - 2 * eps4))
check(G0 - G1 * TH1 ** 2 > 0, "G0 - G1 (0.15)^2 = %.5f > 0, so F > 0 on (0, 0.15]" % (G0 - G1 * TH1 ** 2))
PI4_UP = mp.mpf((iv.pi / 4).b)                    # rigorous upper end: an upper bound of pi/4
ok, nb = prove_positive(zoneF_iv, TH1, PI4_UP)
check(ok and PI4_UP >= mp.pi / 4, "interval arithmetic (%d boxes): F(theta) > 0 on [0.15, u], u = an outward enclosure of pi/4; "
                                  "hence the zone condition holds at theta = pi/(m+2) for every m >= 2 (m = 2 included)" % nb)
gmin = min(R.zone_F(t) / t for t in [math.pi / 4 * i / 5000 for i in range(1, 5001)])
check(abs(gmin - G0) < 1e-3, "numerically min F/theta = %.5f is attained as theta -> 0 (margin 2.9%% of the leading term)" % gmin)
print("   width shift of the optimal trigonometric Robin supersolution, delta_R(m) (zone condition at theta = pi/(m+delta)):")
dr = {m: R.delta_R(m) for m in (2, 3, 4, 5, 6, 8, 10, 16, 32, 64, 128, 1000, 10 ** 5)}
print("     " + "  ".join("m=%d:%.4f" % (m, d) for m, d in dr.items()))
dlim = 2 * c * (1 - c) / (2 * c - 1)
check(all(d < 2 for d in dr.values()) and abs(dr[10 ** 5] - dlim) < 1e-3,
      "delta_R(m) < 2 and increases to 2c(1-c)/(2c-1) = %.5f" % dlim)
# the constant e^(|beta|(m-c)) of Theorem 1
bm = max(abs(R.beta(math.pi / (m + 2))) * (m - c) for m in range(2, 20001))
bound_bm = max((((c * math.pi) ** 2) / 6 + (c * math.pi) ** 4 / (150 * (m + 2) ** 2)) * (m - c) / (m + 2) ** 2 for m in range(2, 20001))
check(bm < 0.042 and bound_bm < 0.063, "|beta(pi/(m+2))|(m-c) <= %.4f numerically (m <= 20000) and <= %.4f by |beta| <= x^2/6 + x^4/150, "
      "x = c pi/(m+2), for all m >= 2" % (bm, bound_bm))
bm1 = max(abs(R.beta(math.pi / (m + 2))) * (m + 1) for m in range(0, 20001))
bound_bm1 = max((((c * math.pi) ** 2) / 6 + (c * math.pi) ** 4 / (150 * (m + 2) ** 2)) * (m + 1) / (m + 2) ** 2 for m in range(0, 20001))
check(bm1 < 0.113 and bound_bm1 < 0.171, "|beta(pi/(m+2))|(m+1) <= %.4f numerically (0 <= m <= 20000) and <= %.4f by the same bound, for all m >= 0 "
      "(the constant of the Dirichlet strata bound)" % (bm1, bound_bm1))
stamp("section 2")

# ------------------------------------------------------------------------------------------------------------------
header("3. Theorem 1: N_m(L) <= 2^(HL) e^(|beta|(m-c)) kappa(pi/(m+2))^L")
LMAX = 400
MX = LMAX + 8
fl = PP.floor_table(MX)
NM = {}
worst = -1e9
where = None
for m in range(2, LMAX + 2):
    NM[m] = PP.refl_counts_all_L(LMAX, m, fl, MX)
    for L in range(1, LMAX + 1):
        if NM[m][L] == 0:
            continue
        r = math.log(NM[m][L]) - R.robin_bound_log(L, m)
        if r > worst:
            worst, where = r, (L, m)
check(worst < 0, "exact DP, 1 <= L <= 400, 2 <= m <= 401: ln(N_m(L)/bound) <= %.4f < 0 (max at (L,m) = %s)" % (worst, where))
for L in (256, 512, 1024, 2048, 4096):
    Mf = L + 8
    flf = PP.floor_table(Mf)
    wr = -1e9
    ms = list(range(2, 61)) + list(range(65, 121, 5))
    for m in ms:
        lr = R.robin_float_log(L, m, flf, Mf)
        wr = max(wr, lr[L] - R.robin_bound_log(L, m))
    check(wr < 0, "float DP, L = %d, 2 <= m <= 60 and m = 65..120 step 5: ln(N_m(L)/bound) <= %.3f < 0" % (L, wr))
print("   comparison at L = 180 and L = 400 (all quantities divided by 2^L):")
print("      L    m    N_m(L)/2^L       Theorem 1 bound   Theorem A bound (pairpeak)   ln(T1)/ln(exact)")
for LL in (180, 400):
    for m in (2, 3, 4, 5, 6, 8, 10, 12, 16):
        ex = NM[m][LL] / 2 ** LL
        t1 = math.exp(R.robin_bound_log(LL, m) - LL * LN2)
        ta = math.exp(R.theorem_A_log(LL, m))
        print("    %3d   %2d   %.4e        %.4e         %.4e                   %.3f" % (LL, m, ex, t1, ta, math.log(t1) / math.log(ex)))
check(all(R.robin_bound_log(LL, m) - LL * LN2 < R.theorem_A_log(LL, m) for LL in (180, 400) for m in range(2, 17)),
      "at L = 180 and L = 400 Theorem 1 is stronger than Theorem A for every 2 <= m <= 16")
# growth rates: Robin rate vs kappa(pi/(m+2)) vs the Dirichlet rate of A_(m+2)
print("   per-step rates at L = 20000 (second half): ln N_m rate, H ln2 + ln kappa(pi/(m+2)), ln A_(m+2) rate:")
Lr = 20000
Mf = Lr + 8
flf = PP.floor_table(Mf)
rows = []
for m in (3, 5, 8, 12, 16):
    lr = R.robin_float_log(Lr, m, flf, Mf)
    la = R.hard_float_log(Lr, m + 2, flf, Mf)
    h = Lr // 2
    gr = (lr[Lr] - lr[h]) / (Lr - h)
    ga = (la[Lr] - la[h]) / (Lr - h)
    gk = HB * LN2 + math.log(R.kappa(math.pi / (m + 2)))
    rows.append((m, gr, gk, ga))
    print("     m=%2d  Robin %.8f <= kappa-rate %.8f <= Dirichlet(m+2) %.8f" % (m, gr, gk, ga))
check(all(gr < gk < ga for (_, gr, gk, ga) in rows),
      "the Robin growth rate is below H ln 2 + ln kappa(pi/(m+2)), which is below the growth rate of A_(m+2) (EMPIRICAL rates)")


def eff_width(rate):
    """W with H ln 2 + ln kappa(pi/W) = rate (the effective width of a measured growth rate)"""
    lo, hi = 1.5, 1e6
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if HB * LN2 + math.log(R.kappa(math.pi / mid)) < rate:
            lo = mid
        else:
            hi = mid
    return mid


print("   effective widths W (H ln 2 + ln kappa(pi/W) = measured rate): Robin W_R(m) - m and Dirichlet W_D(m+2) - (m+2):")
for (m, gr, gk, ga) in rows:
    print("     m=%2d  W_R - m = %.4f   W_D(m+2) - (m+2) = %.4f" % (m, eff_width(gr) - m, eff_width(ga) - m - 2))
check(all(0.5 < eff_width(gr) - m < 1.2 and 0.35 < eff_width(ga) - m - 2 < 0.45 for (m, gr, gk, ga) in rows),
      "EMPIRICAL effective widths: Robin m + (0.5..1.2, increasing in m), Dirichlet M + 0.42; the trigonometric supersolution uses m + 2")
stamp("section 3")

# ------------------------------------------------------------------------------------------------------------------
header("4. Theorem 2: A_M(L) >= 2^(HL) kappa(pi/M)^L / P_2(M)")
# 4a Robbins-based binomial bound
bad = 0
for n in range(2, 401):
    for u in range(1, n):
        p = u / n
        He = -p * math.log(p) - (1 - p) * math.log(1 - p)
        if math.log(math.comb(n, u)) < n * He - math.log(2 * math.sqrt(n)) - 1e-9:
            bad += 1
check(bad == 0, "C(n,u) >= e^(n H_e(u/n))/(2 sqrt n) for all 1 <= u <= n-1 <= 399 (Robbins' Stirling bounds give the factor 0.675 >= 1/2)")
# KL bound
worst = 0
for p in [0.5 + i * (2 * c - 1) / 2000 for i in range(0, 2001)]:
    D = p * math.log(p / c) + (1 - p) * math.log((1 - p) / (1 - c))
    v = min(p * (1 - p), s2) if p <= c else p * (1 - p)
    worst = max(worst, D - (p - c) ** 2 / (2 * v))
check(worst <= 1e-15, "D(p||c) <= (p-c)^2/(2 min_{x between p and c} x(1-x)) on p in [1/2, 2c-1/2] (grid)")
# 4b cycle lemma (exhaustive)


def arrangements(n, u):
    """all distinct arrangements of u steps +(1-c) and n-u steps -c (as tuples)"""
    out = []
    for pos in itertools.combinations(range(n), u):
        a = [-c] * n
        for q in pos:
            a[q] = 1 - c
        out.append(tuple(a))
    return out


viol = 0
tested = 0
for n in range(1, 11):
    for u in range(0, n + 1):
        D = u - c * n
        arrs = arrangements(n, u)
        for thr in (0.0, 0.2):
            if D <= -thr:
                continue
            good = 0
            for a in arrs:
                s = 0.0
                okk = True
                for x in a:
                    s += x
                    if s <= -thr + 1e-12:
                        okk = False
                        break
                good += okk
            tested += 1
            if good * n < len(arrs):
                viol += 1
check(viol == 0, "cycle lemma, exhaustive for n <= 10 (%d multisets x thresholds): if the total exceeds -b (b >= 0), "
                 "at least a fraction 1/n of the arrangements have all partial sums > -b" % tested)
# 4c Hoeffding without replacement (exhaustive small check)
viol = 0
for n in range(2, 13):
    for u in range(1, n):
        arrs = arrangements(n, u)
        mu = (u - c * n) / n
        for k in range(1, n + 1):
            vals = [sum(a[:k]) - k * mu for a in arrs]
            for t in (0.5, 1.0, 1.5, 2.0):
                pu = sum(1 for v in vals if v >= t) / len(vals)
                pd = sum(1 for v in vals if v <= -t) / len(vals)
                if pu > math.exp(-2 * t * t / k) + 1e-12 or pd > math.exp(-2 * t * t / k) + 1e-12:
                    viol += 1
check(viol == 0, "Hoeffding (1963, Thm 4) for sampling without replacement, P(P_k - k mu >= t) <= exp(-2t^2/k): "
                 "confirmed exhaustively for n <= 12")
# 4d bridges and the end lemma at M = 110, 150
for M in (110.0, 150.0):
    n = R.n_of_M(M)
    nmin = R.nmin_of_M(M)
    check(n >= nmin and 4 * n * math.log(2 * n) <= (M - 2) ** 2, "M = %.0f: n(M) = %d >= n_min(M) = %d and 4 n ln(2n) <= (M-2)^2" % (M, n, nmin))
    worstb = float('inf')
    for y in (0.37, 5.1, M / 2 - 0.3, M / 2 + 0.6, M - 3.2, M - 0.2):
        # x1 = the point of y + Z - c n in [M/2, M/2+1)
        k0 = math.ceil(M / 2 - y + c * n)
        x1 = y + k0 - c * n
        assert M / 2 <= x1 < M / 2 + 1
        e = x1 + math.ceil(-x1 + c * n) - c * n
        if e <= 0:
            e += 1
        assert 0 < e <= 1 + 1e-12
        fa, _ = R.bridge_confined_fraction(y, x1, n, 0.0, M)
        fb, _ = R.bridge_confined_fraction(x1, e, n, 0.0, M)
        worstb = min(worstb, fa * 2 * n, fb * 2 * n)
    check(worstb >= 1, "M = %.0f, n = %d: for 6 starting points y the exact bridge fractions of both phases are >= 1/(2n) "
                       "(min ratio to 1/(2n): %.1f)" % (M, n, worstb))
    ys = [0.05 + i * (M - 0.1) / 59 for i in range(60)]
    qmin = R.q_bridge_min(M, n, ys)
    qlog = R.q_end_log(M, n)
    check(math.log(qmin) > qlog, "M = %.0f: min over 60 starting points of Q_y(confined 2n steps, S_2n in (0,1]) = %.3e >= q(M) = %.3e "
                                 "(the proved bound is crude by a factor e^%.0f)" % (M, qmin, math.exp(qlog), math.log(qmin) - qlog))
# 4e the Dirichlet subsolution inequality pointwise, and Theorem 2 against exact counts
rng = random.Random(4242)
viol = 0
for _ in range(20000):
    M = rng.uniform(4, 60)
    th = math.pi / M
    b = R.beta(th)
    k = R.kappa(th)
    s = rng.uniform(1e-9, M - 1e-9)
    psi = lambda x: math.exp(b * x) * math.sin(th * x)
    up = s + 1 - c
    dn = s - c
    val = (c * psi(up) if up < M else 0.0) + ((1 - c) * psi(dn) if dn > 0 else 0.0)
    if val < k * psi(s) - 1e-12:
        viol += 1
check(viol == 0, "Dirichlet subsolution: P_M psi >= kappa(pi/M) psi at 20000 random (M, s), psi = e^(beta s) sin(pi s/M)")
viol = 0
for _ in range(20000):
    W = rng.uniform(1, 60)
    th = math.pi / (W + 1)
    b = R.beta(th)
    k = R.kappa(th)
    s = rng.uniform(0, W - 1e-9)
    f = lambda x: math.exp(b * x) * math.sin(th * (x + c))
    up = s + 1 - c
    dn = s - c
    val = (c * f(up) if up < W else 0.0) + ((1 - c) * f(dn) if dn > 0 else 0.0)
    if val > k * f(s) + 1e-12:
        viol += 1
check(viol == 0, "Dirichlet supersolution: P_W f <= kappa(pi/(W+1)) f at 20000 random (W, s), f = e^(beta s) sin(pi(s+c)/(W+1))")
viol = 0
for _ in range(20000):
    m = rng.randint(2, 60)
    th = math.pi / (m + 2)
    b = R.beta(th)
    k = R.kappa(th)
    s = rng.uniform(0, m - c - 1e-9)
    g = lambda x: math.exp(b * x) * math.sin(th * (x + c))
    if s >= m - 1:
        val = 2 * (1 - c) * g(s - c)
    else:
        val = c * g(s + 1 - c) + ((1 - c) * g(s - c) if s - c > 0 else 0.0)
    if val > k * g(s) + 1e-12:
        viol += 1
check(viol == 0, "Robin supersolution: P_R g <= kappa(pi/(m+2)) g at 20000 random (m, s), g = e^(beta s) sin(pi(s+c)/(m+2))")
worst = 1e9
where = None
for M in range(4, 60):
    A = PP.hard_counts_all_L(LMAX, M, fl, MX)
    th = math.pi / M
    lk = math.log(R.kappa(th))
    p2 = R.P2_log(M, M1)
    for L in range(1, LMAX + 1):
        r = math.log(A[L]) - (HB * L * LN2 + L * lk - p2)
        if r < worst:
            worst, where = r, (L, M)
check(worst > 0, "exact DP, 4 <= M <= 59, 1 <= L <= 400: A_M(L) >= 2^(HL) kappa(pi/M)^L/P_2(M) (min log-margin %.2f at %s)" % (worst, where))
worst = 1e9
for M in range(4, 60):
    A = PP.hard_counts_all_L(LMAX, M, fl, MX)
    lk = math.log(R.kappa(math.pi / M))
    for L in range(M * M // 2, LMAX + 1):
        worst = min(worst, math.log(A[L]) - (HB * L * LN2 + L * lk) + 3 * math.log(M))
check(worst > 0, "EMPIRICAL: A_M(L) >= 2^(HL) kappa(pi/M)^L / M^3 for M^2/2 <= L <= 400 (the proved P_2 is crude)")
print("   P_2(M) (proved):  " + "  ".join("M=%d: e^%.1f" % (M, R.P2_log(M, M1)) for M in (4, 8, 16, 32, 64, 108, 109, 200, 1000, 10 ** 4, 10 ** 6)))
lp = [R.P2_log(M, M1) / math.log(M) for M in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 8)]
check(all(v < 17 for v in lp), "ln P_2(M)/ln M < 17 for M = 1e3..1e8 (asymptotically 16.81 - 7.9 ln(2 ln M)/ln M); P_2 is polynomial")
stamp("section 4")

# ------------------------------------------------------------------------------------------------------------------
header("5. Corollary 3 (Conjecture R with K = 2, up to a polynomial): N_m(L) <= P(m) A_(m+2)(L)")
worst = 0
where = None
w1 = (0, None)
for m in range(2, LMAX + 1):
    A2 = PP.hard_counts_all_L(LMAX, m + 2, fl, MX)
    A1 = PP.hard_counts_all_L(LMAX, m + 1, fl, MX)
    for L in range(1, LMAX + 1):
        if NM[m][L]:
            r = NM[m][L] / A2[L]
            if r > worst:
                worst, where = r, (L, m)
            r1 = NM[m][L] / A1[L]
            if r1 > w1[0]:
                w1 = (r1, (L, m))
check(worst < 1.0001 and all(math.log(worst) < R.P_log(m, M1) for m in range(2, 200)),
      "exact, L <= 400, all m: max N_m(L)/A_(m+2)(L) = %.7f at %s, far below the proved P(m)" % (worst, where))
check(abs(w1[0] - 1.003255) < 1e-6, "for comparison max N_m(L)/A_(m+1)(L) = %.6f at %s (pairpeak note: 1.003255 at (50,7))" % w1)
print("   behaviour in m and L (exact, L <= 400): max_L N_m/A_(m+K), its argmax, and the value at L = 400")
print("      m   K=1: max  (argmax L)  at L=400   |  K=2: max  (argmax L)  at L=400")
beh = {}
for m in (3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30):
    A1 = PP.hard_counts_all_L(LMAX, m + 1, fl, MX)
    A2 = PP.hard_counts_all_L(LMAX, m + 2, fl, MX)
    r1 = [NM[m][L] / A1[L] for L in range(1, LMAX + 1)]
    r2 = [NM[m][L] / A2[L] for L in range(1, LMAX + 1)]
    i1 = max(range(LMAX), key=lambda i: r1[i])
    i2 = max(range(LMAX), key=lambda i: r2[i])
    beh[m] = (r1[i1], i1 + 1, r1[-1], r2[i2], i2 + 1, r2[-1])
    print("     %2d   %.6f (%3d)    %.4e   |  %.7f (%3d)    %.4e" % (m, r1[i1], i1 + 1, r1[-1], r2[i2], i2 + 1, r2[-1]))
check(all(beh[m][2] < 1 and beh[m][5] < 1 for m in (3, 4, 5, 6, 7, 8, 9, 10, 12)) and beh[3][2] < 1e-5,
      "for fixed m the ratios fall below 1 at large L (eigenvalue gap; e.g. N_3/A_4 = %.1e at L = 400), and the maxima sit at L ~ m^2" % beh[3][2])
print("   proved P(m) = e^(|beta|(m-c)) P_2(m+2):  " + "  ".join("m=%d: e^%.1f" % (m, R.P_log(m, M1)) for m in (2, 5, 10, 50, 107, 200, 1000, 10 ** 4)))
Pstar = max(R.P_log(m, M1) - 17 * math.log(m + 2) for m in range(2, 10 ** 5))
check(Pstar < 0, "P(m) <= (m+2)^17 for 2 <= m < 1e5 (max of ln P(m) - 17 ln(m+2) = %.2f), and ln P(m)/ln m -> 16.81" % Pstar)
anal = max(R.P_log(M - 2, M1) - (16.83 * math.log(M) + 1.38) for M in list(range(3400, 20000, 7)) + [10 ** k for k in range(5, 13)])
check(anal < 0, "numerical confirmation of the analytic bound ln P(m) <= 16.83 ln(m+2) + 1.38 for m + 2 >= 3400 (grid to 1e12); "
                "since 16.83 ln M + 1.38 <= 17 ln M for M >= 3400, P(m) <= (m+2)^17 for every m >= 2")
stamp("section 5")

# ------------------------------------------------------------------------------------------------------------------
header("6. Corollaries 4-5: the private price and the sharp second-order term")


def Lambda(L):
    """max over real a >= 2 of -a ln 3 + L ln kappa(pi/a) (golden-section refinement of a grid maximum)"""
    fa = lambda a: -a * LN3 + L * math.log(R.kappa(math.pi / a))
    grid = [2 + 0.01 * i for i in range(0, int(100 * (4 * L ** (1 / 3) + 10)))]
    a0 = max(grid, key=fa)
    lo, hi = max(2, a0 - 0.02), a0 + 0.02
    for _ in range(80):
        m1 = lo + (hi - lo) * 0.382
        m2 = lo + (hi - lo) * 0.618
        if fa(m1) < fa(m2):
            lo = m1
        else:
            hi = m2
    return fa((lo + hi) / 2), (lo + hi) / 2


def upper_peak(L):
    """rigorous upper bound on ln(2^((1-H)L) rho^peak_L): (c/(1-c)) sum_m 3^-m e^(|beta|(m+1)) kappa(pi/(m+2))^(L-1)"""
    terms = []
    for m in range(0, L + 1):
        th = math.pi / (m + 2)
        terms.append(-m * LN3 + abs(R.beta(th)) * (m + 1) + (L - 1) * math.log(R.kappa(th)))
    tail = -(L + 1) * LN3 + math.log(1.5)          # m > L: 3^-m sum, each probability <= 1
    terms.append(tail)
    mx = max(terms)
    return math.log(c / (1 - c)) + mx + math.log(sum(math.exp(t - mx) for t in terms))


def lower_peak(L):
    """rigorous lower bound: max_a [ -a ln 3 + L ln kappa(pi/a) - ln P_2(a) ] over a real grid"""
    best = -1e18
    for i in range(0, int(100 * (4 * L ** (1 / 3) + 10))):
        a = 4 + 0.01 * i
        best = max(best, -a * LN3 + L * math.log(R.kappa(math.pi / a)) - R.P2_log(a, M1))
    return best


def upper_pi(L):
    """rigorous upper bound on ln(2^((1-H)L) pi_L) from Theorem B(2) (m_top = L) with Theorem 1 for every N"""
    terms = []
    lrho = math.log(PK.rho_float(3, L, PP.floor_table(2 * L + 8), 2 * L + 8)) + eta * L * LN2
    terms.append(math.log(2 * L) + (1 - L) * LN3 + lrho)
    for m in range(2, L):
        th = math.pi / (m + 3)
        terms.append(math.log(2 * L) + (1 - m) * LN3 + abs(R.beta(th)) * (m + 1 - c) + L * math.log(R.kappa(th)))
    th = math.pi / 4
    terms.append(math.log(2) + abs(R.beta(th)) * (2 - c) + L * math.log(R.kappa(th)))
    mx = max(terms)
    return mx + math.log(sum(math.exp(t - mx) for t in terms))


print("      L     Lambda(L)  -kappa3 L^(1/3)   lower bound   exact ln(2^(eta L) rho^peak)   upper bound   upper bound for ln(2^(eta L) pi_L)")
rows = []
for L in (100, 200, 400, 1000, 2000):
    lam_, a_ = Lambda(L)
    lo_ = lower_peak(L)
    up_ = upper_peak(L)
    rp = PK.rho_peak_float(3, L, smax=min(70.0, 3 * L ** (1 / 3) + 12))
    ex = math.log(rp["hi"]) + eta * L * LN2
    ex_lo = math.log(rp["lo"]) + eta * L * LN2
    upi = upper_pi(L)
    rows.append((L, lam_, lo_, ex_lo, ex, up_, upi))
    print("   %5d   %9.3f   %9.3f        %9.3f     %9.3f                    %9.3f     %9.3f" % (L, lam_, -R.KAPPA3 * L ** (1 / 3), lo_, ex, up_, upi))
check(all(lo_ <= ex_lo and ex <= up_ for (L, lam_, lo_, ex_lo, ex, up_, upi) in rows),
      "rigorous bracket: lower bound (Theorem 2) <= ln(2^(eta L) rho^peak_L) <= upper bound (Dirichlet supersolution) at L = 100..2000")
check(all(abs(lam_ + R.KAPPA3 * L ** (1 / 3)) < 0.05 for (L, lam_, *_) in rows),
      "Lambda(L) = -kappa_3 L^(1/3) + O(L^(-1/3)): |Lambda + kappa_3 L^(1/3)| < 0.05 on L = 100..2000")
check(all(up_ - lo_ < 17 * math.log(L) for (L, lam_, lo_, ex_lo, ex, up_, upi) in rows),
      "the rigorous bracket has width < 17 ln L (the O(log L) of Corollary 5); at L = 1000 it is [%.1f, %.1f] against the peak note's "
      "elementary [-319.9, -2.09]" % (rows[3][2], rows[3][5]))
check(all(upi - lo_ < 25 * math.log(L) for (L, lam_, lo_, ex_lo, ex, up_, upi) in rows),
      "pi_L <= L^25 rho^peak_L on this range (from the two rigorous bounds; Corollary 4 is the general statement)")
big = []
for L in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    lam_, a_ = Lambda(L)
    big.append((L, (lam_ + R.KAPPA3 * L ** (1 / 3)) * L ** (1 / 3), a_ / L ** (1 / 3)))
    print("     L=%.0e: (Lambda + kappa_3 L^(1/3)) L^(1/3) = %.5f ; argmax a / L^(1/3) = %.5f (theta* = %.5f)" % (L, big[-1][1], big[-1][2], R.THETA_STAR))
check(all(abs(v) < 2 for (_, v, _) in big) and abs(big[-1][2] - R.THETA_STAR) < 1e-2,
      "(Lambda(L) + kappa_3 L^(1/3)) L^(1/3) stays bounded and the maximiser is a ~ theta* L^(1/3)")
Pst = max(R.P_log(m, M1) for m in range(2, 1001))
print("   Corollary 4 at L = 1000 (formal): pi_L <= (9 (27 L + 18) max_(m<=L) P(m) + 4L 2^-L) rho^peak_L = e^%.1f rho^peak_L" % (
    math.log(9 * (27 * 1000 + 18)) + Pst))
stamp("section 6")

# ------------------------------------------------------------------------------------------------------------------
header("7. Theorem 6 (RM reduction): the sharp Conjecture R follows from V_k(y) - V_k(y+1) <= (2 - 2/C) F_k(y)")
KMAX = 400
ts = [t for t in range(1, 2000) if (c * t) % 1 > c][:60]
print("   max over zone points y = z_t - c (z_t = m - {ct}; the first 60 zone times t, t <= %d) and k <= %d of (V_k(y) - V_k(y+1))/F_k(y):" % (ts[-1], KMAX))
res = {}
for K in (1, 2):
    for m in (4, 6, 8, 12, 16, 20):
        M = m + K
        wr = -1e9
        for t in ts:
            z = m - ((c * t) % 1)
            y = z - c
            a = R.counts_from(y, M, KMAX)
            bq = R.counts_from(y + 1, M, KMAX)
            F = R.counts_from(y, m - 1, KMAX)
            for k in range(KMAX + 1):
                if F[k] > 0:
                    wr = max(wr, (a[k] - bq[k]) / F[k])
        res[(K, m)] = wr
        print("     K=%d m=%2d: max ratio %.5f  ->  C_min = 2/(2 - ratio) = %.4f" % (K, m, wr, 2 / (2 - wr)))
check(all(v < 0.4 for (K, m), v in res.items() if K == 1) and all(v < 0.14 for (K, m), v in res.items() if K == 2),
      "RM(1, C) holds with C = 1.25 and RM(2, C) with C = 1.075 at all tested zone points and k <= 400 (FINITE-EXACT on this set; "
      "not a proof of Conjecture R)")
# the F-slack is necessary: plain monotonicity V_k(y) <= V_k(y+1) (which would give K_0 = 1) fails for K = 1 and K = 4
mono = {}
for K, m in ((1, 12), (4, 12)):
    M = m + K
    mn = 1e9
    for t in ts[:40]:
        y = m - ((c * t) % 1) - c
        a = R.counts_from(y, M, 200)
        bq = R.counts_from(y + 1, M, 200)
        mn = min(mn, min(bq[k] / a[k] for k in range(1, 201) if a[k]))
    mono[(K, m)] = mn
    print("     K=%d m=%d: min over 40 zone points and k <= 200 of V_k(y+1)/V_k(y) = %.6f" % (K, m, mn))
check(mono[(1, 12)] < 0.87 and mono[(4, 12)] < 1,
      "plain monotonicity fails (min ratio %.4f for K = 1 and %.6f for K = 4 at m = 12): the F-slack in (RM) is needed" % (mono[(1, 12)], mono[(4, 12)]))
# direct check of the induction's conclusion N_m(L) <= C A_(m+1)(L) - (C-1) F_L(0) with C = 1.25
viol = 0
for m in range(2, 40):
    A1 = PP.hard_counts_all_L(LMAX, m + 1, fl, MX)
    F0 = R.counts_from(0.0, m - 1, LMAX)
    for L in range(1, LMAX + 1):
        if 4 * NM[m][L] > 5 * A1[L] - F0[L]:          # exact integer form of N <= 1.25 A - 0.25 F
            viol += 1
check(viol == 0, "exact, 2 <= m < 40, L <= 400: N_m(L) <= 1.25 A_(m+1)(L) - 0.25 F_L(0) (the form produced by Theorem 6)")
stamp("section 7")
print("\nALL CHECKS PASSED")
