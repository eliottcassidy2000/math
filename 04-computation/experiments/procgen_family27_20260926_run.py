#!/usr/bin/env python3
"""procgen_family27_20260926_run.py -- the single runner of the family-27 lane.

Usage:  python3 -u procgen_family27_20260926_run.py [X]      (default X = 2^32)
Writes the report to stdout (redirect to 05-knowledge/results/procgen_family27_20260926.out).
Every printed claim is a check(...) that raises AssertionError on failure; lines starting with
"  |" are data (tables), not claims.

Sections
  A  the Moran function g of the owner's backward recursion and its constants
  B  the exact rise law (Ville / optional stopping) on words; the window rise spectrum; P(W)
  C  W_k = |Bad_k| exactly; THM-4495 cross-checks
  D  OEIS b-files (A006877/8, A006884/5, A060412/3, A217934): exact re-derivation and record laws
  E  exhaustive scan n <= X (C program, validated against a brute-force Python reference)
"""
import hashlib
import math
import os
import resource
import subprocess
import sys
import time
from fractions import Fraction

import mpmath as mp

EXP = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(EXP))
SCR = os.path.join(REPO, "scratch", "procgen_family27")
sys.path.insert(0, EXP)
import procgen_family27_20260926_theory as th  # noqa: E402
import procgen_family27_20260926_bfiles as bf  # noqa: E402

T0 = time.time()
NCHECK = [0]


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)
    NCHECK[0] += 1
    print("[check %3d] %s" % (NCHECK[0], msg), flush=True)


def data(msg):
    print("  | " + msg, flush=True)


PROD = [True]   # EMPIRICAL tolerances are calibrated for the production size X >= 2^31


def emp_check(cond, msg):
    if PROD[0]:
        check(cond, msg)
    else:
        data("(EMPIRICAL, not checked below X = 2^31) " + msg)


def section(title):
    print("\n==== %s ====  (t=%.1fs)" % (title, time.time() - T0), flush=True)


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def orbitT(n):
    o = [n]
    while o[-1] != 1:
        x = o[-1]
        o.append(x >> 1 if x % 2 == 0 else (3 * x + 1) >> 1)
    return o


# ======================================================================================
def section_A():
    section("A. The Moran function g(s) = 2^-s + (1/3)(3/2)^s of the backward recursion")
    g1 = Fraction(1, 2) + Fraction(1, 3) * Fraction(3, 2)
    g2 = Fraction(1, 4) + Fraction(1, 3) * Fraction(9, 4)
    check(g1 == 1 and g2 == 1, "g(1) = 1/2 + 1/2 = 1 and g(2) = 1/4 + 3/4 = 1 exactly (Fractions)")
    p0 = (Fraction(1) + Fraction(1)) / 2
    p1 = (Fraction(3, 2) + Fraction(1, 2)) / 2
    check(p0 == 1 and p1 == 1, "phi(0) = phi(1) = 1 exactly, phi(t) = ((3/2)^t + (1/2)^t)/2 (AM-fairness 3+1 = 4)")
    worst = max(abs(th.g(s) - th.phi(s - 1)) for s in [mp.mpf(i) / 7 - 3 for i in range(60)])
    check(worst < mp.mpf(10) ** -35, "g(s) = phi(s-1) identically (Lagarias-Weiss M_BP(t) = M_RRW(t+1)); max error %.1e on 60 points" % worst)
    inside = all(th.g(mp.mpf(1) + mp.mpf(i) / 100) < 1 for i in range(1, 100))
    outside = all(th.g(mp.mpf(i) / 50) > 1 for i in range(-100, 50)) and all(th.g(2 + mp.mpf(i) / 50) > 1 for i in range(1, 100))
    check(inside and outside, "g < 1 exactly on (1,2) and g > 1 on [-2,1) and (2,4] (grid of 350 points): the only roots are s = 1, 2")
    c = th.moran_constants()
    check(abs(c["R1"] - c["R1_closed"]) < 1e-25 and abs(c["R1"] - mp.mpf("6.95212")) < 1e-5,
          "residue at s=1: R1 = -1/g'(1) = 1/ln(2/sqrt3) = %s (K-L typical total stopping time 6.95212 ln n)" % mp.nstr(c["R1"], 12))
    check(abs(c["R2"] - c["R2_closed"]) < 1e-25 and abs(c["R2"] - mp.mpf("7.645")) < 1e-3,
          "residue at s=2: R2 = 1/g'(2) = 1/((3/4)ln3 - ln2) = %s (K-L time-to-peak slope 7.645)" % mp.nstr(c["R2"], 12))
    check(abs(c["sstar"] - c["sstar_closed"]) < 1e-25 and abs(c["sstar"] - 1 - mp.mpf("0.488077")) < 1e-6,
          "argmin g = 1 + lambda*, lambda* = log_3(ln2/ln(3/2)) = %s (THM-4487 Chernoff tilt)" % mp.nstr(c["sstar"] - 1, 12))
    check(abs(c["gmin"] - c["gmin_closed"]) < 1e-25 and abs(1 - c["h"] - mp.mpf("0.050044")) < 1e-6,
          "min g = 2^-(1-h), h = h(log_3 2) = %s: 1 - h = %s (glide / dimension exponent)" % (mp.nstr(c["h"], 12), mp.nstr(1 - c["h"], 8)))
    agree = max(abs(c["gamma_BP"] - c["gamma_entropy"]), abs(c["gamma_BP"] - c["gamma_LW"]), abs(c["pstar"] - c["p_at_sopt"]))
    check(agree < 1e-20, "three forms of the delay constant agree to %.1e: 1/max_s(-ln g(s)/s) = 1/(ln2 - p* ln3) with H(p*) = p* ln3 = LW fixed point gamma gLW(1/gamma) = 1" % agree)
    check(abs(c["gamma_BP"] - mp.mpf("41.677647")) < 1e-6 and abs(c["beta_BP"] - mp.mpf("0.02399")) < 1e-5 and abs(c["pstar"] - mp.mpf("0.609091")) < 2e-6,
          "gamma = %s, beta_BP = %s, ones-ratio p* = %s reproduce Lagarias-Weiss / K-L (41.677647, 0.02399, 0.609091)"
          % (mp.nstr(c["gamma_BP"], 12), mp.nstr(c["beta_BP"], 8), mp.nstr(c["pstar"], 8)))
    check(abs(c["zeta"] - c["one_minus_h34"]) < 1e-30 and abs(c["zeta"] - c["tilted_drift_bits"]) < 1e-25,
          "zeta := (3/4)log2 3 - 1 = 1 - h(3/4) = phi'(1)/ln2 = %s (junction of the window rise spectrum)" % mp.nstr(c["zeta"], 12))
    for a in (8, 12, 16, 21.238915, 30, 40):
        e1, e2 = th.delay_spectrum_LW(a), th.delay_spectrum_entropy(a)
        check(abs(e1 - e2) < 1e-25, "delay spectrum 1 - a gLW(1/a) = a(H(p) - p ln3) at a = %g: %s" % (a, mp.nstr(e1, 10)))
    check(abs(th.delay_spectrum_LW(c["R1"]) - 1) < 1e-20 and abs(th.delay_spectrum_LW(c["gamma_BP"])) < 1e-20,
          "delay spectrum equals 1 at the typical ratio R1 and 0 at gamma (records at the zero)")
    check(abs(c["glide_const"] - mp.mpf("19.982227")) < 1e-5 and abs(c["glide_const_std"] - mp.mpf("32.58961")) < 1e-4,
          "glide-record constant 1/(1-h) = %s T-steps per log2 n (standard map (1+log_3 2)/(1-h) = %s)"
          % (mp.nstr(c["glide_const"], 10), mp.nstr(c["glide_const_std"], 10)))

    def Fexc(x):
        if x < 0.3:
            return 0.0
        return 1.0 + sum(2 * (1 - 4 * k * k * x * x) * math.exp(-2 * k * k * x * x) for k in range(1, 60))
    N, hh = 100000, 10.0 / 100000
    m1 = sum((1 - Fexc((i + 0.5) * hh)) * hh for i in range(N))
    m2 = sum(2 * (i + 0.5) * hh * (1 - Fexc((i + 0.5) * hh)) * hh for i in range(N))
    check(abs(m1 - math.sqrt(math.pi / 2)) < 1e-7 and abs(m2 - math.pi ** 2 / 6) < 1e-6,
          "Brownian excursion maximum law P(M<=x) = 1 + 2 sum (1-4k^2x^2)e^(-2k^2x^2): mean %.9f = sqrt(pi/2), E M^2 = %.8f = pi^2/6" % (m1, m2))
    return c, Fexc


# ======================================================================================
def section_B(c):
    section("B. Exact rise law on words (Ville / optional stopping) and the window rise spectrum")
    cases = 0
    for k in range(1, 41):
        for i in range(1, 41):
            N, eps = th.rise_words(k, i)
            if N == 0:
                continue
            P = Fraction(N, 2 ** k)
            one = 1 - eps
            cm, cnt = th.cond_mean_overshoot(k, i)
            ok = (cnt == N and P == one / cm and cm * cm >= 2 ** i and cm * cm < Fraction(9, 4) * 2 ** i
                  and P * P * 2 ** i <= one * one and Fraction(4, 9) * one * one < P * P * 2 ** i)
            if not ok:
                raise AssertionError("rise law fails at k=%d i=%d" % (k, i))
            cases += 1
    check(cases == 920, "Theorem R on words, all %d cases k <= 40, W = 2^(i/2), i <= 40: P(tau<=k) = (1-eps_k)/E[M_tau|tau<=k] EXACTLY, W <= E[M_tau|tau<=k] < 3W/2, hence (2/3)(1-eps)/W < P <= (1-eps)/W" % cases)
    for beta in (1.02, 1.08, 1.15, 1.18, 1.2, 1.25, 1.3, 1.4, 1.5, 1.55):
        e1 = th.window_rise_spectrum(beta)
        e2 = th.window_rise_spectrum_bruteforce(beta, grid=4000)
        check(abs(e1 - e2) < 2e-4, "window rise spectrum at beta=%.2f: closed form %s = max_t[1 - t + t h((t+b)/(t log2 3))] %s" % (beta, mp.nstr(e1, 7), mp.nstr(e2, 7)))
    PW = {}
    for i in range(1, 121):
        lo, up = th.riser_probability(2 ** (i / 2))
        PW[i] = (lo, up)
        W = 2 ** (i / 2)
        if not (2.0 / 3.0 < lo * W and up * W <= 1.0 + 1e-12 and up - lo < 1e-9 * lo):
            raise AssertionError("P(W) bracket fails at i=%d" % i)
    wp = [PW[i][0] * 2 ** (i / 2) for i in range(40, 121)]
    check(True, "P(W) = P(sup_j 3^(o_j)/2^j >= W) bracketed to relative width < 1e-9 for W = 2^(i/2), i = 1..120, and 2/3 < W P(W) <= 1 in all 120 cases")
    check(0.82 < min(wp) and max(wp) < 0.845, "W P(W) for W in [2^20, 2^60] lies in [%.4f, %.4f]: the Cramer-Lundberg constant C = lim W P(W) (non-lattice walk) is about 0.83" % (min(wp), max(wp)))
    for i in (2, 4, 8, 14, 20, 40, 80, 120):
        data("W = 2^%-4g  P(W) = %.10e   W P(W) = %.6f" % (i / 2, PW[i][0], PW[i][0] * 2 ** (i / 2)))
    return PW


# ======================================================================================
def section_C():
    section("C. W_k = |Bad_k| (no-descent words) exactly")
    W = th.W_counts(1100)
    check(W[1:12] == [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128], "W_1..W_11 = 1,1,2,3,4,8,13,19,38,64,128 (THM-4479 table)")
    check([sum(W[1:T]) for T in (12, 16, 20, 24)] == [281, 2903, 31730, 367698],
          "sum_{1<=t<T} W_t = 281, 2903, 31730, 367698 for T = 12,16,20,24 (THM-4495 / THM-4487 brute-force counts)")
    h = 0.949955527188331
    norm = [math.exp(math.log(W[k]) + 1.5 * math.log(k) - h * k * math.log(2)) for k in range(500, 1101)]
    check(9.6 < min(norm) and max(norm) < 11.1, "W_k k^(3/2) 2^(-hk) in [%.4f, %.4f] for 500 <= k <= 1100 (THM-4495 window [9.6648, 11.0517])" % (min(norm), max(norm)))
    return W


# ======================================================================================
def section_D(c, W):
    section("D. OEIS record sequences headed by 27 (b-files, exact re-derivation)")
    cache = os.path.join(SCR, "oeis_cache")
    shas = bf.fetch_all(cache)
    for k in sorted(shas):
        data("%s sha256 %s" % (k, shas[k]))
    A877 = bf.read_bfile(os.path.join(cache, "b006877.txt"))
    A878 = bf.read_bfile(os.path.join(cache, "b006878.txt"))
    A884 = bf.read_bfile(os.path.join(cache, "b006884.txt"))
    A885 = bf.read_bfile(os.path.join(cache, "b006885.txt"))
    A412 = bf.read_bfile(os.path.join(cache, "b060412.txt"))
    A413 = bf.read_json_data(os.path.join(cache, "A060413.json"))
    A934 = bf.read_json_data(os.path.join(cache, "A217934.json"))
    dl = [n for _, n in A877]
    pl = [n for _, n in A884]
    gl = [n for _, n in A412]
    S = {}
    for n in set(dl) | set(pl) | set(gl) | {27, 5649499}:
        S[n] = bf.stats(n)
    check(len(A877) == len(A878) == 148 and all(S[n]["dS"] == v for (_, n), (_, v) in zip(A877, A878)),
          "all 148 A006877 delay records re-derived: standard delay = A006878 (largest n = %d)" % dl[-1])
    check(len(A885) == 97 and all(S[n]["tS"] == v for (_, n), (_, v) in zip(A884, A885)),
          "all 97 A006885 path-record maxima re-derived from A006884 (98 terms; the 98th has no A006885 value in the b-file)")
    check(len(gl) == 35 and all(S[n]["glT"] == a and S[n]["glS"] == b for n, a, b in zip(gl, A413, A934[1:])),
          "all 35 A060412 glide records re-derived: T-glide = A060413, standard glide = A217934")
    s27 = S[27]
    check((s27["dT"], s27["ones"], s27["dS"], s27["tT"], s27["tS"], s27["glT"], s27["glS"], s27["ipeak"]) == (70, 41, 111, 4616, 9232, 59, 96, 45),
          "27: T-delay 70 (41 odd steps, standard 111), T-peak 4616 at T-step 45 (standard 9232), glide 59 T / 96 standard")
    check(dl.index(27) == 8 and pl.index(27) == 5 and gl.index(27) == 3,
          "27 is the 9th delay record, the 6th path record and the 4th glide record (A006877, A006884, A060412)")
    # --- path records: count law 2(H_{X+1} - 1) and shape
    dev = []
    for k, n in enumerate(pl, start=1):
        pred = 2 * (bf.harmonic(n + 1) - 1)
        dev.append((k - pred) / math.sqrt(pred))
    check(max(abs(x) for x in dev) < 2.5, "path records: #records <= n_k versus 2(H_(n_k+1) - 1) (Frechet-scale model, independent of C): max |z| = %.2f over all 98 records; at the last record 98 vs %.2f" % (max(abs(x) for x in dev), 2 * (bf.harmonic(pl[-1] + 1) - 1)))
    for j in (3, 6, 9, 12, 15, 18, 21):
        cnt = sum(1 for n in pl if n <= 10 ** j)
        data("path records <= 10^%-2d: %3d   2(H_(X+1) - 1) = %6.2f" % (j, cnt, 2 * (bf.harmonic(10 ** j + 1) - 1)))
    rr = sorted(S[n]["tT"] / n ** 2 for n in pl if n > 10 ** 6)
    data("path records n > 1e6 (%d): median t(n)/n^2 = %.4f (independent Frechet records: median C/(2 ln 2) = %.2f for C = 0.83)" % (len(rr), rr[len(rr) // 2], 0.83 / (2 * math.log(2))))
    for name, lst in (("delay", dl), ("path", pl), ("glide", gl)):
        data("%s records: %d up to ln X = %.2f, i.e. %.3f per unit of ln X" % (name, len(lst), math.log(lst[-1]), len(lst) / math.log(lst[-1])))
    xs = [math.log(n) for n in pl if n >= 10]
    ys = [k for k, n in enumerate(pl, start=1) if n >= 10]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    slope = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)
    check(1.8 < slope < 2.1, "least-squares slope of #path records against ln n is %.3f (theta + 1 = s2 = 2)" % slope)
    rlist = [(n, S[n]["tT"] / n ** 2) for n in pl if n >= 3 and S[n]["tT"] > n * n]
    check([n for n, _ in rlist] == [27, 319804831, 1410123943, 3716509988199, 9016346070511, 1254251874774375, 10709980568908647, 1980976057694848447],
          "path records with t(n) > n^2 (rho > 2): 8 in all of A006884 (27 has r = %.3f); 10709980568908647 (r = %.3f) is not in K-L Table 3" % (rlist[0][1], rlist[6][1]))
    rhos = [(math.log(S[n]["tT"]) / math.log(n), n) for n in pl if n >= 3]
    check(max(rhos)[1] == 27 and sorted(rhos)[-2][0] < 2.1,
          "rho(n) = ln t(n)/ln n over all A006884 terms n >= 3 is maximal at n = 27 (%.6f); the runner-up is %.4f at n = %d" % (max(rhos)[0], sorted(rhos)[-2][0], sorted(rhos)[-2][1]))
    tp = sorted(S[n]["ipeak"] / math.log(S[n]["tT"] / n) for n in pl if n >= 10 ** 4)
    med = tp[len(tp) // 2]
    mean = sum(tp) / len(tp)
    check(6.5 < med < 7.7 and 6.7 < mean < 7.9, "path records n >= 1e4 (84): time-to-peak / ln(t/n) median %.3f mean %.3f vs R2 = 7.6446 (heuristic)" % (med, mean))
    # --- glide records vs the exact W_k model
    lp = [None] + [math.log2(W[L - 1]) - (L - 1) if W[L - 1] < 2 ** 1000 else (W[L - 1].bit_length() - 1 + math.log2(W[L - 1] / 2 ** (W[L - 1].bit_length() - 1))) - (L - 1) for L in range(1, 1101)]

    def Lstar(n):
        lx, best = math.log2(n), None
        for L in range(1, 1101):
            if lx + lp[L] >= 0:
                best = L
        return best
    ratios = []
    for n, g in zip(gl, A413):
        Ls = Lstar(n)
        data("glide record %-22d glide_T %5d   L*(n) %5d   ratio %.3f   glide/log2 n %.3f" % (n, g, Ls, g / Ls, g / math.log2(n)))
        if n >= 10 ** 4:
            ratios.append(g / Ls)
    check(0.78 < min(ratios) and max(ratios) < 1.18, "glide records n >= 1e4: glide / L*(n) in [%.3f, %.3f], L*(n) = max{L : n W_(L-1) 2^-(L-1) >= 1} (exact W_k model)" % (min(ratios), max(ratios)))
    check(abs(59 / Lstar(27) - 3.6875) < 1e-9, "27's glide 59 is 3.69 times its model value L*(27) = 16 (the extreme outlier of the table)")
    # --- 27's branch among the records
    def inB(n):
        return 3077 in set(S[n]["orbit"])
    fd = sum(1 for n in dl if inB(n))
    fp = sum(1 for n in pl if inB(n))
    fg = sum(1 for n in gl if inB(n))
    data("records in 27's branch B27 = Pred*(3077): delay %d/%d, path %d/%d, glide %d/%d" % (fd, len(dl), fp, len(pl), fg, len(gl)))
    check(fd == 104 and fp == 36 and fg == 13, "records lying in 27's branch: delay 104/148 (70%), path 36/98 (37%), glide 13/35 (37%)")
    check(abs(S[5649499]["dT"] / math.log(5649499) - 24.699176) < 1e-6,
          "gamma(5649499) = 24.699176 < gamma(3732423) = 24.714906: K-L Table 1 row 10 is a ones-ratio record, not a gamma record")
    return dict(dl=dl, pl=pl, gl=gl, A878=A878, A885=A885, A413=A413, A934=A934, S=S)


# ======================================================================================
def parse_scan(fn):
    d = {}
    for line in open(fn):
        t = line.split()
        if t:
            d.setdefault(t[0], []).append(t[1:])
    return d


def compare_reports(a, b):
    bad = 0
    for k in ["ROOTS", "ROOT", "BLK", "HFULL", "HWIN", "HDT", "HDS", "HGL", "HLAB", "HB27", "C27", "X", "HDTB", "HGLEN"]:
        if a.get(k) != b.get(k):
            bad += 1

    def recs(d):
        r = {}
        for t in d.get("REC", []):
            r.setdefault(t[0], []).append((int(t[1]), float(t[2])))
        return r
    ra, rb = recs(a), recs(b)
    for k in set(ra) | set(rb):
        x, y = ra.get(k, []), rb.get(k, [])
        if not (len(x) == len(y) and all(p[0] == q[0] and abs(p[1] - q[1]) < 1e-6 for p, q in zip(x, y))):
            bad += 1
    if sorted(map(tuple, a.get("RHO2", []))) != sorted(map(tuple, b.get("RHO2", []))):
        bad += 1
    ea = {t[0]: t[1] for t in a.get("EXL", [])}
    eb = {t[0]: t[1] for t in b.get("EXL", [])}
    if ea != eb:
        bad += 1
    return bad


def section_E(X, c, W, PW, D, Fexc):
    section("E. Exhaustive scan of 1 <= n <= X = %d (shortcut map T)" % X)
    os.makedirs(os.path.join(SCR, "bin"), exist_ok=True)
    src = os.path.join(EXP, "procgen_family27_20260926_scan.c")
    b_prod = os.path.join(SCR, "bin", "family27_scan")
    b_small = os.path.join(SCR, "bin", "family27_scan_small")
    subprocess.run(["cc", "-O3", "-o", b_prod, src, "-lm"], check=True)
    subprocess.run(["cc", "-O2", "-DLOGN0=12", "-o", b_small, src, "-lm"], check=True)
    # validation against the brute-force reference
    Xv = 2 ** 17
    rv = subprocess.run([b_small, str(Xv)], check=True, capture_output=True, text=True).stdout
    pv = subprocess.run([sys.executable, os.path.join(EXP, "procgen_family27_20260926_reference.py"), str(Xv)], check=True, capture_output=True, text=True).stdout
    fa, fb = os.path.join(SCR, "val_c.txt"), os.path.join(SCR, "val_py.txt")
    open(fa, "w").write(rv)
    open(fb, "w").write(pv)
    check(compare_reports(parse_scan(fa), parse_scan(fb)) == 0,
          "C scanner (N0 = 2^12, stream path exercised) = brute-force Python reference at X = 2^17: every histogram, record list, rho>2 list and excursion count identical")
    # second validation: a build with N0 = 2^20 against the production build (N0 = 2^24) at X = 2^25
    b_mid = os.path.join(SCR, "bin", "family27_scan_n20")
    subprocess.run(["cc", "-O3", "-DLOGN0=20", "-o", b_mid, src, "-lm"], check=True)
    f1, f2 = os.path.join(SCR, "val_prod_25.txt"), os.path.join(SCR, "val_n20_25.txt")
    with open(f1, "w") as f:
        subprocess.run([b_prod, str(2 ** 25)], check=True, stdout=f, stderr=subprocess.DEVNULL)
    with open(f2, "w") as f:
        subprocess.run([b_mid, str(2 ** 25)], check=True, stdout=f, stderr=subprocess.DEVNULL)
    check(compare_reports(parse_scan(f1), parse_scan(f2)) == 0,
          "production build (N0 = 2^24) = build with N0 = 2^20 at X = 2^25: identical reports (memo/stream boundary moved by a factor 16)")
    out = os.path.join(SCR, "scan_%d.txt" % X)
    t1 = time.time()
    with open(out, "w") as f:
        subprocess.run([b_prod, str(X)], check=True, stdout=f, stderr=subprocess.DEVNULL)
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    data("scan wall time %.1f s; peak RSS of child processes %.1f MB (macOS ru_maxrss bytes)" % (time.time() - t1, ru.ru_maxrss / 2 ** 20))
    check(ru.ru_maxrss < 700 * 2 ** 20, "every child process stayed below 700 MB RSS (peak %.1f MB)" % (ru.ru_maxrss / 2 ** 20))
    d = parse_scan(out)
    check(d["END"] == [[]] and int(d["X"][0][0]) == X, "scan report complete (END marker, X = %d)" % X)
    roots = [int(t[1]) for t in d["ROOT"]]
    R = len(roots)
    reach = []
    for i in range(R):
        oi = set(orbitT(roots[i]))
        reach.append([1 if roots[k] in oi else 0 for k in range(R)])
    BLK = {int(t[0]): int(t[1]) for t in d["BLK"]}
    check(sum(BLK.values()) == X, "every n <= X counted exactly once over the dyadic blocks")
    HL = {int(t[0]): list(map(int, t[1:])) for t in d["HLAB"]}
    recs = {}
    for t in d["REC"]:
        recs.setdefault(t[0], []).append((int(t[1]), t[2]))

    # ---------------- E1: records vs OEIS
    dl, pl, gl = D["dl"], D["pl"], D["gl"]
    A878 = dict(D["A878"])
    dS = [(n, int(v)) for n, v in recs["dS"]]
    exp_d = [(n, D["S"][n]["dS"]) for n in dl if n <= X]
    check(dS == exp_d, "standard-delay records n <= X from the scan = A006877 terms <= X (%d records, values = A006878)" % len(dS))
    tS = [(n, int(v)) for n, v in recs["tS"]]
    exp_p = [(n, D["S"][n]["tS"]) for n in pl if n <= X]
    check(tS == exp_p, "path records n <= X from the scan = A006884 terms <= X (%d records, maxima = A006885)" % len(tS))
    gT = [(n, int(v)) for n, v in recs["gT"]]
    gS = [(n, int(v)) for n, v in recs["gS"]]
    exp_g = [(n, D["S"][n]["glT"]) for n in gl if n <= X]
    check(gT == exp_g and [n for n, _ in gS] == [n for n, _ in gT], "glide records n <= X (T-steps and standard steps give the same %d starting values) = A060412 terms <= X, values A060413" % len(gT))
    kl = [(3, 4.551196), (7, 5.652882), (9, 5.916555), (27, 21.238915), (230631, 22.512720), (626331, 23.899366), (837799, 24.122828),
          (1723519, 24.303826), (3732423, 24.714906), (6649279, 26.479917), (8400511, 26.907006), (63728127, 32.943545)]
    gam = [(n, float(v)) for n, v in recs["gam"]]
    check(len(gam) == len([r for r in kl if r[0] <= X]) and all(a[0] == b[0] and abs(a[1] - b[1]) < 1e-6 for a, b in zip(gam, kl)),
          "records of gamma(n) = sigma_T(n)/ln n over odd 3 <= n <= X = K-L Table 1 rows <= X (%d rows, K-L row 5649499 excluded)" % len(gam))
    rho = [(n, float(v)) for n, v in recs["rho"]]
    check([n for n, _ in rho] == [3, 27] and abs(rho[1][1] - math.log(4616) / math.log(27)) < 1e-9,
          "rho(n) = ln t(n)/ln n over odd 3 <= n <= X has only the records 3 and 27: rho(27) = %.9f is the maximum over all 3 <= n <= X (even n never exceed their odd part)" % rho[1][1])
    glr = [n for n, _ in recs["glr"]]
    check(glr[:2] == [3, 27], "glide/log2 n records over odd n <= X: %s (27 held the record from 27 to %d)" % (glr, glr[2] if len(glr) > 2 else X))

    # ---------------- E2: the rho > 2 clusters
    rho2 = sorted((int(t[0]), int(t[1]), int(t[2])) for t in d["RHO2"])
    clusters = {}
    for n, t, lab in rho2:
        clusters.setdefault(t, []).append(n)
    cl = sorted(clusters.items(), key=lambda kv: min(kv[1]))
    for t, ns in cl:
        data("t(n) > n^2 cluster with peak %d: %s" % (t, ns))
    check(len(rho2) == (21 if X >= 1880165257 else len(rho2)) and cl[0][1] == [27, 31, 41, 47, 54, 55, 62, 63],
          "n <= X with t(n) > n^2: %d numbers in %d clusters (one per peak); the first is {27,31,41,47,54,55,62,63} = (27's branch) cap [1, sqrt 4616]" % (len(rho2), len(cl)))
    B27small = [n for n in range(1, 68) if 3077 in set(orbitT(n))]
    check(B27small == [27, 31, 41, 47, 54, 55, 62, 63], "27's branch below sqrt(4616) = 67.9 is exactly {27,31,41,47,54,55,62,63}")

    # ---------------- E3: Theorem R' (window rise, integers) exact on every block
    HW = {int(t[0]): list(map(int, t[1:])) for t in d["HWIN"]}
    HF = {int(t[0]): list(map(int, t[1:])) for t in d["HFULL"]}
    kmax = max(b for b in BLK if BLK[b] == 2 ** b)
    tested = exact = 0
    for k in range(1, kmax + 1):
        cum = [sum(HW[k][i:]) for i in range(len(HW[k]))]
        for i in range(1, int(2 * k * 0.585) + 2):
            lo, _ = th.rise_words(k, i)
            hi, _ = th.rise_words(k, i, shift=True)
            cc = cum[i] if i < len(cum) else 0
            if not (lo <= cc <= hi):
                raise AssertionError("window sandwich fails k=%d i=%d" % (k, i))
            tested += 1
            exact += (lo == cc)
    check(tested > 0, "Theorem R' on integers, all %d (block k <= %d, W = 2^(i/2)) cases: N_k(W) <= #{n in [2^k,2^(k+1)) : max_{j<=k} T^j n >= W n} <= N_k(W - (3/4)^k)" % (tested, kmax))
    check(exact == tested, "the lower bound is attained in all %d cases: block counts equal the exact word counts (no carry effect at these W)" % tested)
    # full orbit vs window vs P(W) at the top complete block
    kk = kmax
    rows = []
    for i in (2, 4, 6, 8, 10, 14, 20, 26, 30, 40):
        cw = sum(HW[kk][i:]) / 2 ** kk
        cf = sum(HF[kk][i:]) / 2 ** kk
        P = PW[i][0]
        rows.append((i, cw, cf, P))
        data("block %d  W = 2^%-4g window %.6e  full orbit %.6e  P(W) %.6e  full/P = %.4f  window/P = %.4f" % (kk, i / 2, cw, cf, P, cf / P, cw / P))
    dev1 = max(abs(cf / P - 1) for i, cw, cf, P in rows if i <= 20)
    dev2 = max(abs(cf / P - 1) for i, cw, cf, P in rows if 20 < i <= 40)
    emp_check(dev1 < 0.02,
          "EMPIRICAL: at block %d the FULL-orbit W-riser density equals P(W) within %.2f%% for W <= 2^10, although the window alone gives only %.0f%% of it at W = 2^10" % (kk, 100 * dev1, 100 * [r for r in rows if r[0] == 20][0][1] / PW[20][0]))
    emp_check(dev2 < 0.12, "EMPIRICAL: full-orbit W-riser density / P(W) within %.1f%% for 2^10 < W <= 2^20 at block %d" % (100 * dev2, kk))

    # ---------------- E4: Theorem G, exact glide identity in the window
    HG = {}
    for t in d["HGLEN"]:
        HG.setdefault(int(t[0]), {})[int(t[1])] = int(t[2])
    tested = 0
    for k in range(2, kmax + 1):
        for L in range(1, int(k * math.log(2) / math.log(3)) + 2):
            cnt = sum(v for LL, v in HG[k].items() if LL >= L)
            pred = 2 ** (k - L + 1) * W[L - 1]
            if cnt != pred:
                raise AssertionError("glide identity fails k=%d L=%d: %d vs %d" % (k, L, cnt, pred))
            tested += 1
    check(tested > 0, "Theorem G on integers, all %d (block k <= %d, L - 1 <= k log_3 2) cases: #{n in [2^k,2^(k+1)) : glide_T(n) >= L} = 2^(k-L+1) W_(L-1) EXACTLY" % (tested, kmax))

    # ---------------- E5: 27's branch and the nested trees along 27's orbit
    def dens(v, bl=None):
        kidx = roots.index(v)
        bl = bl or sorted(HL)
        num = sum(HL[b][l] for b in bl for l in range(R) if reach[l][kidx])
        return num / sum(BLK[b] for b in bl), num
    d3077, n3077 = dens(3077)
    blocks = [b for b in range(24, kmax + 1)]
    dloc = [dens(3077, [b])[0] for b in blocks]
    data("27's branch B27 = Pred*(3077): N(X) = %d, density %.6f; block densities %s" % (n3077, d3077, " ".join("%.6f" % x for x in dloc)))
    check(all(0.3925 < x < 0.3929 for x in dloc), "FINITE-EXACT: B27 has density in (0.3925, 0.3929) in every block 2^24..2^%d (overall %.6f at X): 27's branch is a positive-density family" % (kmax, d3077))
    # independent Monte Carlo at the top of the range (plain Python orbits, fixed seed)
    import random
    rng = random.Random(20260926)
    nmc, hit = 20000, 0
    lo_mc = 2 ** (kmax)
    for _ in range(nmc):
        x = rng.randrange(lo_mc, 2 * lo_mc)
        while x != 1 and x != 3077:
            x = x >> 1 if x % 2 == 0 else (3 * x + 1) >> 1
        hit += (x == 3077)
    pmc = hit / nmc
    se = math.sqrt(dloc[-1] * (1 - dloc[-1]) / nmc)
    check(abs(pmc - dloc[-1]) < 4 * se, "independent Monte Carlo (20000 random n in [2^%d, 2^%d), plain Python orbits): B27 fraction %.4f +- %.4f vs scan %.6f" % (kmax, kmax + 1, pmc, se, dloc[-1]))
    haar = 2 * float(c["R1"]) / 3077
    check(d3077 / haar > 80, "the Haar tree-density law cR/a predicts 2R/3077 = %.5f: the actual branch is %.0f times denser (27's orbit carries the fat trees of 31, 41, 47, ...)" % (haar, d3077 / haar))
    HB = {int(t[0]): list(map(int, t[1:])) for t in d["HB27"]}
    fr = [[x / sum(HB[b]) for x in HB[b]] for b in blocks]
    check(all(abs(x - 1 / 3) < 3e-4 for f in fr for x in f), "class fractions of B27 members mod 3 are 1/3 within 3e-4 in every block 2^24..2^%d (kappa = 1/3)" % kmax)
    bs = list(range(min(24, kmax - 3), kmax))
    s_loc = [math.log2(sum(HB[b + 1]) / sum(HB[b])) for b in bs]
    kap = [HB[b + 1][2] / sum(HB[b + 1]) for b in bs]
    resid = [abs(2 ** (-s) + k_ * 1.5 ** s - 1) for s, k_ in zip(s_loc, kap)]
    check(all(abs(s - 1) < 3e-3 for s in s_loc) and max(resid) < 3e-3,
          "Proposition B on B27: local exponents s in [%.4f, %.4f] and kappa = class-2 fraction satisfy |2^-s + kappa (3/2)^s - 1| < 3e-3" % (min(s_loc), max(s_loc)))
    orbit27 = orbitT(27)
    prof = [(j, v, dens(v)[0]) for j, v in enumerate(orbit27)]
    check(all(prof[j][2] <= prof[j + 1][2] + 1e-15 for j in range(len(prof) - 1)) and prof[-1][2] == 1.0,
          "densities of the nested trees Pred*(T^j 27) increase along the orbit and reach 1 at the root 1 (every n <= X reaches 1)")
    for j, v, dv in prof:
        if j in (0, 1, 3, 4, 10, 12, 31, 41, 43, 44, 45, 47, 57, 60, 64, 66, 67):
            data("orbit index %2d  T^j(27) = %5d  (class %d)  density of Pred*: %.6f   Haar cR/a: %.6f" % (j, v, v % 3, dv, [0, 1, 2][v % 3] * float(c["R1"]) / v))
    # exact ladder decompositions (Proposition 11 of the inverse-tree note) as count identities
    def N(v):
        return dens(v)[1]

    def chain(a):
        return sum(1 for i in range(64) if a * 2 ** i <= X)
    lad3077 = [2051, 8205, 32821, 131285, 525141, 2100565, 8402261]
    miss = N(3077) - chain(3077) - sum(N(v) for v in lad3077)
    check(0 <= miss < 1e-6 * X, "count identity N_3077(X) = #{3077 2^i <= X} + sum_{j<=6} N_(S^j 2051)(X) + (rungs j >= 7, which have S^j 2051 > 2^24): remainder %d >= 0" % miss)
    trunk = [5, 21, 85, 341, 1365, 5461, 21845, 87381, 349525, 1398101, 5592405]
    miss2 = X - 2 - chain(4) - sum(N(v) for v in trunk)
    check(0 <= miss2 < 1e-6 * X, "count identity X = 2 + #{4 2^i <= X} + sum over trunk rungs (4^j-1)/3 <= 5592405 of N_rung(X) + (higher rungs): remainder %d >= 0" % miss2)
    for name, lad in (("3077", lad3077), ("trunk", trunk), ("41", [109, 437, 1749, 6997, 27989, 111957, 447829, 1791317, 7165269])):
        data("ladder of %s: %s" % (name, "  ".join("%d:%.3e" % (v, dens(v)[0]) for v in lad)))
    d5 = dens(5)[0]
    check(abs(d5 - 0.9379) < 1e-3 and abs(dens(2051)[0] / d3077 - 1) < 3e-4,
          "owner's 4x recursion at the root: Pred*(5) carries %.4f of N and Pred*(16) %.4f (Haar model 3/4 : 1/4); along 27's orbit rung 2051 carries %.5f of 3077's tree" % (d5, 1 - d5, dens(2051)[0] / d3077))

    # ---------------- E6: like-27 counters
    C27 = {int(t[0]): list(map(int, t[1:])) for t in d["C27"]}
    P171 = th.riser_probability(4616 / 27)[0]
    top = kmax
    rise_d = C27[top][0] / 2 ** top
    emp_check(abs(rise_d / P171 - 1) < 0.01, "numbers rising like 27 (t(n) >= (4616/27) n): block-%d density %.5e vs P(170.96) = %.5e (EMPIRICAL, ratio %.4f); window-only density %.5e" % (top, rise_d, P171, rise_d / P171, C27[top][1] / 2 ** top))
    Ed = float(th.delay_spectrum_LW(70 / math.log(27)))
    bb = list(range(top - 3, top + 1))
    yy = [math.log2(C27[b][2]) for b in bb]
    mb, my_ = sum(bb) / 4, sum(yy) / 4
    ls = sum((b - mb) * (y - my_) for b, y in zip(bb, yy)) / sum((b - mb) ** 2 for b in bb)
    data("numbers with delays like 27: block counts %s" % [C27[b][2] for b in range(top - 8, top + 1)])
    emp_check(abs(ls - Ed) < 0.06, "numbers with delays like 27 (sigma_T(n) >= 21.2389 ln n): least-squares exponent over blocks %d..%d is %.3f vs Lagarias-Weiss 1 - a gLW(1/a) = %.4f" % (top - 3, top, ls, Ed))
    exp_gl = 2 ** top * W[int(math.ceil(59 * top / math.log2(27))) - 1] / 2 ** (int(math.ceil(59 * top / math.log2(27))) - 1)
    data("numbers with glides like 27 (glide_T >= 12.408 log2 n): block counts %s; model at block %d: %.2f" % ([C27[b][3] for b in range(top - 6, top + 1)], top, exp_gl))
    check(sum(C27[b][4] for b in C27 if b >= 5) == 0, "no n in [32, X] has rho(n) >= rho(27) = 2.55998 (FINITE-EXACT)")

    # ---------------- E7: excursions of long glides (numbers like 27 in the parity-word sense)
    EXL = {int(t[0]): (int(t[1]), float(t[2]), float(t[3])) for t in d["EXL"]}
    offs = []
    big = []
    for lo_, hi_ in ((20, 30), (30, 40), (40, 50), (50, 60), (60, 80), (80, 100), (100, 150), (150, 200)):
        cnt = sum(EXL[L][0] for L in EXL if lo_ <= L < hi_)
        s1 = sum(EXL[L][1] for L in EXL if lo_ <= L < hi_)
        s1L = sum(EXL[L][1] / EXL[L][0] * math.sqrt(L) * EXL[L][0] for L in EXL if lo_ <= L < hi_)
        mu = s1 / cnt
        off = (math.sqrt(math.pi / 2) * sum(math.sqrt(L) * EXL[L][0] for L in EXL if lo_ <= L < hi_) - s1L) / cnt
        offs.append(off)
        if cnt >= 50000:
            big.append(off)
        data("glide L in [%3d,%3d): %10d odd n, mean normalised excursion height u = %.4f, sqrt(L)(sqrt(pi/2) - u) = %.3f" % (lo_, hi_, cnt, mu, off))
    emp_check(len(big) >= 6 and all(0.8 < o < 1.0 for o in big), "EMPIRICAL: in every glide-length bin with >= 50000 samples (%d bins, 20 <= L < 200) the mean u = max_j ln(T^j n/n)/(sigma sqrt L) equals sqrt(pi/2) - (0.9 +- 0.1)/sqrt L: offsets %s (Brownian excursion mean plus a constant discrete offset)" % (len(big), " ".join("%.3f" % o for o in big)))
    EXH = {int(t[0]): list(map(int, t[1:])) for t in d["EXH"]}
    sig = math.sqrt((math.log(2) / math.log(3)) * (1 - math.log(2) / math.log(3))) * math.log(3)
    u27 = math.log(4616 / 27) / (sig * math.sqrt(59))
    H2 = EXH[2]
    below = sum(H2[: int(u27 / 0.02)]) + H2[int(u27 / 0.02)] * ((u27 / 0.02) % 1)
    pct = below / sum(H2)
    emp_check(0.55 < pct < 0.8 and abs(Fexc(u27) - 0.5564) < 1e-3,
          "27's glide excursion: u(27) = %.4f sits at the %.0f-th percentile of u over all n <= X with glide in [40,60), and at the %.1f-th percentile of the Brownian excursion law: a typical long glide" % (u27, 100 * pct, 100 * Fexc(u27)))

    # ---------------- E8: delay enrichment of 27's branch
    HDT = {int(t[0]): list(map(int, t[1:])) for t in d["HDT"]}
    HDTB = {int(t[0]): list(map(int, t[1:])) for t in d["HDTB"]}
    fracs = []
    for a in (0, 8, 12, 16, 20, 24):
        i = 4 * a
        num = sum(sum(HDTB[b][i:]) for b in range(20, top + 1))
        den = sum(sum(HDT[b][i:]) for b in range(20, top + 1))
        fracs.append((a, num / den, den))
        data("n in [2^20, X] with sigma_T(n) >= %2d ln n: %11d numbers, fraction in 27's branch %.4f" % (a, den, num / den))
    emp_check(fracs[-1][1] > 1.5 * fracs[0][1] and all(fracs[i][1] <= fracs[i + 1][1] + 0.02 for i in range(len(fracs) - 1)),
          "EMPIRICAL: the share of 27's branch among long-delay numbers rises from %.3f (all n) to %.3f (sigma_T >= 24 ln n): long delays favour orbits that finish along 27's long tail" % (fracs[0][1], fracs[-1][1]))
    return d


# ======================================================================================
def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 2 ** 32
    PROD[0] = X >= 2 ** 31
    print("procgen_family27_20260926_run.py  X = %d  python %s  mpmath %s" % (X, sys.version.split()[0], mp.__version__))
    for f in ("procgen_family27_20260926_run.py", "procgen_family27_20260926_theory.py", "procgen_family27_20260926_bfiles.py",
              "procgen_family27_20260926_reference.py", "procgen_family27_20260926_scan.c"):
        data("sha256 %s  %s" % (sha(os.path.join(EXP, f)), f))
    c, Fexc = section_A()
    PW = section_B(c)
    W = section_C()
    D = section_D(c, W)
    section_E(X, c, W, PW, D, Fexc)
    ru_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ru_ch = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    section("SUMMARY")
    check(ru_self < 700 * 2 ** 20, "runner peak RSS %.1f MB < 700 MB (children %.1f MB)" % (ru_self / 2 ** 20, ru_ch / 2 ** 20))
    print("ALL %d CHECKS PASSED   wall time %.1f s" % (NCHECK[0], time.time() - T0))


if __name__ == "__main__":
    main()
