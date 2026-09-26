#!/usr/bin/env python3
"""procgen_peak_20260926_run -- the single runner of the peak lane (session collatz-procgen-20260922).

Reproduce:  python3 04-computation/experiments/procgen_peak_20260926_run.py > 05-knowledge/results/procgen_peak_20260926.out
Every claim printed as 'ok:' is a check() that aborts the run on failure.  Lines starting with '[time]'
are the only non-deterministic lines.  Note: 05-knowledge/results/procgen_peak_20260926_peak_discounted_price.md.
"""
import math
import os
import resource
import sys
import time
from fractions import Fraction

import mpmath
import numpy as np
from scipy.stats import binom

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_peak_20260926_lib import (H2, M_L, Mstar, Mstar_closed_bound, block_lower_log2, brute, check,
                                       chernoff_sum_exact, chernoff_sum_log2, cq, elem_lower_norm,
                                       elem_upper_norm, floor_table, kappa, rho_band_exact, rho_exact,
                                       rho_peak_exact, rho_peak_float)
from procgen_peak_20260926_integers import run_integers
from procgen_peak_20260926_pairing import run_pairing

T0 = time.time()
RSS_UNIT = 1 if sys.platform == "darwin" else 1024


def tick(label):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * RSS_UNIT / 2 ** 20
    print(f"[time] {label}: {time.time() - T0:.1f}s elapsed, peak RSS {rss:.0f} MB", flush=True)


def sec(title):
    print("\n" + "=" * 100 + "\n" + title + "\n" + "=" * 100, flush=True)


_CACHE = {}


def RPF(q, L):
    """float layer cake, truncated at s <= 45 (q=3) or 60 (q=5,7); cached."""
    key = (q, L)
    if key not in _CACHE:
        _CACHE[key] = rho_peak_float(q, L, 45 if q == 3 else 60)
    return _CACHE[key]


# =====================================================================================================
sec("0. Constants")
for q in (3, 5, 7, 9):
    c = cq(q)
    k, z, lam, s2 = kappa(q)
    print(f"   q={q}: c=log_q 2={c:.6f}  1-H(c)={1 - H2(c):.6f}  sigma^2=c(1-c)={s2:.6f}  lambda=(1-c)/c={lam:.6f}  "
          f"z=q/max(1,lambda)={z:.6f}  kappa_q={k:.6f}")
check(abs((1 - H2(cq(3))) - 0.050044) < 5e-7, "1 - H(log_3 2) = 0.050044 (THM-4478's eta)")
check(abs((1 - H2(cq(5))) - 0.013911) < 5e-7, "1 - H(log_5 2) = 0.013911")
check(all(cq(q) * (q + 1) >= 1 for q in range(3, 102, 2)),
      "c_q (q+1) >= 1 for every odd 3 <= q <= 101 (Chernoff tilt z = q c/(1-c) >= 1; in general 2^(q+1) >= q)")
check(abs(kappa(3)[0] - 2.107580) < 1e-5 and abs(kappa(5)[0] - 2.435966) < 1e-5 and abs(kappa(7)[0] - 2.410440) < 1e-5,
      "kappa_3 = 2.107580, kappa_5 = 2.435966, kappa_7 = 2.410440 (Mogulskii constants, z_3 = 3, z_5 = 3.7824, z_7 = 3.8731)")
tick("section 0")

# =====================================================================================================
sec("1. Exact small-L validation: brute force = exact layer cake = float layer cake; the digest's numbers")
EX = {}
for q in (3, 5, 7):
    for L in range(1, 17):
        b = brute(q, L)
        e = rho_peak_exact(q, L)
        assert b["rho"] == e["rho"] and b["rho_peak"] == e["rho_peak"] and b["strata"] == e["strata"], (q, L)
        EX[(q, L)] = e
check(True, "q=3,5,7 and 1<=L<=16: enumeration of all 2^L words and the exact layer-cake DP give identical "
            "rho_L, rho^peak_L (exact rationals) and identical stratum counts |S_j|")
b18 = brute(3, 18)
e18 = rho_peak_exact(3, 18)
check(b18["rho_peak"] == e18["rho_peak"] and b18["rho"] == e18["rho"], "q=3, L=18: enumeration = exact layer cake")
worst = 0.0
for q, Ls in ((3, (24, 32, 48, 64)), (5, (24, 32, 48)), (7, (24, 32, 48))):
    for L in Ls:
        ex = rho_peak_exact(q, L)
        EX[(q, L)] = ex
        fl = rho_peak_float(q, L, 80)
        rel = abs(fl["lo"] - float(ex["rho_peak"])) / float(ex["rho_peak"])
        worst = max(worst, rel, (fl["hi"] - fl["lo"]) / fl["lo"])
check(worst < 1e-12, f"float layer cake (truncated at s<=80) agrees with the exact rationals to {worst:.1e} "
                     f"(q=3: L=24..64; q=5,7: L=24..48)")
print("   exact values (q=3):  L  rho_L  rho^peak_L  E[1/w*|bad]")
for L in (4, 8, 12, 16, 24, 32, 48, 64):
    e = EX[(3, L)]
    print(f"     {L:3d}  {float(e['rho']):.6e}  {float(e['rho_peak']):.6e}  {float(e['rho_peak'] / e['rho']):.6f}")
dig = {16: 0.0825, 60: 0.0102, 100: 3.03e-3, 200: 4.03e-4, 400: 2.79e-5}
got = {}
for L in dig:
    if L == 16:
        got[L] = float(EX[(3, 16)]["rho_peak"] / EX[(3, 16)]["rho"])
    else:
        f = RPF(3, L)
        got[L] = f["lo"] / f["rho"]
    print(f"   L={L}: E[1/w* | bad] = {got[L]:.6e}  (digest: {dig[L]:.3g})")
check(all(float(f"{got[L]:.3g}") == dig[L] for L in dig),
      "the digest's numbers E[1/w*|bad] = 0.0825, 0.0102, 3.03e-3, 4.03e-4, 2.79e-5 (L=16,60,100,200,400) are "
      "reproduced to 3 significant figures")
tick("section 1")

# =====================================================================================================
sec("2. T1: the two-sided peak theorem (upper construction on integers; lower capacity constants)")
print("2a. peak-catch modification G = 1 on E = {first peak point of every actual-bad source}, sources <= N = 10^6")
N = 10 ** 6
INT = {}
for q in (3, 5, 7):
    for L in (8, 10, 12, 14, 16):
        M = L + 2
        fl = floor_table(q, M)
        r = run_integers(q, L, N, fl, M)
        e = EX.get((q, L)) or rho_peak_exact(q, L)
        rp = float(e["rho_peak"])
        nbadw = e["nbad"]
        INT[(q, L)] = (r, rp, float(e["rho"]))
        bound = N * rp + nbadw + len(r["exc"])
        assert r["all_descend"] and r["wrong"] == 0 and r["nE"] <= r["nsrc_peak_le_N"] <= bound, (q, L)
        print(f"   q={q} L={L:2d}: bad density {r['nbad'] / N:.5f} (rho_L {float(e['rho']):.5f})  "
              f"#bad sources with peak<=N: {r['nsrc_peak_le_N'] / N:.6f}*N (rho^peak {rp:.6f})  "
              f"|E cap [1,N]| = {r['nE'] / N:.6f}*N = {r['nE'] / N / rp:.3f} rho^peak  exceptional X = {r['exc']}")
check(True, "q=3,5,7, L=8,10,12,14,16: every 2<=n<=10^6 falls below itself within L steps under G; "
            "|E cap [1,N]| <= #{bad n: peak(n) <= N} <= N rho^peak + |Bad_L| + |X| (the finite form of "
            "eps_L <= rho^peak); no lift of a bad word is actually good")
check(all(INT[(3, L)][0]["exc"] == [] for L in (8, 10, 12, 14, 16)) and
      all(INT[(5, L)][0]["exc"] == [13, 17] for L in (8, 10, 12, 14, 16)),
      "exceptional sets in [2,10^6]: empty for q=3; {13, 17} for q=5 (the minima of the 5n+1 cycles "
      "13->33->83->208->104->52->26 and 17->43->108->54->27->68->34)")
check(all(0.45 < INT[k][0]["nE"] / N / INT[k][1] < 0.85 for k in INT),
      "measured edit density = 0.45..0.85 rho^peak (shared peaks: several bad sources on one orbit share its "
      "peak point) -- EMPIRICAL")
check(all(abs(INT[k][0]["nsrc_peak_le_N"] / N / INT[k][1] - 1) < 0.03 for k in INT if k[1] <= 12),
      "for L<=12 (>= 244 lifts per residue class below N): #{bad n : peak(n) <= N} = N rho^peak (1 +- 3%): the "
      "upper-bound count is tight before collisions")

print("2b. single-scale capacity M*_L = sum_{k<L} sum_{e=e_min(k)}^{k} (floor(e/q)+1)")
MS = {}
for q in (3, 5, 7):
    row = []
    for L in (8, 16, 32, 64, 128, 256, 512):
        m = Mstar(q, L)
        MS[(q, L)] = m
        assert m <= Mstar_closed_bound(q, L)
        c = cq(q)
        row.append(f"L={L}:{m} ({m / ((1 - c * c) * L ** 3 / (6 * q)):.3f})")
    print(f"   q={q}: " + "  ".join(row))
check(True, "M*_L <= L(L+1)(2L-2+3q)/(6q) for q=3,5,7 and L = 8..512 (closed-form bound)")
check(all(abs(MS[(q, 512)] / ((1 - cq(q) ** 2) * 512 ** 3 / (6 * q)) - 1) < 0.06 for q in (3, 5, 7)),
      "M*_L ~ (1-c^2) L^3/(6q) (ratio within 6% of 1 at L=512)")

print("2c. lower bounds for eps_L: single-scale rho^peak/M*_L, stratified max_j |S_j| 2^-L/(2^{j+1} M_L(2^{j+1})), "
      "THM-4478 Theorem A max_t rho_L(2^t)/(2^t M_L(2^t))")
for q in (3, 5):
    for L in (16, 32, 64, 128):
        fr = rho_peak_float(q, L, 60, want_strata=True)
        rp = fr["lo"]
        single = rp / Mstar(q, L)
        strat = max(fr["strata"][j] / (2 ** (j + 1) * M_L(q, L, 2 ** (j + 1))) for j in fr["strata"] if fr["strata"][j] > 0)
        thA = 0.0
        tbest = None
        for t in range(0, 40):
            val = rho_band_exact(q, L, 2 ** t, 1) / 2 ** L / (2 ** t * M_L(q, L, 2 ** t))
            if val > thA:
                thA, tbest = val, t
        best = max(single, strat, thA)
        assert best <= rp
        if (q, L) in INT:
            assert best <= INT[(q, L)][0]["nE"] / N
        print(f"   q={q} L={L:3d}: rho^peak={rp:.3e}  single={single:.3e}  stratified={strat:.3e}  "
              f"ThmA={thA:.3e} (K=2^{tbest})  best lower / rho^peak = 1/{rp / best:.0f}")
        if L in (64, 128):
            tot = sum(fr["pstrata"].values())
            top = sorted(fr["pstrata"].items(), key=lambda kv: -kv[1])[:5]
            print("        peak-weighted strata (j: share of rho^peak, 2^j <= w* < 2^(j+1)): " +
                  ", ".join(f"{j}: {v / tot:.3f}" for j, v in sorted(top)))
check(True, "all three rigorous lower bounds lie below the rigorous upper bound rho^peak (and below the measured "
            "E density at L=16); the gap is the polynomial capacity factor")
# the slope condition at j = L
okL = True
for q in (3, 5, 7):
    L = 12
    fl = floor_table(q, L + 1)
    Nw = 1 << L
    idx = np.arange(Nw, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(L)[None, :]) & 1)
    e = np.zeros((Nw, L + 1), dtype=np.int64)
    e[:, 1:] = np.cumsum(bits, axis=1)
    x = e * math.log(q) - np.arange(L + 1)[None, :] * math.log(2)
    lo = np.array([0] + [int(fl[L + 1 + j]) + 1 for j in range(1, L + 1)])
    bad = np.all(e[:, 1:] >= lo[None, 1:], axis=1)
    gap = (x[bad, L] - x[bad, :L].max(axis=1)).max()
    okL &= gap <= math.log(q / 2) + 1e-12
check(okL, "every bad word has w_L <= (q/2) w* (checked by enumeration at L=12, q=3,5,7): a stratum "
           "2^j <= w* < 2^(j+1) lies in THM-4478's band B_L(q 2^j) (upper condition at j=L included) and in "
           "B'_L(2^(j+1)) (upper condition only for j<=L-1, which is all Theorem A's proof uses)")
tick("section 2")

# =====================================================================================================
sec("3. T2: the exponent 1 - H(log_q 2) for every odd q")
worst = 0
for q in (3, 5, 7, 9):
    for L in range(1, 41):
        e = EX.get((q, L)) or rho_peak_exact(q, L)
        ch = chernoff_sum_exact(q, L)
        c = cq(q)
        assert e["rho_peak"] <= ch, (q, L)
        assert float(ch) <= (q / 2) * 2 ** (-(1 - H2(c)) * L) * (1 + 1e-12), (q, L)
check(True, "exact, q=3,5,7,9, 1<=L<=40: rho^peak_L <= (q/2) sum_{e>cL} C(L,e) q^-e <= (q/2) 2^-(1-H(c))L")
for q in (3, 5, 7):
    for L in (100, 200, 400, 800):
        f = RPF(q, L)
        lc = chernoff_sum_log2(q, L)
        c = cq(q)
        assert math.log2(f["hi"]) <= lc <= math.log2(q / 2) - (1 - H2(c)) * L + 1e-9, (q, L)
check(True, "float, q=3,5,7, L=100..800: log2 rho^peak <= log2 Chernoff sum <= log2(q/2) - (1-H(c))L")
for q in (3, 5, 7):
    for L in (100, 400, 1000):
        f = RPF(q, L)
        bl = max(block_lower_log2(q, L, b) for b in range(2, 80))
        assert bl <= math.log2(f["lo"]), (q, L)
        print(f"   q={q} L={L:4d}: log2 rho^peak = {math.log2(f['lo']):10.3f}   block minorant (best b) = {bl:10.3f}   "
              f"-(1-H)L = {-(1 - H2(cq(q))) * L:9.3f}")
check(True, "the cyclic-rotation block minorant (THM-4478 section 4 with c_q, valid for c_q < 1/2) lies below "
            "rho^peak for q=3,5,7")


def catch_high_bound(q, L, W):
    """orchestrator's Proposition 4: 1/W + P[Bin(L-1,1/2): q^e < W 2^(L-1)] (exact integer sum)."""
    tot = 0
    for e in range(L):
        if q ** e < W * 2 ** (L - 1):
            tot += math.comb(L - 1, e)
    return 1 / W + tot / 2 ** (L - 1)


def theta_q(q):
    c = cq(q)
    lo, hi = 0.0, min(1 / c - 1 - 1e-12, 0.5)
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid - (1 - H2(min(0.4999999, c * (1 + mid)))) > 0:
            hi = mid
        else:
            lo = mid
    return lo


th5 = theta_q(5)
print(f"   q=5 catch-high balance exponent theta_5 = {th5:.6f}")
for L in (200, 400, 800):
    f = RPF(5, L)
    chb = catch_high_bound(5, L, 2 ** (th5 * L))
    print(f"   q=5 L={L}: catch-high bound {chb:.3e}   rho^peak {f['lo']:.3e}   (1/L)log2 rho^peak = "
          f"{math.log2(f['lo']) / L:.5f}   ratio {f['lo'] / chb:.2e}")
    assert f["lo"] < 1e-4 * chb
check(True, "q=5: the peak construction beats the catch-high bound by > 10^4 at L=200,400,800")
for q in (3, 5, 7):
    L = 1000
    f = RPF(q, L)
    dev = -math.log2(f["lo"]) / L - (1 - H2(cq(q)))
    norm = (math.log(f["lo"]) + (1 - H2(cq(q))) * L * math.log(2)) / L ** (1 / 3)
    print(f"   q={q} L=1000: -(1/L) log2 rho^peak = {-math.log2(f['lo']) / L:.6f} vs 1-H(c) = {1 - H2(cq(q)):.6f}; "
          f"[ln rho^peak + (1-H)L ln2]/L^(1/3) = {norm:.4f}")
    assert 0 < dev < 0.05 and -3.5 < norm < -1.5
check(True, "at L=1000 the per-step exponent exceeds 1-H(c) by < 0.05 bit and the excess is -Theta(L^(1/3)) "
            "(normalised log between -3.5 and -1.5 times L^(1/3)), q=3,5,7")
tick("section 3")

# =====================================================================================================
sec("4. T3: second order exp(-Theta(L^(1/3))); elementary lemma checks; sharp constant (Mogulskii) vs data")
# 4a tilting identity
for q in (3, 5):
    L = 12
    c = cq(q)
    lam = (1 - c) / c
    Nw = 1 << L
    idx = np.arange(Nw, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(L)[None, :]) & 1)
    e = np.zeros((Nw, L + 1), dtype=np.int64)
    e[:, 1:] = np.cumsum(bits, axis=1)
    S = e - c * np.arange(L + 1)[None, :]
    fl = floor_table(q, L + 1)
    lo = np.array([0] + [int(fl[L + 1 + j]) + 1 for j in range(1, L + 1)])
    bad = np.all(e[:, 1:] >= lo[None, 1:], axis=1)
    Mx = S[:, :L].max(axis=1)
    Qw = c ** e[:, L] * (1 - c) ** (L - e[:, L])
    EQ = float(np.sum((Qw * lam ** S[:, L] * float(q) ** (-Mx))[bad]))
    lhs = float(EX[(q, 12)]["rho_peak"])
    assert abs(2 ** (-(1 - H2(c)) * L) * EQ / lhs - 1) < 1e-12
check(True, "tilting identity rho^peak = 2^-(1-H)L E_Q[lambda^S_L q^-M; bad] (Q: iid Bernoulli(c) letters, "
            "S_j = e_j - jc) verified by enumeration at L=12, q=3,5")
# 4b Lemma C
mpmath.mp.dps = 40
okm = True
for q in (3, 5, 7, 9, 27, 101):
    c = mpmath.log(2) / mpmath.log(q)
    s2 = c * (1 - c)
    for n in range(1, 61):
        m4 = mpmath.fsum(mpmath.binomial(n, k) * c ** k * (1 - c) ** (n - k) * (k - c * n) ** 4 for k in range(n + 1))
        okm &= abs(m4 - (3 * n * n * s2 ** 2 + n * s2 * (1 - 6 * s2))) < mpmath.mpf(10) ** -25 * (1 + m4)
check(okm, "Lemma C(i): E_Q S_n^4 = 3 n^2 sigma^4 + n sigma^2 (1 - 6 sigma^2) (40-digit check, n<=60, "
           "q in {3,5,7,9,27,101})")
mn_sign, mn_exit = 1.0, 1.0
for q in (3, 5, 7):
    c = cq(q)
    s2 = c * (1 - c)
    FM = 20000
    fl = floor_table(q, FM)
    for n in range(1, 3001):
        if n * s2 < 1:
            continue
        f = int(fl[FM + n])
        mn_sign = min(mn_sign, binom.sf(f, n, c), binom.cdf(f, n, c))
    for a in (1, 2, 3, 5, 8, 13, 21, 34):
        n0 = math.ceil(2 * a * a / s2)
        for n in (n0, n0 + 1, 2 * n0, 3 * n0, 4 * n0):
            if n > FM:
                continue
            f = int(fl[FM + n])
            mn_exit = min(mn_exit, binom.sf(f + a, n, c) + binom.cdf(f - a, n, c))
check(mn_sign >= 1 / 16, f"Lemma C(ii): min P(S_n>0), P(S_n<0) over n sigma^2>=1, n<=3000, q=3,5,7 is "
                         f"{mn_sign:.4f} >= 1/16 (proved bound)")
check(mn_exit >= 1 / 14, f"Lemma C(iv): min P(|S_n| > a) over n sigma^2 >= 2a^2 (tested a<=34) is {mn_exit:.4f} "
                         f">= 1/14 (proved bound)")
# 4c elementary brackets
DATA = {}
for q, Ls in ((3, (25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000)),
              (5, (25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000)),
              (7, (25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000))):
    for L in Ls:
        DATA[(q, L)] = RPF(q, L)
    tick(f"layer cake q={q}")
okb = True
for q in (3, 5, 7):
    for L in (200, 400, 1000):
        f = DATA[(q, L)]
        norm = math.log(f["lo"]) + (1 - H2(cq(q))) * L * math.log(2)
        up = math.log(elem_upper_norm(q, L))
        low = elem_lower_norm(q, L)
        okb &= low <= norm <= up
        print(f"   q={q} L={L:4d}: elementary lower {low:10.2f} <= ln[2^(1-H)L rho^peak] = {norm:8.3f} <= elementary upper {up:8.3f}")
check(okb, "the explicit elementary bounds of Lemma C / Theorem 3 bracket the exact values (they are far from "
           "sharp: constants ~31.8 and ~0.36 against the true kappa_3 = 2.108)")
eta = 1 - H2(cq(3))
okr = all(DATA[(3, L)]["rho"] * 2 ** (eta * L) * L * (L + 1) >= 1 / 2.72
          for L in (25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000))
check(okr, "rho_L(3) >= 2^(-eta L)/(2.72 L(L+1)) on the tested range (proved for L>=10 by the cycle lemma)")
print("   q=3 table:   L    rho_L        rho^peak_L    ratio        ln ratio    ln ratio/L^(1/3)")
for L in (25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000):
    f = DATA[(3, L)]
    r = f["lo"] / f["rho"]
    print(f"            {L:5d}  {f['rho']:.4e}  {f['lo']:.4e}  {r:.4e}  {math.log(r):10.4f}  {math.log(r) / L ** (1 / 3):8.4f}"
          f"   (bracket rel. width {(f['hi'] - f['lo']) / f['lo']:.1e})")
for q in (5, 7):
    print(f"   q={q} table:   L    rho_L     rho^peak_L    ln[2^((1-H)L) rho^peak]   /L^(1/3)")
    for L in (25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000):
        f = DATA[(q, L)]
        nm = math.log(f["lo"]) + (1 - H2(cq(q))) * L * math.log(2)
        print(f"            {L:5d}  {f['rho']:.4f}  {f['lo']:.4e}  {nm:12.4f}  {nm / L ** (1 / 3):8.4f}")


def fit(Ls, ys, beta=None, kap=None):
    Ls = np.array(Ls, float)
    ys = np.array(ys, float)
    cols, y = [], ys.copy()
    names = []
    if kap is None:
        cols.append(-Ls ** (1 / 3))
        names.append("kappa")
    else:
        y = y + kap * Ls ** (1 / 3)
    if beta is None:
        cols.append(np.log(Ls))
        names.append("beta")
    else:
        y = y - beta * np.log(Ls)
    cols += [np.ones_like(Ls), Ls ** (-1 / 3)]
    names += ["gamma", "delta"]
    A = np.stack(cols, 1)
    sol = np.linalg.lstsq(A, y, rcond=None)[0]
    return dict(zip(names, sol)), float(np.abs(A @ sol - y).max())


k3 = kappa(3)[0]
L3 = [L for L in (100, 150, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000)]
y3 = [math.log(DATA[(3, L)]["lo"] / DATA[(3, L)]["rho"]) for L in L3]
fb, rb = fit(L3, y3, beta=2 / 3)
ff, rf = fit(L3, y3)
fk, rk = fit(L3, y3, kap=k3)
print(f"   q=3 ln(rho^peak/rho), L>=100: beta=2/3 fixed -> {fb} (max resid {rb:.4f});  free -> {ff} (max resid {rf:.4f}); "
      f" kappa fixed -> {fk} (max resid {rk:.4f})")
check(abs(fb["kappa"] / k3 - 1) < 0.05, f"q=3: with the heuristic prefactor L^(2/3) the fitted kappa = {fb['kappa']:.4f} is "
                                        f"within 5% of kappa_3 = {k3:.4f} -- EMPIRICAL")
check(abs(ff["kappa"] / k3 - 1) < 0.15, f"q=3: a free 4-parameter fit gives kappa = {ff['kappa']:.4f}, within 15% of "
                                        f"kappa_3 -- EMPIRICAL")
check(0.2 < fk["beta"] < 1.2, f"q=3: with kappa fixed at kappa_3 the fitted prefactor exponent is beta = {fk['beta']:.3f} "
                              f"(heuristic 2/3) -- EMPIRICAL")
Lh = [L for L in L3 if L >= 400]
sl = [(y3[L3.index(Lh[i + 1])] - y3[L3.index(Lh[i])]) / (Lh[i + 1] ** (1 / 3) - Lh[i] ** (1 / 3)) for i in range(len(Lh) - 1)]
pred = [-k3 + 2 / ((Lh[i] ** (1 / 3) + Lh[i + 1] ** (1 / 3)) / 2) for i in range(len(Lh) - 1)]
print("   q=3 local slopes d ln ratio / d L^(1/3): " + ", ".join(f"{s:.3f} (pred {p:.3f})" for s, p in zip(sl, pred)))
check(all(abs(s - p) < 0.08 for s, p in zip(sl, pred)),
      "q=3: local slopes between L=400 and 2000 match -kappa_3 + 2/L^(1/3) within 0.08 -- EMPIRICAL")
for q in (5, 7):
    kq = kappa(q)[0]
    Lq = [L for L in (200, 300, 400, 500, 600, 800, 1000)]
    yq = [math.log(DATA[(q, L)]["lo"]) + (1 - H2(cq(q))) * L * math.log(2) for L in Lq]
    ff, rf = fit(Lq, yq)
    fb, rb = fit(Lq, yq, beta=-5 / 6)
    kz = 1.5 * (math.pi ** 2 * cq(q) * (1 - cq(q))) ** (1 / 3) * math.log(q) ** (2 / 3)
    print(f"   q={q} ln[2^((1-H)L) rho^peak], L>=200: free -> {ff} (max resid {rf:.4f}); beta=-5/6 -> {fb} (resid {rb:.4f})")
    check(abs(ff["kappa"] / kq - 1) < 0.10 and abs(fb["kappa"] / kq - 1) < 0.05,
          f"q={q}: fitted kappa {ff['kappa']:.4f} (free) and {fb['kappa']:.4f} (beta=-5/6) agree with kappa_{q} = "
          f"{kq:.4f} (z = q c/(1-c)) within 10% / 5%, and exclude the naive z = q value {kz:.4f} -- EMPIRICAL")
    assert abs(ff["kappa"] - kz) > 0.2
tick("section 4")

# =====================================================================================================
sec("5. T5: the pairing family -- source-level rescue (THM-4475) versus partner-isolated peak catch")
N5 = 10 ** 6
CAP = 150_000_000
PR = {}
for L in (8, 10, 12, 14, 16, 18, 20, 24, 28, 32):
    f = RPF(3, L)
    afb = run_pairing(L, N5, "afb", CAP)
    pk = run_pairing(L, N5, "peak", CAP)
    PR[L] = (afb, pk, f["rho"], f["lo"])
    assert afb["fails"] == 0 and afb["stuck"] == 0 and pk["fails"] == 0 and pk["stuck"] == 0, L
    print(f"   L={L:2d}: THM-4475 rule {afb['dens_half']:.6f}   peak catch {pk['dens_half']:.6f}   2rho_L {2 * f['rho']:.6f}   "
          f"2rho^peak {2 * f['lo']:.6f}   peak/2rho^peak {pk['dens_half'] / (2 * f['lo']):.3f}   peak/2rho_L "
          f"{pk['dens_half'] / (2 * f['rho']):.3f}   (rescues {pk['rescued']}, fallback A/F/B {pk.get('fallback_afb', 0)}, "
          f"multi {pk['multi']}, mean rank of flip {pk['rank_sum'] / max(1, pk['single']):.2f})")
    tick(f"pairing L={L}")
# post-flip freshness lemma: for a prefix word u_0..u_k with u_k = 1, as n runs over its class mod 2^(k+1+l),
# b = (T^k(n) - 1)/2 runs over all residues mod 2^l (so the post-flip parity word of length l is uniform)
def T3(y):
    return y >> 1 if y % 2 == 0 else (3 * y + 1) >> 1


def parity(x, i):
    for _ in range(i):
        x = T3(x)
    return x & 1


okf = True
for k in range(0, 7):
    for pref in range(1 << (k + 1)):
        if not (pref >> k) & 1:
            continue
        # the residue r mod 2^(k+1) whose first k+1 parities are pref (bit i = parity of T^i r)
        r = next(x for x in range(1 << (k + 1)) if all(parity(x, i) == ((pref >> i) & 1) for i in range(k + 1)))
        for l in (1, 3, 6, 8):
            seen = set()
            for t in range(1 << l):
                x = r + (1 << (k + 1)) * t + (1 << (k + 1 + l)) * 7
                for _ in range(k):
                    x = T3(x)
                assert x & 1
                seen.add(((x - 1) >> 1) % (1 << l))
            okf &= len(seen) == 1 << l
check(okf, "post-flip freshness (PROVED: T^k is affine with odd slope 3^e_k on the class, so b = b_0 + 3^e_k t): "
           "for every prefix of length k+1 <= 7 ending in an odd step, b = (T^k n - 1)/2 hits every residue mod 2^l "
           "(l = 1,3,6,8) exactly once as n runs over the class -- FINITE-EXACT")
check(abs(PR[8][0]["dens_quarter"] - 0.06806) < 5e-6 and abs(PR[16][0]["dens_quarter"] - 0.02850) < 5e-6,
      "the re-implemented THM-4475 rule reproduces the audited flip densities 0.06806 (L=8) and 0.02850 (L=16) "
      "on pairs <= N/4, N = 10^6")
check(all(PR[L][1]["fails"] == 0 for L in PR),
      "partner-isolated peak catch: every 3 <= n <= 10^6 falls below itself within L (L = 8..32): a valid "
      "finite section of a member of P_L -- FINITE-EXACT")
check(all(PR[L][1]["dens_half"] < 0.5 * PR[L][0]["dens_half"] for L in PR),
      "the peak catch uses less than half the flips of THM-4475's source-level rule at every tested L")
check(all(1.0 < PR[L][1]["dens_half"] / (2 * PR[L][3]) < 1.6 for L in PR),
      "peak-catch flip density / (2 rho^peak) lies in (1.0, 1.6) for L = 8..32 -- EMPIRICAL")
rs = [PR[L][1]["dens_half"] / (2 * PR[L][2]) for L in sorted(PR)]
check(rs[-1] < 0.06 and rs[0] > 0.15,
      f"peak-catch flip density / (2 rho_L) falls from {rs[0]:.3f} (L=8) to {rs[-1]:.3f} (L=32): the pairing price "
      f"behaves like rho^peak, not like rho_L, on the tested range -- EMPIRICAL")
tick("section 5")

sec("6. Done")
tick("total")
print("ALL CHECKS PASSED")
