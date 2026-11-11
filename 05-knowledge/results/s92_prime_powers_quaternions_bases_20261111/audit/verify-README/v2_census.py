#!/usr/bin/env python3
"""Verification of the corrected S92 note (README.md): the censuses quoted in sections 1 and 5, recomputed with an
independent program (own enumeration bounds, own bookkeeping; nothing imported from the session or the auditors).

Family: ordered triples (h1, h2, s), one-apex states psi(h) = h - q(h), y(h) = 2h + 1, x_b(h) = 1 + b psi(h),
    N_b = x_b(h1) x_b(h2) + b (2^s y1 y2 - [s = 0]) / 2            (b = 2: N = x1 x2 + 2^s y1 y2 - [s = 0]).
Sheet 0 is s = 0, sheet 1 is s >= 1.  Value V = N_b on sheet 1, V = 2 N_b + b on sheet 0.  For an odd modulus m prime
to b a triple with h1 < h2 has type d = gcd(x_b(h1), y(h1), m); the first-order prediction for the class c of V is
sum over d | gcd(m, c) of (d/m) D^(d) (D^(d) = number of ordered triples of type d; the diagonal h1 = h2 is counted
in the classes and not in the types, as in the session's program).

Usage:  python3 v2_census.py            (sizes up to 2^24)
        python3 v2_census.py big        (adds the census of F1 at 2^28 for the eight moduli of section 1, the octaves
                                         2^24 .. 2^28 with the counts of P4, the registered window [2^27, 2^28) by
                                         difference, the conditional flatness at 12 and the second digit of base 3
                                         at 2^26 and 2^28; needs 0.6 GB; this is the run stored in v2_census.out)
        python3 v2_census.py huge       (only the census of F1 at 2^30; needs 2.5 GB; output v2_census_30.out)
Exit status 1 if a number quoted in the note is not reproduced.
"""
import math
import sys
import time

import numpy as np
from numba import njit, prange

NT = 8
OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        OK = False


@njit(cache=True)
def build_q(size):
    q = np.zeros(size, np.int32)
    for h in range(1, size):
        k = h >> 1
        q[h] = q[k] + (((k + 1) >> 1) if (h & 1) else 0)
    return q


@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(parallel=True, cache=True)
def census(q, b, X, sheet, moduli, dindex, cnt, diag, dig2, flags, scr, vraw):
    """cnt[t, j, di, c]: ordered triples with h1 < h2 (weight 2), type index di for moduli[j], class c of V.
    diag[t, j, c]: triples with h1 = h2.  dig2[t, c]: second base-b digit (N_b - 1)/b mod b, all triples.
    flags[t, 0]: type-d triples with d not dividing V; flags[t, 1]: N_b != 1 mod b; flags[t, 2]: table too short."""
    H = q.shape[0]
    nm = moduli.shape[0]
    K = 62
    for t in prange(NT):
        s_lo = 0 if sheet == 0 else 1
        s_hi = 0 if sheet == 0 else K
        for s in range(s_lo, s_hi + 1):
            coef = 0
            if s >= 1:
                if s - 1 > 40:
                    break
                coef = b << (s - 1)
                if coef > X:
                    break
            h1 = np.int64(t)
            while h1 < H:
                y1 = 2 * h1 + 1
                lowx = 1 + b * (h1 // 2)                    # x_b(h) >= 1 + b (h // 2), since psi(h) >= h // 2
                if s == 0:
                    c11 = b * ((y1 * y1 - 1) // 2)
                else:
                    c11 = coef * y1 * y1
                if lowx * lowx + c11 > X:                    # nondecreasing in h1: no later h1 contributes
                    break
                xb1 = 1 + b * (h1 - np.int64(q[h1]))
                g1 = gcd(xb1, y1)
                for j in range(nm):
                    d = gcd(g1, moduli[j])
                    scr[t, j, 0] = dindex[j, d]
                    scr[t, j, 1] = d
                h2 = h1
                while True:
                    if h2 >= H:
                        flags[t, 2] += 1
                        break
                    y2 = 2 * h2 + 1
                    if s == 0:
                        c3 = b * ((y1 * y2 - 1) // 2)
                    else:
                        c3 = coef * y1 * y2
                    if xb1 * (1 + b * (h2 // 2)) + c3 > X:   # nondecreasing in h2
                        break
                    N = xb1 * (1 + b * (h2 - np.int64(q[h2]))) + c3
                    if N <= X:
                        V = N if (sheet == 1 or vraw) else 2 * N + b
                        if N % b != 1:
                            flags[t, 1] += 1
                        if h2 == h1:
                            dig2[t, ((N - 1) // b) % b] += 1
                            for j in range(nm):
                                diag[t, j, V % moduli[j]] += 1
                        else:
                            dig2[t, ((N - 1) // b) % b] += 2
                            for j in range(nm):
                                c = V % moduli[j]
                                cnt[t, j, scr[t, j, 0], c] += 2
                                if c % scr[t, j, 1] != 0:
                                    flags[t, 0] += 1
                    h2 += 1
                h1 += NT


def divisors(m):
    return [d for d in range(1, m + 1) if m % d == 0]


def run(q, b, K, sheet, moduli, vraw=False):
    X = 1 << K
    divs = [divisors(m) for m in moduli]
    nd = max(len(d) for d in divs)
    mm = max(max(moduli), b)
    dindex = np.zeros((len(moduli), max(moduli) + 1), np.int64)
    for j, m in enumerate(moduli):
        for i, d in enumerate(divs[j]):
            dindex[j, d] = i
    cnt = np.zeros((NT, len(moduli), nd, mm), np.int64)
    diag = np.zeros((NT, len(moduli), mm), np.int64)
    dig2 = np.zeros((NT, mm), np.int64)
    flags = np.zeros((NT, 3), np.int64)
    scr = np.zeros((NT, len(moduli), 2), np.int64)
    census(q, b, X, sheet, np.array(moduli, dtype=np.int64), dindex, cnt, diag, dig2, flags, scr, vraw)
    cnt, diag, dig2, flags = cnt.sum(axis=0), diag.sum(axis=0), dig2.sum(axis=0), flags.sum(axis=0)
    if flags[2]:
        raise SystemExit(f"table of states too short for b = {b}, 2^{K}, sheet {sheet}")
    out = dict(total=int(dig2.sum()), viol=int(flags[0]), bad=int(flags[1]),
               digit=float(np.abs(b * dig2[:b] / dig2.sum() - 1).max()), ndiag=int(diag[0].sum()), mod={})
    for j, m in enumerate(moduli):
        c = cnt[j, :len(divs[j]), :m].astype(np.float64)
        dg = diag[j, :m].astype(np.float64)
        tot = c.sum() + dg.sum()
        types = c.sum(axis=1)
        count = c.sum(axis=0) + dg
        pred = np.zeros(m)
        for cl in range(m):
            g0 = math.gcd(m, cl)
            for i, d in enumerate(divs[j]):
                if g0 % d == 0:
                    pred[cl] += d / m * types[i]
        rem = m * np.abs(count - pred).max() / tot
        fac = m * count / tot
        F = float(types[-1] / tot) if len(divs[j]) == 2 else None
        # second convention (the one of every_base.py and bridges_exact.py S10, prime moduli): class 0 against
        # 1 + (p - 1) F, the other classes against 1 - F, the diagonal triples being left in the total
        rem2 = max(abs(fac[0] - (1 + (m - 1) * F)), float(np.abs(fac[1:] - (1 - F)).max())) if F is not None else None
        out["mod"][m] = dict(rem=float(rem), rem2=rem2, fac0=float(fac[0]), F=F, fac=fac, types=types / tot, count=count)
    return out


def huge():
    """the census of F1 at 2^30 (the note quotes the audit's census, proof-5134-AE/summary_tables.out)."""
    t0 = time.time()
    q = build_q((1 << 29) + 2)
    MOD = [9, 27, 81, 25, 125, 49, 45, 75]
    want = {9: (0.01196, 0.02082), 27: (0.01637, 0.02344), 81: (0.01740, 0.02625), 25: (0.00318, 0.00341),
            125: (0.00341, 0.00414), 49: (0.05765, 0.07922), 45: (0.02050, 0.02881), 75: (0.01578, 0.02061)}
    r = [run(q, 2, 30, sh, MOD) for sh in (0, 1)]
    check((r[0]["total"], r[1]["total"]) == (3658114303, 4649378323), f"2^30: triples sheet 0 / sheet 1 = {r[0]['total']} / {r[1]['total']}")
    check(r[0]["viol"] == 0 and r[1]["viol"] == 0, "2^30: every triple of type d has d | N - c_sigma (0 violations)")
    line = ";  ".join(f"{m}: {r[0]['mod'][m]['rem']:.5f} / {r[1]['mod'][m]['rem']:.5f}" for m in MOD)
    print(f"    remainder of the first-order law (times m), sheet 0 / sheet 1, 2^30: {line}")
    check(all(abs(r[sh]["mod"][m]["rem"] - want[m][sh]) < 6e-6 for m in MOD for sh in (0, 1)),
          "2^30: the sixteen remainders of the audit's summary_tables.out are reproduced to five decimals")
    check((round(r[0]["mod"][9]["rem"], 3), round(r[1]["mod"][9]["rem"], 3), round(r[0]["mod"][49]["rem"], 3), round(r[1]["mod"][49]["rem"], 3))
          == (0.012, 0.021, 0.058, 0.079), "note: 0.012 / 0.021 at 9 and 0.058 / 0.079 at 49")
    # P3-type comparison at 49: measured factor against the first-order value and the limit, groups g0 = 7 and 49
    for sh in (0, 1):
        o = r[sh]["mod"][49]
        F1, F7, F49 = (float(v) for v in o["types"])
        fac = o["fac"]
        g7 = float(np.mean([fac[c] for c in range(49) if c % 7 == 0 and c != 0]))
        g49 = float(fac[0])
        phi7, phi49 = F1 + 7 * F7, F1 + 7 * F7 + 49 * F49
        lim7, lim49 = (1 - 1 / 49) * (1 + 1 / 7), 1 + 1 / 7 - 1 / 343
        closer = [abs(g7 - phi7) < abs(g7 - lim7), abs(g49 - phi49) < abs(g49 - lim49)]
        print(f"    49, sheet {sh}: g0 = 7 measured {g7:.4f}, first-order {phi7:.4f}, limit {lim7:.4f}; g0 = 49 measured {g49:.4f}, "
              f"first-order {phi49:.4f}, limit {lim49:.4f}; closer to the first-order value: {closer}")
        check(closer == [False, False], f"2^30, sheet {sh}: at 49 both groups with g0 > 1 are closer to the limit than to the first-order value")
    print(f"({time.time() - t0:.0f} s)")
    print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
    sys.exit(0 if OK else 1)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "huge":
        huge()
    big = len(sys.argv) > 1 and sys.argv[1] == "big"
    t0 = time.time()
    q = build_q((1 << 27) + 2 if big else (1 << 23) + 2)
    print(f"table of q(h) for h < {q.shape[0]}  ({time.time() - t0:.1f} s)")

    print("1. F1 (b = 2) at odd moduli, section 1 of the note (stored outputs cum_22.out, cum_24.out of the session)")
    MOD = [9, 27, 81, 25, 125, 49, 45, 75]
    want_tot = {22: (10532579, 13218857), 24: (45886662, 57818155)}
    want_rem24 = {9: (0.02075, 0.02420), 27: (0.02229, 0.02904), 81: (0.02410, 0.04160), 25: (0.00373, 0.00478),
                  125: (0.00761, 0.00840), 49: (0.16339, 0.14140), 45: (0.03974, 0.05634), 75: (0.05046, 0.05598)}
    for K in (22, 24):
        r = [run(q, 2, K, sh, MOD) for sh in (0, 1)]
        check((r[0]["total"], r[1]["total"]) == want_tot[K], f"2^{K}: triples sheet 0 / sheet 1 = {r[0]['total']} / {r[1]['total']}")
        check(r[0]["viol"] == 0 and r[1]["viol"] == 0, f"2^{K}: every triple of type d has d | N - c_sigma (0 violations)")
        line = ";  ".join(f"{m}: {r[0]['mod'][m]['rem']:.5f} / {r[1]['mod'][m]['rem']:.5f}" for m in MOD)
        print(f"    remainder of the first-order law (times m), sheet 0 / sheet 1, 2^{K}: {line}")
        if K == 24:
            good = all(abs(r[sh]["mod"][m]["rem"] - want_rem24[m][sh]) < 6e-6 for m in MOD for sh in (0, 1))
            check(good, "2^24: the sixteen remainders of cum_24.out are reproduced to five decimals")
            check([round(r[0]["mod"][m]["rem"], 3) for m in (9, 27, 81)] == [0.021, 0.022, 0.024]
                  and (round(r[0]["mod"][49]["rem"], 3), round(r[1]["mod"][49]["rem"], 3)) == (0.163, 0.141),
                  "note: sheet 0 at 2^24: 0.021, 0.022, 0.024 at 9, 27, 81; 49: 0.163 / 0.141")
            g9 = [(round(float(r[sh]["mod"][9]["fac"][0]), 4)) for sh in (0, 1)]
            print(f"    factor of the distinguished class modulo 9 at 2^24: sheet 0 {g9[0]}, sheet 1 {g9[1]} (cum_24.out: 1.4532, 1.4561)")
            check(g9 == [1.4532, 1.4561], "the distinguished class modulo 9 at 2^24 as in cum_24.out")

    if big:
        want_rem28 = {9: (0.02383, 0.01854), 27: (0.02921, 0.02708), 81: (0.03138, 0.03081), 25: (0.00354, 0.00374),
                      125: (0.00477, 0.00483), 49: (0.09316, 0.10282), 45: (0.03779, 0.03549), 75: (0.02744, 0.02939)}
        r = [run(q, 2, 28, sh, MOD) for sh in (0, 1)]
        check((r[0]["total"], r[1]["total"]) == (854415478, 1083261684), f"2^28: triples sheet 0 / sheet 1 = {r[0]['total']} / {r[1]['total']}")
        check(r[0]["viol"] == 0 and r[1]["viol"] == 0, "2^28: every triple of type d has d | N - c_sigma (0 violations)")
        line = ";  ".join(f"{m}: {r[0]['mod'][m]['rem']:.5f} / {r[1]['mod'][m]['rem']:.5f}" for m in MOD)
        print(f"    remainder of the first-order law (times m), sheet 0 / sheet 1, 2^28: {line}")
        check(all(abs(r[sh]["mod"][m]["rem"] - want_rem28[m][sh]) < 6e-6 for m in MOD for sh in (0, 1)),
              "2^28: the sixteen remainders of cum_28.out are reproduced to five decimals")
        g9 = [round(float(r[sh]["mod"][9]["fac"][0]), 4) for sh in (0, 1)]
        check(g9 == [1.4511, 1.4619], f"2^28: factor of the distinguished class modulo 9: {g9[0]} (sheet 0), {g9[1]} (sheet 1); limit 35/27 = 1.2963")
        g25 = [[round(float(r[sh]["mod"][25]["fac"][c]), 4) for c in (0,)] for sh in (0, 1)]
        check(g25 == [[1.2392], [1.2480]], f"2^28: factor of the distinguished class modulo 25: {g25[0][0]}, {g25[1][0]} (cum_28.out: 1.2392, 1.2480)")
        noise = [math.sqrt(m / r[1]["total"]) for m in MOD]
        check(all(min(r[0]["mod"][m]["rem"], r[1]["mod"][m]["rem"]) > 10 * nz for m, nz in zip(MOD, noise)),
              "2^28: every remainder is more than ten times the sampling scale sqrt(m/#T)")

    if big:
        print("1a. the octaves 2^24 .. 2^28 (section 1: oscillation at the powers of 3, the count of P4)")
        PP = [9, 27, 81, 25, 125, 49]
        rem = {}
        for K in (24, 25, 26, 27, 28):
            rr = [run(q, 2, K, sh, PP) for sh in (0, 1)]
            for m in PP:
                for sh in (0, 1):
                    rem[(K, m, sh)] = rr[sh]["mod"][m]["rem"]
            print(f"    2^{K}: " + ";  ".join(f"{m}: {rem[(K, m, 0)]:.5f} / {rem[(K, m, 1)]:.5f}" for m in PP))
        p3 = [rem[(K, m, sh)] for K in (24, 25, 26, 27, 28) for m in (9, 27, 81) for sh in (0, 1)]
        check(round(min(p3), 3) == 0.018 and round(max(p3), 3) == 0.042,
              f"at 9, 27, 81 the remainder stays between {min(p3):.5f} and {max(p3):.5f} over the octaves 2^24 .. 2^28 (note: 0.018 and 0.042)")
        check(all(rem[(28, m, 0)] > rem[(24, m, 0)] for m in (9, 27, 81)) and all(rem[(28, m, 1)] < rem[(24, m, 1)] for m in (9, 27, 81)),
              "from 2^24 to 2^28 it grows at 9, 27, 81 on sheet 0 (no net decrease) and decreases on sheet 1")
        check(all(rem[(28, m, sh)] < rem[(24, m, sh)] for m in (25, 125, 49) for sh in (0, 1)), "it decreases at 25, 125, 49 from 2^24 to 2^28 on both sheets")

        def smaller(K0, K1):
            return sum(1 for m in PP for sh in (0, 1) if rem[(K1, m, sh)] < rem[(K0, m, sh)])

        check((smaller(26, 28), smaller(26, 27), smaller(27, 28)) == (8, 4, 10),
              f"P4-type counts: {smaller(26, 28)} of 12 from 2^26 to 2^28, {smaller(26, 27)} of 12 from 2^26 to 2^27, {smaller(27, 28)} of 12 from 2^27 to 2^28")
        grew = sorted((m, sh) for m in PP for sh in (0, 1) if rem[(28, m, sh)] >= rem[(26, m, sh)])
        check(grew == [(9, 0), (27, 1), (81, 1), (125, 0)], f"from 2^26 to 2^28 the remainder grew at (modulus, sheet) = {grew}")
        dec = sorted(round(rem[(26, m, sh)] - rem[(28, m, sh)], 5) for m in PP for sh in (0, 1) if rem[(28, m, sh)] < rem[(26, m, sh)])
        check(dec[:2] == [0.00013, 0.00013], f"the two smallest of the eight decreases are {dec[0]:.5f} and {dec[1]:.5f} (noise scale sqrt(m/#T) = 0.00015 to 0.00018 there)")

        print("1a'. the registered window [2^27, 2^28) by difference of two cumulative censuses (window_prereg_28_27.out)")
        WM = [9, 25, 49]
        cw = {}
        for sh in (0, 1):
            hi = run(q, 2, 28, sh, WM)
            lo = run(q, 2, 27, sh, WM)
            for m in WM:
                cw[(m, sh)] = hi["mod"][m]["count"] - lo["mod"][m]["count"]
            if sh == 0:
                check(hi["total"] - lo["total"] == 442236323, f"window: {hi['total'] - lo['total']} triples on sheet 0")
            else:
                check(hi["total"] - lo["total"] == 561401507, f"window: {hi['total'] - lo['total']} triples on sheet 1")
        # registered a-priori predictions (PREREGISTRATION.md, read only) and limits, groups g0 > 1
        reg = {(9, 3): (1.2143, 1.2156, 32 / 27), (9, 9): (1.4249, 1.4377, 35 / 27),
               (25, 5): (1.2308, 1.2412, 144 / 125), (25, 25): (1.2409, 1.2496, 149 / 125),
               (49, 7): (1.1098, 1.0967, 384 / 343), (49, 49): (1.1188, 1.1051, 391 / 343)}
        want_meas = {(9, 3): (1.2214, 1.2272), (9, 9): (1.4483, 1.4645), (25, 5): (1.2299, 1.2398), (25, 25): (1.2381, 1.2460),
                     (49, 7): (1.1691, 1.1744), (49, 49): (1.1949, 1.1959)}
        good_meas, closer = True, {}
        for (m, g0), (p0, p1, lim) in reg.items():
            for sh, pred in ((0, p0), (1, p1)):
                c = cw[(m, sh)]
                tot = c.sum()
                cls = [x for x in range(m) if math.gcd(m, x) == g0]
                meas = float(np.mean([m * c[x] / tot for x in cls]))
                good_meas &= round(meas, 4) == want_meas[(m, g0)][sh]
                closer[(m, g0, sh)] = abs(meas - pred) < abs(meas - lim)
        check(good_meas, "window: the twelve group means at 9, 25, 49 (g0 > 1) of window_prereg_28_27.out are reproduced to four decimals")
        check(all(closer[k] for k in closer if k[0] in (9, 25)), "window: at 9 and 25 every group with g0 > 1 is closer to the registered prediction than to the limit")
        check(not any(closer[k] for k in closer if k[0] == 49), "window: at 49 none of the four groups with g0 > 1 is (the comparison of P3 fails there)")

    print("1b. conditional flatness at 12 (section 11: 'flat to 0.0002 at 12 at 2^28'; numbers of the audit, HYP-9665 (A))")
    want_flat = {22: (0.0021, 0.0014), 24: (0.0009, 0.0005), 26: (0.0005, 0.0001), 28: (0.0002, 0.0002)}
    for K in (22, 24) + ((26, 28) if big else ()):
        fl = []
        for sh in (0, 1):
            cn = run(q, 2, K, sh, [12], vraw=True)["mod"][12]["count"]       # classes of N itself modulo 12
            worst = 0.0
            for e in range(3):
                tot_e = sum(cn[c] for c in range(12) if c % 3 == e)
                for a in (1, 3):
                    cls = [c for c in range(12) if c % 3 == e and c % 4 == a][0]
                    worst = max(worst, abs(2 * cn[cls] / tot_e - 1))
            fl.append(worst)
        print(f"    2^{K}: R_flat at 12 = {fl[0]:.4f} / {fl[1]:.4f}")
        check((round(fl[0], 4), round(fl[1], 4)) == want_flat[K], f"2^{K}: {want_flat[K][0]:.4f} / {want_flat[K][1]:.4f}")

    print("2. base b, sheet 1, N_b <= 2^24 (every_base.out (4))")
    want1 = {2: 57818155, 3: 35564915, 4: 24947862, 5: 18838778, 6: 14919260, 10: 7634213}
    want_rem7 = {2: 0.133, 3: 0.067, 4: 0.125, 5: 0.042, 6: 0.105, 10: 0.060}
    want_dig = {2: 0.000, 3: 0.040, 4: 0.001, 5: 0.002, 6: 0.052, 10: 0.003}
    res1 = {}
    for b in (2, 3, 4, 5, 6, 10):
        primes = [p for p in (3, 5, 7, 11, 13) if b % p]
        r = run(q, b, 24, 1, primes)
        res1[b] = r
        cells = "; ".join(f"p={p}: class 0 {r['mod'][p]['fac0']:.3f}, F p^2 = {r['mod'][p]['F'] * p * p:.2f}, rem {r['mod'][p]['rem']:.3f}" for p in primes)
        print(f"    b = {b:2d}: {r['total']} triples; second digit flat to {r['digit']:.4f}; {cells}")
        check(r["total"] == want1[b] and r["viol"] == 0 and r["bad"] == 0,
              f"b = {b}: {want1[b]} triples, N_b = 1 (mod b), every p-divisible smaller state in the class 0")
        check(round(r["mod"][7]["rem2"], 3) == want_rem7[b] and abs(r["mod"][7]["rem"] - r["mod"][7]["rem2"]) < 1e-3,
              f"b = {b}: remainder at 7 is {want_rem7[b]:.3f} ({r['mod'][7]['rem2']:.5f}; {r['mod'][7]['rem']:.5f} with the other treatment of the diagonal)")
        check(abs(r["digit"] - want_dig[b]) < 6e-4, f"b = {b}: second digit flat to {want_dig[b]:.3f}")
    check(max(r["mod"][p]["rem2"] for b, r in res1.items() for p in r["mod"]) < 0.1335
          and max(r["mod"][p]["rem2"] for b, r in res1.items() for p in r["mod"] if p != 7) <= 0.0515,
          "sheet 1: remainder at most 0.133 (at 7), at most 0.051 at the other primes")
    fp2 = (res1[3]["mod"][7]["F"] * 49, res1[10]["mod"][7]["F"] * 49, res1[5]["mod"][13]["F"] * 169)
    check([round(v, 2) for v in fp2] == [3.37, 3.28, 3.38], f"early boosts: F p^2 = {fp2[0]:.2f} (b = 3, p = 7), {fp2[1]:.2f} (b = 10, p = 7), {fp2[2]:.2f} (b = 5, p = 13)")
    check(round(res1[10]["digit"], 4) == 0.0034, f"b = 10: second digit flat to {res1[10]['digit']:.4f} (the corrected 0.0034 of HYP-9665 (C))")

    print("3. base b, sheet 0, N_b <= 2^22 (bridges_exact.out S10)")
    want0 = {2: 10532579, 3: 6030660, 5: 2896626, 10: 1042874}
    want0_rem7 = {2: 0.167, 3: 0.101, 5: 0.078, 10: 0.118}
    res0 = {}
    for b in (2, 3, 5, 10):
        primes = [p for p in (3, 5, 7, 11, 13) if b % p]
        r = run(q, b, 22, 0, primes)
        res0[b] = r
        cells = "; ".join(f"p={p}: class {(-b * pow(2, -1, p)) % p} factor {r['mod'][p]['fac0']:.3f}, F p^2 = {r['mod'][p]['F'] * p * p:.2f}, rem {r['mod'][p]['rem']:.3f}" for p in primes)
        print(f"    b = {b:2d}: {r['total']} triples; {cells}")
        check(r["total"] == want0[b] and r["viol"] == 0, f"b = {b}: {want0[b]} triples; p-divisible smaller states give 2 N_b + b = 0 (mod p): 0 violations")
        check(round(r["mod"][7]["rem2"], 3) == want0_rem7[b] and abs(r["mod"][7]["rem"] - r["mod"][7]["rem2"]) < 1e-3,
              f"b = {b}: remainder at 7 is {want0_rem7[b]:.3f} in the convention of the session "
              f"({r['mod'][7]['rem2']:.5f}; {r['mod'][7]['rem']:.5f} when the diagonal triples are not put in the type 1; {r['ndiag']} diagonal triples)")
    check(max(r["mod"][p]["rem2"] for b, r in res0.items() for p in r["mod"]) < 0.1675
          and max(r["mod"][p]["rem2"] for b, r in res0.items() for p in r["mod"] if p != 7) <= 0.0575,
          "sheet 0: remainder at most 0.167 (at 7), at most 0.057 at the other primes")

    print("4. second base-3 digit by value (section 5; numbers of the audit, hypotheses/second_digit_b3.out)")
    want_d = {16: (0.131, 0.127), 20: (0.071, 0.073), 24: (0.041, 0.040), 26: (0.031, 0.030), 28: (0.024, 0.023)}
    devs = {}
    for K in (16, 18, 20, 22, 24) + ((26, 28) if big else ()):
        d = tuple(run(q, 3, K, sh, [5])["digit"] for sh in (0, 1))
        devs[K] = d
        print(f"    2^{K}: max_c |3 share - 1| = {d[0]:.4f} / {d[1]:.4f}")
        if K in want_d:
            check((round(d[0], 3), round(d[1], 3)) == want_d[K], f"2^{K}: {want_d[K][0]:.3f} / {want_d[K][1]:.3f}")
    ks = sorted(devs)
    rate = [(devs[ks[-1]][sh] / devs[ks[0]][sh]) ** (1.0 / (ks[-1] - ks[0])) for sh in (0, 1)]
    print(f"    decay per binary digit from 2^{ks[0]} to 2^{ks[-1]}: {rate[0]:.4f} / {rate[1]:.4f}; sqrt(3)/2 = {math.sqrt(3) / 2:.4f}")
    check(all(abs(v - math.sqrt(3) / 2) < 0.012 for v in rate), "the deviation decays at about sqrt(3)/2 per binary digit")
    # the affine form: alpha_1 = 2^-1 x2 + (2^(s-1) y1) y2 + const (mod b) on sheet 1, 2^-1 x2 + 2^-1 y1 y2 + const on sheet 0
    good = True
    qq = q[:4000].astype(np.int64)
    rng = np.random.default_rng(1)
    for b in (3, 5, 7, 9, 15):
        i2 = pow(2, -1, b)
        for _ in range(4000):
            h1, h2 = (int(v) for v in rng.integers(0, 4000, 2))
            s = int(rng.integers(0, 6))
            p1, p2, y1, y2 = h1 - int(qq[h1]), h2 - int(qq[h2]), 2 * h1 + 1, 2 * h2 + 1
            a1 = p1 + p2 + ((y1 * y2 << s) - (1 if s == 0 else 0)) // 2
            x2 = 2 * p2 + 1
            if s >= 1:
                form = i2 * x2 + (2 ** (s - 1) * y1) * y2 + (p1 - i2)
            else:
                form = i2 * x2 + i2 * y1 * y2 + (p1 - i2 - i2)
            good &= (a1 - form) % b == 0
    check(good, "modulo an odd b, alpha_1 is an affine form in (x2, y2) with the unit coefficient 1/2 at x2 (20,000 random joins, b = 3, 5, 7, 9, 15)")
    print(f"({time.time() - t0:.0f} s)")
    print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
    sys.exit(0 if OK else 1)


if __name__ == "__main__":
    main()
