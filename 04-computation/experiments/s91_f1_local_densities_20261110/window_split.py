#!/usr/bin/env python3
"""S91 (THM-5128 / HYP-9660): residue classes of the F1 family on a window of N, split by p-divisible rooted states,
against an a-priori prediction.

F1:  N = x1 x2 + 2^s y1 y2 - [s = 0],  x(h) = 2 psi(h) + 1 = y(h) - 2 q(h),  y(h) = 2h + 1,  ordered triples (h1, h2, s).
Sheets: s = 0 (distinguished class N = -1 mod p) and s >= 1 (distinguished class N = 0 mod p).
A triple is p-divisible when p | gcd(x(h1), y(h1)) or p | gcd(x(h2), y(h2)); then N lies in the distinguished class.

MEASURED on the odd N of [2^LO, 2^K) (exact enumeration of all triples):
  T          number of triples on the sheet,
  E          (mean of r on the distinguished class - mean over the other classes) / that mean,
  f_div      proportion of p-divisible triples,
  E_nd       the same excess after removing the p-divisible triples,
  ripple     max relative deviation of a non-distinguished class from the mean of the non-distinguished classes.
PREDICTED without the window (inputs: the list of p-divisible states, the uniform law of rho = x/y on [1/2, 1] and
the first-order class law (*) of THM-5128; no window data):
  W(h1, s) = (1/x1) int_{f0}^{f1} w(f) df,  w(f) = (B - A)/f for f < A,  B/f for A <= f < B,  0 for f >= B,
             f0 = x1/2 + 2^s y1,  f1 = x1 + 2^s y1   (A = 2^LO, B = 2^K; the number of odd y2 with A <= y2 f < B),
  T_pred = sum_{h1, s} W,   f = (sum over p-divisible h1 of W) / T_pred,
  Delta = 2 f - 1/p^2,   E_pred = p Delta / (1 - Delta),   f_div_pred = 2f - f^2,   E_nd_pred = p (f^2 - 1/p^2) / (1 - Delta).

Usage:
  python3 window_split.py K LO                        enumerate, print measured against predicted, and evaluate the
                                                      tolerance rules below on this window (a calibration report)
  python3 window_split.py K LO write-prereg FILE      write the predictions and the test intervals for the window to
                                                      FILE (json) WITHOUT enumerating it
  python3 window_split.py K LO FILE.json              enumerate and evaluate the intervals stored in FILE.json

Tolerance rules (fixed on the windows [2^18, 2^21) and [2^21, 2^23) before the pre-registered run):
  se(p) = sqrt(0.75 p / #odd N)   (r has relative variance about 0.75 on a sheet)
  slow primes (rate rho_p >= 0.85, where the second-order term is visible): 3, 7, 31, 73, 127
  A  total mass: all triples within 0.2 %, each sheet within 2 % of T_pred
  B  weights: #triples with h1 = h within 0.5 % of W(h) for h in 0, 1, 2, 3, 8, 16, 64, 128, 1024, each sheet
  C  divisible fraction: relative error <= 3 % (first divisible state < 2^12) or 15 % (later), plus 4 / sqrt(count)
  D  excess, fast primes and the Fermat / Wagstaff primes: |E - E_pred| <= 0.008 + 0.08 |E_pred| + 4 se
  E  excess, slow primes: |E - E_pred| <= 0.08 (p = 3), 0.15 (p = 7), 0.05 (p = 31, 73, 127)
  F  controls next to the Fermat / Wagstaff primes: |E - E_pred| <= 0.008 + 4 se
  G  non-divisible deficiency, fast primes <= 131: |E_nd - E_nd_pred| <= 0.012 + 4 se
  H  signs: E >= 0.10 for p in 3, 5, 7, 11, 17, 257, 65537; -0.05 <= E <= 0.02 for the sixteen primes listed in SIGN_LOW
"""
import json
import math
import sys

import numpy as np
from numba import njit, prange

PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
          109, 113, 127, 131, 241, 251, 257, 263, 677, 683, 691, 2729, 2731, 2741, 8191, 43691, 65521, 65537, 65539]
SLOW = {3: 0.08, 7: 0.15, 31: 0.05, 73: 0.05, 127: 0.05}
SPECIAL = [241, 257, 683, 2731, 43691, 65537]            # prime factors of 2^m + 1 (m = 12, 8, 11, 13, 17, 16)
CONTROLS = [251, 263, 677, 691, 2729, 2741, 8191, 65521, 65539]
SIGN_HIGH = [3, 5, 7, 11, 17, 257, 65537]
SIGN_LOW = [29, 31, 37, 61, 67, 71, 73, 79, 83, 97, 101, 103, 109, 113, 127, 131]
WEIGHT_H = [0, 1, 2, 3, 8, 16, 64, 128, 1024]


@njit(cache=False)
def build_q(size):
    q = np.zeros(size, np.int32)
    for h in range(1, size):
        half = h >> 1
        if h & 1:
            q[h] = q[half] + ((half + 1) >> 1)
        else:
            q[h] = q[half]
    return q


@njit(cache=False)
def build_mask(q, primes):
    size = q.shape[0]
    mask = np.zeros(size, np.int64)
    np_ = primes.shape[0]
    first = np.full(np_, -1, np.int64)
    for h in range(size):
        y = 2 * h + 1
        x = y - 2 * q[h]
        m = 0
        for j in range(np_):
            p = primes[j]
            if y % p == 0 and x % p == 0:
                m |= (1 << j)
                if first[j] < 0:
                    first[j] = h
        mask[h] = m
    return mask, first


@njit(cache=False)
def w_closed(x1, Y, A, B):
    f0 = 0.5 * x1 + Y
    f1 = 1.0 * x1 + Y
    w = 0.0
    if f0 < A:
        top = f1 if f1 < A else A
        w += (B - A) * math.log(top / f0)
    lo = f0 if f0 > A else A
    hi = f1 if f1 < B else B
    if lo < hi:
        # A <= f < B: all odd y2 in [1, B/f); their number is ceil((B/f - 1)/2), which is B/(2f) on average
        w += B * math.log(hi / lo)
    return w / x1


@njit(parallel=True, cache=False)
def predicted(q, A, B, K, mask, nprimes, nchunks):
    size = q.shape[0]
    T = np.zeros((nchunks, 2))
    F = np.zeros((nchunks, nprimes, 2))
    step = size // nchunks + 1
    for c in prange(nchunks):
        for h1 in range(c * step, min(size, (c + 1) * step)):
            y1 = 2 * h1 + 1
            x1 = y1 - 2 * q[h1]
            m = mask[h1]
            for s in range(K):
                Y = y1 << s
                if Y >= B:
                    break
                w = w_closed(x1, Y, A, B)
                sheet = 0 if s == 0 else 1
                T[c, sheet] += w
                if m != 0:
                    mm = m
                    j = 0
                    while mm != 0:
                        if mm & 1:
                            F[c, j, sheet] += w
                        mm >>= 1
                        j += 1
    return T.sum(axis=0), F.sum(axis=0)


@njit(cache=False)
def enumerate_window(q, A, B, K, mask, nprimes, wlim):
    size = (B - A) // 2
    r0 = np.zeros(size, np.int32)
    r1 = np.zeros(size, np.int32)
    W = np.zeros((wlim, 2), np.int64)
    cdiv = np.zeros((nprimes, 2), np.int64)
    tot = np.zeros(2, np.int64)
    H = q.shape[0]
    for s in range(K):
        sheet = 0 if s == 0 else 1
        d0 = 1 if s == 0 else 0
        for h1 in range(H):
            y1 = 2 * h1 + 1
            Y = y1 << s
            if Y >= B:
                break
            x1 = y1 - 2 * q[h1]
            y2max = (2 * (B - 1 + d0)) // (2 * Y + x1)
            if y2max < 1:
                continue
            h2max = (y2max - 1) // 2 + 1
            if h2max >= H:
                h2max = H - 1
            h2min = ((A // (Y + x1)) - 1) // 2 - 1
            if h2min < 0:
                h2min = 0
            m1 = mask[h1]
            for h2 in range(h2min, h2max + 1):
                y2 = 2 * h2 + 1
                x2 = y2 - 2 * q[h2]
                N = x1 * x2 + Y * y2 - d0
                if N >= A and N < B:
                    idx = (N - A) >> 1
                    if sheet == 0:
                        r0[idx] += 1
                    else:
                        r1[idx] += 1
                    tot[sheet] += 1
                    if h1 < wlim:
                        W[h1, sheet] += 1
                    m = m1 | mask[h2]
                    if m != 0:
                        j = 0
                        while m != 0:
                            if m & 1:
                                cdiv[j, sheet] += 1
                            m >>= 1
                            j += 1
    return r0, r1, W, cdiv, tot


def predictions(q, A, B, K, mask):
    Tp, Fp = predicted(q, A, B, K, mask, len(PRIMES), 48)
    out = {}
    for j, p in enumerate(PRIMES):
        for sheet in (0, 1):
            f = Fp[j, sheet] / Tp[sheet]
            delta = 2 * f - 1 / p ** 2
            out[(p, sheet)] = dict(f=float(f), E=float(p * delta / (1 - delta)), fdiv=float(2 * f - f * f),
                                   End=float(p * (f * f - 1 / p ** 2) / (1 - delta)))
    return Tp, out


def weight_pred(q, h1, sheet, A, B, K):
    y1 = 2 * h1 + 1
    x1 = y1 - 2 * int(q[h1])
    srange = (0,) if sheet == 0 else range(1, K)
    return sum(w_closed(float(x1), float(y1 << s), float(A), float(B)) for s in srange if (y1 << s) < B)


def build_tests(Tp, pred, first, nodd, q, A, B, K):
    """the test intervals of the rules in the module docstring, as a json-serialisable dict."""
    def se(p):
        return math.sqrt(0.75 * p / nodd)

    fast_small = [p for p in PRIMES if p <= 131 and p not in SLOW]
    tests = {}
    tot_pred = float(Tp.sum())
    tests["A_total_mass"] = {"items": [
        {"key": "total:2", "lo": tot_pred * 0.998 / nodd, "hi": tot_pred * 1.002 / nodd},
        {"key": "total:0", "lo": float(Tp[0]) * 0.98 / nodd, "hi": float(Tp[0]) * 1.02 / nodd},
        {"key": "total:1", "lo": float(Tp[1]) * 0.98 / nodd, "hi": float(Tp[1]) * 1.02 / nodd}]}
    tests["B_weights"] = {"items": [
        {"key": f"weight:{h}:{sheet}", "lo": 0.995 * weight_pred(q, h, sheet, A, B, K), "hi": 1.005 * weight_pred(q, h, sheet, A, B, K)}
        for h in WEIGHT_H for sheet in (0, 1)]}
    items = []
    for j, p in enumerate(PRIMES):
        if p <= 131 or p in (241, 257, 683):
            for sheet in (0, 1):
                fd = pred[(p, sheet)]["fdiv"]
                count = fd * float(Tp[sheet])
                rel = (0.03 if 0 <= first[j] < 4096 else 0.15) + 4 / math.sqrt(max(count, 1.0))
                items.append({"key": f"class:{p}:{sheet}:fdiv", "lo": fd * (1 - rel), "hi": fd * (1 + rel)})
    tests["C_divisible_fraction"] = {"items": items}
    items = []
    for p in fast_small + SPECIAL:
        for sheet in (0, 1):
            e = pred[(p, sheet)]["E"]
            tol = 0.008 + 0.08 * abs(e) + 4 * se(p)
            items.append({"key": f"class:{p}:{sheet}:E", "lo": e - tol, "hi": e + tol})
    tests["D_excess_first_order"] = {"items": items}
    tests["E_excess_slow_primes"] = {"items": [
        {"key": f"class:{p}:{sheet}:E", "lo": pred[(p, sheet)]["E"] - tol, "hi": pred[(p, sheet)]["E"] + tol}
        for p, tol in SLOW.items() for sheet in (0, 1)]}
    tests["F_controls"] = {"items": [
        {"key": f"class:{p}:{sheet}:E", "lo": pred[(p, sheet)]["E"] - 0.008 - 4 * se(p), "hi": pred[(p, sheet)]["E"] + 0.008 + 4 * se(p)}
        for p in CONTROLS for sheet in (0, 1)]}
    tests["G_nondivisible_deficiency"] = {"items": [
        {"key": f"class:{p}:{sheet}:End", "lo": pred[(p, sheet)]["End"] - 0.012 - 4 * se(p),
         "hi": pred[(p, sheet)]["End"] + 0.012 + 4 * se(p)} for p in fast_small for sheet in (0, 1)]}
    tests["H_signs"] = {"items": [{"key": f"class:{p}:{sheet}:E", "lo": 0.10, "hi": 10.0} for p in SIGN_HIGH for sheet in (0, 1)]
                        + [{"key": f"class:{p}:{sheet}:E", "lo": -0.05, "hi": 0.02} for p in SIGN_LOW for sheet in (0, 1)]}
    return tests


def evaluate(tests, results, label):
    allok = True
    for name, test in tests.items():
        bad = []
        for item in test["items"]:
            val = results[item["key"]]
            if not (item["lo"] <= val <= item["hi"]):
                bad.append(f"{item['key']} = {val:.6g} not in [{item['lo']:.6g}, {item['hi']:.6g}]")
        allok &= not bad
        print(f"{label} {name}: {'CONFIRMED' if not bad else 'REFUTED'}  ({len(test['items'])} items"
              + (")" if not bad else f", {len(bad)} outside)  " + "; ".join(bad[:10])))
    return allok


def main():
    K = int(sys.argv[1])
    LO = int(sys.argv[2])
    mode = sys.argv[3] if len(sys.argv) > 3 else ""
    A, B = 1 << LO, 1 << K
    q = build_q(B // 2 + 2)
    primes = np.array(PRIMES, dtype=np.int64)
    mask, first = build_mask(q, primes)
    Tp, pred = predictions(q, A, B, K, mask)
    nodd = (B - A) // 2
    print(f"F1 on odd N in [2^{LO}, 2^{K}):  predicted triples per odd N: sheet s = 0: {Tp[0] / nodd:.4f}, "
          f"sheet s >= 1: {Tp[1] / nodd:.4f}, total {Tp.sum() / nodd:.4f}")
    tests = build_tests(Tp, pred, first, nodd, q, A, B, K)
    if mode == "write-prereg":
        spec = {"window": f"odd N in [2^{LO}, 2^{K})", "script": "window_split.py", "rules": "see the module docstring",
                "predictions": {f"{p}:{sheet}": pred[(p, sheet)] for p in PRIMES for sheet in (0, 1)},
                "predicted_triples_per_odd_N": {"sheet0": float(Tp[0] / nodd), "sheet1": float(Tp[1] / nodd)},
                "first_divisible_state": {str(p): int(first[j]) for j, p in enumerate(PRIMES)},
                "tests": tests}
        json.dump(spec, open(sys.argv[4], "w"), indent=1)
        print("pre-registration written to", sys.argv[4], "(no enumeration was run)")
        for p in PRIMES:
            a, b = pred[(p, 0)], pred[(p, 1)]
            print(f"  p = {p:6d} first divisible h = {int(first[PRIMES.index(p)]):9d}   s = 0: E {a['E']:+.5f}  f_div {a['fdiv']:.6f}"
                  f"  E_nd {a['End']:+.5f}    s >= 1: E {b['E']:+.5f}  f_div {b['fdiv']:.6f}  E_nd {b['End']:+.5f}")
        return
    r0, r1, W, cdiv, tot = enumerate_window(q, A, B, K, mask, len(PRIMES), 1 << 12)
    results = {"total:0": tot[0] / nodd, "total:1": tot[1] / nodd, "total:2": tot.sum() / nodd}
    print(f"measured triples per odd N: sheet s = 0: {tot[0] / nodd:.4f}, sheet s >= 1: {tot[1] / nodd:.4f}, "
          f"total {tot.sum() / nodd:.4f};  misses (r = 0): {int(((r0 + r1) == 0).sum())}")
    ratios = {}
    for h in WEIGHT_H:
        for sheet in (0, 1):
            results[f"weight:{h}:{sheet}"] = float(W[h, sheet])
            ratios[(h, sheet)] = round(float(W[h, sheet]) / weight_pred(q, h, sheet, A, B, K), 5)
    print("weights of small h1, measured count / W(h1):", ratios)
    Nodd = np.arange(A + 1, B, 2, dtype=np.int64)
    print(f"{'p':>6} {'sheet':>5} {'E meas':>9} {'E pred':>9} {'f_div meas':>11} {'f_div pred':>11} {'E_nd meas':>10} "
          f"{'E_nd pred':>10} {'ripple':>8} {'se(E)':>8}")
    for j, p in enumerate(PRIMES):
        res = Nodd % p
        for sheet, r in ((0, r0), (1, r1)):
            cs = (p - 1) if sheet == 0 else 0
            counts = np.bincount(res, weights=r, minlength=p)
            nclass = np.bincount(res, minlength=p)
            means = counts / nclass
            others = np.delete(means, cs)
            mo = others.mean()
            E = (means[cs] - mo) / mo
            fdiv = cdiv[j, sheet] / tot[sheet]
            End = (means[cs] - cdiv[j, sheet] / nclass[cs] - mo) / mo
            ripple = float(np.max(np.abs(others - mo)) / mo)
            var = float((r.astype(np.float64) ** 2).mean() - r.mean() ** 2)
            se = math.sqrt(var / nclass[cs]) / mo
            pr = pred[(p, sheet)]
            results[f"class:{p}:{sheet}:E"] = float(E)
            results[f"class:{p}:{sheet}:fdiv"] = float(fdiv)
            results[f"class:{p}:{sheet}:End"] = float(End)
            print(f"{p:6d} {('s=0' if sheet == 0 else 's>=1'):>5} {E:+9.5f} {pr['E']:+9.5f} {fdiv:11.6f} {pr['fdiv']:11.6f} "
                  f"{End:+10.5f} {pr['End']:+10.5f} {ripple:8.5f} {se:8.5f}")
    if mode.endswith(".json"):
        spec = json.load(open(mode))
        ok = evaluate(spec["tests"], results, "PRE-REGISTERED")
        print("ALL PRE-REGISTERED TESTS CONFIRMED" if ok else "SOME PRE-REGISTERED TESTS REFUTED")
    else:
        ok = evaluate(tests, results, "CALIBRATION (in-sample rules)")
        print("calibration: all rules hold on this window" if ok else "calibration: some rules fail on this window")


if __name__ == "__main__":
    main()
