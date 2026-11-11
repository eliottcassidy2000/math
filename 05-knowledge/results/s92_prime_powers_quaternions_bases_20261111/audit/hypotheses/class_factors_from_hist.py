#!/usr/bin/env python3
# Auditor "hypotheses": class factors read from the saved joint histograms (even_moduli_hist_E.json,
# [sheet, type w.r.t. 315, N mod 2520]) of the own census.
#  (a) THM-5134 "Numerical support": at 9, X = 2^28, factors of the classes of valuation 0, 1, >= 2 of N - c_sigma;
#  (b) HYP-9665 (A): factors (2^j m) * share at 12, 20, 36 by (N mod 2^j, valuation class mod m) at 2^22..2^28,
#      against 2 L_m, and the spread over the odd classes mod 2^j inside each class mod m.
import numpy as np, math, sys, json

LCM = 2520
csig = {0: -1, 1: 0}


def Lpk(p, k, e):
    # L_(p^k)(e)
    q = p ** k
    e %= q
    if e == 0:
        v = k
    else:
        v = 0
        while e % p == 0:
            e //= p
            v += 1
    if v < k:
        return (1 + 1 / p) * (1 - p ** (-v - 1))
    return 1 + 1 / p - p ** (-k - 1)


def main():
    print("# class_factors_from_hist.py (auditor hypotheses)")
    for E in (22, 24, 26, 28):
        hist = np.array(json.load(open("even_moduli_hist_%d.json" % E))["hist"], dtype=np.int64).sum(axis=1)  # [sheet, N mod 2520]
        print("X = 2^%d" % E)
        for sh in (0, 1):
            T = hist[sh].sum()
            # (a) modulo 9
            cnt9 = np.zeros(9)
            for c in range(LCM):
                cnt9[c % 9] += hist[sh, c]
            g = {0: [], 1: [], 2: []}
            for c in range(9):
                e = (c - csig[sh]) % 9
                v = 2 if e == 0 else (1 if e % 3 == 0 else 0)
                g[v].append(9 * cnt9[c] / T)
            print("  sheet %d, modulo 9: factors by valuation of N - c_s: v=0 %.3f, v=1 %.3f, v>=2 %.3f   (limits 0.889, 1.185, 1.296)" % (
                sh, np.mean(g[0]), np.mean(g[1]), np.mean(g[2])))
        # (b) even moduli
        for (j, p, k) in ((2, 3, 1), (2, 5, 1), (2, 3, 2), (3, 3, 2)):
            m = p ** k
            Mod = (1 << j) * m
            for sh in (0, 1):
                T = hist[sh].sum()
                cnt = np.zeros(Mod)
                for c in range(LCM):
                    cnt[c % Mod] += hist[sh, c]
                rows = []
                for vv in range(k + 1):
                    cls = [c for c in range(1, Mod, 2) if (lambda e: (k if e == 0 else min(k, next(i for i in range(k + 1) if (e // p ** i) % p != 0))))((c - csig[sh]) % m) == vv]
                    f = [Mod * cnt[c] / T for c in cls]
                    lim = 2 * Lpk(p, k, p ** vv if vv < k else 0)
                    rows.append("v=%d: mean %.4f [min %.4f, max %.4f] over %d classes, limit 2L = %.4f" % (vv, np.mean(f), min(f), max(f), len(cls), lim))
                print("  sheet %d, modulo %d = 2^%d * %d: (2^j m) * share on the odd classes by valuation of N - c_s mod %d:  %s" % (sh, Mod, j, m, m, "; ".join(rows)))


if __name__ == "__main__":
    main()
