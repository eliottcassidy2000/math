#!/usr/bin/env python3
"""The Fourier coefficients of the 3-adic Syracuse law at the powers of two, to level 120, without the law.

The exact frequency recursion  mu_hat_n(t) = sum_(a>=1) 2^-a e((t 2^-a mod 3^n)/3^n) mu_hat_(n-1)(t 2^-a mod 3^(n-1))
closes on the family t = 2^j (j in Z, the inverse powers included): t = 2^j gives t 2^-a = 2^(j-a).  With the
valuation sum truncated at a <= AMAX (error <= 2^-AMAX per level) and mu_hat_0 = 1, the values m_n(j) = mu_hat_n(2^j mod 3^n)
follow for every j and n by a recursion over ~5000 exponents per level.  For h <= 18 the maximal primitive
coefficient M(h) sits at +-2^s (S20 note, FFT of the full law), so max_j |m_h(j)| reproduces M(h) there and gives a
LOWER bound for M(h) at every level (an upper bound only if the argmax stays on the powers of two, OBSERVED to 18).
Output: max_j |m_h(j)|, its argmax, ratios, and fits (geometric vs shifted power law) over h <= HMAX.
Run: python 04-computation/experiments/collatz_five_mirrors_powers_of_two_20260929.py [HMAX]   (about a minute)
"""
from __future__ import annotations

import cmath
import math
import sys
from fractions import Fraction

if __name__ == "__main__":
    HMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    AMAX = 40
    JMAX = 2 * HMAX + 20
    JMIN = -AMAX * HMAX - AMAX
    weights = [2.0 ** (-a) for a in range(1, AMAX + 1)]
    # level 0: mu_hat_0 = 1 at every frequency
    prev = {j: 1.0 + 0j for j in range(JMIN, JMAX + 1)}
    print("h, max_j |mu_hat_h(2^j)|, argmax j (s - h), ratio, FFT value M(h) where known")
    known = {1: 0.577350, 2: 0.377924, 3: 0.252237, 4: 0.176999, 5: 0.129274, 6: 0.096106, 7: 0.075870, 8: 0.060891,
             9: 0.048026, 10: 0.038278, 11: 0.031944, 12: 0.026458, 13: 0.022052, 14: 0.019128, 15: 0.016284,
             16: 0.014409, 17: 0.012511, 18: 0.011187}
    hist = []
    prevmax = None
    for n in range(1, HMAX + 1):
        mod = 3 ** n
        lo = -AMAX * (HMAX - n)  # level n feeds levels above it: exponent j at level HMAX needs j - a_1 - ... at level n
        pw = {j: pow(2, j, mod) for j in range(lo - AMAX, JMAX + 1)}
        cur = {}
        for j in range(lo, JMAX + 1):
            s = 0j
            for a in range(1, AMAX + 1):
                r = pw[j - a]  # 2^(j-a) mod 3^n
                ph = cmath.exp(2j * math.pi * float(Fraction(r, mod)))
                s += weights[a - 1] * ph * prev.get(j - a, 0j)
            cur[j] = s
        # the coefficient at 2^j mod 3^n depends on j only through 2^j mod 3^n; exponents differing by the order
        # L = 2 3^(n-1) coincide, which the recursion respects automatically.
        best = max(cur, key=lambda j: abs(cur[j]))
        mx = abs(cur[best])
        ratio = f"{mx / prevmax:.4f}" if prevmax else "-"
        kn = f"{known[n]:.6f}" if n in known else ""
        print(f"   h={n:3d}: {mx:.6e}  argmax j={best} (s-h={best - n:+d})  ratio {ratio}  {kn}")
        hist.append((n, mx))
        prevmax = mx
        prev = cur
    # fits over the upper half of the range
    hs = [h for h, m in hist if h >= HMAX // 3 and m > 0]
    ms = [m for h, m in hist if h >= HMAX // 3 and m > 0]
    # geometric: log m = a + b h ; shifted power law: log m = a - alpha log(h + c) with c in a small grid
    import statistics
    xbar, ybar = statistics.mean(hs), statistics.mean([math.log(m) for m in ms])
    b = sum((h - xbar) * (math.log(m) - ybar) for h, m in zip(hs, ms)) / sum((h - xbar) ** 2 for h in hs)
    a = ybar - b * xbar
    res_geo = max(abs(math.log(m) - (a + b * h)) for h, m in zip(hs, ms))
    best_pl = None
    for c10 in range(0, 200):
        c = c10 / 10
        xs = [math.log(h + c) for h in hs]
        xb = statistics.mean(xs)
        bb = sum((x - xb) * (math.log(m) - ybar) for x, m in zip(xs, ms)) / sum((x - xb) ** 2 for x in xs)
        aa = ybar - bb * xb
        res = max(abs(math.log(m) - (aa + bb * x)) for x, m in zip(xs, ms))
        if best_pl is None or res < best_pl[0]:
            best_pl = (res, c, -bb, aa)
    print(f"   fits over h = {hs[0]}..{hs[-1]}: geometric ratio {math.exp(b):.5f} (max log-residual {res_geo:.4f}); "
          f"shifted power law C (h + {best_pl[1]:.1f})^(-{best_pl[2]:.3f}) (max log-residual {best_pl[0]:.4f})")
    # the no-descent probability: P_h = P(a_1 + ... + a_j < j log_2 3 for all j <= h) for i.i.d. geometric(1/2) valuations,
    # by dynamic programming over the prefix sum; its rate is e^(-I(log_2 3)) = 3^(h* - 1), h* = h(log_3 2) = 0.94996 (THM-4476)
    L23 = math.log2(3)
    hstar = -(math.log2(3) ** -1) * math.log2(math.log2(3) ** -1) - (1 - 1 / L23) * math.log2(1 - 1 / L23)
    print(f"   h* = h(log_3 2) = {hstar:.5f}; 3^(h* - 1) = {3 ** (hstar - 1):.5f}; e^(-I(log_2 3)) with I(a) = a H(1/a) - a ln 2: {math.exp(L23 * (hstar * math.log(2)) - L23 * math.log(2)):.5f}")
    dist = {0: 1.0}  # prefix sum -> probability, conditioned on no descent so far (unnormalised)
    Pnd = []
    for j in range(1, HMAX + 1):
        nd = {}
        for P, pr in dist.items():
            for a in range(1, 60):
                Q = P + a
                if Q < j * L23:
                    nd[Q] = nd.get(Q, 0.0) + pr * 2.0 ** (-a)
        dist = nd
        Pnd.append(sum(dist.values()))
    print("   no-descent probability P_h and the ratio max_j |mu_hat_h(2^j)| / P_h:")
    for h, m in hist:
        if h in (5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120):
            print(f"      h={h:3d}: P_h = {Pnd[h-1]:.4e} (ratio to previous level {Pnd[h-1]/Pnd[h-2]:.4f}); M/P_h = {m / Pnd[h-1]:.4f}; (M/P_h)^(1/h) = {(m / Pnd[h-1]) ** (1 / h):.4f}")
    print("   local exponents from doublings h -> 2h: " + ", ".join(f"{h}->{2*h}: {math.log(dict(hist)[h] / dict(hist)[2*h], 2):.3f}" for h in (10, 15, 20, 25, 30, 40, 50, 60) if 2 * h <= HMAX))
    print("DONE")
