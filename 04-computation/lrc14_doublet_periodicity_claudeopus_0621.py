#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: DOES THE GENUINE-WIDE MAXIMIZER CLOSE BY THM-563-STYLE PERIODICITY?

DISCOVERY (lrc14_genuine_wide_maximizer): the genuine-wide p0 maximizer is the
stable family  E_w = consec_{k-1} U {w, w+1}  (consec base + a TIGHT far DOUBLET).
It is the [m-2, 2] coherent-block partition (the largest genuine-wide-compatible
one; [m-1,1] = single-far is EXCLUDED from genuine-wide).

THM-563 closes SINGLE-far by proving w*Delta_w is EXACTLY PERIODIC in w
(period 7*lcm(base)), so sup_w Delta_w*w is a finite period-max.

PREDICTION (claude-opus): the DOUBLET deviation closes the SAME way. Adding the
locked pair {w, w+1} contributes phases frac(w*t) and frac((w+1)*t)=frac(w*t + t),
BOTH periodic in integer w with period = denominator(t) | 7*lcm(base). So
    w * Delta_w^{doublet}   should ALSO be exactly periodic in w, period 7*lcm(base).
If TRUE: the genuine-wide MAXIMIZER closes by a finite period-max EXACTLY like
single-far, and the entire wide region unifies under THM-563's mechanism.

This script TESTS the prediction empirically (exact rationals):
  - base = consec_{k-1} = {0,1,...,k-2}
  - plateau Phi = lim_{w->inf} p0(E_w)  (estimated by stabilization at large coprime w)
  - check whether w * (p0(E_w) - Phi) is periodic in w with period P = 7*lcm(base)
    by comparing the multiset/sequence on [w0, w0+P) vs [w0+P, w0+2P).
Also reports the empirical period-max and compares to 15*margin (the THM-563 closure test).
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL


def lcm_list(xs):
    return reduce(lambda a, b: a * b // gcd(a, b), xs, 1)


def doublet_config(k, w):
    base = list(range(k - 1))  # 0..k-2
    return tuple(sorted(set(base + [w, w + 1])))


def test_periodicity(k, w0=15, verbose=True):
    base = list(range(1, k - 1))  # nonzero base speeds 1..k-2
    P = 7 * lcm_list(base)
    # plateau: p0 at a large position; doublet far away -> decorrelated limit.
    # Use several large w to confirm stabilization of the PERIODIC pattern's mean.
    # We test periodicity of the integer sequence a_w = w*(p0(E_w) - Phi).
    # Phi unknown exactly; but PERIODICITY of w*Delta_w is equivalent to:
    #   w*p0(E_w) - (w+P)*p0(E_{w+P}) being an ARITHMETIC/linear-in-Phi shift.
    # Cleaner test: a_w - a_{w+P} = w*p0_w-(w+P)*p0_{w+P} + Phi*P. If w*p0_w-(w)*Phi periodic,
    #   then g_w := w*p0(E_w) satisfies g_{w+P} - g_w = Phi*P + (periodic_w difference).
    # If w*Delta_w is periodic, then h_w := w*p0(E_w) - w*Phi is periodic =>
    #   D_w := p0(E_w) - p0(E_{w+P}) should equal Phi*(1/w - 1/(w+P))-ish... messy with unknown Phi.
    #
    # DIRECT test independent of Phi: w*Delta_w periodic  <=>  the second difference
    #   w*p0_w - 2*(w+P)*p0_{w+P} + (w+2P)*p0_{w+2P}  is periodic-consistent, i.e. equals
    #   the same for all w (the Phi*P linear term cancels in the second difference).
    # Define S_w = w*p0(E_w). If S_w = w*Phi + periodic(w), then
    #   S_w - 2 S_{w+P} + S_{w+2P} = (w -2(w+P)+(w+2P))*Phi + [per - 2per + per] = 0.
    # So the SECOND DIFFERENCE of S over step P must VANISH identically iff w*Delta_w periodic.
    secdiff = {}
    ok = True
    worst = F(0)
    for w in range(w0, w0 + P):
        Sw = w * p0_fast(doublet_config(k, w))
        SwP = (w + P) * p0_fast(doublet_config(k, w + P))
        Sw2P = (w + 2 * P) * p0_fast(doublet_config(k, w + 2 * P))
        sd = Sw - 2 * SwP + Sw2P
        secdiff[w] = sd
        if sd != 0:
            ok = False
            if abs(sd) > worst:
                worst = abs(sd)
    return P, ok, worst, secdiff


def estimate_plateau_and_periodmax(k, w0=15):
    """If periodic, Phi = (S_{w+P}-S_w)/P is constant; period-max = max Delta_w*w over a period."""
    base = list(range(1, k - 1))
    P = 7 * lcm_list(base)
    # Phi from linear slope of S_w = w*p0 over one period step
    S = {}
    for w in range(w0, w0 + 2 * P + 1):
        S[w] = w * p0_fast(doublet_config(k, w))
    # slope estimates (should be constant = Phi if periodic)
    slopes = set((S[w + P] - S[w]) for w in range(w0, w0 + P) if (w + P) in S)
    Phi = None
    if len(slopes) == 1:
        Phi = next(iter(slopes)) / P
    # period-max of Delta_w * w = S_w - w*Phi
    pmax = None
    if Phi is not None:
        vals = [S[w] - w * Phi for w in range(w0, w0 + P)]
        pmax = max(vals)
    return P, Phi, pmax


def main():
    print("=" * 78)
    print("DOUBLET PERIODICITY TEST  (does the genuine-wide maximizer close like THM-563?)")
    print("claude-opus 2026-06-21   E_w = consec_{k-1} U {w, w+1}")
    print("=" * 78)
    for k in range(6, 10):  # keep periods tractable: lcm(1..k-2)
        base = list(range(1, k - 1))
        P = 7 * lcm_list(base)
        print(f"\nk={k}  base=consec_{k-1}={tuple(range(k-1))}  predicted period P=7*lcm{tuple(base)}={P}")
        if P > 30000:
            print(f"   (period {P} large; testing anyway, may be slow)")
        Pv, ok, worst, _ = test_periodicity(k, w0=15)
        verdict = "PERIODIC (second diff vanishes identically)" if ok else f"NOT periodic (worst 2nd-diff={worst})"
        print(f"   PERIODICITY VERDICT: {verdict}")
        if ok:
            _, Phi, pmax = estimate_plateau_and_periodmax(k, w0=15)
            cap = CAP.get(k)
            q = QVAL.get(k)
            print(f"   plateau Phi = {Phi} = {float(Phi):.6f}" if Phi is not None else "   plateau: slope not constant")
            if pmax is not None:
                print(f"   doublet period-max(Delta_w * w) = {pmax} = {float(pmax):.5f}")
                if cap is not None and Phi is not None:
                    margin = cap - Phi
                    closes = pmax < 15 * margin
                    print(f"   cap_{k}={cap}  margin=cap-Phi={margin}={float(margin):.5f}  "
                          f"15*margin={float(15*margin):.4f}  -> Delta_w<margin for w>=15? {closes}")
    print("\n" + "=" * 78)
    print("If PERIODIC at every k: the genuine-wide doublet maximizer closes by a FINITE")
    print("period-max EXACTLY like THM-563 single-far. Leg (C) binding case unifies under THM-563.")


if __name__ == "__main__":
    main()
