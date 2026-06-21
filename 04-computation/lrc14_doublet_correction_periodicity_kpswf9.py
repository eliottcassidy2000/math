#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD C (kps-wf9): IS THE CORRECTION c(f)=p0(E_f)-p0_inf EXACTLY-PERIODIC-TIMES-1/f?

THM-563 closes single-far because f*(p0(S_f)-Phi_1) is EXACTLY periodic in f.
For the doublet, opus showed f*(p0(E_f)-Phi) is NOT exactly periodic (2nd-diff decays).
BUT we now have the EXACT plateau p0_inf (frozen-phase). Two structural tests:

  TEST 1 (the THM-563 analogue for the doublet correction):
     is  g(f) := f * c(f) = f*(p0(E_f) - p0_inf)  EXACTLY periodic in f (period P=7*lcm(base))?
     <=> second difference of S(f):=f*p0(E_f) over step P vanishes identically
         AND the slope (S(f+P)-S(f))/P equals p0_inf exactly.
     (opus tested 2nd-diff of f*p0; here we ALSO pin the slope to the EXACT p0_inf.)

  TEST 2 (almost-periodic = periodic + decaying): if NOT exactly periodic, measure how the
     2nd-difference decays in f: sup over [f0,f0+P) of |S(f)-2S(f+P)+S(f+2P)| as f0 grows.
     If it -> 0 like 1/f, the correction is periodic-part + a 1/f-summable tail (BV in 1/f),
     so sup_f f*|c(f)| is finite (rigorous if the decay is provable).

  TEST 3 (the BvK / two-scale split): decompose c(f) into
     - a PERIODIC part  c_per(f)  (period P) = the f->inf-removed oscillation, and
     - a RESIDUAL  c_res(f) = c(f) - c_per(f)  that -> 0.
     We extract c_per by averaging c(f) over f in a residue class mod P at LARGE f (where
     residual ~0), then measure sup|c_res(f)| over [15,...] and its decay. This pins the
     EXACT periodic envelope whose max governs the f>=15 sup of p0(E_f).

Exact rationals. The point: pin the EXACT structure that bounds sup_{f>=15} p0(E_f).
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce, lru_cache
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL
from lrc14_doublet_frozen_phase_kpswf9 import frozen_p0


def lcm_list(xs):
    return reduce(lambda a, b: a * b // gcd(a, b), xs, 1)


def base_nz(k):
    return list(range(1, k - 1))


def doublet(k, f, g=1):
    return tuple(sorted(set(list(range(k - 1)) + [f, f + g])))


@lru_cache(maxsize=None)
def p0c(E):
    return p0_fast(E)


def test_period(k, f0=15):
    base = base_nz(k)
    P = 7 * lcm_list(base)
    pinf = frozen_p0(k, 1)
    # TEST 1: second difference of S(f)=f*p0(E_f) over step P; and slope vs p0_inf
    worst_sd = F(0)
    slope_set = set()
    for f in range(f0, f0 + P):
        Sf = f * p0c(doublet(k, f))
        SfP = (f + P) * p0c(doublet(k, f + P))
        Sf2P = (f + 2 * P) * p0c(doublet(k, f + 2 * P))
        sd = Sf - 2 * SfP + Sf2P
        if abs(sd) > worst_sd:
            worst_sd = abs(sd)
        slope_set.add(SfP - Sf)  # = P*p0_inf + periodic-diff if structure holds
    exact_periodic = (worst_sd == 0)
    slope_matches = (len(slope_set) == 1 and next(iter(slope_set)) == P * pinf)
    return P, pinf, exact_periodic, worst_sd, slope_matches, len(slope_set)


def test_decay(k, f0=15, n_blocks=6):
    """2nd-difference sup over [f0+jP, ...) as j grows -> does it decay like 1/f?"""
    base = base_nz(k)
    P = 7 * lcm_list(base)
    out = []
    for j in range(n_blocks):
        base_f = f0 + j * 50  # shift the window forward (not by P, just to grow f)
        worst = F(0)
        for f in range(base_f, base_f + min(P, 60)):
            Sf = f * p0c(doublet(k, f))
            SfP = (f + P) * p0c(doublet(k, f + P))
            Sf2P = (f + 2 * P) * p0c(doublet(k, f + 2 * P))
            sd = Sf - 2 * SfP + Sf2P
            if abs(sd) > worst:
                worst = abs(sd)
        out.append((base_f, worst, base_f * worst))
    return P, out


def extract_periodic_envelope(k, f_large=2000, span_periods=1):
    """At large f, c(f)=p0(E_f)-p0_inf is ~ periodic (residual decayed). Extract the periodic
    envelope c_per[r] for residue r mod P by sampling at large f, and report its MAX (the
    asymptotic oscillation amplitude). Then compare sup_{f>=15} p0(E_f) to p0_inf + max c_per."""
    base = base_nz(k)
    P = 7 * lcm_list(base)
    if P > 3000:
        return P, None, None  # keep tractable
    pinf = frozen_p0(k, 1)
    # sample one period at large f
    cper = {}
    for r in range(P):
        f = f_large + r
        cper[r] = p0c(doublet(k, f)) - pinf
    max_cper = max(cper.values())
    min_cper = min(cper.values())
    return P, max_cper, min_cper


def main():
    print("=" * 92)
    print("THREAD C: structure of the doublet correction c(f)=p0(E_f)-p0_inf  (kps-wf9)")
    print("=" * 92)
    for k in range(8, 11):  # k=8,9,10 (P=420,2940,2520 tractable; k>=11 P large)
        print(f"\nk={k}")
        P, pinf, exper, wsd, slmatch, nsl = test_period(k)
        print(f"  P=7*lcm={P}  p0_inf={pinf}={float(pinf):.6f}")
        print(f"  TEST1 f*c(f) exactly periodic? {exper}  (worst 2nd-diff of f*p0 = {wsd}={float(wsd):.6f})")
        print(f"        slope == P*p0_inf exactly? {slmatch}  (#distinct slopes={nsl})")
        Pd, decay = test_decay(k)
        print(f"  TEST2 2nd-diff DECAY (window base_f, worst|2nd-diff|, base_f*worst):")
        for bf, w, prod in decay:
            print(f"        f~{bf}: |2nd-diff|={float(w):.6f}  f*that={float(prod):.4f}")
        Pe, mx, mn = extract_periodic_envelope(k)
        if mx is not None:
            cap = CAP[k]
            print(f"  TEST3 asymptotic periodic envelope of c(f): max={float(mx):+.6f} min={float(mn):+.6f}")
            print(f"        => asymptotic sup p0(E_f) ~ p0_inf+max_cper = {float(pinf+mx):.6f}; "
                  f"cap={float(cap):.6f}; margin={float(cap-(pinf+mx)):+.6f}")
    print("\n" + "=" * 92)
    print("If TEST1 false but TEST2 shows clean 1/f decay: c(f) = periodic_envelope + O(1/f),")
    print("the f>=15 sup is the small-f peak (f=21), handled by exact finite check; the tail")
    print("f>=F0 is bounded by p0_inf + envelope-max + decay/F0 < cap.")


if __name__ == "__main__":
    main()
