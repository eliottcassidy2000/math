#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD C (kps-wf9): the FINITE-f CORRECTION  c(f) = p0(E_f) - p0_inf  for the binding
adjacent doublet  E_f = consec_{k-1} U {f, f+1}.

We have (frozen_phase script, VERIFIED) the EXACT plateau
    p0_inf(k) = lim_{f->inf} p0(E_f)   [g=1 adjacent = worst].
The doublet sup over f>=15 is governed by  p0(E_f) = p0_inf + c(f).  We need:
  (A) the SIGN of c(f): is p0(E_f) <= p0_inf for all f>=15 (correction <= 0)?  If so the
      plateau p0_inf is an UPPER bound for the doublet and  cap - p0_inf  is the margin.
  (B) the SIZE / DECAY of c(f): tabulate c(f), test f*|c(f)| and f^2*|c(f)| bounded,
      report sup over the window and the trend (is it the THM-563 1/f decay or faster?).
  (C) the actual sup_f p0(E_f) over [15, Fmax] EXACT, and  cap - sup  (the real margin),
      vs the >=0.16 target.
  (D) compare to plateau-based bound: max(p0_inf, sup over a finite window) -- which one
      is the binding upper bound, and does either give margin >= 0.16?

ALSO: the boundary/peak region. Although f->inf gives the plateau, the EMPIRICAL sup of
p0(E_f) was at f=20/21 (tight). We separately tabulate the small-f peak region 15..40 to
locate the global sup over all f and confirm it is < cap with the target margin.

Exact rationals throughout.
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

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL, MARGIN
from lrc14_doublet_frozen_phase_kpswf9 import frozen_p0


def doublet(k, f, g=1):
    return tuple(sorted(set(list(range(k - 1)) + [f, f + g])))


@lru_cache(maxsize=None)
def p0c(E):
    return p0_fast(E)


def analyze(k, Fmax=900, peak_hi=45):
    pinf = frozen_p0(k, 1)
    cap = CAP[k]
    # full scan
    vals = {}
    for f in range(15, Fmax + 1):
        vals[f] = p0c(doublet(k, f, 1))
    # (A) sign of correction
    c = {f: vals[f] - pinf for f in vals}
    n_pos = sum(1 for f in c if c[f] > 0)
    n_zero = sum(1 for f in c if c[f] == 0)
    n_neg = sum(1 for f in c if c[f] < 0)
    max_pos = max((c[f], f) for f in c)  # most positive correction (exceeds plateau)
    # (B) decay
    supfc = max((abs(c[f]) * f, f) for f in c)
    supf2c = max((abs(c[f]) * f * f, f) for f in c)
    # tail sup (f large): sup |c(f)| for f >= thr
    tail = {}
    for thr in (15, 50, 100, 200, 400, 700):
        sub = [abs(c[f]) for f in c if f >= thr]
        tail[thr] = max(sub) if sub else None
    # (C) actual sup over window
    fsup = max(vals, key=lambda f: vals[f])
    psup = vals[fsup]
    # (D) peak region detail
    peak = {f: vals[f] for f in range(15, peak_hi + 1)}
    return dict(k=k, pinf=pinf, cap=cap, n_pos=n_pos, n_zero=n_zero, n_neg=n_neg,
                max_pos=max_pos, supfc=supfc, supf2c=supf2c, tail=tail,
                fsup=fsup, psup=psup, peak=peak, vals=vals, c=c)


def main():
    print("=" * 92)
    print("THREAD C: FINITE-f CORRECTION c(f)=p0(E_f)-p0_inf for binding doublet (kps-wf9)")
    print("E_f = consec_{k-1} U {f,f+1}.  Target: sup_f p0(E_f) < cap_k with margin >= 0.16")
    print("=" * 92)
    Fmax = 900
    summary = []
    for k in range(8, 13):
        r = analyze(k, Fmax)
        pinf = r['pinf']; cap = r['cap']
        print(f"\nk={k}  cap={cap}={float(cap):.6f}  p0_inf={pinf}={float(pinf):.6f}  "
              f"cap-p0_inf={float(cap-pinf):+.6f}")
        print(f"  (A) SIGN of c(f)=p0(E_f)-p0_inf over [15,{Fmax}]: "
              f"pos={r['n_pos']} zero={r['n_zero']} neg={r['n_neg']}")
        mp, mpf = r['max_pos']
        print(f"      most POSITIVE correction = {float(mp):+.6f} at f={mpf}  "
              f"(if <=0 everywhere, p0_inf is an UPPER bound)")
        print(f"  (B) DECAY: sup_f f*|c(f)| = {float(r['supfc'][0]):.4f} (f={r['supfc'][1]}); "
              f"sup_f f^2*|c(f)| = {float(r['supf2c'][0]):.2f} (f={r['supf2c'][1]})")
        print(f"      tail sup|c(f)| for f>=: " +
              "  ".join(f"{thr}:{float(v):.5f}" for thr, v in r['tail'].items() if v is not None))
        print(f"  (C) ACTUAL sup_f p0(E_f) over [15,{Fmax}] = {float(r['psup']):.6f} at f={r['fsup']}")
        margin_actual = cap - r['psup']
        margin_plateau = cap - max(pinf, r['psup'])
        print(f"      cap - sup_f = {float(margin_actual):+.6f}   (>= 0.16? {margin_actual >= F(16,100)})")
        print(f"      UPPER BOUND = max(p0_inf, window-sup) = {float(max(pinf, r['psup'])):.6f}; "
              f"cap - that = {float(margin_plateau):+.6f}  (>=0.16? {margin_plateau >= F(16,100)})")
        # peak region
        peakvals = sorted(r['peak'].items(), key=lambda t: t[1], reverse=True)[:5]
        print(f"  (D) top-5 peak-region p0(E_f), f in 15..45: " +
              ", ".join(f"f={f}:{float(v):.5f}" for f, v in peakvals))
        summary.append((k, cap, pinf, r['psup'], margin_actual, margin_plateau))
    print("\n" + "=" * 92)
    print("SUMMARY (the doublet upper bound and margin)")
    print("=" * 92)
    print(f"  {'k':>2} {'cap':>10} {'p0_inf':>10} {'win-sup':>10} {'cap-sup':>10} {'cap-UB':>10} {'>=0.16?':>8}")
    for k, cap, pinf, psup, ma, mp in summary:
        ub_margin = min(ma, mp)
        print(f"  {k:>2} {float(cap):>10.6f} {float(pinf):>10.6f} {float(psup):>10.6f} "
              f"{float(ma):>10.6f} {float(mp):>10.6f} {str(mp>=F(16,100)):>8}")
    print("\nREADING: if c(f) <= 0 (p0_inf is an upper bound) OR window-sup+plateau both clear")
    print("cap by >= 0.16, the binding doublet closes. The tight k is the binding margin.")


if __name__ == "__main__":
    main()
