#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: THE BINDING GENUINE-WIDE DOUBLET — SIGNED BOUND + FINITE-CHECK CUTOFF.

THE ONE REMAINING GAP (mac-mini-S21): the genuine-wide finite-M error bound.
Actual error ~0.01, GAP (margin) 0.12-0.26, but the THM-557 BV bound 7*C(m,2)/M is
~30 at M=15 (USELESS) — the SAME signed-cancellation overcounting THM-563 fixed for
single-far. HYP-2794: the binding genuine-wide config is the doublet family
    E_M = consec_{k-1} U {M, M+1}.
The doublet's w*Delta_w is NOT exactly periodic (HYP-2794) but is ALMOST periodic
(small decaying 2nd-difference). KEY: if  M*error(M)  is BOUNDED (periodic + decaying
tail), then error(M) <= C/M <= C/15, closing if C < 15*margin.

This script (exact rationals) for k=8..12:
  (1) sup_{M in [15, Mmax]} p0(E_M) and where — confirm < cap_k with the actual margin.
  (2) estimate plateau Phi_2 = mean of p0(E_M) over a late window; M*error(M) = M*(p0-Phi_2);
      report sup_M |M*error(M)| = the SIGNED 'doublet period-max' analogue (vs BV 7*C(m,2)).
  (3) THM-557-style rigorous finite cutoff M_BV = ceil(7*C(m,2)/margin_2): beyond it,
      error <= 7*C(m,2)/M < margin by the BV bound -> p0 < cap PROVED; [15, M_BV] is a
      finite exact check. Report M_BV and the exact max over [15, min(M_BV,Mmax)].
  (4) the OVERCOUNT ratio  (7*C(m,2)) / (sup_M M*|error|)  — how lossy the BV bound is
      (the analogue of THM-563's 125x for single-far). A small signed bound => tiny cutoff.

If sup_M M*error is small (~0.2) and < 15*margin (~2.4): the binding genuine-wide doublet
closes by a THM-563-style SIGNED bounded-deviation argument, collapsing the cutoff from
~1500 to ~tens. That is the concrete close of the ONE remaining gap (for this family).
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from math import comb, ceil
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL


def doublet(k, M):
    return tuple(sorted(set(list(range(k - 1)) + [M, M + 1])))


def analyze(k, Mmax=600, late_lo=400):
    m = k - 1
    cap = CAP[k]
    # (1) sup over a wide window
    vals = {}
    for M in range(15, Mmax + 1):
        vals[M] = p0_fast(doublet(k, M))
    Msup = max(vals, key=lambda M: vals[M])
    psup = vals[Msup]
    # (2) plateau estimate: average over a late window (doublet far away ~ decorrelated)
    late = [vals[M] for M in range(late_lo, Mmax + 1)]
    Phi2 = sum(late) / len(late)
    # signed M*error
    signed = {M: M * (vals[M] - Phi2) for M in vals}
    Msig = max(signed, key=lambda M: abs(signed[M]))
    sup_Merr = abs(signed[Msig])
    margin2 = cap - Phi2
    # (3) BV cutoff (THM-557 style, single-block constant 7*C(m,2) as a safe over-estimate)
    BV = 7 * comb(m, 2)
    M_BV = ceil(BV / margin2) if margin2 > 0 else None
    # exact max over [15, min(M_BV, Mmax)] (the finite check window we actually cover)
    hi = min(M_BV, Mmax) if M_BV else Mmax
    maxwin = max(vals[M] for M in range(15, hi + 1))
    # (4) overcount ratio
    overcount = float(BV) / float(sup_Merr) if sup_Merr > 0 else None
    return dict(k=k, m=m, cap=cap, Msup=Msup, psup=psup, Phi2=Phi2, margin2=margin2,
                sup_Merr=sup_Merr, Msig=Msig, BV=BV, M_BV=M_BV, hi=hi, maxwin=maxwin,
                overcount=overcount, vals_at_21=vals.get(21))


def main():
    print("=" * 80)
    print("BINDING GENUINE-WIDE DOUBLET: SIGNED BOUND + FINITE CUTOFF  (the ONE remaining gap)")
    print("claude-opus 2026-06-21   E_M = consec_{k-1} U {M, M+1}")
    print("=" * 80)
    for k in range(8, 13):
        r = analyze(k)
        print(f"\nk={k}  m={r['m']}  cap={r['cap']}={float(r['cap']):.5f}")
        print(f"  (1) sup_M p0 = {float(r['psup']):.6f} at M={r['Msup']}  "
              f"(p0 at M=21: {float(r['vals_at_21']):.6f});  cap-sup = {float(r['cap']-r['psup']):+.6f}")
        print(f"  (2) plateau Phi_2 ~ {float(r['Phi2']):.6f}  margin_2=cap-Phi_2 ~ {float(r['margin2']):.6f}")
        print(f"      SIGNED doublet bound: sup_M |M*error(M)| = {float(r['sup_Merr']):.4f} at M={r['Msig']}")
        print(f"      closes if < 15*margin_2 = {float(15*r['margin2']):.4f}  -> "
              f"{'YES' if r['sup_Merr'] < 15*r['margin2'] else 'NO'}")
        print(f"  (3) BV cutoff M_BV = ceil(7*C({r['m']},2)/margin_2) = {r['M_BV']}  "
              f"(BV const 7*C(m,2)={r['BV']})")
        print(f"      exact max p0 over [15,{r['hi']}] = {float(r['maxwin']):.6f} < cap? "
              f"{'YES' if r['maxwin'] < r['cap'] else 'NO'}")
        print(f"  (4) OVERCOUNT ratio BV / sup_M(M*err) = {r['overcount']:.1f}x  "
              f"(THM-563 single-far analogue was ~125x)")
    print("\n" + "=" * 80)
    print("READING: if (2) 'YES' everywhere, the binding doublet closes by a SIGNED bounded-")
    print("deviation bound (THM-563 analogue) with cutoff collapsed from ~1500 to ~tens.")
    print("Combined with 'maximizer = consec+doublet' (HYP-2794), this closes leg-C binding case.")


if __name__ == "__main__":
    main()
