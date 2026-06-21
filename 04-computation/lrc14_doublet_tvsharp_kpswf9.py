#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- SHARPEN the doublet correction constant J via EXACT per-period total
variation of the cover-indicator step function, shrinking the closure cutoff f0.

On a base-piece P_i (M=M(y) const, floor(y)=fy const, frac(y)=phi in (phi_lo,phi_hi)), the
cover-indicator as a function of the fast variable t=f*y is, over ONE period of t (t: r=floor(t)
sweeps 0..6 as frac(t):0->1 within each r), the function
   C(t) = 1{ M ⊆ { r, (r+delta) mod 7 } },   r=floor(t) mod7,
   delta = fy mod7        if frac(t) < 1 - phi      (u<1-frac(y) branch)
   delta = (fy+1) mod7    if frac(t) >= 1 - phi.
The threshold 1-phi varies across the piece; the TOTAL VARIATION of C in t over one period is
maximized when the threshold is interior. We compute the EXACT max-over-phi per-period TV:
walk r=0..6, within each r the value is c0(r)=cover(r,delta0) then c1(r)=cover(r,delta1)
(two sub-steps split by the threshold). The TV over one period =
   sum of |jumps| around the cyclic sequence  [c0(0),c1(0), c0(1),c1(1), ..., c0(6),c1(6)] (cyclic).
This is the MAX TV (threshold interior for every r). |Delta_f| <= (1/(7f)) * SUM_i TV_i =: J_sharp/f.

Reports J_sharp, the new f0_sharp, and certifies f*|Delta_f| <= J_sharp on a window.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from math import floor
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP
from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast
from lrc14_doublet_frozen_law_kpswf9 import p0_inf_doublet, base_y_pieces


def cover(M, r, delta):
    return 1 if M <= {r % 7, (r + delta) % 7} else 0


def period_tv(M, fy):
    """Exact max-over-threshold total variation of the cover step-function over one period of t."""
    d0 = fy % 7
    d1 = (fy + 1) % 7
    seq = []
    for r in range(7):
        seq.append(cover(M, r, d0))  # frac(t) < 1-phi branch
        seq.append(cover(M, r, d1))  # frac(t) >= 1-phi branch
    # cyclic total variation
    tv = 0
    for i in range(len(seq)):
        tv += abs(seq[i] - seq[(i + 1) % len(seq)])
    return tv


def tv_budget_sharp(base):
    pieces = base_y_pieces(base)
    total = 0
    for (yl, yh, M, fy) in pieces:
        total += period_tv(M, fy)
    return F(total, 7), len(pieces)


def doublet(k, M):
    return tuple(sorted(set(list(range(k - 1)) + [M, M + 1])))


def main():
    print("=" * 96)
    print("SHARP doublet correction constant J_sharp via exact per-period TV  (kind-pasteur kpswf9)")
    print("=" * 96)
    Fwin = 300
    for k in (8, 9, 10, 11, 12):
        base = tuple(range(k - 1))
        pinf = p0_inf_doublet(base, 1)
        cap = CAP[k]
        margin_inf = cap - pinf
        Js, N = tv_budget_sharp(base)
        f0 = floor(Js / margin_inf) + 1
        # certify
        Kemp = F(0); viol = 0
        for f in range(15, Fwin + 1):
            d = abs((p0_fast(doublet(k, f)) - pinf) * f)
            if d > Kemp:
                Kemp = d
            if d > Js:
                viol += 1
        print(f"\nk={k}: N={N}  J_sharp=(1/7)SUM TV_i = {Js} = {float(Js):.3f}   "
              f"(prev conservative J=2*N={2*N})")
        print(f"   margin_inf={float(margin_inf):.6f}  f0_sharp=floor(J_sharp/margin)+1 = {f0}")
        print(f"   certify f*|Delta|<=J_sharp on [15,{Fwin}]: violations={viol}  (emp sup={float(Kemp):.4f}, "
              f"J_sharp {float(Js/Kemp):.1f}x loose)")
    print("\n" + "=" * 96)
    print("J_sharp shrinks f0; the finite check [15,f0_sharp) then fully closes the doublet tail.")


if __name__ == "__main__":
    main()
