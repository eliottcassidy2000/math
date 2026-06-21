#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- FULL corrected-cardinality closure of the binding genuine-wide doublet.

THE BINDING GENUINE-WIDE CONFIG (HYP-2794), CORRECT CARDINALITY:
    E_f = consec_{k-2} U {f, f+1}   (|E| = k runners),  compared to cap_k.
  At f=21:  k=10 -> {0..7,21,22}, p0 = 1301/2940 = 0.442517  (matches the task statement).

FULL ARGUMENT (exact rationals), per k=8..12:

  (1) FROZEN LAW (proved):  p0(E_f) -> p0_inf(consec_{k-2})  as f->inf.
      sector s_e(x)=floor(7 e x) mod 7 = floor(e y) mod 7, y=7x in [0,7).
      EXACT increment identity:  s_{f+1}(y) = s_f(y) + floor(frac(f y) + y)  (mod 7),
      so the doublet's two far sectors are {s, s+delta}, delta = floor(y)+[u>=1-frac(y)] mod7,
      u=frac(f y) -> Unif[0,1) and s=s_f -> Unif(Z/7) jointly (Weyl, single fast phase t=f y).
      p0_inf = (1/7) INT_0^7 [ (1-frac y) A(M(y),floor y mod7) + (frac y) A(M(y),(floor y+1) mod7) ] dy,
      A(M,d) = (1/7)#{s: M ⊆ {s,(s+d) mod7}},  M(y)={1..6}\{floor(e y)mod7 : e in base}.

  (2) CORRECTION RATE (proved):  on each base-piece P_i, y->t=f y is linear, so
      |INT_{P_i} phi(f y) dy - |P_i| <phi>| <= TV_i / f  (Koksma; one incomplete period, range<=1),
      where TV_i = exact per-period total variation of the cover-indicator in t (period_tv).
      Hence  |p0(E_f) - p0_inf| <= J_sharp / f,  J_sharp = (1/7) SUM_i TV_i.  CERTIFIED on a window.

  (3) CLOSURE:  for f >= f0 := floor(J_sharp/(cap-p0_inf))+1,  p0 <= p0_inf + J_sharp/f < cap.
      EXACT finite check on [15, f0]:  p0(E_f) < cap (max reported).  => p0 < cap for ALL f>=15.

ROBUST MARGIN:  cap_k - max_{f>=15} p0(E_f) >= 0.16 at every k (the task's robust target).
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
from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, QVAL
from lrc14_doublet_frozen_law_kpswf9 import p0_inf_doublet, base_y_pieces
from lrc14_doublet_tvsharp_kpswf9 import period_tv


def base_of(k):
    """k-runner genuine-wide doublet base = consec_{k-2}."""
    return tuple(range(k - 2))


def cfg(k, f):
    return tuple(sorted(set(list(base_of(k)) + [f, f + 1])))


def J_sharp_of(base):
    return F(sum(period_tv(M, fy) for (yl, yh, M, fy) in base_y_pieces(base)), 7)


def main():
    print("=" * 100)
    print("FULL CLOSURE: binding genuine-wide doublet  E_f = consec_{k-2} U {f,f+1}  (|E|=k)  < cap_k")
    print("kind-pasteur kpswf9   (corrected cardinality; exact rationals)")
    print("=" * 100)
    Mrobust = 3000  # range to assert the robust margin (peak is at small f; tail -> p0_inf)
    overall = True
    rows = []
    for k in (8, 9, 10, 11, 12):
        base = base_of(k)
        pinf = p0_inf_doublet(base, 1)
        cap = CAP[k]
        margin_inf = cap - pinf
        Js = J_sharp_of(base)
        f0 = floor(Js / margin_inf) + 1
        # exact finite check on [15, f0] AND certify f*|Delta| <= Js
        maxp0 = F(-1); fmax = None; maxfd = F(0); ffd = None; bad = 0
        for f in range(15, f0 + 1):
            pv = p0_fast(cfg(k, f))
            if pv >= cap:
                bad += 1
            fd = abs((pv - pinf) * f)
            if fd > Js:
                bad += 1
            if pv > maxp0:
                maxp0, fmax = pv, f
            if fd > maxfd:
                maxfd, ffd = fd, f
        tail_bound = pinf + Js / f0
        robust_margin = cap - maxp0  # max over [15,f0] is the global max (tail<=p0_inf+Js/f0<that)
        closed = (bad == 0) and (tail_bound < cap)
        overall = overall and closed
        rows.append((k, base, pinf, cap, Js, f0, maxp0, fmax, robust_margin, tail_bound, closed))
        print(f"\nk={k}  base=consec_{k-2}={base}  cap={cap}={float(cap):.6f}")
        print(f"  (1) p0_inf = {pinf} = {float(pinf):.7f}   cap - p0_inf = {float(margin_inf):.6f}")
        print(f"  (2) J_sharp = (1/7)SUM TV_i = {Js} = {float(Js):.3f};  "
              f"certified f*|Delta_f| <= J_sharp on [15,{f0}] (max {float(maxfd):.4f} at f={ffd}); bad={bad}")
        print(f"  (3) f0 = floor(J_sharp/(cap-p0_inf))+1 = {f0}")
        print(f"      finite check [15,{f0}]: max p0 = {maxp0} = {float(maxp0):.6f} at f={fmax}")
        print(f"      tail f>=f0: p0 <= p0_inf+J_sharp/f0 = {float(tail_bound):.6f} < cap  ({tail_bound<cap})")
        print(f"  ROBUST MARGIN cap - max_p0 = {float(robust_margin):.6f}  "
              f"({'>= 0.16 OK' if robust_margin >= F(16,100) else '< 0.16'})")
        print(f"  ==> k={k}: {'CLOSED (< cap for all f>=15)' if closed else 'NOT CLOSED'}")

    print("\n" + "=" * 100)
    print("SUMMARY TABLE  (binding genuine-wide doublet)")
    print(f"  {'k':>3} {'p0_inf':>10} {'cap':>10} {'max_p0':>10} {'@f':>4} {'robust_margin':>13} {'J_sharp':>8} {'f0':>5} {'closed':>7}")
    for (k, base, pinf, cap, Js, f0, maxp0, fmax, rm, tb, closed) in rows:
        print(f"  {k:>3} {float(pinf):>10.6f} {float(cap):>10.6f} {float(maxp0):>10.6f} {fmax:>4} "
              f"{float(rm):>13.6f} {float(Js):>8.3f} {f0:>5} {str(closed):>7}")
    print(f"\n  OVERALL doublet genuine-wide tail CLOSED (all f>=15, k=8..12): {overall}")
    print(f"  ROBUST MARGIN >= 0.16 at every k: {all(rm >= F(16,100) for (*_, rm, tb, c) in [(r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],r[10]) for r in rows])}")
    print("=" * 100)


if __name__ == "__main__":
    main()
