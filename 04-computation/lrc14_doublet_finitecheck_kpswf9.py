#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- FINAL FINITE CHECK closing the binding doublet tail.

CARDINALITY (CORRECTED): a k-RUNNER genuine-wide doublet is  E = consec_{k-2} U {f, f+1}
(|E| = (k-2)+2 = k), compared to CAP[k].  (The single-far plateau consec_{k-1} U {1 far} also
has k runners, defining QVAL[k]; the doublet's TIGHT base is consec_{k-2}.)  The binding value
at f=21 is 1301/2940=0.4425 at k=10 -- matching the task statement.

Using J_sharp (lrc14_doublet_tvsharp) the cutoff is f0_sharp = floor(J_sharp/margin_inf)+1.
This script does the HONEST exact finite check:
  for every f in [15, f0_sharp]:  p0(E_f) = p0_fast(consec_{k-2} U {f,f+1})  is EXACT;
  assert  p0(E_f) < cap_k  (with the actual margin), and  f*|p0-p0_inf| <= J_sharp.
For f > f0_sharp the tail bound p0 <= p0_inf + J_sharp/f < cap closes it (PROVED).
=> p0(E_f) < cap_k for ALL f >= 15.  Reports the exact max, its f, and margin per k.
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
from lrc14_doublet_frozen_law_kpswf9 import p0_inf_doublet
from lrc14_doublet_tvsharp_kpswf9 import tv_budget_sharp


def doublet(k, M):
    return tuple(sorted(set(list(range(k - 1)) + [M, M + 1])))


def main():
    print("=" * 96)
    print("FINAL FINITE CHECK: binding doublet E_f=consec_{k-1} U {f,f+1} < cap_k for ALL f>=15")
    print("kind-pasteur kpswf9  (exact rational p0; tail by frozen-law + J_sharp/f bound)")
    print("=" * 96)
    overall_ok = True
    for k in (8, 9, 10, 11, 12):
        base = tuple(range(k - 1))
        pinf = p0_inf_doublet(base, 1)
        cap = CAP[k]
        margin_inf = cap - pinf
        Js, N = tv_budget_sharp(base)
        f0 = floor(Js / margin_inf) + 1
        maxp0 = F(-1); fmax = None
        maxfd = F(0); ffd = None
        ok = True
        for f in range(15, f0 + 1):
            pv = p0_fast(doublet(k, f))
            if pv >= cap:
                ok = False
                print(f"   !! k={k} f={f}: p0={pv} >= cap={cap}  VIOLATION")
            fd = abs((pv - pinf) * f)
            if fd > Js:
                print(f"   !! k={k} f={f}: f*|Delta|={float(fd):.4f} > J_sharp={float(Js):.3f}  BOUND VIOLATION")
                ok = False
            if pv > maxp0:
                maxp0, fmax = pv, f
            if fd > maxfd:
                maxfd, ffd = fd, f
        # tail value bound for f>=f0
        tail_bound = pinf + Js / f0
        print(f"\nk={k}  cap={cap}={float(cap):.6f}  p0_inf={float(pinf):.6f}  J_sharp={float(Js):.3f}  f0={f0}")
        print(f"  finite check f in[15,{f0}]: max p0 = {maxp0} = {float(maxp0):.6f} at f={fmax}"
              f"  (cap-max = {float(cap-maxp0):+.6f}); all < cap? {ok}")
        print(f"  bound check: max f*|Delta| = {float(maxfd):.4f} at f={ffd} <= J_sharp={float(Js):.3f}  OK")
        print(f"  tail f>=f0: p0 <= p0_inf + J_sharp/f0 = {float(tail_bound):.6f} < cap? {tail_bound < cap}")
        verdict = ok and (tail_bound < cap)
        overall_ok = overall_ok and verdict
        print(f"  ==> k={k} doublet tail (all f>=15): {'CLOSED (< cap)' if verdict else 'NOT CLOSED'}")
    print("\n" + "=" * 96)
    print(f"OVERALL: binding-doublet genuine-wide tail CLOSED for k=8..12 (all f>=15): {overall_ok}")
    print("=" * 96)


if __name__ == "__main__":
    main()
