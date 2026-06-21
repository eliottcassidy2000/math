#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- THREAD B: the EXACT frozen-phase f->inf law for the binding doublet.

SECTOR CONVENTION (matches p0_fast / missed_distribution EXACTLY, verified):
  point coordinate x in [0,1); SLOW coordinate y = 7*x in [0,7).
  sector of speed e:  s_e = floor(7*e*x) mod 7 = floor(e*y) mod 7   (== codex sector(e,y)).
  p0(E) = measure_x{ {s_e : e in E} ⊇ {1..6} } = (1/7) * measure_y{ ... },  y in [0,7).

DOUBLET far pair {f, f+1}, g=1. As f->inf the FAST phase u = frac(7 f x) = frac(f y) equidistributes
in [0,1), INDEPENDENT (in the limit) of s_f = floor(f y) mod 7 ~ Unif(Z/7). EXACT increment law:
  s_{f+1} = floor((f+1) y) mod 7 = ( s_f + floor( u + frac(f... )) ) ...  cleanly:
  floor((f+1) y) = floor(f y + y) = floor(f y) + floor( frac(f y) + y ).
  With u = frac(f y) in [0,1) and y in [0,7): floor(u + y) = floor(y) + [ u >= 1 - frac(y) ].
  => increment  delta = s_{f+1} - s_f  ==  ( floor(y) + [u >= 1-frac(y)] )  mod 7.
So in the f->inf limit, at slow y:
  delta = floor(y) mod 7         with prob 1 - frac(y),
  delta = (floor(y)+1) mod 7     with prob frac(y);
  s = s_f ~ Unif(Z/7), independent of delta.
The far pair contributes the sector SET {s, (s+delta) mod 7}.

THE EXACT FROZEN LIMITING LAW.  Base(y) = {floor(e*y) mod 7 : e in base}, M(y)={1..6}\Base(y).
For each y, average over s~Unif(Z/7) and the two-branch delta, the indicator that {s,s+delta} ⊇ M(y):
  Pcov(y) = (1-frac(y)) * A(M(y), floor(y) mod 7)  +  frac(y) * A(M(y), (floor(y)+1) mod 7),
  A(M, d) = (1/7) * #{ s in Z/7 : M ⊆ {s, (s+d) mod 7} }.
  p0_inf(base, g=1) = (1/7) * INT_0^7 Pcov(y) dy.
Base(y) is piecewise-constant in y (breakpoints = where some floor(e*y) jumps), so each integral is
exact: on a piece [yl,yh) with constant M and constant floor(y) (refine breakpoints to integers too,
since floor(y) and frac(y) appear), INT (1-frac(y)) and INT frac(y) are exact rationals.

EXACT rationals throughout. (A) verifies vs p0_fast at large f; (B) computes p0_inf vs cap_k.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive
from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, QVAL


def lcm_list(xs):
    return reduce(lambda a, b: a // gcd(a, b) * b, xs, 1)


def A_cover(M, d):
    """(1/7) * #{ s in Z/7 : M ⊆ {s, (s+d) mod 7} }.  M = frozenset of inner sectors (1..6)."""
    cnt = 0
    for s in range(7):
        pair = {s, (s + d) % 7}
        if M <= pair:
            cnt += 1
    return F(cnt, 7)


def base_y_pieces(base):
    """Pieces over y in [0,7) where Base(y)={floor(e*y) mod7 : e in base} is constant AND floor(y)
    is constant. Returns list of (yl, yh, M, floor_y) with EXACT Fraction endpoints.
    Breakpoints: integer points 0..7 (for floor(y)) and y=j/e for e in base (where floor(e*y) jumps)."""
    nz = [int(e) for e in base if e]
    bps = set()
    for j in range(0, 8):
        bps.add(F(j))
    for e in nz:
        # floor(e*y) jumps at y = j/e, j=1..7e-1 within (0,7)
        for j in range(1, 7 * e):
            bps.add(F(j, e))
    bps = sorted(b for b in bps if F(0) <= b <= F(7))
    pieces = []
    for yl, yh in zip(bps, bps[1:]):
        if yh <= yl:
            continue
        ymid = (yl + yh) / 2
        hit = set()
        for e in nz:
            hit.add(int(e * ymid) % 7)   # floor(e*ymid) mod 7
        M = frozenset(s for s in range(1, 7) if s not in hit)
        fy = int(ymid)  # floor(y), constant on the piece since we split at integers
        pieces.append((yl, yh, M, fy))
    return pieces


def p0_inf_doublet(base, g=1):
    """EXACT f->inf frozen-phase limit of p0(base U {f, f+g}) for g=1 (adjacent doublet)."""
    assert g == 1
    total = F(0)
    for yl, yh, M, fy in base_y_pieces(base):
        d0 = fy % 7
        d1 = (fy + 1) % 7
        A0 = A_cover(M, d0)
        A1 = A_cover(M, d1)
        L = yh - yl
        # frac(y) = y - fy on this piece. INT_yl^yh frac(y) dy = INT (y-fy) dy
        int_frac = (yh * yh - yl * yl) / 2 - fy * L
        int_1mfrac = L - int_frac
        # Pcov(y) = (1-frac)*A0 + frac*A1 ; INT Pcov dy:
        total += A0 * int_1mfrac + A1 * int_frac
    return total / 7  # the (1/7) from dx = dy/7


def refined_profile(base):
    """Base measures by missing-sector count, for reporting (in x-measure, = (1/7)*y-length)."""
    from collections import Counter
    acc = Counter()
    for yl, yh, M, fy in base_y_pieces(base):
        acc[len(M)] += (yh - yl) / 7
    return acc


def verify_frozen(base, fs=(101, 503, 1009, 10007, 100003, 1000003)):
    pinf = p0_inf_doublet(base, 1)
    rows = []
    for f in fs:
        E = tuple(sorted(set(list(base) + [f, f + 1])))
        pv = p0_fast(E)
        rows.append((f, pv, float(pv - pinf)))
    return pinf, rows


def main():
    print("=" * 94)
    print("THREAD B: EXACT FROZEN-PHASE f->inf LAW for the binding doublet  (kind-pasteur kpswf9)")
    print("  sector s_e = floor(7 e x) mod 7 = floor(e y) mod 7, y=7x in [0,7); doublet increment delta")
    print("  delta = floor(y)+[u>=1-frac(y)] mod 7 (u=frac(f y) ~ Unif), s_f ~ Unif Z/7, independent.")
    print("=" * 94)

    # ---------- (A) VERIFY ----------
    print("\n(A) VERIFY frozen law: p0_fast(consec_{k-1} U {f,f+1}) -> p0_inf as f->inf")
    for k in (8, 9, 10):
        base = tuple(range(k - 1))
        pinf, rows = verify_frozen(base)
        print(f"\n  k={k}  base=consec_{k-1}={base}")
        print(f"    p0_inf (exact) = {pinf} = {float(pinf):.9f}")
        print(f"    {'f':>9} {'p0_fast(f)':>14} {'p0_fast-p0_inf':>16} {'f*(diff)':>12}")
        for f, pv, diff in rows:
            print(f"    {f:>9} {float(pv):>14.9f} {diff:>16.9f} {f*diff:>12.4f}")

    # ---------- (B) EXACT p0_inf vs cap ----------
    print("\n" + "=" * 94)
    print("(B) EXACT p0_inf for base=consec_{k-1}, doublet g=1  vs cap_k")
    print("=" * 94)
    print(f"  {'k':>3} {'p0_inf (exact)':>26} {'p0_inf~':>11} {'cap_k':>10} {'cap-p0_inf':>12} {'<cap?':>6}")
    results = {}
    for k in (8, 9, 10, 11, 12):
        base = tuple(range(k - 1))
        pinf = p0_inf_doublet(base, 1)
        cap = CAP[k]
        results[k] = (pinf, cap)
        print(f"  {k:>3} {str(pinf):>26} {float(pinf):>11.7f} {float(cap):>10.6f} "
              f"{float(cap - pinf):>12.6f} {str(pinf < cap):>6}")
    print("\n  margin cap_k - p0_inf (the robust target is >= 0.16):")
    for k in (8, 9, 10, 11, 12):
        pinf, cap = results[k]
        print(f"    k={k}: cap - p0_inf = {float(cap - pinf):.6f}  ({'>=0.16' if cap-pinf>=F(16,100) else '<0.16'})")
    print("\n" + "=" * 94)


if __name__ == "__main__":
    main()
