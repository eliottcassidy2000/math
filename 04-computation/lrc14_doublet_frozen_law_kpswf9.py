#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""kind-pasteur kpswf9 -- THREAD B: the EXACT frozen-phase f->inf law for the binding doublet.

SETUP (matches p0_fast EXACTLY):
  sector of speed e at point x in [0,1):  s_e(x) = floor(e*x) mod 7.
  p0(E) = measure{ x in [0,1) : { s_e(x) : e in E } ⊇ {1,2,3,4,5,6} }   (the 6 inner sectors).

DOUBLET far pair {f, f+g}, g=1 (binding/adjacent). As f->inf, the FAST phase frac(f*x)
equidistributes and the relevant SLOW coordinate is x itself. Two exact facts (for x in [0,1)):
  s_f(x)     = floor(f*x) mod 7               -- equidistributes over Z/7 as f->inf (fixed slow x)
  s_{f+1}(x) = floor((f+1)x) mod 7 = (s_f(x) + [frac(f*x) >= 1-x]) mod 7
So the doublet sectors are { s , s + b } with s = s_f(x) uniform in Z/7 (limit) and
  b = [frac(f*x) >= 1-x] in {0,1}, P(b=1) = x  (frac(f*x) uniform, INDEPENDENT of s in the limit).

THE FROZEN LIMITING LAW (g=1). Let M(x) = {1..6} \ Base(x) = sectors the BASE misses
(Base(x) = {floor(e*x) mod 7 : e in base}). The far pair must supply ALL of M(x). With
t=|M(x)|, in the f->inf limit (s uniform Z/7, b Bernoulli(x), independent):
  t=0: cover prob 1.
  t=1 (missing sector m): {s} covers iff s=m (prob 1/7 over the b=0 branch, weight 1-x);
       {s,s+1} covers iff s in {m-1,m} (prob 2/7 over the b=1 branch, weight x).
       => Pcov = (1-x)*(1/7) + x*(2/7) = (1+x)/7.
  t=2 missing {m1,m2}: a single sector covers <=1 of 2 => need b=1 AND {s,s+1}={m1,m2}
       i.e. the two missing sectors are CYCLICALLY CONSECUTIVE and s = the lower one.
       => Pcov = x*(1/7)  if {m1,m2} consecutive (cyclically in Z/7), else 0.
  t>=3: 0.
So
  p0_inf(base, g=1) = INT_0^1 [ 1*1{t=0} + ((1+x)/7)*1{t=1} + (x/7)*1{t=2 & consec} ] dx.

The base Base(x) is piecewise-constant in x (breakpoints = base-alone wall breakpoints), so each
integral is a sum over pieces of c0*(length) + c1*INT x + c2*INT x with the right coefficients.
EXACT rationals throughout.

This script:
  (A) VERIFY the frozen law against p0_fast(base U {f,f+1}) for large f (should -> p0_inf).
  (B) compute p0_inf EXACTLY for base=consec_{k-1}, k=8..12 (and the dilated/other bases).
  (C) report p0_inf vs cap_k -- is p0_inf < cap with margin?
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

ALL_INNER = 0b1111110


def lcm_list(xs):
    return reduce(lambda a, b: a // gcd(a, b) * b, xs, 1)


def base_pieces(base):
    """Return list of (lo, hi, missing_set) over x in [0,1) as EXACT Fractions, where
    missing_set = inner sectors {1..6} NOT hit by the base on the open interval (lo,hi).
    Breakpoints are the base-alone wall breakpoints (where some floor(e*x) jumps)."""
    nz = [int(e) for e in base if e]
    if not nz:
        # base = {0}: floor(0*x)=0, sector 0, misses all of 1..6 everywhere
        return [(F(0), F(1), frozenset(range(1, 7)))]
    l = lcm_list(nz)
    d = 7 * l  # work with x = num/d
    bps = {0, d}
    for e in nz:
        step = l // e  # in d-units, breakpoints of floor(e*x) at x = j/(7e) => num = j*step
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x)
            x += step
    bps = sorted(b for b in bps if 0 <= b <= d)
    pieces = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = lo + hi  # 2*midpoint in d-units; x_mid = mid/(2d)
        hit = set()
        for e in nz:
            # floor(e * x_mid) mod 7, x_mid = mid/(2d), e*x_mid = e*mid/(2d)
            hit.add((e * mid // (2 * d)) % 7)
        missing = frozenset(s for s in range(1, 7) if s not in hit)
        pieces.append((F(lo, d), F(hi, d), missing))
    return pieces


def is_cyclic_consecutive(pair):
    """Two sectors in Z/7 (values in 1..6) are CYCLICALLY consecutive iff difference is 1 mod 7."""
    a, b = sorted(pair)
    return (b - a) == 1 or (a + 7 - b) == 1  # adjacent on the 7-cycle


def p0_inf_doublet(base, g=1):
    """EXACT f->inf frozen-phase limit of p0(base U {f, f+g}) for g=1 (adjacent doublet)."""
    assert g == 1, "frozen law below derived for g=1 (binding/adjacent)"
    total = F(0)
    for lo, hi, missing in base_pieces(base):
        t = len(missing)
        L = hi - lo
        intx = (hi * hi - lo * lo) / 2  # INT_lo^hi x dx
        if t == 0:
            total += L                       # coeff 1
        elif t == 1:
            total += (L + intx) / 7          # INT (1+x)/7
        elif t == 2:
            if is_cyclic_consecutive(missing):
                total += intx / 7            # INT x/7
            # else 0
        # t>=3: 0
    return total


def refined_profile(base):
    """For reporting: (q0, q1, q2_consec, q2_nonconsec, q_ge3) measures of the base."""
    q0 = q1 = q2c = q2n = qg3 = F(0)
    for lo, hi, missing in base_pieces(base):
        L = hi - lo
        t = len(missing)
        if t == 0:
            q0 += L
        elif t == 1:
            q1 += L
        elif t == 2:
            if is_cyclic_consecutive(missing):
                q2c += L
            else:
                q2n += L
        else:
            qg3 += L
    return q0, q1, q2c, q2n, qg3


def verify_frozen(base, fs=(101, 211, 503, 1009, 2003, 5003, 10007)):
    """Compare p0_fast(base U {f,f+1}) to p0_inf as f grows."""
    pinf = p0_inf_doublet(base, 1)
    rows = []
    for f in fs:
        E = tuple(sorted(set(list(base) + [f, f + 1])))
        if reduce(gcd, [e for e in E if e]) != 1:
            # make primitive by ensuring base has a 1 or gcd 1; consec bases already have 1
            pass
        pv = p0_fast(E)
        rows.append((f, pv, float(pv - pinf)))
    return pinf, rows


def main():
    print("=" * 92)
    print("THREAD B: EXACT FROZEN-PHASE f->inf LAW for the binding doublet  (kind-pasteur kpswf9)")
    print("  sector s_e(x)=floor(e*x) mod 7;  far pair {f,f+1} -> {s, s+[frac(f x)>=1-x]}, s~Unif Z/7")
    print("=" * 92)

    # ---------- (A) VERIFY the frozen law against p0_fast at large f ----------
    print("\n(A) VERIFY frozen law: p0_fast(consec_{k-1} U {f,f+1}) -> p0_inf  as f->inf")
    for k in (8, 9, 10):
        base = tuple(range(k - 1))
        pinf, rows = verify_frozen(base)
        print(f"\n  k={k}  base=consec_{k-1}={base}")
        print(f"    p0_inf (exact)   = {pinf} = {float(pinf):.8f}")
        print(f"    {'f':>7} {'p0_fast(f)':>14} {'p0_fast - p0_inf':>20} {'f*(p0-p0_inf)':>16}")
        for f, pv, diff in rows:
            print(f"    {f:>7} {float(pv):>14.8f} {diff:>20.8f} {f*diff:>16.5f}")

    # ---------- (B) EXACT p0_inf for the binding consec doublet, k=8..12 ----------
    print("\n" + "=" * 92)
    print("(B) EXACT p0_inf for base=consec_{k-1}, doublet g=1  vs cap_k")
    print("=" * 92)
    print(f"  {'k':>3} {'p0_inf (exact)':>22} {'p0_inf~':>10} {'cap_k':>10} {'cap-p0_inf':>12} {'Q(k-1)':>9} {'<cap?':>6}")
    results = {}
    for k in (8, 9, 10, 11, 12):
        base = tuple(range(k - 1))
        pinf = p0_inf_doublet(base, 1)
        cap = CAP[k]
        q = QVAL[k]
        results[k] = (pinf, cap)
        print(f"  {k:>3} {str(pinf):>22} {float(pinf):>10.6f} {float(cap):>10.6f} "
              f"{float(cap - pinf):>12.6f} {float(q):>9.6f} {str(pinf < cap):>6}")
        rp = refined_profile(base)
        print(f"        refined base profile (q0,q1,q2consec,q2nonconsec,q>=3) = "
              f"({float(rp[0]):.4f},{float(rp[1]):.4f},{float(rp[2]):.4f},{float(rp[3]):.4f},{float(rp[4]):.4f})")

    print("\n" + "=" * 92)
    print("READING: if p0_fast(f) -> p0_inf in (A) and p0_inf < cap with margin in (B), the frozen")
    print("law is the correct f->inf target. The finite-f correction is handled in a separate script.")
    print("=" * 92)


if __name__ == "__main__":
    main()
