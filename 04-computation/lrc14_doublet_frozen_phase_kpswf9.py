#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD C (kps-wf9): the FROZEN-PHASE limit law for the binding doublet, and the
finite-f correction. This is the CORRECT decomposition (the naive 'peel M+1 vs S_M'
is loose: D(M) grows because S_M keeps the huge runner M -- HYP-2776).

SETUP. Position x in [0,1). Runner e sits in inner-sector  s_e(x) = floor(7*e*x) mod 7.
p0(E) = measure{ x : {s_e(x) : e in E, e>0} covers all 6 inner sectors 1..6 }.
(sector 0 = home; "covers all 6 inner" = the all-6 measure used repo-wide.)

The binding doublet: E_f = consec_{k-1} U {f, f+1},  consec_{k-1} = {0,1,...,k-2}.
Far pair {f, f+1}. Write the two far sectors via the SLOW/FAST split:
    7*f*x       = 7*u           where u = frac(f*x)   (FAST phase, equidistributes)
    7*(f+1)*x   = 7*f*x + 7*x   = 7*u + y   (mod 7),   y = frac(7*x) in [0,1)*7 -> use 7x mod 7
So  s_f(x)   = floor(7u) mod 7,
    s_{f+1}(x) = floor((7u + 7x) mod 7) ... carefully reduce mod 7.
The BASE sectors s_0..s_{k-2}(x) are SLOW (functions of x directly). The two far
sectors depend on the fast phase u AND the slow gap (7x mod 7).

FROZEN LAW (codex HYP-2796). As f->inf, for a.e. slow x the fast phase u
equidistributes uniformly on [0,1). So
    p0_inf := lim_{f->inf} p0(E_f)
            = INT_0^1 dx  [ over u in [0,1) uniform: P( base sectors(x) together with
                            {floor(7u), floor((7u+7x) mod 7)} cover all 6 inner ) ].
This is a DOUBLE integral (slow x, fast u), computable EXACTLY in rationals because
everything is piecewise-constant with breakpoints at rationals with denominator
dividing 7 (for u) and 7*lcm(base) (for x).

THIS SCRIPT:
  (1) compute p0_inf EXACTLY for k=8..12 via the frozen double integral.
  (2) compare to the empirical large-f plateau Phi_2 of the doublet (should match).
  (3) the FINITE-f correction c(f) = p0(E_f) - p0_inf: tabulate, test |c(f)| <= C/f
      (DECAY rate), report sup_f f*|c(f)|.
  (4) GENERAL GAP g: do the same frozen law for far pair {f, f+g}; confirm g=1 (adjacent)
      gives the LARGEST p0_inf (worst doublet), i.e. larger g decorrelates -> lower plateau.

If p0_inf < cap_k with margin >= 0.16 AND f*|c(f)| is bounded by C with
p0_inf + C/15 < cap_k, the doublet closes rigorously.
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

ALL_INNER = 0b1111110  # sectors 1..6


def lcm_list(xs):
    return reduce(lambda a, b: a * b // gcd(a, b), xs, 1)


def base_speeds(k):
    return list(range(1, k - 1))  # nonzero base speeds 1..k-2 (0 contributes sector 0 only)


def doublet(k, M, g=1):
    return tuple(sorted(set(list(range(k - 1)) + [M, M + g])))


@lru_cache(maxsize=None)
def p0c(E):
    return p0_fast(E)


def frozen_p0(k, g=1):
    """EXACT frozen-phase limit p0_inf for far pair {f, f+g}, f->inf.

    Slow x in [0,1): base sector of speed e is floor(7*e*x) mod 7. Far pair arguments:
        7*f*x        = 7u            (u = frac(f*x), the FAST phase)
        7*(f+g)*x    = 7u + 7*g*x    (mod 7)
    so  s_f   = floor(7u) mod 7,   s_{f+g} = floor((7u + 7*g*x) mod 7) mod 7.
    p0_inf = INT_x [ INT_u  1[base(x) U {s_f,s_{f+g}} covers 1..6] du ].

    EXACT and EFFICIENT: x is piecewise constant with denom dividing 7*lcm(base U {g}).
    For each x-cell (constant base mask + constant offset off=frac(7*g*x_mid)), the inner
    u-integral is exact: as u runs [0,1), the PAIR (s_f, s_{f+g}) is piecewise-constant with
    breakpoints at u = j/7 (j=1..6) and u = (j - off_frac)/7 i.e. where 7u or 7u+off crosses
    an integer. We enumerate those <=14 u-breakpoints, take each sub-interval's midpoint
    (exact rational), evaluate the pair, and accumulate the sub-interval length when covered.
    """
    base = base_speeds(k)
    Lx = lcm_list(base + [g])
    Dx = 7 * Lx  # x-cell count; midpoint of cell ix is (2ix+1)/(2Dx)
    total = F(0)
    for ix in range(Dx):
        xmid = F(2 * ix + 1, 2 * Dx)
        bmask = 0
        for e in base:
            bmask |= 1 << (int(7 * e * xmid) % 7)
        if (bmask & ALL_INNER).bit_count() == 6:
            total += 1  # base alone already covers -> full u-measure 1 contributes
            continue
        off = (7 * g * xmid) % 7  # in [0,7); the additive shift of the partner's argument
        # u-breakpoints in (0,1): j/7 for j=1..6 ; and where 7u+off crosses integer:
        #   7u = n - off  => u = (n-off)/7 for n with 0 < n-off < 7
        bps = {F(0), F(1)}
        for j in range(1, 7):
            bps.add(F(j, 7))
        # off has denom dividing 7*Lx/7 = Lx ... compute crossings
        n_lo = int(off) + 1
        for n in range(n_lo, n_lo + 8):
            u = (F(n) - off) / 7
            if 0 < u < 1:
                bps.add(u)
        bps = sorted(bps)
        covered_len = F(0)
        for lo, hi in zip(bps, bps[1:]):
            if hi <= lo:
                continue
            umid = (lo + hi) / 2
            sf = int(7 * umid) % 7
            sfg = int((7 * umid + off) % 7) % 7
            mask = bmask | (1 << sf) | (1 << sfg)
            if (mask & ALL_INNER).bit_count() == 6:
                covered_len += hi - lo
        total += covered_len  # inner u-integral over [0,1)
    return total / Dx  # average over x (each x-cell weight 1/Dx)


def main():
    print("=" * 90)
    print("THREAD C: FROZEN-PHASE limit law for the binding doublet  (kps-wf9)")
    print("E_f = consec_{k-1} U {f, f+g};  p0_inf = lim_{f->inf} p0(E_f) via slow(x)/fast(u) split")
    print("=" * 90)
    for k in range(8, 13):
        cap = CAP[k]; q = QVAL[k]
        print(f"\nk={k}  cap={cap}={float(cap):.6f}  Q(k-1)={float(q):.6f}")
        # frozen limit for g=1 (adjacent, worst) and g=2,3 to confirm decorrelation
        pinf1 = frozen_p0(k, g=1)
        print(f"  FROZEN p0_inf (g=1, adjacent) = {pinf1} = {float(pinf1):.6f}")
        print(f"     cap - p0_inf = {float(cap - pinf1):+.6f}   (margin if this is the plateau)")
        print(f"     Q(k-1) - p0_inf = {float(q - pinf1):+.6f}")
        # empirical large-f plateau check: p0 at a few large coprime f
        emp = [(f, float(p0c(doublet(k, f, 1)))) for f in (301, 401, 503, 601)]
        print(f"     empirical p0(E_f) at large f (g=1): {emp}")
        # decorrelation: g=2,3
        for g in (2, 3):
            pg = frozen_p0(k, g=g)
            print(f"  FROZEN p0_inf (g={g}) = {float(pg):.6f}   (should be <= g=1 if adjacent is worst)")
    print("\n" + "=" * 90)
    print("If p0_inf(g=1) < cap with margin and is the plateau, the doublet sup is governed")
    print("by p0_inf + finite-f correction. Next: bound the finite-f correction (decay).")


if __name__ == "__main__":
    main()
