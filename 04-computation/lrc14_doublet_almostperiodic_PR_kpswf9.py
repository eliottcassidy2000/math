#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A (kps-2026-06-21-Swf9): the ALMOST-PERIODIC P/R DECOMPOSITION of the binding
genuine-wide DOUBLET signed error.

THE ONE OPEN LEG (OPEN-Q-108, leg C / genuine-wide). The binding maximizer is
    E_M = consec_{k-1} U {M, M+1}   (consec base + a TIGHT far doublet, opus HYP-2794).
THM-563 closes SINGLE-far by exact periodicity of w*Delta_w. The DOUBLET is NOT exactly
periodic (refuted: 2nd-diff ~0.01-0.03 nonzero, decaying) -> ALMOST-periodic.

DECOMPOSITION (this script, exact rationals):
    g(M) := M * Delta_M    where Delta_M = p0(E_M) - Phi      (Phi = exact plateau, below)
    g(M) = P(M) + R(M)
  - Phi  = the EXACT plateau lim_{M->inf} p0(E_M). We compute it EXACTLY (not a late-window
           average) as the FROZEN-doublet boundary value: the measure where the base hits its
           sectors AND the far doublet {M,M+1}, in the M->inf equidistributed limit, hits the
           remaining sectors with the frozen two-far joint law over the slow phase y in [0,7).
           Equivalently Phi = mean over y in[0,7) of [base covers via its arcs] x [two frozen
           far runners at speeds (1,1) offset cover the rest]. We compute it via the
           period-exact integral of the FROZEN model (period 7*lcm(base) in the base phase,
           and an explicit average over the frozen far phase).
  - P(M) = the exactly-periodic part: the contribution of the BASE arcs (breakpoints j/(7e),
           e in base) to M*Delta_M. By the THM-563 sawtooth mechanism this is periodic in M
           with period 7*lcm(base). We isolate it as the part of g(M) supported on base
           breakpoints (the far runners enter only through which sectors are "already covered").
  - R(M) = g(M) - P(M): the correction from the w-DEPENDENT breakpoints j/(7M), j/(7(M+1))
           where the two far runners interact with each other and with the base (HYP-2794 part2).

This script:
  (1) compute Phi EXACTLY (frozen-doublet integral) for k=8..12 and compare to opus late-window.
  (2) compute g(M)=M*Delta_M exactly; build P(M) by period-averaging g over the base period
      P_per=7*lcm(base) (the unique periodic component); R(M)=g(M)-P(M).
  (3) DECAY TEST: tabulate M*|R(M)| and M^2*|R(M)|. Is R(M)=O(1/M)? sign-controlled?
  (4) period-max(P): exact max of the periodic part over one period (finite computation).
  (5) tail bound on sup_{M>=M0}|R(M)|.
  (6) ASSEMBLE: p0(E_M) = Phi + g(M)/M = Phi + (P(M)+R(M))/M. Show < cap_k for M>=15 via
      [period-max(P) + sup_{M>=M0}|R|]/M < margin (target 0.16); finite window [15,M0] exact.

Exact rationals throughout. p0 engine = repo p0_fast (matches measS7 / p0_repo exactly).
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL, MARGIN

ALL_INNER = 0b1111110


def lcm_list(xs):
    return reduce(lambda a, b: a * b // gcd(a, b), xs, 1)


def doublet(k, M):
    return tuple(sorted(set(list(range(k - 1)) + [M, M + 1])))


# ----------------------------------------------------------------------------------
# EXACT PLATEAU Phi = lim_{M->inf} p0(E_M), via the FROZEN-DOUBLET model.
#
# p0(E) = measure{x in [0,1): for every inner sector s in {1..6}, some e in E has
#          floor(7*e*x) ~ in sector s}, i.e. {e*x mod 1 : e in E} touches all 6 inner
#          1/7-arcs.  (Sector of e at x = floor(7*(e*x mod 1)) but the engine uses
#          (e*midnum//den2)%7 == floor(7*e*x) mod 7; inner sectors are bits 1..6.)
#
# As M->inf, the two far runners M, M+1 have phases (M*x mod1, (M+1)*x mod1). Write
# the slow variable y = x (the base sees e*x for small e). On a fine base-cell the base
# sectors are fixed; the far phases (M*x, (M+1)*x) = (theta, theta + x) with theta = M*x
# mod1 equidistributing INDEPENDENTLY of x within the cell as M->inf, but the OFFSET
# between the two far runners is (M+1)*x - M*x = x = y (the slow base coordinate). So the
# frozen far law at base-coordinate y is: pick theta uniform in [0,1); the two far runners
# sit at sectors floor(7*theta) and floor(7*(theta+y)). They cover sector pair
#   {floor(7 theta) , floor(7(theta+y)) mod ... } -- a y-dependent FROZEN doublet law.
#
# Phi = INTEGRAL over y in[0,1) of  Pr_theta[ base_sectors(y) U {s_far1, s_far2} >= all 6 inner ].
# This is exact: for fixed y, base sectors are piecewise-constant in y (breakpoints j/(7e),
# e in base), and the theta-probability that the two frozen far runners supply EXACTLY the
# missing inner sectors is a rational (sum of arc lengths). We compute it exactly.
# ----------------------------------------------------------------------------------

def base_sectors_at(base_nz, y_lo, y_hi):
    """Inner-sector bitmask covered by the BASE (nonzero speeds) on the open base-cell
    (y_lo,y_hi). Uses the cell midpoint. Returns mask over bits 0..6 (bit s set if some
    base runner is in sector s); inner = bits 1..6."""
    # midpoint as exact fraction
    ymid = (y_lo + y_hi) / 2
    mask = 0
    for e in base_nz:
        # sector = floor(7 * frac(e*ymid)) ; engine convention floor(7 e y) mod 7
        val = e * ymid  # fraction
        sec = int((val * 7)) % 7  # floor(7 e y) mod 7  (val*7 fraction -> floor via int)
        # int() of a Fraction truncates toward zero; e*ymid>=0 so == floor
        mask |= 1 << sec
    return mask


def frozen_far_cover_prob(missing_mask, y):
    """Pr_{theta in[0,1)} that two frozen far runners at sectors a=floor(7 theta),
    b=floor(7(theta+y)) cover the set `missing_mask` (subset of inner bits 1..6).
    Exact rational. theta ranges [0,1); a,b are 1/7-sector indices.
    For the missing set to be covered by just {a,b}: need missing_mask subset {a,b}.
    => |missing|<=2 and {a,b} contains all missing sectors.
    We integrate over theta in 1/7-cells, with b determined by theta and y."""
    miss = [s for s in range(7) if (missing_mask >> s) & 1]
    if len(miss) > 2:
        return F(0)
    miss_set = set(miss)
    # theta in [0,1): break into cells where a=floor(7 theta) constant AND
    # b=floor(7(theta+y)) constant. b-cell breakpoints at theta = (j - 7y)/7 mod ...
    # Simplest exact approach: a in 0..6 over theta in [a/7,(a+1)/7). Within, theta+y has
    # floor(7(theta+y)) = floor(7 theta + 7y). Let fy = 7*y (a Fraction). For theta in
    # [a/7,(a+1)/7), 7 theta in [a, a+1). 7 theta + 7y in [a+7y, a+1+7y). b jumps once at
    # the integer crossing inside. Handle by subdividing.
    fy = 7 * y  # Fraction
    total = F(0)
    for a in range(7):
        # theta-cell [a/7,(a+1)/7); parameter u=7 theta in [a,a+1).
        # b = floor(u + fy). On u in [a,a+1), u+fy in [a+fy, a+1+fy); integer breakpoints:
        lo_u = F(a)
        hi_u = F(a + 1)
        # breakpoints where floor(u+fy) jumps: u = n - fy for integer n in (a+fy, a+1+fy)
        import math
        nlo = math.floor(a + fy) + 1
        nhi = math.ceil(a + 1 + fy) - 1
        bps_u = [lo_u]
        for n in range(nlo, nhi + 1):
            ub = F(n) - fy
            if lo_u < ub < hi_u:
                bps_u.append(ub)
        bps_u.append(hi_u)
        bps_u = sorted(set(bps_u))
        for ulo, uhi in zip(bps_u, bps_u[1:]):
            if uhi <= ulo:
                continue
            umid = (ulo + uhi) / 2
            b = int(umid + fy) % 7  # floor(umid+fy) mod7
            cover = {a % 7, b}
            if miss_set <= cover:
                total += (uhi - ulo) / 7  # d theta = du/7
    return total


def phi_exact(k):
    """EXACT plateau Phi = lim p0(E_M) via frozen-doublet integral over y in[0,1)."""
    base_nz = list(range(1, k - 1))  # 1..k-2
    # base-cell breakpoints in y: j/(7e) for e in base_nz, plus far-OFFSET breakpoints
    # for the frozen law these come from frozen_far_cover_prob's internal handling; but
    # the base mask changes at j/(7e). Also frozen prob changes with y (via fy=7y crossing).
    l = lcm_list(base_nz)
    den = 7 * l
    bps = {F(0), F(1)}
    for e in base_nz:
        for j in range(7 * e + 1):
            bps.add(F(j, 7 * e))
    # the frozen far prob has y-breakpoints where fy=7y hits the cell structure; those are
    # at y = n/7 (i.e. 7y integer) and at sector reorganizations; include y=n/7.
    for n in range(8):
        bps.add(F(n, 7))
    bps = sorted(b for b in bps if F(0) <= b <= F(1))
    total = F(0)
    for ylo, yhi in zip(bps, bps[1:]):
        if yhi <= ylo:
            continue
        bmask = base_sectors_at(base_nz, ylo, yhi) & ALL_INNER
        if (bmask).bit_count() == 6:
            # base alone covers all 6 inner -> far irrelevant, full cell counts
            total += (yhi - ylo)
            continue
        missing = (~bmask) & ALL_INNER
        ymid = (ylo + yhi) / 2
        p = frozen_far_cover_prob(missing, ymid)
        total += (yhi - ylo) * p
    return total


def main():
    print("=" * 84)
    print("DOUBLET ALMOST-PERIODIC P/R DECOMPOSITION  (THREAD A, kps-Swf9)")
    print("E_M = consec_{k-1} U {M, M+1}    g(M)=M*(p0(E_M)-Phi) = P(M)+R(M)")
    print("=" * 84)

    # ---- STEP 1: exact plateau Phi vs opus late-window estimate ----
    print("\n[STEP 1] EXACT plateau Phi (frozen-doublet integral) vs opus late-window avg")
    print("-" * 84)
    opus_phi2 = {8: 0.261700, 9: 0.411009, 10: 0.490553, 11: 0.565827, 12: 0.629083}
    PHI = {}
    for k in range(8, 13):
        ph = phi_exact(k)
        PHI[k] = ph
        cap = CAP[k]
        margin = cap - ph
        print(f"  k={k}: Phi_exact = {ph} = {float(ph):.6f}   (opus late-win {opus_phi2[k]:.6f}, "
              f"diff {float(ph)-opus_phi2[k]:+.6f})")
        print(f"        cap_{k}={float(cap):.6f}  margin=cap-Phi={float(margin):.6f}  target 0.16")
    return PHI


if __name__ == "__main__":
    main()
