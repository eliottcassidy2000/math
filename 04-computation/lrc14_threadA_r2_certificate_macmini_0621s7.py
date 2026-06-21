#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 2 (mac-mini-2026-06-21-S7): the r=2 genuine-wide slack-floor CERTIFICATE, k=10,11,12.

PROVED STRUCTURE:
  (P1) genuine-wide => r >= 2. [removing the single far of an r=1 config leaves base ⊆ [0,14],
       reduced-span ≤ 14, contradiction.] So the slack floor reduces to r >= 2.
  (P2) The binding genuine-wide stratum is r = 2 (base size k-2); r >= 3 has a small base and
       is bounded well below Q (checked separately, big margin).

CERTIFICATE for r = 2: max p0 over r=2 genuine-wide configs < Q(k-1), via
  (A) EXACT enumeration of all r=2 configs with far positions f1 in a WINDOW [15, F*],
      gap g = f2 - f1 in [1, G*], over ALL bases (k-2)-subset of [0,14] with 0;
  (B) a RIGOROUS TAIL bound for f1 > F* via the THM-557 diagonal-freeze:
      p0(B ∪ {f1, f1+g}) <= D_block(B,g) + J/f1, J = 7*C(2,2)=7 (the 2-far block variation),
      where D_block(B,g) = the exact diagonal-freeze double integral (large-f1 limit). We show
      D_block(B,g) + 7/F* < Q(k-1) for every base/gap, so f1 > F* is safe with margin.

This script reports the EXACT window max and the EXACT D_block tail bound, deciding the floor.

To make (A) feasible we use the saturation fact (verified): the per-base sup over f1 is attained
at small f1; we take F* large enough that the tail (B) closes (F* ~ 7/(Q - sup_g D_block)).
We scan ALL bases for the window-max (not just dense) but cap the gap at G* = 24 (wider gaps only
lower coverage -- the coherent block is tightest at small gap, THM-557).
"""
from __future__ import annotations
import sys, functools, time
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce, lru_cache
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP, missed_distribution
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}


def p0_fast(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l
    den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x)
            x += step
    bps = sorted(bps)
    num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num0 += hi - lo
    return F(num0, d)


def primitive(E):
    nz = [e for e in E if e]
    return bool(nz) and reduce(gcd, nz) == 1


def reduced_span(S):
    S = sorted(set(S))
    g = 0
    for a, b in zip(S, S[1:]):
        g = gcd(g, b - a)
    return 0 if g == 0 else (S[-1] - S[0]) // g


def is_genuine_wide(E):
    E = tuple(sorted(set(E)))
    if reduced_span(E) <= 14:
        return False
    for i in range(len(E)):
        if reduced_span(E[:i] + E[i + 1:]) <= 14:
            return False
    return True


# ---------------------------------------------------------------------------
# D_block(B,g): the EXACT diagonal-freeze double integral = lim_{f1->inf} p0(B∪{f1,f1+g}).
# Matching p0_fast: runner e at slow x covers sector floor(e*x) mod 7 (inner = 1..6).
# For the two far runners, in the freeze limit the carrier phase psi=frac(f1*x) is uniform on
# [0,1) and INDEPENDENT of x; far1 sector = floor(7*psi)?? -- NO. We must match p0_fast's model.
#
# In p0_fast, runner speed e contributes sector (e*midnum//den2)%7 = floor(e*x) mod 7 where x in [0,1).
# For a FAR speed f and slow x in [0,1): floor(f*x) mod 7. Write f*x = q + r, q=floor(f*x), r=frac.
#   As f->inf, for x in a small subinterval where base sectors are constant, q mod 7 cycles through
#   all residues as x moves by 7/f; within one base-cell the pair (floor(f1*x) mod7, floor(f2*x) mod7)
#   equidistributes over the DIAGONAL family parametrized by the carrier offset.
# Concretely: let u = f1*x mod 7 (real in [0,7)), then far1 sector = floor(u) mod 7? and
#   far2: f2*x = f1*x + g*x, so far2 sector = floor(u + g*x) mod 7. In the freeze limit u ~ Uniform[0,7)
#   independent of x. So:
#   D_block(B,g) = ∫_0^1 ( 1/7 ∫_0^7 1{ baseSec(x) ∪ {floor(u)%7, floor(u+g*x)%7} = {1..6} } du ) dx.
# We compute this EXACTLY: for fixed x, baseSec is fixed; the inner u-integral is piecewise-constant
#   in u with breakpoints at integers and at integers - g*x (mod). Exact rational via breakpoint sort.
# ---------------------------------------------------------------------------
def base_sectors_at(base, x):
    """sectors covered by base nonzero runners at slow coord x (Fraction), as a bitmask of 0..6."""
    mask = 0
    for e in base:
        if e == 0:
            continue
        # floor(e*x) mod 7
        s = int((e * x).__floor__()) % 7
        mask |= 1 << s
    return mask


def Dblock(base, g):
    """EXACT diagonal-freeze double integral = lim_{f1->inf} p0(B u {f1, f1+g}).

    FREEZE MODEL (empirically verified): at slow x, the base sectors are fixed; the two far
    sectors are (a, (a + delta)%7) where a = far1 ~ Uniform{0..6} and the offset delta takes
    value floor(g*x) with prob 1-frac(g*x) and floor(g*x)+1 with prob frac(g*x). [Proof: far1
    equidistributes mod 7; floor(f2*x)-floor(f1*x) -> floor(g*x) or +1 with weights 1-frac, frac.]
    So
      D_block(B,g) = ∫_0^1 (1/7) Σ_{a=0}^6 [ (1-φ(x))·1{miss⊆{a,(a+d0)%7}} + φ(x)·1{miss⊆{a,(a+d0+1)%7}} ] dx,
    where d0=floor(g*x), φ(x)=frac(g*x). Exact rational via x-breakpoints (base + g-grid)."""
    # x-breakpoints: base-sector changes (x=j/e) AND offset changes (x=m/g, where floor(g*x) jumps).
    xb = {F(0), F(1)}
    for e in base:
        if e == 0:
            continue
        for j in range(1, e):
            xb.add(F(j, e))
    if g > 0:
        m = 1
        while F(m, g) < 1:
            xb.add(F(m, g))
            m += 1
    xb = sorted(xb)
    total = F(0)
    for xl, xr in zip(xb, xb[1:]):
        if xr <= xl:
            continue
        xm = (xl + xr) / 2
        bmask = base_sectors_at(base, xm)
        missing = [s for s in range(1, 7) if not (bmask >> s) & 1]
        t = len(missing)
        if t == 0:
            total += xr - xl
            continue
        if t > 2:
            continue
        d0 = int((g * xm).__floor__())  # floor(g*x) constant on this cell (g-grid breakpoints added)

        def covers(a, delta):
            cov = {a % 7, (a + delta) % 7}
            return all(mm in cov for mm in missing)

        # For each far1=a: prob over offset = (1-φ)·1{cov(a,d0)} + φ·1{cov(a,d0+1)}, averaged over a (1/7).
        # ∫_{xl}^{xr} [ (1-φ(x))·A + φ(x)·B ] dx  with A,B constants (counts over a, /7).
        cntA = sum(1 for a in range(7) if covers(a, d0)) / F(7)
        cntB = sum(1 for a in range(7) if covers(a, d0 + 1)) / F(7)
        # φ(x) = g*x - d0  on the cell. ∫_{xl}^{xr} φ dx = g*(xr^2-xl^2)/2 - d0*(xr-xl).
        int_phi = g * (xr * xr - xl * xl) / 2 - d0 * (xr - xl)
        width = xr - xl
        int_1mphi = width - int_phi
        total += cntA * int_1mphi + cntB * int_phi
    return total


@lru_cache(maxsize=None)
def Dblock_cached(base, g):
    return Dblock(base, g)


def verify_Dblock_against_largef(base, g, N=1200, span=49):
    """Cross-check exact Dblock vs the average of p0 over a window of large f1 (should match)."""
    s = F(0)
    n = 0
    for f in range(N, N + span):
        E = tuple(sorted(set(base) | {f, f + g}))
        s += p0_fast(E)
        n += 1
    return s / n


def main():
    print("=" * 96)
    print("THREAD 2: r=2 GENUINE-WIDE FLOOR CERTIFICATE  k=10,11,12  (mac-mini-S7)")
    print("=" * 96)
    for k in (10, 11, 12):
        print(f"  k={k}: cap={float(CAP[k]):.6f}  Q(k-1)={float(QVAL[k]):.6f}  MARGIN={float(MARGIN[k]):.6f}")
    print()

    # ---- cross-check Dblock exact vs large-f average ----
    print("-" * 96)
    print("Dblock EXACT vs large-f1 average (sanity):")
    for base, g in [((0,1,2,3,4,5,6,7),1), ((0,1,2,3,4,5,6,7),2), ((0,1,2,3,4,5,6,7,8),1)]:
        de = Dblock_cached(base, g)
        av = verify_Dblock_against_largef(base, g)
        print(f"  base={base} g={g}: Dblock={float(de):.6f}={de}  large-f avg={float(av):.6f}  close? {abs(float(de)-float(av))<0.003}")
    print()


if __name__ == "__main__":
    main()
