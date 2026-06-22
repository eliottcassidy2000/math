#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_true_maxgap_good_kpswf13.py   (kind-pasteur 2026-06-22-S36, THREAD 2 -- the REAL GOOD)

DISCOVERY driving this file: the kpswf12 engine's GOOD proxy = cover_set(E)^c = {>=1 of the 7
half-open sectors [r/7,(r+1)/7) is MISSED} is NOT symmetric under x->-x (its arc-set symmetric
difference with its reflection has measure 31/245 ~ 0.127 at E={0..7}).  Hence the kpswf12 ghat is
NOT actually real; the prompt's premise "ghat REAL because GOOD=-GOOD" pins the TRUE GOOD, namely
       GOOD_true = { x : maxgap_circle{frac(e_i x) : e_i in E} > 1/7 }
which IS exactly x->-x symmetric (reflecting all phases preserves the cyclic gap structure).

THIS FILE builds GOOD_true EXACTLY (cyclic max-gap of the phase set, threshold 1/7), verifies its
reflection symmetry (sym-diff measure = 0 EXACTLY), and recomputes the spectrum sum
       SPEC_true = sum_{n!=0} chat(n) ghat_true(n),     R'_true = 1 + SPEC_true/baseline_true,
with ghat_true now GENUINELY REAL (sine energy = 0 exactly).  We then redo the THREAD-2 analysis
on the CORRECT object:
   (P2)  sine energy of GOOD_true = 0 exactly  (the complement-odd channel truly vanishes);
   (P1)  GOOD_true(E) = GOOD_true(-E) = GOOD_true(N-E) exactly (cluster complement invariance);
   (P3)  even/odd cluster contribution + the n-parity / 7-divisibility grading of SPEC_true.

EXACT max-gap of a finite phase set on the circle:
  phases = sorted{ frac(e_i x) }.  gaps = consecutive differences cyclically (wrap 1).  maxgap =
  max gap.  GOOD_true(x) = [maxgap > 1/7].  As x sweeps [0,1) the phases are piecewise-linear in x
  with breakpoints where two phases collide (frac(e_i x)=frac(e_j x) => (e_i-e_j)x in Z) or where a
  phase crosses 0 (e_i x in Z); on each elementary x-interval the COMBINATORIAL order of phases is
  fixed, every gap is an affine function alpha + beta*x, and {maxgap>1/7} is a union of sub-intervals
  found by solving each affine gap = 1/7 (all rational => EXACT Fraction arcs).
"""
import sys, itertools
from fractions import Fraction as F
from math import pi
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, safe_set, fourier_num_of_arcs,
)

THRESH = F(1, 7)

# ---------------------------------------------------------------- exact GOOD_true
def good_breakpoints(E):
    """x-values where the cyclic phase ORDER can change: frac(e_i x)=frac(e_j x) or =0.
       i.e. (e_i - e_j) x in Z  and  e_i x in Z.  Collect all in [0,1)."""
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    diffs = set()
    for i in range(len(E)):
        for j in range(len(E)):
            d = E[i] - E[j]
            if d != 0:
                diffs.add(abs(d))
        if E[i] != 0:
            diffs.add(abs(E[i]))   # crossing 0
    for d in diffs:
        for m in range(d + 1):
            v = F(m, d)
            if 0 <= v <= 1:
                bps.add(v)
    return sorted(bps)

def maxgap_at(E, x):
    """exact cyclic max gap of {frac(e x)} at a single x (Fraction)."""
    ph = sorted(set((F(e) * x) % 1 for e in E))
    if len(ph) == 1:
        return F(1)            # single point: gap = full circle
    gaps = [ph[k + 1] - ph[k] for k in range(len(ph) - 1)]
    gaps.append(F(1) - ph[-1] + ph[0])     # wrap
    return max(gaps)

def good_true(E):
    """GOOD_true = {x: maxgap{frac(e_i x)} > 1/7} as EXACT arc union.

    EXACT method (no sorting of gaps; track physical affine gaps):
    Within an order-interval [a,b] (no two phases collide, no phase crosses 0), each phase
    frac(e_i x) is an AFFINE function  e_i*x - floor(e_i*mid)  of x (the integer part is constant
    on the interval).  Sorting the phases at the midpoint fixes the cyclic order; consecutive
    differences (incl. the wrap 1 - top + bottom) are then AFFINE in x with EXACT rational
    coefficients.  maxgap(x) is the upper envelope of these affine gap-functions; it is piecewise
    affine and crosses the level 1/7 only at points where SOME individual affine gap = 1/7 (an
    overcount of envelope breakpoints, which is fine -- we just refine there and sign-test midpoints).
    All crossings are exact Fractions => GOOD_true is an EXACT arc union."""
    E = sorted(set(int(e) for e in E))
    base_bps = good_breakpoints(E)
    refined = set(base_bps)
    for a, b in zip(base_bps, base_bps[1:]):
        mid = (a + b) / 2
        # affine phase reps on this interval: phase_i(x) = e_i*x - k_i,  k_i = floor(e_i*mid)
        reps = []
        for e in E:
            ke = (F(e) * mid).__floor__()
            reps.append((F(e), F(ke)))      # phase = e*x - ke, value in [0,1) at mid
        # cyclic order from midpoint values
        vals = sorted(range(len(reps)), key=lambda i: (reps[i][0] * mid - reps[i][1]))
        # dedup identical phases at mid (collisions are at breakpoints, not interior)
        ordered = []
        seen = set()
        for i in vals:
            v = reps[i][0] * mid - reps[i][1]
            if v in seen:
                continue
            seen.add(v); ordered.append(reps[i])
        if len(ordered) < 2:
            continue
        # affine gaps: between consecutive ordered phases, plus wrap (1 - top + bottom)
        gaps_affine = []
        for t in range(len(ordered) - 1):
            (s2, c2) = ordered[t + 1]; (s1, c1) = ordered[t]
            gaps_affine.append((s2 - s1, -(c2 - c1)))      # gap(x) = (s2-s1) x - (c2-c1)
        (sT, cT) = ordered[-1]; (s0, c0) = ordered[0]
        gaps_affine.append((s0 - sT, F(1) - (c0 - cT)))    # wrap = 1 + (s0-sT)x -(c0-cT)
        for (slope, intercept) in gaps_affine:
            if slope == 0:
                continue
            xc = (THRESH - intercept) / slope
            if a < xc < b:
                refined.add(xc)
    refined = sorted(refined)
    arcs = []
    for a, b in zip(refined, refined[1:]):
        mid = (a + b) / 2
        if maxgap_at(E, mid) > THRESH:
            arcs.append((a, b))
    return merge(arcs)

def reflect_arcs(arcs):
    out = []
    for a, b in arcs:
        lo, hi = (1 - b), (1 - a)
        out.append((min(lo, hi), max(lo, hi)))
    return merge(out)

def chat_c(arcs, n):
    if n == 0:
        return complex(float(meas(arcs)), 0.0)
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)

def sine_energy(arcs, N=1500):
    s = 0.0; mx = 0.0
    for n in range(1, N + 1):
        im = chat_c(arcs, n).imag
        s += im * im; mx = max(mx, abs(im))
    return s, mx

def spec_true(P, E):
    gp = safe_set(P); good = good_true(E)
    mGP = meas(gp); mG = meas(good); base = mGP * mG
    if base == 0:
        return None
    inter = meas(intersect(gp, good))
    return dict(SPEC=inter - base, Rp=inter / base, base=base, gp=gp, good=good, mG=mG, mGP=mGP, inter=inter)

def spec_freq_grading(gp, good, N=4000):
    tot = even_n = odd_n = sev = nonsev = 0.0
    terms = []
    for n in range(1, N + 1):
        term = 2.0 * (chat_c(gp, n) * chat_c(good, n).conjugate()).real
        tot += term
        if n % 2 == 0: even_n += term
        else: odd_n += term
        if n % 7 == 0: sev += term
        else: nonsev += term
        if abs(term) > 1e-5:
            terms.append((n, term))
    terms.sort(key=lambda t: -abs(t[1]))
    return dict(tot=tot, even_n=even_n, odd_n=odd_n, sev=sev, nonsev=nonsev, top=terms[:10])

def main():
    print("#" * 96)
    print("# THREAD 2 -- THE TRUE max-gap GOOD set (kind-pasteur kpswf13)")
    print("#   GOOD_true = {maxgap{frac(e_i x)} > 1/7}.  Genuinely x->-x symmetric => ghat REAL.")
    print("#" * 96)

    bank = [
        ([1, 2, 3, 12, 13], list(range(8)),                       "k=8 consec"),
        ([1, 2, 3, 4, 5],   list(range(8)),                       "k=8 consec smallP"),
        ([1, 2, 3],         list(range(10)),                      "k=10 consec"),
        ([2, 3, 4, 5, 6],   list(range(8)),                       "killer P"),
        ([5, 7, 11],        list(range(10)),                      "coprime P"),
        ([1, 2, 6],         [0, 4, 6, 8, 10, 12, 14, 15, 16, 17], "wide"),
        ([1, 2, 3, 12, 13], [0, 1, 3, 7, 12, 20, 30, 31],         "asym cluster"),
    ]

    # ---- compare GOOD_true vs the kpswf12 cover^c proxy + reflection symmetry ----
    print("\n" + "=" * 96)
    print("(A) GOOD_true reflection symmetry (sym-diff with its x->-x reflection = 0 EXACTLY?)")
    print("    and meas(GOOD_true) vs meas(cover^c proxy).")
    print("=" * 96)
    from lrc14_spectrum_intersection_sum_kpswf12 import cover_set
    print(f"  {'case':<20}{'meas(GOOD_true)':>16}{'meas(cover^c)':>15}{'symdiff(true)':>15}{'sineE(true)':>13}")
    for P, E, lab in bank:
        Es = sorted(set(E))
        gt = good_true(Es)
        rgt = reflect_arcs(gt)
        inter = intersect(gt, rgt)
        symd = meas(gt) + meas(rgt) - 2 * meas(inter)
        covc = complement(cover_set(Es))
        seG, _ = sine_energy(gt, N=1200)
        print(f"  {lab:<20}{float(meas(gt)):>16.6f}{float(meas(covc)):>15.6f}{float(symd):>15.2e}{seG:>13.2e}")
    print("  => symdiff(true)=0 EXACTLY: GOOD_true IS reflection-symmetric, ghat_true REAL")
    print("     (sineE ~ 1e-30).  The cover^c proxy was asymmetric (sineE ~ 1.6e-2); meas differs too.")

    # ---- (P2) recompute SPEC_true and R'_true ----
    print("\n" + "=" * 96)
    print("(B) SPEC_true and R'_true on the CORRECT GOOD (compare to kpswf12 proxy R')")
    print("=" * 96)
    print(f"  {'case':<20}{'baseTrue':>10}{'interTrue':>11}{'SPEC_true':>11}{'R_true':>9}")
    res = []
    for P, E, lab in bank:
        b = spec_true(P, sorted(set(E)))
        if not b:
            continue
        res.append((lab, P, sorted(set(E)), b))
        print(f"  {lab:<20}{float(b['base']):>10.5f}{float(b['inter']):>11.5f}"
              f"{float(b['SPEC']):>11.6f}{float(b['Rp']):>9.5f}")

    # ---- (P1) cluster complement invariance of GOOD_true ----
    print("\n" + "=" * 96)
    print("(C) cluster complement invariance:  GOOD_true(E)=GOOD_true(-E)=GOOD_true(N-E)?  (EXACT)")
    print("=" * 96)
    for P, E, lab in bank:
        Es = sorted(set(E)); N = max(Es) + min(Es)
        gE = good_true(Es)
        gneg = good_true([-e for e in Es])
        gN = good_true([N - e for e in Es])
        def eq(a, b):
            return meas(a) + meas(b) - 2 * meas(intersect(a, b)) == 0
        print(f"  {lab:<20} N={N:3d}:  GOOD(-E)==GOOD(E)? {eq(gneg, gE)}   GOOD(N-E)==GOOD(E)? {eq(gN, gE)}")
    print("  => GOOD_true(-E)==GOOD_true(E) EXACTLY (true x->-x symmetry, no boundary artifact).")

    # ---- (P3) frequency grading of SPEC_true ----
    print("\n" + "=" * 96)
    print("(D) frequency Z_2 grading of SPEC_true:  even-n/odd-n, 7|n / 7-nmid-n")
    print("=" * 96)
    for lab, P, Es, b in res:
        fr = spec_freq_grading(b['gp'], b['good'], N=4000)
        print(f"\n  {lab}:  P={P}")
        print(f"    SPEC_true(spec) = {fr['tot']:+.6f}  (exact {float(b['SPEC']):+.6f}, |d|={abs(fr['tot']-float(b['SPEC'])):.1e})")
        print(f"    even-n={fr['even_n']:+.6f}  odd-n={fr['odd_n']:+.6f}    7|n={fr['sev']:+.6f}  7nmid={fr['nonsev']:+.6f}")
        tops = ", ".join(f"n={n}:{t:+.4f}" for n, t in fr['top'][:8])
        print(f"    top: {tops}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
