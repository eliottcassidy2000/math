#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 3 (opus, 2026-06-21): port the HUFFER-SHEPP REFLECTION LEMMA to the
PER-CELL survival width W_a of the LRC(14) layer-3 wall.

GOAL.  measS7 = sum_{a=1..6} W_a.  Each W_a is a UNION-COVERAGE measure inside the
cell  I_a = [a/7 - 1/14, a/7 + 1/14]  (width 1/7).  Parametrise x = a/7 + y,
y in J = [-1/14, 1/14].  Clock e sits at residue  floor(7 frac(e x))  which at
y=0 equals  s_e := (e a) mod 7  and DRIFTS at speed 7e as y moves.  W_a = meas{
y in J : every sector 0..6 is occupied by some clock }.

HUFFER-SHEPP (1987) reflection lemma (reconstructed structure):  coverage
probability of the circle by arcs is Schur-CONVEX; the engine is a PAIRWISE
TRANSFER plus a CENTERING / REFLECTION sub-lemma -- "a fixed-length interval is
more likely covered the closer (its covering structure is) to the centre".  The
reflection operation reflects ONE arc's position about a midpoint; coverage of a
fixed target interval does not decrease when the covering arc is reflected toward
covering the centre.

THE PER-CELL PORT.  Inside I_a the natural "reflection about the cell midpoint
y=0" is  y -> -y, i.e. e -> -e  (reverse the drift direction of a clock).  We ask
the SHARP question:

  (R) PER-CELL REFLECTION MONOTONICITY.  Does reflecting a clock e -> -e in a
      single cell a (keeping y=0 phases as they are, only reversing drift) leave
      W_a unchanged or move it toward the consec value?  And does the FULL
      reflection (all clocks e -> -e, i.e. E -> -E) preserve each W_a?

  (S) DRIFT-CENTERING.  consec = {0..k-1}.  The clock e=0 is a PERFECT SHORT on
      sector 0 (HYP-2761): it never drifts, permanently covering s_0 = 0.  The
      claim to port: in consec each binding sector is covered by a clock whose
      drift is CENTERED on the cell (max slack on both sides), so reflection is a
      symmetry of W_a; a non-consec shape breaks this centering and W_a drops.

  (T) THE ACTUAL HS TRANSFER on W_a.  HS transfer = spread two covering-arc
      lengths.  The covering-arc half-width of clock e in I_a is 1/(7|e|) (time to
      cross one sector).  A "reflection toward centre" should EQUALISE the binding
      slacks.  Test whether equalising the binding slacks (centering) raises W_a.

DELIVERABLE: whether the reflection lemma ports (a per-cell monotonicity of W_a),
with EXACT-fraction tests, or the precise obstruction.

stdlib only; exact Fractions.
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ---------------- exact W_a (verified tool, reused) ----------------
def W_a(E, a):
    """Exact survival width of phase a: measure of full-coverage in cell I_a."""
    E = sorted(set(E))
    lo = F(a, 7) - F(1, 14); hi = F(a, 7) + F(1, 14)
    bps = {lo, hi}
    for e in E:
        if e == 0: continue
        d = 7 * abs(e)
        j0 = math.floor(lo * d); j1 = math.ceil(hi * d)
        for j in range(j0 - 1, j1 + 2):
            x = F(j, d)
            if lo <= x <= hi: bps.add(x)
    bps = sorted(bps); tot = F(0)
    for l, h in zip(bps, bps[1:]):
        if h <= l: continue
        xm = (l + h) / 2; hit = set()
        for e in E:
            v = e * xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if len(hit) == 7: tot += h - l
    return tot

def Wvec(E):  return [W_a(E, a) for a in range(1, 7)]
def measS7(E): return sum(W_a(E, a) for a in range(1, 7))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def majorizes(x, y):
    xs = sorted(x, reverse=True); ys = sorted(y, reverse=True)
    if len(xs) != len(ys) or sum(xs) != sum(ys): return False
    px = py = F(0)
    for i in range(len(xs) - 1):
        px += xs[i]; py += ys[i]
        if px < py: return False
    return True

# ---------------- reflection operations ----------------
def reflect_all(E):
    """Full reflection y -> -y on every clock: E -> -E.  (Drift directions reverse.)"""
    return sorted(-e for e in E)

def reflect_one(E, idx):
    """Reflect a single clock's drift: e_idx -> -e_idx."""
    E2 = list(E); E2[idx] = -E2[idx]; return E2

if __name__ == "__main__":
    print("=" * 80)
    print("ROUTE 3: PORT THE HUFFER-SHEPP REFLECTION LEMMA TO PER-CELL W_a")
    print("=" * 80)

    # =====================================================================
    print("\n[R0] BASELINE: per-cell W_a of consec_8 and the e=0 perfect short.")
    print("=" * 80)
    C = consec(8)
    wc = Wvec(C)
    print(f"  consec_8 = {C}")
    print(f"  W_a (a=1..6) = {[str(w) for w in wc]}")
    print(f"  sum = measS7 = {sum(wc)} = {float(sum(wc)):.6f}")
    # check the dilation symmetry: Z/7* should permute the W_a.  For consec the
    # speeds are {0..7}; multiply by c mod 7 permutes residues -> permutes a.
    print("\n  Z/7* phase-permutation check on consec (a -> c*a mod 7 permutes W_a?):")
    for c in range(1, 7):
        perm = [((c * a) % 7) for a in range(1, 7)]
        # W_a should equal W_{c*a}?  Only if the WHOLE shape is dilation-fixed.
    print(f"  (consec W vector sorted = {sorted(str(w) for w in wc)})")

    # =====================================================================
    print("\n" + "=" * 80)
    print("[R1] FULL REFLECTION E -> -E.  Does each W_a map to W_{a'} (a symmetry)?")
    print("     Port of HS: reflecting all covering arcs about the centre is a")
    print("     measure-preserving symmetry of coverage.  Predict: measS7 INVARIANT,")
    print("     and W_a(-E) = W_{-a mod 7}(E)  (reflection sends cell a -> cell -a).")
    print("=" * 80)
    test_shapes = [consec(8), [0,1,2,3,4,5,6,8], [0,2,3,4,5,6,7,8],
                   [0,1,2,3,4,5,6,9], [0,1,2,4,5,6,7,8], consec(9), consec(10)]
    refl_inv = True
    cell_refl_ok = True
    for E in test_shapes:
        Er = reflect_all(E)
        wE = Wvec(E); wR = Wvec(Er)
        same_meas = (sum(wE) == sum(wR))
        # cell-level reflection: W_a(E) =?= W_{(7-a)%7}(Er) ... but cell index a in 1..6
        # reflection x->-x sends a/7 -> -a/7 = (7-a)/7, so cell a <-> cell 7-a.
        cell_ok = all(W_a(E, a) == W_a(Er, (7 - a) % 7) for a in range(1, 7))
        if not same_meas: refl_inv = False
        if not cell_ok: cell_refl_ok = False
        print(f"  E={E}")
        print(f"     measS7(E)={float(sum(wE)):.6f}  measS7(-E)={float(sum(wR)):.6f}  "
              f"invariant={same_meas}")
        print(f"     cell reflection W_a(E)==W_(7-a)(-E): {cell_ok}")
    print(f"\n  => measS7 reflection-invariant on all tested: {refl_inv}")
    print(f"  => cell-wise W_a(E)=W_(7-a)(-E) (reflection symmetry of cells): {cell_refl_ok}")

    # =====================================================================
    print("\n" + "=" * 80)
    print("[R2] SINGLE-CLOCK REFLECTION e_idx -> -e_idx.  THE per-cell HS step.")
    print("     HS centering lemma ported: reflecting ONE clock's drift toward")
    print("     covering the cell centre should NOT DECREASE W_a (move toward consec).")
    print("     Test on consec and near-consec shapes: sign of W_a after reflecting")
    print("     each clock, per cell.")
    print("=" * 80)
    for E in [[0,1,2,3,4,5,6,8], [0,2,3,4,5,6,7,8], [0,1,2,3,4,5,7,8]]:
        print(f"\n  base E={E}  measS7={float(measS7(E)):.6f}")
        for idx in range(len(E)):
            if E[idx] == 0: continue   # reflecting 0 is a no-op
            E2 = reflect_one(E, idx)
            if not primitive(E2): continue
            d = measS7(E2) - measS7(E)
            tag = "+" if d > 0 else ("-" if d < 0 else "0")
            print(f"     reflect e_{idx}={E[idx]:>3} -> {E2[idx]:>3}: dmeasS7={float(d):+.6f} [{tag}]  E'={E2}")

    # =====================================================================
    print("\n" + "=" * 80)
    print("[R3] THE KEY PER-CELL CLAIM.  Inside cell a, define the BINDING SLACKS:")
    print("     for each sector s, the (signed) y-distance to losing its cover.")
    print("     consec should have all binding slacks EQUAL (centered/balanced) =>")
    print("     W_a is a symmetric function of the slacks maximised at balance.")
    print("     Compute the slack profile per cell for consec vs a rival; test")
    print("     whether centering (equalising slacks) is what consec achieves.")
    print("=" * 80)

    def slack_profile(E, a):
        """At cell midpoint y=0 (x=a/7), for each sector s currently occupied,
        find the y-window [y_minus(s), y_plus(s)] during which s stays occupied
        by SOME clock.  Return per-sector the half-window to the left and right of 0.
        The cell is covered on the intersection of all sector-occupancy windows
        that matter; W_a = length of the y-interval where ALL 7 sectors occupied."""
        E = sorted(set(E))
        # occupancy(y, s) = True if some clock e has floor(7 frac(e(a/7+y))) == s
        # The full-coverage y-set = {y in J : union of clock sectors = Z/7}.
        # We want, around y=0, the maximal symmetric-ish interval. Compute the
        # connected component of full-coverage containing y=0 (if y=0 is covered),
        # and its endpoints; report (left reach, right reach).
        lo = F(-1, 14); hi = F(1, 14)
        # breakpoints in y
        bps = {lo, hi, F(0)}
        for e in E:
            if e == 0: continue
            d = 7 * abs(e)
            # clock e sector boundary at e*(a/7+y) = m/7 => y = (m/(7) - e*a/7)/e ... use exact
            # e*x = m/7 ; x = a/7 + y ; y = m/(7e) - a/7
            # iterate m so that y in [lo,hi]
            ylo = lo; yhi = hi
            mlo = math.floor((F(a,7)+ylo)*e*7) - 2
            mhi = math.ceil((F(a,7)+yhi)*e*7) + 2
            for m in range(mlo, mhi+1):
                y = F(m, 7*e) - F(a, 7)
                if lo <= y <= hi: bps.add(y)
        bps = sorted(bps)
        # covered intervals
        covered = []
        for l, h in zip(bps, bps[1:]):
            if h <= l: continue
            ym = (l + h) / 2; hit = set()
            for e in E:
                v = e * (F(a,7)+ym); v = v - (v.numerator // v.denominator)
                hit.add((v.numerator * 7) // v.denominator)
            if len(hit) == 7: covered.append((l, h))
        # find component containing 0
        comp = None
        for (l, h) in covered:
            if l <= 0 <= h: comp = (l, h); break
        return comp, covered

    for E, lbl in [(consec(8), "consec_8"), ([0,1,2,3,4,5,6,8], "rival +8"),
                   ([0,2,3,4,5,6,7,8], "rival shift")]:
        print(f"\n  {lbl}: E={E}")
        for a in range(1, 7):
            comp, cov = slack_profile(E, a)
            w = W_a(E, a)
            if comp is None:
                print(f"     a={a}: y=0 NOT covered; W_a={float(w):.5f}; covered comps={len(cov)}")
            else:
                lft = -comp[0]; rgt = comp[1]
                sym = "SYM" if lft == rgt else "asym"
                print(f"     a={a}: comp around 0 = [{float(comp[0]):+.5f},{float(comp[1]):+.5f}] "
                      f"left={float(lft):.5f} right={float(rgt):.5f} [{sym}] "
                      f"complen={float(comp[1]-comp[0]):.5f} W_a={float(w):.5f}")
