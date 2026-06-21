#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 2 -- JOINT coupling, WHOLE-CIRCLE version.

ATTACK 2/3 (conductance majorization & monotonicity) both REFUTED in
lrc14_joint_coupling_kpswf7.py: measS7 is NOT monotone/Schur-convex in the
conductance vector c_r (rivals can have HIGHER total conductance yet LOWER measS7,
because doubling the smallest magnitude makes c_1=9/8>1 but wastes coverage).
=> the conductance vector is the WRONG order parameter for the wall.

NEW ATTACK -- the genuine Huffer-Shepp WHOLE-CIRCLE coupling.

Setup. For a fixed offset y in [0,1), the LRC "all-sectors-hit at y" event is
  HIT(E,y) = { the 7 residues {round(7 * frac(e*y))} over e in E cover all of Z/7 }.
measS7(E) = meas{ y in [0,1) : HIT(E,y) and y not in cell_0 } (the 6 active cells).

The HS theorem (circle coverage by random arcs, Schur-convex in arc lengths) is a
WHOLE-CIRCLE statement. The right analogue: think of HIT as a function of y; the
UNCOVERED set U(E) = {y : not HIT} is a union of arcs. HS-style extremality would
say: the SORTED gap (uncovered-arc) length vector of consec is MAJORIZED by every
rival's => consec has the smallest "spread of holes" => largest covered measure.

We test, exactly:
  (J1) The uncovered-arc LENGTH MULTISET of E (within the 6 active cells).
       Does consec's uncovered-length vector get MAJORIZED by every rival's?
       (uncovered = 1-measS7 over active region; if majorization holds AND the
        relevant functional is Schur-CONVEX in hole-lengths, coverage is maximized.)
  (J2) DIRECT MONOTONE COUPLING on y. Is there a fixed measure-preserving
       bijection T:[0,1)->[0,1) with HIT(E, T(y)) => HIT(consec, y) for all y?
       Sufficient test (necessary condition): the SORTED indicator (decreasing
       rearrangement) of HIT(consec,.) pointwise-dominates that of HIT(E,.).
       Since indicators are 0/1, decreasing rearrangement domination
       <=> meas HIT(consec) >= meas HIT(E). That is just the wall restated -- so
       J2 reduces to the wall. The CONTENT is whether a STRUCTURED (not abstract)
       coupling exists: we test the BINDING-SECTOR coupling.
  (J3) BINDING-SECTOR / max-min coupling. At each y, let bind(E,y)= the residue
       class that is HARDEST to cover (the missing one as y leaves the covered comp).
       Conjecture: consec's binding profile is the "flattest" (max-min optimal).
       Test the per-sector survival function and whether consec's vector of
       per-sector marginal coverages is majorized.

OUTPUT: which coupling (if any) gives a clean majorization certifying consec-max.
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def frac(x):  # fractional part of a Fraction
    return x - (x.numerator // x.denominator)

def sector(x):  # which of 7 sectors does frac in [0,1) land in: floor(7*frac)
    f = frac(x)
    return (f.numerator * 7) // f.denominator

def hit_at(E, y):
    """True if {sector(e*y): e in E} == all of Z/7."""
    s = set()
    for e in E:
        s.add(sector(e * y))
        if len(s) == 7: return True
    return len(s) == 7

def breakpoints_full(E):
    """All breakpoints of HIT(E,.) on [0,1): where any sector(e*y) jumps,
       i.e. e*y = j/7 => y = j/(7|e|)."""
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        d = 7 * abs(e)
        for j in range(0, d+1):
            bps.add(F(j, d))
    return sorted(bps)

def covered_intervals(E):
    """Return list of (lo,hi) maximal intervals on [0,1) where HIT(E,.) is True."""
    bps = breakpoints_full(E)
    ivs = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo + hi) / 2
        if hit_at(E, xm):
            ivs.append((lo, hi))
    # merge adjacent
    merged = []
    for lo, hi in ivs:
        if merged and merged[-1][1] == lo:
            merged[-1] = (merged[-1][0], hi)
        else:
            merged.append((lo, hi))
    return merged

def cell_bounds(a):
    return (F(a,7)-F(1,14), F(a,7)+F(1,14))

def measS7_active(E):
    """Coverage within the 6 active cells a=1..6 (excluding cell_0)."""
    ivs = covered_intervals(E)
    tot = F(0)
    for a in range(1,7):
        lo, hi = cell_bounds(a)
        for (l,h) in ivs:
            ll = max(l, lo); hh = min(h, hi)
            if hh > ll: tot += hh - ll
    return tot

def uncovered_lengths_active(E):
    """Sorted (desc) multiset of uncovered-arc lengths WITHIN the active region
       [1/7-1/14, 6/7+1/14] (the union of the 6 active cells, which is contiguous
       except for the gaps at cell boundaries -- but cells a=1..6 are contiguous:
       cell1=[1/14,3/14], cell2=[3/14,5/14],... cell6=[11/14,13/14]. They TILE
       [1/14, 13/14]!). So active region = [1/14, 13/14], length 12/14=6/7."""
    lo_act, hi_act = F(1,14), F(13,14)
    ivs = covered_intervals(E)
    # restrict covered to active region
    cov = []
    for (l,h) in ivs:
        ll = max(l, lo_act); hh = min(h, hi_act)
        if hh > ll: cov.append((ll,hh))
    cov.sort()
    # complement within [lo_act, hi_act]
    holes = []
    cur = lo_act
    for (l,h) in cov:
        if l > cur: holes.append(l - cur)
        cur = max(cur, h)
    if cur < hi_act: holes.append(hi_act - cur)
    return sorted(holes, reverse=True)

def full_residue(E): return set(e % 7 for e in E) == set(range(7))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return tuple(range(k))

def majorizes(u, v):
    su = sorted(u, reverse=True); sv = sorted(v, reverse=True)
    n = max(len(su), len(sv))
    su += [F(0)]*(n-len(su)); sv += [F(0)]*(n-len(sv))
    eq = sum(su) == sum(sv)
    pu = pv = F(0); ok = True
    for i in range(n):
        pu += su[i]; pv += sv[i]
        if pu < pv: ok = False
    return ok, eq

def stratum(k, box):
    out=[]
    for combo in itertools.combinations(range(1,box+1), k-1):
        E=(0,)+combo
        if full_residue(E) and primitive(E): out.append(E)
    return out

if __name__ == "__main__":
    print("="*88)
    print("WHOLE-CIRCLE coupling: active region [1/14,13/14] tiled by 6 cells")
    print("="*88)
    # sanity: active region length and consec coverage match measS7 from route3
    for k, box in [(8,14),(9,14),(10,13)]:
        C = consec(k)
        mA = measS7_active(C)
        print(f"k={k}: measS7_active(consec)={mA}={float(mA):.6f}")

    for k, box in [(8,14),(9,14),(10,13)]:
        print(f"\n{'#'*84}\n# k={k}, box=[0,{box}]\n{'#'*84}")
        C = consec(k)
        hC = uncovered_lengths_active(C)   # consec's holes (desc)
        mC = measS7_active(C)
        print(f"consec measS7_active={float(mC):.6f}, total hole length={float(sum(hC)):.6f}")
        print(f"  consec uncovered-arc lengths (desc, floats) = {[round(float(x),4) for x in hC]}")

        S = stratum(k, box)
        print(f"  stratum size = {len(S)}")

        # (J1) does EVERY rival's hole-vector MAJORIZE consec's?
        # (If holes majorize => more concentrated holes; but total hole = 1-meas differs.
        #  We need: meas(consec)>=meas(E) <=> sum holes(consec) <= sum holes(E).
        #  Schur route only clean if total hole is FIXED -- check that first.)
        totalholes = set(sum(uncovered_lengths_active(E)) for E in S)
        print(f"  distinct total-hole lengths over stratum = {len(totalholes)} (total hole NOT fixed => pure Schur on holes can't apply directly)")

        # The RIGHT test: is consec's total hole the MINIMUM? (= the wall)
        minhole = min(sum(uncovered_lengths_active(E)) for E in S)
        print(f"  min total hole over stratum = {float(minhole):.6f}; consec total hole = {float(sum(hC)):.6f}; consec is min? {sum(hC)==minhole}")

        # (J3) per-SECTOR marginal coverage vector. For each residue r in 1..6,
        #   m_r(E) = meas{y in active : sector NOT-missing is r-related}...
        # Cleaner: per-cell coverage vector W_a, and test if consec's W-vector,
        # SORTED, is majorized by rivals (Schur on the 6 cell coverages).
        Wc = [F(0)]*6
        for a in range(1,7):
            lo,hi = cell_bounds(a);
            tot=F(0)
            for (l,h) in covered_intervals(C):
                ll=max(l,lo); hh=min(h,hi)
                if hh>ll: tot+=hh-ll
            Wc[a-1]=tot
        print(f"  consec W-cell vector (sorted desc) = {[round(float(x),4) for x in sorted(Wc,reverse=True)]}")
        # is consec's W-vector majorized by every rival? (Schur-convex sum => consec MIN sum -- wrong direction)
        wmaj_into=0; wmaj_from=0
        for E in S:
            if E==C: continue
            WE=[F(0)]*6
            for a in range(1,7):
                lo,hi=cell_bounds(a); tot=F(0)
                for (l,h) in covered_intervals(E):
                    ll=max(l,lo); hh=min(h,hi)
                    if hh>ll: tot+=hh-ll
                WE[a-1]=tot
            mi,_=majorizes(Wc,WE)   # consec majorizes rival?
            mo,_=majorizes(WE,Wc)   # rival majorizes consec?
            if mi: wmaj_into+=1
            if mo: wmaj_from+=1
        print(f"  consec W-vec majorizes rival: {wmaj_into}/{len(S)-1}; rival majorizes consec: {wmaj_from}/{len(S)-1}")
