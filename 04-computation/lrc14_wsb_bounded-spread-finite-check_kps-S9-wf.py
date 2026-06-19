#!/usr/bin/env python3
"""
lrc14_wsb_bounded-spread-finite-check_kps-S9-wf.py   (kind-pasteur-2026-06-19-S9-wf)

ANGLE: bounded-spread finite check, made RIGOROUS, dovetailing with the wide-spread bound.

THE OPEN STATEMENT (k=8,9,10 the dangerous rows):
    meas(S7(E)) <= cap_k   for ALL primitive E = {0=e_1<...<e_k}, |E|=k.
    cap_8 = 2243/5880, cap_9 = 1979/4004, cap_10 = 55/91.
    meas(S7(E)) = meas{ x in [0,1) : {floor(7 e_i x) mod 7 : i} = Z/7 }  (every sector hit).

WHY METRIC SPREAD IS THE WRONG NOTION (the subtlety the prompt flags):
    A shape {0,1,...,6, BIG} has arbitrarily large metric span yet behaves like the small
    consecutive cluster {0..6} plus an isolated far point.  So a finite check over "span<=B"
    is NOT exhaustive: the residual span>B is infinite.

THE RIGOROUS FINITE REDUCTION (this script's contribution):
    We combine three PROVED tools to make the bounded part a GENUINELY FINITE exact check.

    (L0) SCALE INVARIANCE (PROVED, canon THM-532/HYP-2606): meas(S7(dE)) = meas(S7(E)).
         => WLOG gcd(E)=1.

    (L1) MONOTONICITY UNDER ADDING POINTS (PROVED here, trivial): if E subset E' then
         S7(E) subset S7(E') pointwise (more indices can only hit more sectors), so
         meas(S7(E)) <= meas(S7(E')).  [used for the domination direction]

    (L2) SUBSET DOMINATION = THM-536-B2 (PROVED): if E subset {0,1,...,N} then
         meas(S7(E)) <= meas(S7({0,1,...,N})) = meas(S7(AP_{N+1})).

    These give the certificate for SMALL-SPAN E (span <= N*(k), where
    meas(S7(AP_{N*+1})) <= cap_k).  This is the EASY part.

    The HARD part is large span.  We make it finite by the GAP-STRUCTURE characterization,
    NOT metric span.  Concretely we prove (verify exactly) the GAP-CAPPING LEMMA:

    (L3) GAP-CAPPING LEMMA (the new finite reduction).  Sort E={0=e_1<...<e_k}.  Let the
         consecutive gaps be g_i = e_{i+1}-e_i (i=1..k-1).  CLAIM (verified exactly below):
            meas(S7(E)) is determined, AND maximized at fixed gap-MULTISET, by the
            arrangement; and CAPPING each gap g_i at 7 does NOT decrease meas(S7).
         More precisely we test the operational statement we actually need:
            Replacing any gap g_i >= 7 by g_i' = g_i - 7 (i.e. pulling the top block in by 7)
            leaves meas(S7) UNCHANGED when that gap is "saturated" -- equivalently a gap of
            size >=7 contributes like a gap of size (g_i mod 7) shifted, because residues
            floor(7 e x) mod 7 only see e mod (denominators).  We TEST whether
            meas(S7) depends only on the gaps reduced mod 7 (NO -- must verify), and if not,
            we find the true finite invariant.

    We DO NOT assume L3; we TEST it and report the exact truth.  Whatever the true finite
    invariant is, we then enumerate its (finite) class set and run the exact cap check.

DELIVERABLE: an exact, exhaustive certificate for the bounded part for k=8,9,10, over the
CORRECT finite family, with margins; plus confirmation consec is the max on that family.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------------------------------
# EXACT meas(S7).  Breakpoints of x->floor(7 e x) mod 7 are at x=m/(7e).  Exact rational.
# ----------------------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(int(e) for e in E))
    Enz = [e for e in E if e != 0]
    bps = set([F(0), F(1)])
    for e in Enz:
        # breakpoints where floor(7 e x) jumps: x = m/(7e), m=0..7e
        for m in range(0, 7*e + 1):
            bps.add(F(m, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        res = set(int(7*e*xm) % 7 for e in E)   # e=0 gives residue 0
        if len(res) == 7:
            total += (x1 - x0)
    return total

def measS7_AP(m):
    return measS7(tuple(range(m)))

# canon cap_k = min_{|P|=13-k} meas(G_P)
cap = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91),
       11: F(66,91), 12: F(6,7), 13: F(1)}

def primitive(E):
    g = 0
    for e in E:
        g = gcd(g, e)
    return g == 1

print("="*96)
print("STEP 0: anchor values (reproduce canon).")
print("="*96)
apv = {m: measS7_AP(m) for m in range(1, 16)}
print("  meas(S7(AP_m)), m=1..15:")
for m in range(1,16):
    print(f"    AP_{m:2d}: {apv[m]} = {float(apv[m]):.6f}")
print()
print("  consec_k anchors (should match canon meas(S7)):")
for k in [8,9,10,11,12,13]:
    v = measS7_AP(k)
    print(f"    k={k}: meas(S7(consec))={v} = {float(v):.6f}   cap_k={cap[k]}={float(cap[k]):.6f}   "
          f"<=cap? {v<=cap[k]}  slack={cap[k]-v}={float(cap[k]-v):+.6f}")
print()

print("="*96)
print("STEP 1: N*(k) -- max span certified by subset-domination (THM-536-B2).")
print("  E subset {0..N}, |E|=k => meas(S7(E)) <= meas(S7(AP_{N+1})).  Certified iff <=cap_k.")
print("="*96)
Nstar = {}
for k in sorted(cap):
    ck = cap[k]
    best = None
    for N in range(k-1, 40):
        m = N+1
        v = apv.get(m)
        if v is None:
            v = measS7_AP(m); apv[m] = v
        if v <= ck:
            best = N
        else:
            break
    Nstar[k] = best
    print(f"  k={k}: cap_k={float(ck):.5f}  N*={best}  "
          f"(meas(S7(AP_{best+1}))={float(apv[best+1]):.5f} <= cap_k; "
          f"meas(S7(AP_{best+2}))={float(apv.get(best+2, measS7_AP(best+2))):.5f} > cap_k)")
print()
print("  => Subset-domination CERTIFIES every primitive E with span <= N*(k).  RESIDUAL: span > N*(k).")
print()
