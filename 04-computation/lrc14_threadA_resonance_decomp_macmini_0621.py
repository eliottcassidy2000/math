#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- THE RESONANCE DECOMPOSITION OF measS7.

GOAL: test whether the EXACT structure measS7(E) = P2(E) + sum_{p/q} R_{p/q}(E),
with each resonance R now exact, proves consec maximizes measS7 over the
FULL-RESIDUE stratum (k=8, the e_i mod 7 = Z/7 shapes), and -- crucially --
whether the maximization follows from a CLEAN per-resonance / aggregate-Markoff
property of consec's residue structure, NOT a 256-shape finite check.

SETUP.  measS7(E) = meas{ x in [0,1) : {floor(7 frac(e_i x)) : i} = Z/7 }.
The "all 7 sectors covered" event.  measS7 = P(N=0) where N = # missed sectors.

RESONANCE STRUCTURE.  Near a rational x = a/7 (a coprime to 7), each frac(e_i x)
sits near (e_i a mod 7)/7 = (e_i mod 7)*a / 7.  So on a small arc around a/7,
the covered sectors are approximately a*{e_i mod 7} mod 7.  This = Z/7 iff
{e_i mod 7} = Z/7 (full-residue).  These x=a/7 arcs are the DOMINANT resonances
(period-7 / apex-prime resonances).  But measS7 also has cover on arcs where the
*spread* of the e_i (not just residues mod 7) creates full cover -- the
"non-resonant" / background P2(E).

This script:
  (1) Verify the full-residue stratum & cap-stratum separation (HYP-2749) at k=8.
  (2) Build the EXACT measS7 and its decomposition by base-point structure.
  (3) Identify WHICH arcs carry the measS7 cover (resonant a/7 vs background)
      and compute per-shape the resonant vs background split.
  (4) Test consec-max on the full-residue stratum and report whether the
      resonant part alone, the background alone, or only their sum is maximized
      by consec.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
def covered_sectors(E, xm):
    """set of sectors floor(7 frac(e*xm)) over e in E."""
    secs = set()
    for e in E:
        v = e * xm
        v = v - (v.numerator // v.denominator)   # frac in [0,1)
        secs.add((v.numerator * 7) // v.denominator)
    return secs

def measS7_arcs(E):
    """Exact measS7; return total AND list of (lo,hi,xm) covering arcs."""
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1):
            bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0); arcs = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo + hi) / 2
        if len(covered_sectors(E, xm)) == 7:
            total += hi - lo
            arcs.append((lo, hi, xm))
    return total, arcs

def measS7(E):
    return measS7_arcs(E)[0]

def residues(E):
    return frozenset(e % 7 for e in E)

def is_full_residue(E):
    return residues(E) == frozenset(range(7))

def primitive(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1

def consec(k):
    return list(range(k))

# Classify each covering arc by its nearest resonance a/7.
def nearest_resonance_7(xm):
    """Return (a, dist) where a/7 is the nearest multiple of 1/7 to xm."""
    best = None
    for a in range(8):
        d = abs(xm - F(a, 7))
        if best is None or d < best[1]:
            best = (a, d)
    return best

# ---------------------------------------------------------------------------
CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91)}

def split_resonant(E, eps=F(1, 28)):
    """Split measS7 into (resonant near a/7 within eps) + (background)."""
    tot, arcs = measS7_arcs(E)
    res = F(0); bg = F(0)
    res_by_a = defaultdict(lambda: F(0))
    for lo, hi, xm in arcs:
        a, d = nearest_resonance_7(xm)
        if d < eps and a not in (0, 7):  # interior resonance
            res += hi - lo
            res_by_a[a] += hi - lo
        else:
            bg += hi - lo
    return tot, res, bg, dict(res_by_a)

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("#" * 78)
    print("# THREAD A: RESONANCE DECOMPOSITION of measS7 -- consec-max test")
    print("#" * 78)

    k = 8; W = 11
    C = consec(k); mC = measS7(C)
    print(f"\n=== k={k}, consec={C}, measS7(consec)={mC} = {float(mC):.6f}, cap_{k}={CAP[k]}={float(CAP[k]):.6f} ===")
    print(f"    consec is full-residue? {is_full_residue(C)}  (residues {sorted(residues(C))})")

    # (1) Stratum separation: full vs non-full at k=8 span<=W
    bank = [(0,) + r for r in itertools.combinations(range(1, W + 1), k - 1)]
    bank = [E for E in bank if primitive(E)]
    full = [E for E in bank if is_full_residue(E)]
    nonfull = [E for E in bank if not is_full_residue(E)]
    mful = max(measS7(list(E)) for E in full)
    mnon = max(measS7(list(E)) for E in nonfull)
    print(f"\n(1) STRATUM (k={k}, span<= {W}): {len(bank)} shapes; "
          f"full-res {len(full)}, non-full {len(nonfull)}")
    print(f"    max measS7 over FULL-residue   = {float(mful):.6f}")
    print(f"    max measS7 over NON-full       = {float(mnon):.6f}  "
          f"(< cap by {float(CAP[k]-mnon):.4f}; trivially safe: {mnon < CAP[k]})")

    # (2) consec-max over the full-residue stratum exactly
    beaters = [(E, measS7(list(E))) for E in full if measS7(list(E)) > mC + F(1, 10**12)]
    ties = [E for E in full if measS7(list(E)) == mC and tuple(E) != tuple(C)]
    print(f"\n(2) consec-max over FULL-residue stratum (k={k}, span<= {W}):")
    print(f"    consec measS7 = {float(mC):.6f}; beaters = {len(beaters)}; ties (non-consec) = {len(ties)}")
    if beaters:
        for E, v in sorted(beaters, key=lambda t: -t[1])[:5]:
            print(f"      BEATER {list(E)}: {float(v):.6f}")
    if ties:
        for E in ties[:5]:
            print(f"      TIE {list(E)}")

    # (3) resonant / background split for consec and top competitors
    print(f"\n(3) RESONANT (near a/7) vs BACKGROUND split of measS7:")
    tot, res, bg, rba = split_resonant(C)
    print(f"    consec={C}: tot={float(tot):.6f} resonant={float(res):.6f} background={float(bg):.6f}")
    print(f"      resonant-by-a: {{ {', '.join(f'{a}:{float(v):.4f}' for a,v in sorted(rba.items()))} }}")
    # top 5 full-residue shapes by measS7
    ranked = sorted(full, key=lambda E: -measS7(list(E)))[:6]
    for E in ranked:
        if tuple(E) == tuple(C): continue
        tot, res, bg, rba = split_resonant(list(E))
        print(f"    {list(E)}: tot={float(tot):.6f} resonant={float(res):.6f} background={float(bg):.6f}")

    # Does consec maximize the RESONANT part alone? the BACKGROUND alone?
    print(f"\n(4) Does consec maximize RESONANT alone / BACKGROUND alone over full-res?")
    res_C = split_resonant(C)[1]; bg_C = split_resonant(C)[2]
    res_beat = sum(1 for E in full if split_resonant(list(E))[1] > res_C + F(1,10**12))
    bg_beat = sum(1 for E in full if split_resonant(list(E))[2] > bg_C + F(1,10**12))
    print(f"    consec resonant={float(res_C):.6f}: #shapes with HIGHER resonant = {res_beat} (consec-max resonant? {res_beat==0})")
    print(f"    consec background={float(bg_C):.6f}: #shapes with HIGHER background = {bg_beat} (consec-max bg? {bg_beat==0})")
