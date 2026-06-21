#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- THE PER-RESONANCE EXACT WIDTH near x=a/7.

Finding from lrc14_threadA_resonance_decomp: consec MAXIMIZES the resonant part
of measS7 (the cover near x=a/7, a=1..6) with 0 beaters over the full-residue
stratum, while it does NOT maximize the background. So the action is in the
resonances. This script computes the EXACT per-resonance width and proves the
local structure of each resonance arc.

LOCAL MODEL near x = a/7 (a coprime to 7).  Write x = a/7 + y, y small.
frac(e x) = frac(e a / 7 + e y) = (r_e + 7 e y)/7 + ... where r_e = e a mod 7.
The sector floor(7 frac(e x)) = r_e  EXACTLY as long as 7 e y stays within the
sector, i.e. the offset 7 e y has not crossed into the next sector. The covered
set is {r_e} = a*{e mod 7} mod 7 = Z/7 (since full-residue) AT y=0.

As y moves away from 0, sectors start to FAIL (an e crosses a boundary) one at a
time. The resonance WIDTH = the y-interval on which ALL 7 sectors remain covered.

This is EXACTLY a covering-by-arcs problem: each e in E covers a sub-interval of
y around 0, and measS7-resonance-a = length of the intersection-of-coverage where
the UNION of sectors = Z/7. We compute it exactly and relate it to the residues.

KEY QUESTION: is the per-resonance width a CLEAN function of the residue multiset
(e mod 7) and the "spread" (which e realizes each residue with what magnitude)?
And does consec (residues {0..6}, smallest magnitudes 0..7) MAXIMIZE each
per-resonance width? If YES per-resonance, then consec-max-resonant is a clean
local theorem (no finite check).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def covered_sectors(E, xm):
    secs = set()
    for e in E:
        v = e * xm; v = v - (v.numerator // v.denominator)
        secs.add((v.numerator * 7) // v.denominator)
    return secs

def measS7_arcs(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0); arcs = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo + hi) / 2
        if len(covered_sectors(E, xm)) == 7:
            total += hi - lo; arcs.append((lo, hi))
    return total, arcs

def resonance_width(E, a):
    """Exact total cover length on the arc strictly between a/7-1/14 and a/7+1/14
       (the half-cell around a/7 where the a-resonance lives)."""
    tot, arcs = measS7_arcs(E)
    lo_b, hi_b = F(2*a-1, 14), F(2*a+1, 14)   # the cell of width 1/7 centered at a/7
    w = F(0)
    for lo, hi in arcs:
        l = max(lo, lo_b); h = min(hi, hi_b)
        if h > l: w += h - l
    return w

def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

if __name__ == "__main__":
    print("#"*78); print("# PER-RESONANCE EXACT WIDTH near x=a/7  (THREAD A)"); print("#"*78)
    k = 8; W = 11
    C = consec(k)
    print(f"\nconsec={C}, per-resonance widths (a=1..6):")
    for a in range(1, 7):
        print(f"   a={a}: width={resonance_width(C,a)} = {float(resonance_width(C,a)):.6f}")
    tot = sum(resonance_width(C, a) for a in range(1, 7))
    print(f"   sum a=1..6: {tot} = {float(tot):.6f}")

    # Per-resonance consec-max over full-residue stratum
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
    full = [E for E in bank if primitive(E) and is_full_residue(E)]
    print(f"\nFull-residue stratum: {len(full)} shapes. Per-resonance consec-max test:")
    for a in range(1, 7):
        wC = resonance_width(C, a)
        beat = [(E, resonance_width(list(E), a)) for E in full
                if resonance_width(list(E), a) > wC + F(1, 10**12)]
        bestE = max(full, key=lambda E: resonance_width(list(E), a))
        print(f"   a={a}: consec width={float(wC):.6f}  beaters={len(beat)}  "
              f"argmax={'consec' if tuple(bestE)==tuple(C) else list(bestE)}"
              + (f" (={float(resonance_width(list(bestE),a)):.6f})" if tuple(bestE)!=tuple(C) else ""))

    # AGGREGATE resonant consec-max
    def total_res(E): return sum(resonance_width(E, a) for a in range(1, 7))
    tC = total_res(C)
    tbeat = [(E, total_res(list(E))) for E in full if total_res(list(E)) > tC + F(1,10**12)]
    print(f"\nAGGREGATE resonant (sum a=1..6): consec={float(tC):.6f}  beaters={len(tbeat)}")
    if tbeat:
        for E,v in sorted(tbeat, key=lambda t:-t[1])[:5]: print(f"   BEATER {list(E)}: {float(v):.6f}")

    # Per-resonance structure: width depends on which e realizes each residue.
    # The resonance at a/7 only "sees" residues r_e = e*a mod 7 and the LOCAL
    # slopes 7*e. Print, for consec, the (e, e mod 7, e*a mod 7) table to see the
    # residue->magnitude assignment that consec uses.
    print(f"\nConsec residue assignment (e, e mod 7) and the local slope 7e:")
    for e in C:
        print(f"   e={e}: e mod 7={e%7}, local slope 7e={7*e}")
