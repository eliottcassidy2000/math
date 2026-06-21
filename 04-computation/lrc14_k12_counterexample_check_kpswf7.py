#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 2 -- VERIFY the apparent k=12 wall counterexample.

lrc14_topcell_split found: at k=12 (box 13), the shape
   E* = (0,1,2,3,4,5,6,7,8,9,10,12)
has measS7(E*) = 0.6333 > measS7(consec_12) = 0.6176.
i.e. consec is BEATEN. Verify EXACTLY and determine scope:
 (a) Is E* full-residue mod 7?  primitive?  min=0?  span<=14?
 (b) EXACT measS7 of both (rational).
 (c) Is the SAME phenomenon present at k=9..11 with a wider box? Re-scan FULL
     stratum span<=14 (the relevant LRC window: LRC(14) => differences live in
     a 14-window, so span<=13 is the natural cap; span 14 means 0..14 i.e. 15
     residues -- check the exact span cap that matters).
 (d) What is the GLOBAL argmax of measS7 over the full-residue span<=B stratum
     for each k and B? Is it consec, or E*-type (consec with a single top gap)?
 (e) IMPORTANT scope question: does the LRC(14) wall even claim consec-max at
     k=12?  The relevant k for LRC(14) -- print which k the proof needs.
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def W_a(E, a):
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

def measS7(E): return sum(W_a(E,a) for a in range(1,7))
def full_residue(E): return set(e % 7 for e in E) == set(range(7))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return tuple(range(k))
def span(E): return max(E)-min(E)

if __name__=="__main__":
    Es = (0,1,2,3,4,5,6,7,8,9,10,12)
    C = consec(12)
    print("=== (a)(b) k=12 candidate counterexample ===")
    print(f"E* = {Es}")
    print(f"  full-residue mod7? {full_residue(Es)}  primitive? {primitive(Es)}  min={min(Es)} span={span(Es)}")
    print(f"  consec_12 = {C}  span={span(C)}")
    mE=measS7(Es); mC=measS7(C)
    print(f"  measS7(E*)     = {mE} = {float(mE):.8f}")
    print(f"  measS7(consec) = {mC} = {float(mC):.8f}")
    print(f"  E* beats consec? {mE>mC}   delta={float(mE-mC):+.8f}")
    print(f"  residue multiset E*: {sorted(e%7 for e in Es)}")
    print(f"  residue multiset C : {sorted(e%7 for e in C)}")

    print("\n=== (d) GLOBAL argmax of measS7 over full-residue, primitive, min=0, span<=B ===")
    for k in [8,9,10,11,12]:
        for B in [k-1, 13, 14]:
            if B < k-1: continue
            best=None; bestval=F(-1); nbetter=0
            cval=measS7(consec(k))
            cnt=0
            for combo in itertools.combinations(range(1,B+1), k-1):
                E=(0,)+combo
                if not full_residue(E): continue
                if not primitive(E): continue
                cnt+=1
                v=measS7(E)
                if v>bestval: bestval=v; best=E
                if v>cval+F(1,10**9): nbetter+=1
            is_consec = (best==consec(k))
            print(f"  k={k} span<=B={B}: stratum={cnt:5d}  argmax={'CONSEC' if is_consec else best}  max={float(bestval):.6f}  consec={float(cval):.6f}  #rivals>consec={nbetter}")
