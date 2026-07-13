#!/usr/bin/env python3
"""
lrc14_rainbow_arcs_klein_S279.py
================================
klein-2026-07-13-S279 (owner: prove the shared Beurling-Selberg / multi-linear estimate).

The density swing-endpoint bound: |U_s^{e'}| <= #cond_s-arcs, cond_s = "the k-1 OTHER offsets cover
exactly {0..6}\{s}". For the k=8 row, k-1=6 others => cond_s = 6 offsets in a RAINBOW (6 distinct
sectors missing s) -- an irreducibly >=6-way (multi-linear) condition.

DECISIVE QUESTION: is #cond_s-arcs (rainbow arcs) BOUNDED as the 6 others SPREAD (diameter grows),
or does it grow with diameter?
 - BOUNDED => density crux is EASY (high-order rarity makes T small, U bounded) => provable.
 - GROWS   => multi-linear HARD (needs Gowers cancellation), same as covering.

Count, for 6 offsets, the maximal x-arcs where they occupy 6 distinct sectors (any 6 of 7), and
the arcs covering exactly {0..6}\s for each s. Vary spread.
"""
import math
from math import gcd
from functools import reduce
NG=2000000  # fine grid to resolve arcs
def secs(offs,x):
    o=0
    for e in offs: o|=1<<(int((e*x%1.0)*7.0)%7)
    return o
def count_arcs(offs):
    """#maximal arcs where offs occupy 6 DISTINCT sectors (rainbow, |occ|==6);
       and per-s #arcs covering exactly {0..6}\s."""
    rainbow=0; in_r=False
    per_s=[0]*7; in_s=[False]*7
    for k in range(1,NG):
        o=secs(offs,k/NG); nf=bin(o).count("1")
        # rainbow = exactly 6 distinct sectors occupied (|occ|==6)
        r=(nf==6)
        if r and not in_r: rainbow+=1
        in_r=r
        for s in range(7):
            cs = (nf==6) and not ((o>>s)&1)  # cover exactly {0..6}\s
            if cs and not in_s[s]: per_s[s]+=1
            in_s[s]=cs
    return rainbow, per_s

print("6 others, growing spread: #rainbow-arcs (|occ|=6) and max_s #cond_s-arcs vs diameter")
print("(BOUNDED => density crux EASY/provable; GROWS => multi-linear hard, like covering)")
print("="*74)
print("  {:32s} {:>5} {:>10} {:>12}".format("6 others","diam","#rainbow","max_s cond_s"))
fams=[
  ("{0,1,2,3,4,5}",   [0,1,2,3,4,5]),
  ("{0,1,2,3,4,10}",  [0,1,2,3,4,10]),
  ("{0,1,2,3,4,20}",  [0,1,2,3,4,20]),
  ("{0,2,4,8,16,32}", [0,2,4,8,16,32]),
  ("{0,1,4,9,16,25}", [0,1,4,9,16,25]),
  ("{0,7,17,31,50,73}",[0,7,17,31,50,73]),
  ("{0,10,23,41,66,99}",[0,10,23,41,66,99]),
  ("{0,20,50,90,140,199}",[0,20,50,90,140,199]),
]
for name,offs in fams:
    rb,per_s=count_arcs(offs)
    print("  {:32s} {:5d} {:10d} {:12d}".format(name,max(offs),rb,max(per_s)))
print("-"*74)
print("  Watch: does #rainbow / max_s cond_s stay bounded (~O(1)) or grow ~diam as offsets spread?")
print("\ndone.")
