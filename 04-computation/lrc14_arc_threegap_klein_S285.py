#!/usr/bin/env python3
"""
lrc14_arc_threegap_klein_S285.py
================================
klein-2026-07-13-S285 (owner: attack the density Weyl bound). NEW ANGLE: three-gap / arithmetic
structure of the R_s arcs. Each offset e' puts "in sector s" = e' equal arcs (width 1/(7e'), period
1/e') -- a perfect AP. R_s = intersection/union of these. Do the R_s arc WIDTHS take FEW distinct
values (three-gap-like => exploitable structured cancellation of Sum_i w_i e(-ell w c_i)), or MANY
(generic => genuine equidistribution, needs external input)?

Also: are the arc MIDPOINTS c_i arithmetically structured (few residue classes under some modulus)?
"""
import math
from collections import Counter
NG=3000000
def sec(e,x): return int((e*x%1.0)*7.0)%7
def occ(E,x):
    o=0
    for e in E: o|=1<<sec(e,x)
    return o
def Rs_arcs(E,s):
    arcs=[]; inR=False; st=0.0
    for k in range(1,NG):
        x=k/NG; o=occ(E,x); cur=(7-bin(o).count("1")==1) and not((o>>s)&1)
        if cur and not inR: st=(k-0.5)/NG; inR=True
        if (not cur) and inR: arcs.append((st,(k-0.5)/NG)); inR=False
    if inR: arcs.append((st,1.0))
    return arcs

print("three-gap check: # DISTINCT R_s-arc widths (rounded) vs # arcs; are widths few?")
print("="*72)
print("  {:28s} {:>2} {:>5} {:>14} {:>14}".format("E'","s","#arcs","#distinct-w","top-3 widths x1e4"))
for E in [[0,1,2,3,4,5,6],[0,3,7,15,30,55,90],[0,10,27,55,99,150,199]]:
    for s in [1,3]:
        arcs=Rs_arcs(E,s); r=len(arcs)
        if r<2: continue
        wid=[round((b-a)*1e6) for a,b in arcs]  # widths in units of 1e-6
        c=Counter(wid); nd=len(c)
        top=sorted(c.items(),key=lambda kv:-kv[1])[:3]
        tops=" ".join("{:.1f}x{}".format(w/100,cnt) for w,cnt in top)  # w/100 => 1e-4 units
        print("  {:28s} {:>2} {:>5} {:>14} {:>14}".format(str(E),s,r,nd,tops))
print("-"*72)
print("  If #distinct-w is SMALL (three-gap-like): the width-weighted Weyl sum groups into few classes,")
print("  each an AP of midpoints => structured cancellation (NEW, beyond large-sieve). If MANY: generic.")
print()
# midpoint structure: for a fixed cluster, check if midpoints c_i are near multiples of some 1/Q
E=[0,3,7,15,30,55,90]; arcs=Rs_arcs(E,1)
mids=sorted(((a+b)/2) for a,b in arcs)
gaps=sorted(round((mids[i+1]-mids[i])*1e6) for i in range(len(mids)-1))
gc=Counter(gaps)
print("  midpoint GAPS for {} s=1: #distinct-gap={}, top-3: {}".format(E,len(gc),
      sorted(gc.items(),key=lambda kv:-kv[1])[:3]))
print("  (three-gap theorem: an AP of points has <=3 distinct gaps. Few gaps => midpoints structured.)")
print("\ndone.")
