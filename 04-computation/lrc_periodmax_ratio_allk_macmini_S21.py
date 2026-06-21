#!/usr/bin/env python3
"""
lrc_periodmax_ratio_allk_macmini_S21.py  (mac-mini-2026-06-21-S21)

THM-563 general closure: max(period-max(B)/margin(B)) over DANGEROUS bounded bases
B subset [0,14], 0 in B, k=|B|+1 = 9,10,11,12.  If max ratio < 15 at every k, the
single-far case is CLOSED window-free for all w>=15.  Reuses the VALIDATED sawtooth
functions (Aj_arcs/Dw_w/measS7/p1f) from lrc_period_max_macmini_0621s6.py.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

# --- exec the validated function definitions from the original script (everything before its main loop) ---
src = open('04-computation/lrc_period_max_macmini_0621s6.py').read()
prefix = src.split('\nfor k in (8,9,10):')[0]
ns = {}
exec(prefix, ns)
Aj_arcs, Dw_w, measS7, p1f = ns['Aj_arcs'], ns['Dw_w'], ns['measS7'], ns['p1f']
def lcm(a,b): return a*b//gcd(a,b)
def prim(E): return reduce(gcd,[e for e in E if e],0)==1

# --- validate against THM-563 binding table ---
print("VALIDATION vs THM-563 table (consec bases):")
for k,exp in [(8,F(1)),(9,F(43,49)),(10,F(1007,980))]:
    B=list(range(k-1)); P=7*reduce(lcm,[e for e in B if e],1); arcs=Aj_arcs(B)
    mx=max(Dw_w(arcs,w) for w in range(15,15+P))
    print(f"  k={k} consec: period-max={mx}={float(mx):.5f} expect {exp}  {'OK' if mx==exp else 'MISMATCH'}")

caps={9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}
PMAX=80000
print("\nGENERAL CHECK -- max period-max/margin over dangerous (margin<0.20) bounded bases:")
for k in [9,10,11,12]:
    worst=F(0); worstB=None; checked=0; skipped=0; danger=0
    for combo in itertools.combinations(range(1,15),k-2):
        B=[0]+list(combo)
        if not prim(B): continue
        plat=measS7(B)+p1f(B)*F(1,7); margin=caps[k]-plat
        if margin<=0:
            print(f"  !! margin<=0 at {B}"); continue
        if margin>=F(1,5): continue
        danger+=1
        P=7*reduce(lcm,[e for e in B if e],1)
        if P>PMAX: skipped+=1; continue
        checked+=1; arcs=Aj_arcs(B); mx=F(0); mw=None
        for w in range(15,15+P):
            v=Dw_w(arcs,w)
            if v>mx: mx=v; mw=w
        ratio=mx/margin
        if ratio>worst: worst=ratio; worstB=(tuple(B),float(mx),float(margin),mw,P)
    print(f"k={k}: dangerous(margin<0.2)={danger} checked={checked} skipped(P>{PMAX})={skipped}")
    if worstB:
        B,mx,mg,mw,P=worstB
        print(f"   WORST ratio={float(worst):.3f} at B={B} (pmax={mx:.4f} margin={mg:.4f} w={mw} P={P}) "
              f"{'CLOSES <15' if worst<15 else '!!! >=15'}")
