#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
w AS A LOG-BARRIER / CENTRAL-PATH PARAMETER: is peeling a far speed an exact
mean-field step p0(C u {w}) = Plat(C) + (w*Delta_w)/w, converging to the analytic
center Plat(C) at rate O(1/w) with the residual bounded by the THM-563 period-max?

kind-pasteur-2026-06-30. Follows the Riesz-depth-ladder + signed/unsigned survey.
The reframing the owner asked for: use the far speed w as a smoothing / log-barrier
parameter inside a Beurling-Selberg / Fourier-positivity SIGNED certificate. The
technical claim under it:

  L(C u {w}) = (6/7)*L_inf(C) + Delta_w      (mean-field: far runner -> its mean 6/7)
  |Delta_w| <= (period-max)/w                 (THM-563 single-far Dedekind bound)

so 1/w is a CENTRAL-PATH coordinate: w=inf is the analytic center (the mean-field core
L_inf(C), provable > 0 by LRC<=13 + a Beurling-Selberg minorant on the BOUNDED core),
and the signed content is confined to the O(1/w) barrier residual Delta_w, which is a
periodic (bounded) Dedekind object -- NOT the divergent absolute tail.

We test the SECTOR proxy p0 (engine validated in lrc14_riesz_depth_ladder_kps.py):
  p0(C u {w}) =?= Plat(C) + wD(cells_C, w)/w,   Plat(C) = p0(C) + p1(C)/7,
and that p0(C u {w}) -> Plat(C) as 1/w -> 0 with residual <= supwD/w.
"""
import sys
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---- VERIFIED ENGINE (verbatim) ----
def G0(y):
    f=y-(y.numerator//y.denominator)
    return Fraction(6,7)*f if f<=Fraction(1,7) else Fraction(6,49)-(f-Fraction(1,7))/7
def cells_of(Ep):
    Ep=sorted(set(Ep)); bps={Fraction(0),Fraction(1)}
    for e in Ep:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); out=[]
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in Ep:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        m=[j for j in range(1,7) if j not in hit]
        if len(m)==1:
            if out and out[-1][2]==m[0] and out[-1][1]==lo: out[-1]=(out[-1][0],hi,m[0])
            else: out.append((lo,hi,m[0]))
    return out
def wD(cells,w):
    return sum(G0(w*b-Fraction(s,7))-G0(w*a-Fraction(s,7)) for (a,b,s) in cells)
def p0(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=hi-lo
    return tot
def p1meas(E): return sum(hi-lo for (lo,hi,s) in cells_of(E))
def plat(C): return p0(C)+p1meas(C)/7
# ------------------------------------

def supwD(C, wmax):
    cells=cells_of(C); best=Fraction(0)
    for w in range(2,wmax):
        if reduce(gcd,tuple(C)+(w,))!=1: continue
        v=abs(wD(cells,w))
        if v>best: best=v
    return best

for C in ([0,1,2,3,4,5,6,7], [0,1,2,3,5,7,9,11]):
    print("="*82)
    print(f"CORE C = {C}   (a bounded near-cover; the 'analytic center' is Plat(C))")
    print("="*82)
    Pl=plat(C); cells=cells_of(C); pm=supwD(C, 6*max(C)+400)
    print(f"  Plat(C) = p0 + p1/7 = {float(Pl):.6f} = {Pl}     [center, w=inf]")
    print(f"  period-max sup_w|w*Delta_w| (barrier-residual constant) ~ {float(pm):.4f}")
    print(f"  {'w':>5} {'p0(C+{w}) direct':>18} {'Plat+wD/w (decomp)':>20} {'match?':>7} {'|resid|=|p0-Plat|':>18} {'bound pm/w':>12}")
    for w in (15,20,30,50,100,200,400,800):
        if reduce(gcd,tuple(C)+(w,))!=1: continue
        direct=p0(sorted(C+[w]))
        decomp=Pl + wD(cells,w)/w
        match = (direct==decomp)
        resid=abs(direct-Pl); bound=pm/w
        print(f"  {w:5d} {float(direct):18.6f} {float(decomp):20.6f} {str(match):>7} {float(resid):18.6f} {float(bound):12.6f}   {'resid<=bound OK' if resid<=bound+Fraction(1,10**9) else '*** resid>bound ***'}")
    print(f"  => as 1/w -> 0, p0(C+{{w}}) -> Plat(C). The signed content is the O(1/w)")
    print(f"     Dedekind residual (bounded, periodic), NOT a divergent absolute tail.\n")

print("INTERPRETATION")
print("  1/w is a CENTRAL-PATH coordinate. w=inf = analytic center Plat(C).")
print("  Peeling far runners one at a time = following the central path inward.")
print("  A Beurling-Selberg minorant proves the CENTER (bounded core, finite) is positive;")
print("  the far runner is the barrier whose 1/w residual is a SIGNED but BOUNDED Dedekind")
print("  sum -- exactly the signed content that positivity-only methods miss, now isolated")
print("  to a finite periodic object instead of the divergent |T|>=3 absolute tail.")
print("DONE.")
