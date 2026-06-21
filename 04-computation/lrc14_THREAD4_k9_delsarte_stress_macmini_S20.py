#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD4_k9_delsarte_stress_macmini_S20.py  (THREAD 4, mac-mini-2026-06-21-S20)

THE RAZOR-THIN ROW.  The Delsarte leg (max_E L_y(E) <= cap_k) is the load-bearing obligation,
and its TIGHTEST row is k=9: cap_9 - max L_y = +0.00138 only.  This is where a cap-violating
shape -- a GENUINE problem for the Delsarte proof leg -- would hide.

This script:
  (1) prints the EXACT k=9 margin (cap_9 - max L_y) as a Fraction;
  (2) STRESSES max_E L_y(E) over k=9 with WIDE span (does the L_y maximizer escape above cap_9
      as span grows?), reporting per-span best L_y and whether ANY shape has L_y > cap_9
      (= a Delsarte-leg cap violation; measS7 would still be < cap, but the *majorant* leg fails);
  (3) reports the argmax shape of L_y at k=9 and how far measS7 of that shape sits below cap.

If any L_y(E) > cap_9 appears, the deg-3 dual is TOO WEAK for k=9 and the Delsarte leg needs a
sharper certificate there (a real, namable residual).  measS7 itself stays safe (it is the cap
that measS7 must clear, and measS7 << L_y).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass

def factorial_moments(E, R=3):
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*abs(e)))
    bps=sorted(bps); S=[F(0)]*(R+1)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; L=x1-x0
        hit=set(int(7*abs(e)*xm)%7 for e in Enz)
        free=6-len([s for s in hit if s!=0])
        for r in range(R+1): S[r]+=L*comb(free,r)
    return S

# k=9 deg-3 dual (R=3): L=1 -(13/18)S_1 +(4/9)S_2 -(1/6)S_3
def Ly9(S): return F(1)-F(13,18)*S[1]+F(4,9)*S[2]-F(1,6)*S[3]
def measS7_from_moments(S): return sum((-1)**r*S[r] for r in range(len(S)))
CAP9=F(1979,4004)

def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1

print("#"*84)
print("# THREAD 4 -- k=9 DELSARTE STRESS (the razor-thin Delsarte-leg row, margin +0.00138)")
print("#"*84)
print(f"\n  cap_9 = {CAP9} = {float(CAP9):.8f}")

best_Ly=F(-10); best_LyE=None; best_span=None
violations=[]
SMAX=20
print(f"\n  per-span max L_y over k=9 primitive shapes, span up to {SMAX}:")
for span in range(8, SMAX+1):
    span_best=F(-10); span_bestE=None
    for rest in itertools.combinations(range(1,span+1),8):
        if rest[-1]!=span: continue   # new shapes (max==span)
        E=(0,)+rest
        if not primitive(E): continue
        S=factorial_moments(list(E),3)
        ly=Ly9(S)
        if ly>span_best: span_best,span_bestE=ly,E
        if ly>best_Ly: best_Ly,best_LyE,best_span=ly,E,span
        if ly>CAP9: violations.append((E,ly))
    if span_bestE is not None:
        over="  <-- > cap_9 !!" if span_best>CAP9 else ""
        print(f"   span={span:2d}: max L_y at span = {float(span_best):.6f}  running max L_y = {float(best_Ly):.6f}{over}")

margin = CAP9 - best_Ly
S_arg=factorial_moments(list(best_LyE),3)
meas_arg=measS7_from_moments(S_arg)
print(f"\n  RUNNING MAX L_y(k=9) = {best_Ly} = {float(best_Ly):.8f}")
print(f"     argmax E = {list(best_LyE)} (span {best_span})")
print(f"     measS7 of that shape = {meas_arg} = {float(meas_arg):.6f}  (cap-measS7 = {float(CAP9-meas_arg):+.6f})")
print(f"  cap_9 - max L_y = {margin} = {float(margin):+.8f}")
print(f"  Delsarte-leg cap clearance at k=9: {'HOLDS (>0)' if margin>0 else '*** FAILS (L_y exceeds cap_9) ***'}")
print(f"\n  L_y > cap_9 violations found: {len(violations)}")
for E,ly in violations[:5]:
    print(f"     VIOLATION E={list(E)} L_y={float(ly):.6f} > cap_9={float(CAP9):.6f}")
if not violations:
    print("     => NO shape has L_y > cap_9 in span<=%d. Delsarte deg-3 dual clears k=9 (barely)." % SMAX)
