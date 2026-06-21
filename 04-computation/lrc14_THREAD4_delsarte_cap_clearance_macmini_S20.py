#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD4_delsarte_cap_clearance_macmini_S20.py  (THREAD 4, mac-mini-2026-06-21-S20)

ADVERSARIAL audit of the DELSARTE leg of the corrected architecture.

Claim (reframe): measS7 <= L_y (a Bonferroni/Delsarte dual majorant per k), and the proof
spends rigor proving L_y <= cap at the binding rows.  L_y is shape-INDEPENDENT-bound: a single
Krawtchouk-nonnegative certificate g_k(t) with g_k(0)>=1, g_k(t)>=0 (t=1..6) gives
   measS7(E) = E[1{N=0}] <= E[g_k(N)] = sum_r coeff_r S_r(E) =: L_y(E)   for ALL E,
where N = #missed inner sectors, S_r = E[C(N,r)] the factorial moments.

THE LIVE QUESTIONS this script answers:
  (Q1) For each dual g_k, verify g_k(0)>=1 and g_k(t)>=0, t=1..6 (validity of the majorant).
  (Q2) Does L_y itself clear the cap?  i.e. is  max_E L_y(E) <= cap_k ?  This is what the
       proof needs (it bounds ALL shapes at once).  If max_E L_y > cap_k, the Delsarte leg
       does NOT close that row by itself -- a GENUINE residual.  Compute max_E L_y exhaustively
       (bounded span) and compare to cap.
  (Q3) Compare measS7(consec) vs L_y(consec) vs cap: how much slack does the majorant cost?
  (Q4) WHICH row is actually binding for the cap (smallest cap - max measS7)?  The reframe says
       k=10; the cap audit found k=8.  Settle it and report what the Delsarte leg must clear.

NOTE on the duals (THM-534 / sector_thm534check):
  k=8  (R=4): L=1 -S_1 +S_2 -(9/10)S_3 +(3/5)S_4,  g=(t-1)(t-2)(t-4)(t-5)/40
  k=9,10(R=3):L=1 -(13/18)S_1 +(4/9)S_2 -(1/6)S_3, g=-(t-2)(t-3)(t-6)/36
  k=11,12,13(R=2):L=1 -(1/2)S_1 +(1/6)S_2,         g=(t-3)(t-4)/12
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass

def factorial_moments(E, R=6):
    """S_r = E[C(N,r)], N = #missed inner sectors. EXACT via x-breakpoints."""
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*abs(e)))
    bps=sorted(bps); S=[F(0)]*(R+1)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; L=x1-x0
        hit=set(int(7*abs(e)*xm)%7 for e in Enz)  # nonzero inner-sector hits
        # sector 0 is pinned; count free INNER sectors among {1..6}
        free=6-len([s for s in hit if s!=0])
        for r in range(R+1): S[r]+=L*comb(free,r)
    return S

def measS7_from_moments(S):
    return sum((-1)**r*S[r] for r in range(len(S)))

duals = {
  8:  (lambda S: F(1)-S[1]+S[2]-F(9,10)*S[3]+F(3,5)*S[4],  lambda t: F((t-1)*(t-2)*(t-4)*(t-5),40)),
  9:  (lambda S: F(1)-F(13,18)*S[1]+F(4,9)*S[2]-F(1,6)*S[3], lambda t: F(-(t-2)*(t-3)*(t-6),36)),
  10: (lambda S: F(1)-F(13,18)*S[1]+F(4,9)*S[2]-F(1,6)*S[3], lambda t: F(-(t-2)*(t-3)*(t-6),36)),
  11: (lambda S: F(1)-F(1,2)*S[1]+F(1,6)*S[2], lambda t: F((t-3)*(t-4),12)),
  12: (lambda S: F(1)-F(1,2)*S[1]+F(1,6)*S[2], lambda t: F((t-3)*(t-4),12)),
  13: (lambda S: F(1)-F(1,2)*S[1]+F(1,6)*S[2], lambda t: F((t-3)*(t-4),12)),
}
CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}

def primitive(E): return reduce(gcd,[abs(e) for e in E if e!=0],0)==1

print("#"*84)
print("# THREAD 4 -- DELSARTE leg audit: does L_y clear the cap?")
print("#"*84)

# (Q1) validity of each dual
print("\n(Q1) majorant validity: g_k(0)>=1 and g_k(t)>=0 for t=1..6")
for k in [8,9,10,11,12,13]:
    _,g=duals[k]; gs=[g(t) for t in range(7)]
    ok=(gs[0]>=1) and all(gs[t]>=0 for t in range(1,7))
    print(f"  k={k:2d}: g(t)={[str(x) for x in gs]}  valid={'YES' if ok else 'NO -- INVALID MAJORANT'}")

# (Q2)+(Q3)+(Q4) exhaustive
print("\n(Q2/Q3/Q4) exhaustive bounded-span: max measS7, max L_y, vs cap")
SPAN={8:14,9:14,10:14,11:14,12:15,13:16}
rows=[]
for k in range(8,14):
    Ly,_=duals[k]; cap=CAP[k]; span=SPAN[k]
    maxmeas=F(0); maxmeasE=None
    maxLy=F(-10); maxLyE=None
    # also check Ly>=meas pointwise (Bonferroni) and Ly>=0
    bonf_fail=0; tested=0
    for rest in itertools.combinations(range(1,span+1),k-1):
        E=(0,)+rest
        if not primitive(E): continue
        tested+=1
        S=factorial_moments(list(E),6)
        meas=measS7_from_moments(S)
        ly=Ly(S)
        if ly+F(1,10**12) < meas: bonf_fail+=1
        if meas>maxmeas: maxmeas,maxmeasE=meas,E
        if ly>maxLy: maxLy,maxLyE=ly,E
    consec=tuple(range(k)); Sc=factorial_moments(list(consec),6)
    meas_c=measS7_from_moments(Sc); ly_c=Ly(Sc)
    cap_margin_meas = cap-maxmeas
    cap_margin_Ly = cap-maxLy
    rows.append((k,maxmeas,maxmeasE,maxLy,maxLyE,cap,cap_margin_meas,cap_margin_Ly,meas_c,ly_c,bonf_fail,tested))
    print(f"\n  k={k:2d} (span<={span}, {tested} shapes, Bonferroni L_y>=measS7 fails: {bonf_fail})")
    print(f"     max measS7 = {float(maxmeas):.6f}  at {list(maxmeasE)}")
    print(f"     max L_y    = {float(maxLy):.6f}  at {list(maxLyE)}   (majorant)")
    print(f"     cap_{k}     = {float(cap):.6f}")
    print(f"     cap - max measS7 = {float(cap_margin_meas):+.6f}   [TRUE cap margin]")
    print(f"     cap - max L_y    = {float(cap_margin_Ly):+.6f}   "
          f"[{'Delsarte leg CLEARS cap' if cap_margin_Ly>0 else '*** Delsarte leg does NOT clear (genuine residual) ***'}]")
    print(f"     consec: measS7={float(meas_c):.6f}  L_y={float(ly_c):.6f}  (majorant slack {float(ly_c-meas_c):+.6f})")

# (Q4) binding row
print("\n" + "="*84)
print("(Q4) BINDING ROW (smallest TRUE cap margin = cap - max measS7):")
print("="*84)
srt=sorted(rows,key=lambda r:r[6])
for r in srt:
    print(f"  k={r[0]:2d}: cap-maxmeasS7={float(r[6]):+.6f}   cap-maxL_y={float(r[7]):+.6f}")
print(f"\n  => TRUE binding row (cap) = k={srt[0][0]} (margin {float(srt[0][6]):+.6f})")
srtL=sorted(rows,key=lambda r:r[7])
print(f"  => binding row for the DELSARTE leg (cap-maxL_y smallest) = k={srtL[0][0]} (margin {float(srtL[0][7]):+.6f})")
print("\n  REFRAME said binding k=10 (margin 0.0999). Check which is right above.")
