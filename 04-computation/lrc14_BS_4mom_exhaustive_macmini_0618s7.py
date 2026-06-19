#!/usr/bin/env python3
"""
lrc14_BS_4mom_exhaustive_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A closure test

CRITICAL TEST: does the 4-wise moment certificate U4(E) := max{p_0 : binomial moments S_1..S_4
matched} <= cap_8 for EVERY primitive k=8 cluster shape in a box?  If yes (exhaustively), the
4-marginal moment certificate CLOSES LRC(14)-S3 at k=8 (the most dangerous row) -- a rigorous,
exact-rational, single-inequality certificate with NO minorant and NO band-limiting.

We scan all primitive 0=e1<...<e8<=B and record max U4 and any U4>cap_8 (the residual).
Also cross-check: U4 >= meas(S7) always (it's a valid upper bound), and report the worst slack
cap_8 - U4 (>0 => closed).
"""
import sys, itertools, math
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def meas_empty(E, J):
    Jset=set(J); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; empty=True
        for e in E:
            if e==0: continue
            if int(((e*xm)%1)*7) in Jset: empty=False; break
        if empty: total+=x1-x0
    return total

def measS7_geom(E):
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            secs.add(int(((e*xm)%1)*7))
        if len(secs)==7: total+=x1-x0
    return total

def U4(E):
    secs=list(range(1,7))
    S=[F(1)]
    for t in range(1,5):
        st=F(0)
        for J in itertools.combinations(secs,t):
            st+=meas_empty(E,list(J))
        S.append(st)
    n=7  # counts 0..6
    c=np.zeros(n); c[0]=-1.0
    A_eq=np.array([[float(math.comb(i,t)) for i in range(n)] for t in range(5)])
    b_eq=np.array([float(S[t]) for t in range(5)])
    res=linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0,None)]*n, method='highs')
    if not res.success: return None
    return res.x[0]

cap8=F(2243,5880); cap8f=float(cap8)

if __name__=="__main__":
    B=int(sys.argv[1]) if len(sys.argv)>1 else 13
    print(f"Exhaustive k=8, 0=e1<...<e8<=B={B}, primitive. cap_8={cap8f:.6f}")
    print("Testing U4(E) <= cap_8 for all shapes (4-moment certificate closure).")
    over=[]; maxU4=(0.0,None); worst_slack=(1e9,None); cnt=0; viol_s7=0
    for combo in itertools.combinations(range(1,B+1),7):
        E=[0]+list(combo)
        g=0
        for e in E: g=math.gcd(g,e)
        if g!=1: continue
        cnt+=1
        u4=U4(E)
        if u4 is None:
            print("LP fail at",E); continue
        s7=float(measS7_geom(E))
        if u4 < s7-1e-9: viol_s7+=1   # validity check
        if u4>maxU4[0]: maxU4=(u4,E)
        sl=cap8f-u4
        if sl<worst_slack[0]: worst_slack=(sl,E)
        if u4>cap8f:
            over.append((u4,E,s7))
    over.sort(reverse=True)
    print(f"\nScanned {cnt} primitive shapes.")
    print(f"Validity (U4 >= meas(S7)) violations: {viol_s7}  (should be 0)")
    print(f"MAX U4 = {maxU4[0]:.6f} at {maxU4[1]}")
    print(f"Worst slack cap_8 - U4 = {worst_slack[0]:+.6f} at {worst_slack[1]}")
    print(f"Shapes with U4 > cap_8 (residual NOT closed by 4-moment cert): {len(over)}")
    for u4,E,s7 in over[:25]:
        print(f"   U4={u4:.5f}  meas(S7)={s7:.5f}  {E}")
    if not over:
        print("\n*** ALL shapes in box have U4 <= cap_8: the 4-MOMENT CERTIFICATE CLOSES k=8 in this box. ***")
    else:
        print(f"\n{len(over)} shapes still overshoot; characterize them (likely AP/dilated-AP residual).")
