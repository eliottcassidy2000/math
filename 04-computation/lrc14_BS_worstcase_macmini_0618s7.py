#!/usr/bin/env python3
"""
lrc14_BS_worstcase_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A frontier

Find the WORST-CASE shapes for the pairwise second-moment certificate
   U(E) = 1 - max(Chung-Erdos, Dawson-Sankoff)   (rigorous upper bd on meas(S7)),
by exhaustive search over primitive k=8 shapes in a bounded box (0=e_1<e_2<...<e_8<=B).
Goal: identify EXACTLY which shapes have U(E) > cap_8 (the certificate's "uncovered residual"),
and confirm whether the consecutive AP {0..7} is the unique / extremal obstruction, OR whether
non-AP near-consecutive shapes also overshoot.  This is what determines if the certificate
reduces LRC(14)-S3(k=8) to a FINITE checkable residual.

Also report: among shapes with U>cap, are they all "AP-like" (small spread, low relation height)?
And the MAX of meas(S7) itself over the box (to confirm consec is the true extremiser, cross-check
THM-532 part C exhaustively at this k in the box).
"""
import sys, itertools, math
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

def cert_U(E):
    secs=list(range(1,7))
    P={j: meas_empty(E,[j]) for j in secs}
    S1=sum(P.values(),F(0))
    if S1==0: return F(1)
    den=F(0)
    for j in secs: den+=P[j]
    Pjl={}
    for a,b in itertools.combinations(secs,2):
        Pjl[(a,b)]=meas_empty(E,[a,b]); den+=2*Pjl[(a,b)]
    uCS=S1*S1/den if den!=0 else F(0)
    S2=sum(Pjl[(a,b)] for a,b in itertools.combinations(secs,2))
    uDS=F(0)
    for t in range(1,7):
        val=F(2)*S1/F(t+1)-F(2)*S2/F(t*(t+1))
        if val>uDS: uDS=val
    return F(1)-max(uCS,uDS)

cap8=F(2243,5880)

if __name__=="__main__":
    import argparse
    B=12  # box bound on max(E); primitives only (gcd-free not required; dilations give same U by scale-inv)
    print(f"Exhaustive k=8 search, 0=e1<...<e8<=B={B}. cap_8={float(cap8):.5f}")
    print("Looking for shapes with U(E) > cap_8 (certificate's uncovered residual).")
    over=[]; maxS7=(F(0),None); maxU=(F(0),None)
    cnt=0
    for combo in itertools.combinations(range(1,B+1),7):
        E=[0]+list(combo); cnt+=1
        # primitivity: skip if gcd>1 (scale-invariant duplicate)
        g=0
        for e in E: g=math.gcd(g,e)
        if g!=1: continue
        U=cert_U(E)
        s7=measS7_geom(E)
        if s7>maxS7[0]: maxS7=(s7,E)
        if U>maxU[0]: maxU=(U,E)
        if U>cap8:
            over.append((U,E,s7))
    over.sort(reverse=True)
    print(f"\nScanned {cnt} combos. Shapes with U>cap_8: {len(over)}")
    print(f"{'U(E)':>9}{'meas(S7)':>10}  E (primitive, max<= {B})")
    for U,E,s7 in over[:30]:
        spread=max(E)-min(E)
        ap = (E==list(range(8)))
        tag=" <-- AP" if ap else ""
        print(f"{float(U):>9.4f}{float(s7):>10.4f}  {E}  spread={spread}{tag}")
    print(f"\nMAX meas(S7) over box: {float(maxS7[0]):.5f} at {maxS7[1]}  (is it the AP {list(range(8))}? {maxS7[1]==list(range(8))})")
    print(f"MAX U(E)      over box: {float(maxU[0]):.5f} at {maxU[1]}  (is it the AP? {maxU[1]==list(range(8))})")
    # characterize the residual: spreads of overshooting shapes
    if over:
        spreads=sorted(set(max(E)-min(E) for _,E,_ in over))
        print(f"\nSpreads among overshooting shapes: {spreads}  (small spread = AP-like, finite residual)")
        # are ALL overshooting shapes 'covering' the small part well? report how many distinct
        print(f"All overshooting shapes have meas(S7) in [{float(min(s for _,_,s in over)):.4f}, {float(max(s for _,_,s in over)):.4f}]")
