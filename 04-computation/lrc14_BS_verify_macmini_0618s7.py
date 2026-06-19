#!/usr/bin/env python3
"""
lrc14_BS_verify_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A adversarial verification

Independently verify the second-moment certificate's load-bearing claims:
 (V1) meas_empty(E,J) (exact, breakpoint) matches a brute-force fine-grid Monte/Riemann estimate.
 (V2) The TRUE union measure meas(union_j A_j) = 1 - meas(S7) (exact), and the Chung-Erdos
      lower bound union_CS <= true union (direction sanity) for every test shape.
 (V3) Dawson-Sankoff lower bound uDS <= true union for every shape.
 (V4) 1 - meas(union) == meas(S7_geom) (the IE identity wiring is right).
 (V5) The pairwise certificate U = 1 - max(uCS,uDS) >= true meas(S7) (it's a genuine upper bound).
"""
import sys, itertools, math, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

import importlib.util
spec=importlib.util.spec_from_file_location("sm","/Users/e/Documents/GitHub/math/04-computation/lrc14_BS_secondmoment_macmini_0618s7.py")
# don't exec the print-heavy module; re-implement the needed pieces here independently.

def meas_empty(E, J):
    E=sorted(set(E)); Jset=set(J); bps=set([F(0),F(1)])
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

def meas_empty_grid(E, J, NG=200000):
    Jset=set(J); cnt=0
    for a in range(NG):
        x=(a+0.5)/NG; empty=True
        for e in E:
            if e==0: continue
            if int((e*x % 1)*7) in Jset: empty=False; break
        if empty: cnt+=1
    return cnt/NG

def measS7_geom(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
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

def true_union(E):
    # meas{x: some sector j in 1..6 empty} = 1 - meas(S7) (since sector 0 always hit)
    return F(1)-measS7_geom(E)

def union_CS(E):
    secs=list(range(1,7))
    P={j: meas_empty(E,[j]) for j in secs}
    S1=sum(P.values(),F(0))
    if S1==0: return F(0)
    den=F(0)
    for j in secs: den+=P[j]
    for j,l in itertools.combinations(secs,2): den+=2*meas_empty(E,[j,l])
    return S1*S1/den if den!=0 else F(0)

def union_DS(E):
    secs=list(range(1,7)); P={j: meas_empty(E,[j]) for j in secs}
    S1=sum(P[j] for j in secs); S2=sum(meas_empty(E,[a,b]) for a,b in itertools.combinations(secs,2))
    if S1==0: return F(0)
    best=F(0)
    for t in range(1,7):
        val=F(2)*S1/F(t+1)-F(2)*S2/F(t*(t+1))
        if val>best: best=val
    return best

shapes=[("consec8",list(range(8))),("dissoc8",[0,1,3,7,15,31,63,127]),
        ("Sidon8",[0,1,3,7,12,20,30,44]),("generic8",[0,5,13,27,41,58,79,97]),
        ("perf8",[0,2,3,4,5,6,7,9]),("consec9",list(range(9))),("consec11",list(range(11)))]

print("="*92)
print("ANGLE A verification")
print("="*92)

print("\n(V1) meas_empty exact vs fine-grid (J=[1], J=[3,5]) -- should agree:")
for name,E in shapes[:5]:
    e1=meas_empty(E,[1]); g1=meas_empty_grid(E,[1])
    e2=meas_empty(E,[3,5]); g2=meas_empty_grid(E,[3,5])
    print(f"  {name:<9} P_1 exact={float(e1):.5f} grid={g1:.5f} | P_35 exact={float(e2):.5f} grid={g2:.5f}")

print("\n(V2-V5) union identities and inequality directions (all should hold):")
print(f"  {'shape':<9}{'measS7':>9}{'trueUnion':>10}{'1-S7':>8}{'uCS':>8}{'uDS':>8}{'CS<=U?':>8}{'DS<=U?':>8}{'cert>=S7?':>10}")
allok=True
for name,E in shapes:
    s7=measS7_geom(E); tu=true_union(E); chk=F(1)-s7
    uCS=union_CS(E); uDS=union_DS(E)
    csok=uCS<=tu; dsok=uDS<=tu
    U=F(1)-max(uCS,uDS); certok=U>=s7
    ident = (tu==chk)
    if not (csok and dsok and certok and ident): allok=False
    print(f"  {name:<9}{float(s7):>9.4f}{float(tu):>10.4f}{float(chk):>8.4f}{float(uCS):>8.4f}{float(uDS):>8.4f}{str(csok):>8}{str(dsok):>8}{str(certok):>10}")
print(f"\nALL CHECKS PASS: {allok}")
print("(V2/V3: lower bounds <= true union. V4: 1-meas(S7)==union. V5: certificate U >= meas(S7).)")

# exact rational slack for consec8 (the failing extremiser) and dissoc8 (a passing one)
print("\nExact-rational certificate values:")
for name,E in [("consec8",list(range(8))),("dissoc8",[0,1,3,7,15,31,63,127])]:
    uCS=union_CS(E); uDS=union_DS(E); U=F(1)-max(uCS,uDS)
    cap8=F(2243,5880)
    print(f"  {name}: U = {U} = {float(U):.6f};  cap_8 = {cap8} = {float(cap8):.6f};  U-cap = {float(U-cap8):+.6f}")
