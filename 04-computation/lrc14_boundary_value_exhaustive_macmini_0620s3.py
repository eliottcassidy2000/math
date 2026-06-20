#!/usr/bin/env python3
"""
lrc14_boundary_value_exhaustive_macmini_0620s3.py  (mac-mini-2026-06-20-S3)  -> THM-548 part (i)+(ii)

(i) EXHAUSTIVE max P_r(B) over bounded B subset {0..14}, |B|=k-r (0 in B), for the dangerous
    rows k=8,9,10 and r=2..(k-1).  P_r(B)=sum_t prof_t(B) c_t(r) is the Fatou boundary value.
    Check max_B P_r(B) <= cap_k.  This is the main-term safety of the true-wide decomposition.

(ii) DECOMPOSE real true-wide rows E=B u F (B=elts<=14, F=elts>14): p0(E) vs P_r(B), and the
     resonance correction R=p0(E)-P_r(B), with the resonant far-pairs flagged (small mf+nf').
"""
import sys, math, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if len(set(sector_of(e*xm) for e in E))==7: tot+=x1-x0
    return tot
def miss_profile(B):
    B=sorted(set(B)); bps=set([F(0),F(1)])
    for e in B:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); prof=[F(0)]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; t=7-len(set(sector_of(e*xm) for e in B))
        if t<=6: prof[t]+=x1-x0
    return prof
def c_t(t,r): return sum((-1)**i*math.comb(t,i)*(1-F(i,7))**r for i in range(t+1))  # exact
def P_r_exact(prof,r): return sum(prof[t]*c_t(t,r) for t in range(7))

caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7)}
print("="*78)
print("(i) EXHAUSTIVE max P_r(B) over bounded B subset{0..14}, |B|=k-r, vs cap_k")
print("="*78)
print(f"{'k':>3}{'r':>3}{'|B|':>4}{'#B':>7}{'max P_r(B)':>13}{'argmax B':>26}{'cap_k':>9}{'margin':>9}")
for k in (8,9,10):
    for r in range(2, k):   # r>=2 true-wide; |B|=k-r>=1
        b=k-r
        if b<1: continue
        base=list(range(1,15)); best=F(-1); bestB=None
        for combo in itertools.combinations(base, b-1):
            B=(0,)+combo; prof=miss_profile(B); Pr=P_r_exact(prof,r)
            if Pr>best: best=Pr; bestB=B
        margin=caps[k]-best
        flag='' if margin>0 else '  <-- OVER!'
        print(f"{k:>3}{r:>3}{b:>4}{math.comb(14,b-1):>7}{float(best):>13.5f}{str(bestB):>26}{float(caps[k]):>9.5f}{float(margin):>+9.5f}{flag}")
print("max_B P_r(B) over bounded B is the worst-case Fatou boundary value; if <cap_k, main term is safe.")

print("\n"+"="*78)
print("(ii) DECOMPOSE real true-wide rows: p0(E) = P_r(B) + R,  R=resonance correction")
print("="*78)
def small_rel(f,g,H=6):
    best=None
    for m in range(-H,H+1):
        for n in range(-H,H+1):
            if m==0 and n==0: continue
            val=abs(m*f+n*g)
            if best is None or val<best[0]: best=(val,(m,n))
    return best
rows=[(0,4,6,8,10,12,14,15,16),(0,1,2,4,8,12,16,20),(0,2,4,6,8,10,12,14,15),(0,4,6,8,10,12,14,15,16,17)]
for E in rows:
    B=tuple(e for e in E if e<=14); Fr=tuple(e for e in E if e>14); r=len(Fr)
    if r==0: continue
    prof=miss_profile(B); Pr=P_r_exact(prof,r); p0=measS7(E); R=p0-Pr
    cap=caps.get(len(E),F(4,7))
    print(f"\nE={E}  (k={len(E)})")
    print(f"   B={B} (bounded), F={Fr} (r={r} far)")
    print(f"   p0(E)={float(p0):.5f}  P_r(B)={float(Pr):.5f}  R=p0-P_r={float(R):+.5f}  cap={float(cap):.5f}  margin={float(cap-p0):+.5f}")
    if r>=2:
        print(f"   far-pair resonances (min|mf+nf'|, (m,n)): "+", ".join(
            f"{a}-{b}:{small_rel(a,b)[0]}{small_rel(a,b)[1]}" for a,b in itertools.combinations(Fr,2)))
print("\nDONE.")
