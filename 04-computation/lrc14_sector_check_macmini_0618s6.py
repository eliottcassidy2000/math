#!/usr/bin/env python3
"""
lrc14_sector_check_macmini_0618s6.py  (mac-mini-2026-06-18-S6)

CRITICAL CHECK of codex HYP-2603 (seven-sector net cap):  meas(S7(E)) <= cap_k ?
S7(E) = {x in [0,1): every fixed sector [j/7,(j+1)/7), j=0..6, is hit by some frac(e x), e in E}.
N(E)  = {x: cluster phases form a 1/7-net (all gaps <= 1/7)} = {maxgap <= 1/7}; N subset S7.
The CRUX (HYP-2602) needs meas(N(E)) <= cap_k = m_k.  HYP-2603 strengthens to meas(S7)<=cap_k.

CONCERN: meas(S7) has 1-dim MAIN term M7(k)=Sum_{T subset {1..6}}(-1)^|T| (1-|T|/7)^{k-1}
(the e=0 point auto-fills sector 0; the k-1 nonzero offsets average over B_T). For k=8,
M7(8) ~ 0.486 > cap_8 = 2243/5880 ~ 0.381. So if the singular-series corrections vanish for
DISSOCIATED (high-relation-height) E, meas(S7) -> 0.486 and HYP-2603 FAILS. Test it.

OUTPUT per shape: exact meas(S7), meas(N)=1-mu_{1/7}, vs cap_8. Trend over increasingly
dissociated/spread shapes. Also the main term M7(k).
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)

# ---- meas(S7): all 7 sectors hit ----
def measS7(E):
    E=sorted(set(E))
    # breakpoints where some sector membership changes: x=m/(7e), e!=0
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps)
    total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        secs=set()
        for e in E:
            y=(e*xm)%1
            secs.add(int(y*7))   # sector index floor(7 y) in 0..6
        if len(secs)==7:
            total+=x1-x0
    return total

# ---- meas(N) = 1 - mu_{1/7}: all gaps <= 1/7 ----
def good_set(E, thr=SEV):
    E=sorted(set(E)); k=len(E); diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1); good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((e*xm)%1,e) for e in E); order=[e for _,e in pts]
        fl=[int((e*xm)//1) for e in order]
        for idx in range(k):
            ec=order[idx]; fc=fl[idx]
            if idx<k-1: en=order[idx+1]; fn=fl[idx+1]; wrap=F(0)
            else: en=order[0]; fn=fl[0]; wrap=F(1)
            A=F(en-ec); Cc=F(fc-fn)+wrap
            if A==0:
                if Cc>thr: good.append((x0,x1));
                continue
            xb=(thr-Cc)/A
            lo,hi=(max(x0,xb),x1) if A>0 else (x0,min(x1,xb))
            if lo<hi: good.append((lo,hi))
    good=sorted(good); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def measN(E): return F(1)-sum((b-a for a,b in good_set(E)),F(0))

cap8 = F(2243,5880)
# main term M7(k)
def M7(k):
    s=F(0)
    for t in range(0,7):
        from math import comb
        s += F((-1)**t * comb(6,t)) * F(7-t,7)**(k-1)
    return s
print(f"main term M7(8) = {M7(8)} = {float(M7(8)):.5f}   cap_8 = {cap8} = {float(cap8):.5f}")
print("="*86)
print(f"{'shape (k=8)':<34}{'meas(S7)':>12}{'meas(N)':>12}{'S7<=cap8?':>12}{'spread':>8}")
print("="*86)
shapes = [
  ("consec {0..7}", list(range(8))),
  ("perf {0,2,3,4,5,6,7,9}", [0,2,3,4,5,6,7,9]),
  ("2-dilate {0,2,..,14}", [2*i for i in range(8)]),
  ("dissoc {0,1,3,7,15,31,63,127}", [0,1,3,7,15,31,63,127]),
  ("Sidon {0,1,3,7,12,20,30,44}", [0,1,3,7,12,20,30,44]),
  ("spread {0,1,2,3,40,41,42,43}", [0,1,2,3,40,41,42,43]),
  ("generic {0,5,13,27,41,58,79,97}", [0,5,13,27,41,58,79,97]),
  ("bigspread {0,1,2,3,100,200,300,400}", [0,1,2,3,100,200,300,400]),
]
for name,E in shapes:
    s7=measS7(E); n=measN(E); sp=max(E)-min(E)
    flag = "OK" if s7<=cap8 else "*** VIOLATES ***"
    print(f"{name:<34}{float(s7):>12.5f}{float(n):>12.5f}{flag:>16}{sp:>8}")
print("\nIf dissoc/spread shapes give meas(S7) ~ 0.486 > cap_8 = 0.381, HYP-2603 is REFUTED")
print("(S7 is too lossy a relaxation of N for high-relation-height clusters).")
print("\nDONE.")
