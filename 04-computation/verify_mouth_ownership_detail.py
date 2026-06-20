#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Detail checks of the exact numerical claims in the mouth-ownership angle."""
from fractions import Fraction
import importlib.util, sys, os

spec = importlib.util.spec_from_file_location(
    "vmo", os.path.join(os.path.dirname(__file__), "verify_mouth_ownership_independent.py"))
# Avoid re-running its __main__; load functions only by reading source manually.

THETA = Fraction(1,14); THR2 = Fraction(426,35035)

def danger_arcs_tagged(d):
    out=[]; den=14*d
    for m in range(0,d+1):
        out.append((Fraction(14*m-1,den), Fraction(14*m+1,den), (d,m,'L'), (d,m,'R')))
    return out

def safe_components(core):
    intervals=[]; wall_at={}
    for d in core:
        if d==0: continue
        for lo,hi,lt,rt in danger_arcs_tagged(d):
            a=max(lo,Fraction(0)); b=min(hi,Fraction(1))
            if a<b: intervals.append((a,b))
            if 0<=lo<=1: wall_at.setdefault(lo,set()).add(lt)
            if 0<=hi<=1: wall_at.setdefault(hi,set()).add(rt)
    intervals.sort(); merged=[]
    for a,b in intervals:
        if merged and a<=merged[-1][1]:
            if b>merged[-1][1]: merged[-1]=(merged[-1][0],b)
        else: merged.append((a,b))
    comps=[]; cur=Fraction(0)
    for a,b in merged:
        if a>cur: comps.append((cur,a))
        if b>cur: cur=b
    if cur<1: comps.append((cur,Fraction(1)))
    return [(lo,hi,wall_at.get(lo,set()),wall_at.get(hi,set())) for lo,hi in comps]

def meas(core): return sum((hi-lo for lo,hi,_,_ in safe_components(core)),Fraction(0))
def core_from(holes,tails): return tuple(sorted([d for d in range(1,14) if d not in set(holes)]+list(tails)))

print("=== CLAIM (1): exact mouth endpoints/lengths match stated ===")
drop6=core_from({6},[])
claimed=[("[29/182,9/56]",Fraction(29,182),Fraction(9,56),Fraction(1,728),3),
         ("[29/168,27/154]",Fraction(29,168),Fraction(27,154),Fraction(5,1848),5),
         ("[127/154,139/168]",Fraction(127,154),Fraction(139,168),Fraction(5,1848),5),
         ("[47/56,153/182]",Fraction(47,56),Fraction(153,182),Fraction(1,728),3)]
actual=[(lo,hi,hi-lo) for lo,hi,_,_ in safe_components(drop6)]
for (lbl,clo,chi,clen,cdet),(alo,ahi,alen) in zip(claimed,actual):
    ok = (clo==alo and chi==ahi and clen==alen)
    print(f"  {lbl}: claimed=({clo},{chi},len{clen})  actual=({alo},{ahi},len{alen})  {'MATCH' if ok else 'MISMATCH'}")

print("\n=== CLAIM (2): per-tower-bit deletion measures ===")
claimed_del={1:Fraction(239,3003),2:Fraction(461,12012),4:Fraction(389,12012),8:Fraction(2243,42042)}
for s,cv in claimed_del.items():
    av=meas(tuple(sorted(set(drop6)-{s})))
    print(f"  del {s}: claimed={cv}  actual={av}  {'MATCH' if cv==av else 'MISMATCH'}")

print("\n=== CLAIM (3): two-tail surviving det=3 mass = 1/364 ; new mass = 20667/1611610 ===")
two=core_from({4,6,10},[20,46])
# surviving mouths are the two det=3 outer mouths: [29/182,9/56] and [47/56,153/182], each 1/728
surv = Fraction(1,728)+Fraction(1,728)
print(f"  surviving det=3 mass = 2*(1/728) = {surv}  claimed 1/364={Fraction(1,364)}  {'MATCH' if surv==Fraction(1,364) else 'MISMATCH'}")
mt=meas(two)
print(f"  meas(two-tail) = {mt} = {float(mt):.9f}")
# 'New mass' = mass NOT in the surviving outer mouths = total - 1/364
new_mass = mt - surv
print(f"  meas - surviving = {new_mass} = {float(new_mass):.9f}")
print(f"  claimed new mass = 20667/1611610 = {float(Fraction(20667,1611610)):.9f}")
print(f"  match? {new_mass==Fraction(20667,1611610)}")
print(f"  claim: new mass alone ({float(new_mass):.6f}) exceeds THR2 ({float(THR2):.6f})? {new_mass>THR2}")

print("\n=== CLAIM (4): det pattern [3,5,5,3] outer=3 inner=5 ===")
def detn(lo,hi,ol,oh):
    rw=[(v,a) for (v,a,sd) in ol if sd=='R']; lw=[(w,b) for (w,b,sd) in oh if sd=='L']
    v,a=min(rw); w,b=min(lw); return v*(14*b-1)-w*(14*a+1)
dets=[detn(lo,hi,ol,oh) for lo,hi,ol,oh in safe_components(drop6)]
print(f"  determinant numerators in order = {dets}  claimed [3,5,5,3]  {'MATCH' if dets==[3,5,5,3] else 'MISMATCH'}")
