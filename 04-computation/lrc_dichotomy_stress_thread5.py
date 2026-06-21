#!/usr/bin/env python3
"""THREAD 5: stress-test the HYP-2788 dichotomy with configs OUTSIDE the heuristic bank.
The dichotomy claims: every wide primitive E (k=8) with p0>Q(7)=0.19660 is single-perturbation-
bounded (remove ONE element -> reduced-span<=14). Counterexample search: find a wide primitive
k=8 config with p0>Q AND remove-one-fails (genuine-wide) AND p0 near cap. If found -> dichotomy
has a hole. Search far positions UP TO 200 (beyond the bank's far_start_hi=34) and multi-cluster
shapes, to probe whether the bank's truncation missed anything.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
ALL_INNER=0b1111110
def sector_of(p): return int((p%1)*7)
def p0(E):
    E=sorted(set(int(x) for x in E)); nz=[e for e in E if e]
    if not nz: return F(0)
    l=reduce(lambda a,b:a//gcd(a,b)*b,nz); d=7*l; den2=2*l
    bps={0,d}
    for e in nz:
        step=l//e; x=0
        for _ in range(7*e+1): bps.add(x); x+=step
    bps=sorted(bps); num=0
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        midnum=lo+hi; mask=0
        for e in nz: mask|=1<<((e*midnum//den2)%7)
        if (mask&ALL_INNER).bit_count()==6: num+=hi-lo
    return F(num,d)
def primitive(E):
    nz=[e for e in E if e]; return reduce(gcd,nz)==1 if nz else False
def reduced_span(S):
    S=sorted(S); g=0
    for a,b in zip(S,S[1:]): g=gcd(g,b-a)
    return 0 if g==0 else (S[-1]-S[0])//g
def remove_one_bounded(E):
    E=tuple(sorted(E))
    for i in range(len(E)):
        if reduced_span(E[:i]+E[i+1:])<=14: return True
    return False
Q=F(38557,196070)  # Q(7); will recompute
# recompute Q(7)=Plat(consec_7) approx = use the known 0.19660
cap=F(2243,5880); k=8
Qf=0.19660
print(f"k={k}: cap={float(cap):.5f} Q(7)~{Qf}. Hunting genuine-wide (remove-one-fails) with p0>Q.")
rng=random.Random(12345)
worst=F(-1); wE=None; n_genuine_above_Q=0; n_total=0
nearcap_genuine=[]
# Structured multi-cluster families with LARGE far positions (beyond bank truncation)
cands=[]
for M in (20,30,40,50,70,100,150,200):
    for csz in (2,3,4):
        cl=[]; c=0
        while len(cl)<k:
            for t in range(csz): cl.append(c*M+t)
            c+=1
        cands.append(tuple(sorted(set(cl))[:k]))
# dilated APs and two-scale
for d in (2,3,5,7):
    for M in (30,50,80,120):
        for s1 in range(2,k-1):
            blk=tuple(d*i for i in range(s1))+tuple(M+i for i in range(k-s1))
            cands.append(tuple(sorted(set(blk))[:k]))
# random wide with far up to 200
for _ in range(120000):
    E=tuple(sorted(set([0]+rng.sample(range(1,200),k-1))))
    if len(E)==k: cands.append(E)
for E in cands:
    if len(E)!=k or max(E)<=14 or not primitive(E): continue
    n_total+=1
    v=p0(E)
    if float(v)>Qf:
        n_genuine_above_Q  # count all above Q
        if not remove_one_bounded(E):
            n_genuine_above_Q+=1
            if v>worst: worst=v; wE=E
            if float(v)>Qf: nearcap_genuine.append((v,E))
nearcap_genuine.sort(reverse=True)
print(f"  scanned {n_total} wide primitive configs.")
print(f"  GENUINE-WIDE (remove-one-fails) with p0>Q(7): {n_genuine_above_Q}")
if wE is not None:
    print(f"  WORST genuine-wide-above-Q: p0={float(worst):.5f} (cap-p0={float(cap-worst):.5f}) E={wE}")
    print(f"  -> dichotomy claim 'near-cap=>single-perturbation' would be VIOLATED if this is near cap")
    for v,E in nearcap_genuine[:8]:
        print(f"     p0={float(v):.5f} cap-p0={float(cap-v):.5f} E={E}")
else:
    print(f"  NO genuine-wide config with p0>Q found -> dichotomy holds on this (much wider) bank.")
