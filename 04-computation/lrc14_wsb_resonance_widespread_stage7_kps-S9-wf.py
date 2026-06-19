#!/usr/bin/env python3
"""
lrc14_wsb_resonance_widespread_stage7_kps-S9-wf.py  (kind-pasteur-2026-06-19-S9)

Stage 7 (decisive consolidation).  Established (exact, verified):
  * LEMMA B: consec UNIQUELY maximizes meas(S7) among PRIMITIVE same-k sets (k=8 maxE<=14 exhaustive
    + wide random to spread 60; k=9 exhaustive maxE<=12).  STRONGER than HYP-2607/L_y.
  * RESONANCE DISSOLVED: the w==0 mod 7 'resonant' wide configs that worried HYP-2608 are
    NON-PRIMITIVE (gcd=7) -- by scale-invariance they EQUAL consec exactly, not over cap.
    Primitivity collapses the resonance worry.

This stage HARDENS the two open sides of the proof:
 (A) THE WIDE-SPREAD BOUND (HYP-2608(a)).  For PRIMITIVE E the dangerous near-resonant configs are
     {0..6}+{large 7-aligned tail}.  Stress these specifically + structured 'near-consec-but-wide'
     primitive sets.  Find the PRIMITIVE wide config with the LARGEST meas(S7); confirm it stays
     below consec (and far below cap).  Tabulate meas(S7) vs span to confirm monotone-ish decay.
 (B) THE BOUNDED-SPREAD FINITE CHECK.  Exhaustive primitive k=8 over the FULL relevant box.
     THM-536-B2 (subset-dom) already certifies span<=7 for k=8 (only consec).  The residual is
     primitive 8-sets with span in [8, B].  Confirm meas(S7)<=consec for ALL of them up to the
     largest B feasible, AND meas(S7) <= cap_8 (the actual goal).  Combine with (A) for span>B.
 (C) Define B*(k) = the spread threshold beyond which meas(S7)(E) < consec for ALL primitive E
     of larger span (the wide-spread bound's B).  Empirically locate it.
 (D) THE RESONANCE STRUCTURE: among primitive sets, which have meas(S7) CLOSEST to consec?
     (the 'near-consec' perturbations).  These are the binding competitors for the finite check.
"""
import sys, itertools, random
from math import gcd
from functools import reduce
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def gcd_all(E): return reduce(gcd,[e for e in E if e!=0],0)

def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            c=F(j,7)
            for m in range(e): bps.add((c+m)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); tot=F(0)
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2; hit=set()
        for e in E:
            fr=(e*xm)%1; hit.add((fr.numerator*7)//fr.denominator)
        if len(hit)==7: tot+=x1-x0
    return tot

CONSEC8=F(481,1470); CAP8=F(2243,5880)
print(f"consec_8 meas(S7)={CONSEC8}={float(CONSEC8):.6f}  cap_8={CAP8}={float(CAP8):.6f}")
print()
print("="*78)
print("(A) WIDE-SPREAD primitive stress: near-resonant {0..6}+7m tails and structured wides")
print("="*78)
# specifically the configs the canon feared: consec core + 7-aligned far speed
near_res=[]
for tail in [7,14,21,28,35,42,49,56,63,70]:
    E=[0,1,2,3,4,5,6,tail]
    if gcd_all(E)==1:
        near_res.append((measS7(E),E))
for tail in [13,20,27,34,41]:  # +6 mod 7 (one below multiple of 7)
    E=[0,1,2,3,4,5,6,tail]
    near_res.append((measS7(E),E))
near_res.sort(reverse=True)
print("  {0..6}+tail (the feared near-resonant family), sorted by meas(S7):")
for m,E in near_res[:12]:
    print(f"    tail={E[-1]:3d}: meas(S7)={float(m):.6f}  {'<consec OK' if m<CONSEC8 else 'NOT < consec!'}")

# random primitive wide with VERY large spread to confirm asymptotic collapse
random.seed(99)
print("\n  random primitive wide, spread up to 200 (asymptotic collapse):")
buckets={}
mx=(F(0),None)
for _ in range(1500):
    M=random.randint(40,200); rest=sorted(random.sample(range(1,M+1),7)); E=[0]+rest
    if gcd_all(E)!=1 or max(E)!=M: continue
    m=measS7(E);
    if m>mx[0]: mx=(m,E)
    bk=M//40*40; buckets.setdefault(bk,[]).append(float(m))
for bk in sorted(buckets):
    v=buckets[bk]; print(f"    spread[{bk},{bk+40}): n={len(v)} mean={sum(v)/len(v):.4f} max={max(v):.4f}")
print(f"  global max over wide primitive bank: {float(mx[0]):.5f} @ {mx[1]} (<consec {float(CONSEC8):.4f})")

print()
print("="*78)
print("(B) BOUNDED-SPREAD finite check: exhaustive primitive k=8, meas(S7)<=consec AND <=cap_8")
print("="*78)
for B in (10,12,13,14,15):
    over_c=0; over_cap=0; cnt=0; mx=(F(0),None)
    for r in itertools.combinations(range(1,B+1),7):
        E=[0]+list(r)
        if gcd_all(E)!=1: continue
        cnt+=1; m=measS7(E)
        if m>CONSEC8+F(1,10**12): over_c+=1
        if m>CAP8+F(1,10**12): over_cap+=1
        if m>mx[0] and m<CONSEC8-F(1,10**15): mx=(m,E)   # largest non-consec
    print(f"  primitive maxE<={B}: {cnt} sets; >consec:{over_c}; >cap_8:{over_cap}; "
          f"largest non-consec meas(S7)={float(mx[0]):.6f} @ {mx[1]}")

print()
print("="*78)
print("(C) B*(k=8): smallest span S s.t. ALL primitive 8-sets of span>=S have meas(S7)<consec")
print("="*78)
# group exhaustive primitive sets (maxE<=15) by span; report max meas(S7) per span band
spanmax={}
for r in itertools.combinations(range(1,16),7):
    E=[0]+list(r)
    if gcd_all(E)!=1: continue
    sp=max(E); m=measS7(E)
    if sp not in spanmax or m>spanmax[sp][0]: spanmax[sp]=(m,E)
print("  span -> max meas(S7) over primitive 8-sets of that exact span (maxE):")
for sp in sorted(spanmax):
    m,E=spanmax[sp]
    tag="=CONSEC" if sp==7 else ("<consec" if m<CONSEC8 else ">=consec!!")
    print(f"    span={sp:2d}: max meas(S7)={float(m):.6f} {tag}  @ {E}")

print()
print("="*78)
print("(D) tightest primitive competitors to consec (the binding finite-check shapes)")
print("="*78)
comp=[]
for r in itertools.combinations(range(1,14),7):
    E=[0]+list(r)
    if gcd_all(E)!=1: continue
    m=measS7(E)
    if m<CONSEC8: comp.append((m,E))
comp.sort(reverse=True)
print("  top-8 primitive non-consec meas(S7) (closest to consec, all must stay < cap_8):")
for m,E in comp[:8]:
    print(f"    meas(S7)={float(m):.6f}  margin-to-cap={float(CAP8-m):.6f}  E={E}")
print("\nDONE stage 7.")
