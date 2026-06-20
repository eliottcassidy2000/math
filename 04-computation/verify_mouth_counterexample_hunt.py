#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL counterexample hunt for the mouth-ownership claim.

Claim components to attack:
  A. The drop-6 mouths are walled EXCLUSIVELY by {11,12,13}; tower owns NO mouth wall.
     -> Already confirmed exactly. Probe robustness: are 11,12,13 the ONLY speeds whose
        deletion destroys a mouth? Deleting a non-walling speed should keep all 4 mouths.
  B. Deleting a tower bit keeps 4/4 mouths AND pays THR2.
     -> Confirmed for the 4 single deletions. But the STRONGER global claim is:
        you cannot get sub-THR2 with a tower bit deleted (HYP-2661 carry law).
        Hunt: any 12-core AP-tail (holes subset of 1..13, tails>=14) with {1,2,4,8} NOT
        all present, meas < THR2 ?  Should be NONE.
  C. Probe whether tower bits could EVER wall a mouth in a NON-drop6 core (to test the
     'tower owns no mouth' as a drop-6-local statement, not universal).
"""
from fractions import Fraction
from functools import reduce
from math import gcd
from itertools import combinations

THETA=Fraction(1,14); THR2=Fraction(426,35035)

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

drop6=core_from({6},[])
orig_mouths={(lo,hi) for lo,hi,_,_ in safe_components(drop6)}
tower={1,2,4,8}

print("=== A. Delete each single speed from drop-6: which deletions destroy a mouth? ===")
for s in sorted(drop6):
    c=tuple(sorted(set(drop6)-{s}))
    nm={(lo,hi) for lo,hi,_,_ in safe_components(c)}
    surv=len(orig_mouths & nm)
    note="walls a mouth" if surv<4 else "no mouth touched"
    print(f"  del {s:2d}: 4 mouths surviving={surv}/4  ({note})")
print("  -> mouths destroyed ONLY by deletions of {11,12,13}? checking...")
destroyers=[s for s in drop6 if len(orig_mouths & {(lo,hi) for lo,hi,_,_ in safe_components(tuple(sorted(set(drop6)-{s})))})<4]
print(f"     speeds whose deletion destroys >=1 mouth: {sorted(destroyers)}")
print(f"     == {{11,12,13}}? {set(destroyers)=={11,12,13}}")

print("\n=== B. COUNTEREXAMPLE HUNT: AP-tail 12-core, tower NOT intact, meas<THR2 ? ===")
print("    Scan: holes subset of {1..13} (size 1..4), tails in [14..50], |core|=12, primitive.")
print("    Looking for any sub-THR2 core MISSING at least one of {1,2,4,8}.")
found=[]
cnt=0
# 12-core from AP window {1..13}: drop |holes| and add |tails| with |holes| = |tails|+1.
# Scope: nh=1..4 (drop up to 4 AP positions). Tail range grows modestly with nh.
TAIL_MAX={1:60,2:60,3:48,4:30}
for nh in range(1,5):
    tmax=TAIL_MAX[nh]; ntails=nh-1
    for holes in combinations(range(1,14), nh):
        for tails in combinations(range(14,tmax+1), ntails):
            c=core_from(set(holes), tails)
            if len(c)!=12: continue
            if reduce(gcd,c)!=1: continue
            cnt+=1
            if tower <= set(c):
                continue
            m=meas(c)
            if m<THR2:
                found.append((m,holes,tails))
    print(f"    [nh={nh} done, tail<= {tmax}, running cnt={cnt}]")
found.sort()
print(f"    scanned {cnt} cores; sub-THR2 with broken tower: {len(found)}")
for m,h,t in found[:20]:
    print(f"      *** COUNTEREXAMPLE meas={m}={float(m):.7f} holes={h} tails={t}")
if not found:
    print("    NONE — no sub-THR2 AP-tail 12-core with a tower bit missing (HYP-2661 holds in scan).")

print("\n=== B2. For completeness: ALL sub-THR2 AP-tail 12-cores in the scan (tower intact or not) ===")
allbelow=[]
cnt2=0
for nh in range(1,4):
    tmax=TAIL_MAX[nh]; ntails=nh-1
    for holes in combinations(range(1,14), nh):
        for tails in combinations(range(14,tmax+1), ntails):
            c=core_from(set(holes), tails)
            if len(c)!=12: continue
            if reduce(gcd,c)!=1: continue
            cnt2+=1
            m=meas(c)
            if m<THR2:
                allbelow.append((m,holes,tails,tower<=set(c)))
allbelow.sort()
print(f"    scanned {cnt2}; total sub-THR2: {len(allbelow)}")
for m,h,t,ti in allbelow:
    print(f"      meas={m}={float(m):.7f} holes={h} tails={t} tower_intact={ti}")

print("\n=== C. Does any tower speed EVER wall a drop-6 mouth indirectly? ===")
print("    The 4 mouths are in [29/182,9/56],[29/168,27/154],[127/154,139/168],[47/56,153/182].")
print("    Check: for each tower speed, is any tower danger-wall point inside/adjacent a mouth?")
for s in sorted(tower):
    pts=[]
    for lo,hi,lt,rt in danger_arcs_tagged(s):
        for p in (lo,hi):
            if 0<=p<=1: pts.append(p)
    for (lo,hi,_,_) in safe_components(drop6):
        near=[p for p in pts if lo<=p<=hi]
        if near:
            print(f"    speed {s}: wall(s) {near} INSIDE mouth [{lo},{hi}]  (would contradict claim!)")
print("    (no output above = no tower wall lands inside any mouth)")
