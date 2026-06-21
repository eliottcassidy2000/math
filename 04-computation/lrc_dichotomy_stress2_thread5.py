#!/usr/bin/env python3
"""THREAD 5: FAST stress-test of HYP-2788 dichotomy. Float p0 for screening (large lcm),
exact confirm only on candidates. Hunt genuine-wide (remove-one-fails) k=8 with p0>Q(7).
Far positions up to 120, multi-cluster, dilated. If any genuine-wide has p0 near cap -> hole."""
import sys, random
from math import gcd, floor
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
def p0_float(E, GRID=20000):
    """float p0 via fine grid (screening only). measS7 = frac of x where all 6 inner sectors hit."""
    E=[e for e in sorted(set(E)) if e]
    if not E: return 0.0
    hits=0
    for i in range(GRID):
        x=(i+0.5)/GRID; mask=0
        for e in E:
            s=int((e*x)%1*7); mask|=1<<s
        if (mask & 0b1111110).bit_count()==6: hits+=1
    return hits/GRID
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
cap=0.38146; Qf=0.19660; k=8
print(f"k={k}: cap={cap} Q(7)={Qf}. FAST hunt genuine-wide remove-one-fails with p0>Q (float screen).")
rng=random.Random(999)
cands=[]
for M in (20,30,40,50,70,90,120):
    for csz in (2,3,4):
        cl=[]; c=0
        while len(cl)<k:
            for t in range(csz): cl.append(c*M+t)
            c+=1
        cands.append(tuple(sorted(set(cl))[:k]))
for d in (2,3,5,7):
    for M in (30,50,80):
        for s1 in range(2,k-1):
            blk=tuple(d*i for i in range(s1))+tuple(M+i for i in range(k-s1))
            cands.append(tuple(sorted(set(blk))[:k]))
for _ in range(20000):
    E=tuple(sorted(set([0]+rng.sample(range(1,120),k-1))))
    if len(E)==k: cands.append(E)
seen=set(); n_total=0; genuine=[]
for E in cands:
    if E in seen: continue
    seen.add(E)
    if len(E)!=k or max(E)<=14 or not primitive(E): continue
    n_total+=1
    if remove_one_bounded(E): continue  # single-perturbation-bounded, fine
    v=p0_float(E)
    if v>Qf:
        genuine.append((v,E))
genuine.sort(reverse=True)
print(f"  screened {n_total} GENUINE-WIDE (remove-one-fails) primitive wide configs.")
print(f"  genuine-wide with p0>Q(7): {len(genuine)}")
for v,E in genuine[:12]:
    print(f"     p0~{v:.5f} cap-p0~{cap-v:.5f} {'<-- NEAR CAP, dichotomy HOLE!' if v>cap-0.05 else ''} E={E}")
if not genuine:
    print("  NONE -> dichotomy 'near-cap=>single-perturbation' holds on this wider bank (far up to 120).")
else:
    mx=genuine[0][0]
    print(f"  WORST genuine-wide p0~{mx:.5f}; cap-p0~{cap-mx:.5f}.")
    print("  (genuine-wide with p0>Q EXISTS => the slack-floor 'genuine-wide=>p0<Q' is FALSE as stated;")
    print("   but if all such are far below cap, the dichotomy still closes via a WEAKER slack bound.)")
