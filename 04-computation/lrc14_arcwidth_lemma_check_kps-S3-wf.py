#!/usr/bin/env python3
"""
Independent re-derivation & verification of the ARC-WIDTH LEMMA (THM-526) and the
implication C(S) => M(S) >= 1/14, which is the LOAD-BEARING step.

ARC-WIDTH LEMMA claim:
  Runner v's danger set {tau: ||v tau|| < 1/14} consists of v "teeth", each an open arc
  of width 1/(7v) centered at j/v (j=0..v-1), so consecutive tooth CENTERS are 1/v apart
  and each tooth has half-width 1/(14v) => full width 1/(7v). The SAFE gaps between teeth
  have width 1/v - 1/(7v) = 6/(7v). The longest contiguous run of v-DANGER is a single tooth
  = width 1/(7v).
  CONSEQUENCE: if A's widest safe arc W(A) > 1/(7v), then that arc cannot be entirely covered
  by a single v-tooth; since v-teeth are spaced with safe gaps, the arc must contain a point
  where ||v tau|| >= 1/14, i.e. a point safe for A AND for v => M(A u {v}) >= 1/14.

We verify:
  (1) tooth width = 1/(7v) exactly (the danger arc of a single runner).
  (2) the longest danger RUN of runner v is exactly 1/(7v) (no two teeth merge; gap 6/(7v)>0).
  (3) the IMPLICATION end-to-end: for many sets A and v with W(A)>1/(7v), verify a common
      safe point for A u {v} actually exists (so M >= 1/14), by EXACT computation of
      safe_components(A u {v}) being nonempty.
  (4) ADVERSARIAL: is the implication ever FALSE? i.e. W(A)>1/(7v) but A u {v} has NO safe arc?
      (would refute the lemma.)  Also test the boundary W(A) = 1/(7v) exactly.
"""
from fractions import Fraction as F
from math import gcd
import random

H=F(1,14)
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def safe_components(A,h=H):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def Wwidth_and_arc(A):
    sc=safe_components(A)
    if not sc: return F(0),None
    best=F(0); barc=None
    for a,b in sc:
        if b-a>best: best=b-a; barc=(a,b)
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1:
        wrapw=sc[0][1]+(1-sc[-1][0])
        if wrapw>best: best=wrapw; barc=("wrap",sc[-1][0],sc[0][1])
    return best,barc

# (1)(2) tooth geometry
print("Tooth-geometry check (single runner v):")
for v in [1,2,3,5,7,14,53,100]:
    # danger of runner v alone: complement of safe. safe = safe_components([v]).
    sc=safe_components([v])
    # safe arcs each width 6/(7v); danger arcs (gaps) width 1/(7v)
    safe_widths=set(b-a for a,b in sc)
    # account wrap: arc touching 0 and 1
    print(f"  v={v}: #safe arcs={len(sc)}  safe widths={sorted(set(float(w) for w in safe_widths))[:3]} expected 6/(7v)={float(F(6,7*v))}")
    # the danger run = 1 - total safe? compute max danger run = 1/(7v)
    # total safe measure:
    tot=sum(b-a for a,b in sc)
    danger_meas=1-tot
    print(f"      total safe measure={danger_meas if False else tot} ; danger measure={danger_meas}={float(danger_meas)} ; v*1/(7v)={float(v*F(1,7*v))} (should be 1/7 total danger)")

print()
print("Longest single danger run of v = 1/(7v)?  (verify gap between teeth = 6/(7v) > 0):")
for v in [2,3,7,53,100]:
    # tooth centered at 0: (-1/(14v),1/(14v)); next center 1/v. gap=(1/(14v), 1/v - 1/(14v)) width=6/(7v)
    gap = F(1,v)-F(1,14*v)-F(1,14*v)
    print(f"  v={v}: inter-tooth safe gap width={gap}={float(gap)} =? 6/(7v)={float(F(6,7*v))}  tooth width 1/(7v)={float(F(1,7*v))}")

# (3)(4) the IMPLICATION C => M>=1/14, exact end-to-end.
print()
print("IMPLICATION test: W(A) > 1/(7v)  =>  A u {v} has a safe point (M>=1/14):")
rng=random.Random(11)
viol=0; tested=0; boundary_cases=[]
for trial in range(200000):
    nA=rng.randint(3,12)
    A=sorted(set(rng.randint(1,200) for _ in range(nA)))
    if len(A)<2: continue
    v=rng.randint(1,200)
    if v in A: continue
    W,_=Wwidth_and_arc(A)
    thr=F(1,7*v)
    if W>thr:
        tested+=1
        sc2=safe_components(A+[v])
        if not sc2:
            viol+=1
            if viol<=10: print("  !!! IMPLICATION VIOLATION:",A,"v=",v,"W=",float(W),"thr=",float(thr))
    elif W==thr:
        sc2=safe_components(A+[v])
        boundary_cases.append((A[:],v,bool(sc2)))
print(f"  tested {tested} cases with W>1/(7v); implication violations: {viol}")
print(f"  boundary W==1/(7v) cases: {len(boundary_cases)}; of these, A u v safe nonempty in {sum(1 for _,_,s in boundary_cases if s)}")
if boundary_cases[:5]:
    for A,v,s in boundary_cases[:5]:
        print(f"     boundary A={A} v={v} -> AU{{v}} safe? {s}")
print()
print("CONCLUSION on arc-width lemma:", "PASS (0 violations)" if viol==0 else f"FAIL ({viol} violations)")
