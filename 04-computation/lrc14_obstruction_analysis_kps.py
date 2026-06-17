#!/usr/bin/env python3
"""
lrc14_obstruction_analysis_kps  (kind-pasteur 2026-06-17, idea fuel for PROVE side)

WHY can't L be pushed below 1/1260, and WHY does max-min floor at 1/14?
Analyze the structure of the minimizer {1..11,13,36} (L=1/1260) vs tight {1..11,13,24}.
 - locate the single lonely arc of the minimizer (where is the leftover measure?)
 - which speed's danger arcs FAIL to cover it, and why 24 succeeds but 36 fails.
 - the role of the AP covering [0,1) by danger arcs with NO gap (the '12' is critical).
Also: characterize the gap above 1/14 in max-min (the second-smallest value 2/27).
"""
import sys, io
try: sys.stdout=io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
except Exception: pass
from fractions import Fraction as Fr
from functools import reduce
from math import gcd

def arcs_of(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A

def union(S):
    A=[]
    for v in set(S): A.extend(arcs_of(v))
    A=sorted((a,b) for a,b in A if b>a)
    out=[]; cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch: ch=max(ch,b)
        else: out.append((cl,ch)); cl,ch=a,b
    if ch is not None: out.append((cl,ch))
    return out

def lonely_gaps(S):
    cov=union(S); gaps=[]; prev=Fr(0)
    for a,b in cov:
        if a>prev: gaps.append((prev,a))
        prev=max(prev,b)
    if prev<1: gaps.append((prev,Fr(1)))
    return gaps

def L_exact(S):
    return sum((b-a for a,b in lonely_gaps(S)), Fr(0))

def covers_point(v,t):
    # does ||v t|| <= 1/14 ?
    x=v*t; f=x-(x.numerator//x.denominator); d=min(f,1-f)
    return d<=Fr(1,14), d

print("="*78)
print("MINIMIZER {1..11,13,36}  vs  TIGHT {1..11,13,24}")
print("="*78)
M=[1,2,3,4,5,6,7,8,9,10,11,13,36]
Tt=[1,2,3,4,5,6,7,8,9,10,11,13,24]
AP=list(range(1,14))

for name,S in [("AP {1..13}",AP),("tight 12->24",Tt),("min 12->36",M)]:
    g=lonely_gaps(S); L=sum((b-a for a,b in g),Fr(0))
    print("\n%s :  L=%s=%.9g   #gaps=%d"%(name,L,float(L),len(g)))
    for a,b in g:
        mid=(a+b)/2
        print("   gap [%s, %s] len=%s  center tau=%s"%(a,b,b-a,mid))
        # which speeds are CLOSEST to covering this gap center?
        dl=sorted(((covers_point(v,mid)[1],v) for v in S))
        print("       closest speeds at center (||v tau||): "+
              ", ".join("v=%d:%s"%(v,d) for d,v in dl[:5]))

print("\n"+"="*78)
print("WHY 12 matters: remove 12 from AP, what gap opens? (the role of speed 12)")
print("="*78)
S0=[x for x in AP if x!=12]  # {1..11,13}, only 12 entries
g=lonely_gaps(S0); print("{1..11,13} (12 speeds) L=%s gaps=%s"%(sum((b-a for a,b in g),Fr(0)),
      [(str(a),str(b)) for a,b in g]))
# the missing 13th speed must cover this gap. 24 covers it (tight), 36 does not fully.
for w in [12,24,36,48,60,72,84,18,30,42]:
    Sg=S0+[w]; L=L_exact(Sg)
    print("   add %d: L=%s=%.9g"%(w,L,float(L)))

print("\n"+"="*78)
print("max-min floor: the second value above 1/14 is 2/27 (single 13->26).")
print("The AP touches 1/14 at tau=13/14 (a single critical tau). Perturbing 12->36")
print("RAISES max-min to 3/41>1/14 but OPENS a lonely gap of measure 1/1260.")
print("Trade-off: you cannot simultaneously keep max-min AT 1/14 AND keep L=0 unless")
print("the config is exactly tight (AP or sporadic).")
print("="*78)
# show: the AP's critical tau and the danger structure there
for S,nm in [(AP,"AP"),(Tt,"sporadic")]:
    # find tau achieving max-min and show it sits on a gap boundary collapse
    print("%s: union of danger covers ALL of [0,1) (L=0): %s"%(nm, L_exact(S)==0))
