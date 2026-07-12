#!/usr/bin/env python3
"""cont.40: RESOLVE the framing. Two candidate k=9 base objects:
(J)  sector-emptiness  J = E[N(7-N)]  -- SHIFT-DEPENDENT (positions e_i)
(nu) density floor     nu = mu(goodSet E) -- DIFFERENCE-based => SHIFT-INVARIANT
Compute both for {0..8} vs {1..9} (translates). If nu is equal but J differs, the
base check must be phrased on nu (or J with 0 forced). Which has minimizer matching
THM-661's bar_9 = m_P + 1 - cap_10?"""
from fractions import Fraction as F
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def inter(A,B):
    out=[]
    for a1,b1 in A:
        for a2,b2 in B:
            lo,hi=max(a1,a2),min(b1,b2)
            if lo<hi: out.append((lo,hi))
    return merge(out)
def pair_good(d):
    dd=abs(d)
    if dd==0: return [(F(0),F(1))]
    if d>0: return merge([(F(7*j+1,7*dd),F(7*j+7,7*dd)) for j in range(dd)])
    return merge([(F(7*j,7*dd),F(7*j+6,7*dd)) for j in range(dd)])
def nu(E):
    teeth=[]
    for a in E:
        cur=[(F(0),F(1))]
        for b in E:
            if b!=a: cur=inter(cur,pair_good(b-a))
        teeth.append(cur)
    allint=merge([iv for t in teeth for iv in t])
    return sum(b-a for a,b in allint)
def J(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return sum(p[n]*n*(7-n) for n in range(8))
WMP=F(14249,252252); cap10=F(55,91); bar9=WMP+1-cap10
print(f"bar_9 (THM-661 density-floor target) = m_P + 1 - cap_10 = {bar9} = {float(bar9):.4f}")
print(f"{'family':16s} {'nu=mu(goodSet)':>16s} {'J=E[N(7-N)]':>14s}")
for nm,E in [("{0..8}",list(range(9))),("{1..9}",list(range(1,10))),
             ("{2..10}",list(range(2,11))),("{0,1,..7,9}",[0,1,2,3,4,5,6,7,9])]:
    print(f"  {nm:14s} {float(nu(E)):16.5f} {float(J(E)):14.5f}")
print()
print("SHIFT-INVARIANCE check: nu({0..8}) == nu({1..9})?", nu(list(range(9)))==nu(list(range(1,10))))
print("                        J({0..8}) == J({1..9})? ", J(list(range(9)))==J(list(range(1,10))))
print()
print("=> the DENSITY FLOOR nu is shift-invariant (difference-based, THM-661 object);")
print("   the SECTOR object J is shift-dependent. The base check that feeds the reach")
print("   is nu >= bar_9. Is nu(consec) the min, and does it clear bar_9?")
mn=min((nu(list(range(a,a+9))),a) for a in range(0,6))
print(f"   min nu over consec-shifts = {float(mn[0]):.5f} at shift {mn[1]}, clears bar_9={float(bar9):.4f}: {mn[0]>=bar9}")
