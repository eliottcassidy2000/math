#!/usr/bin/env python3
"""
lrc14_residual_formula_kps  (kind-pasteur, PROVE side)

CLOSED-FORM of the champion residual arc, and a GENERAL two-speed-clash formula.

CLAIM (verified below): the 1/1260 residual lonely arc is a "two-speed clash"
between speed a=5 and speed b=36 with NO other speed covering between them.
The residual arc is
    [ p/a + 1/(14a) ,  q/b - 1/(14b) ]
when p/a + 1/(14a) < q/b - 1/(14b) and consecutive (no center between).
Its length is
    q/b - p/a - 1/(14a) - 1/(14b)
For a=5,b=36,p=2,q=15:  15/36 - 2/5 - 1/70 - 1/504.

GENERAL two-speed clash length between adjacent danger arcs of speeds a,b at
fractions p/a, q/b (with the b-arc to the right of the a-arc):
    Delta = (q/b - p/a) - 1/(14a) - 1/(14b)
The lonely residual is max(0, Delta) per clash. Summed over all uncovered clashes
gives L. We tabulate the smallest achievable POSITIVE Delta combinatorially to see
why 1/2520 (and hence 2x = 1/1260) is the floor for AP single-perturbations.
"""
from fractions import Fraction as F
import math

# verify champion residual length
a,b,p,q = 5,36,2,15
Delta = (F(q,b)-F(p,a)) - F(1,14*a) - F(1,14*b)
print("Champion clash a=5,b=36,p=2,q=15:")
print(f"  q/b - p/a = {F(q,b)-F(p,a)}")
print(f"  minus 1/(14a)+1/(14b) = {F(1,14*a)+F(1,14*b)}")
print(f"  residual Delta = {Delta} = {float(Delta):.6e}")
print(f"  2*Delta = {2*Delta} (two symmetric clashes) ; 1/1260={F(1,1260)}")
print(f"  match: {2*Delta==F(1,1260)}")

print("\n--- General principle ---")
print("Residual length of a clash between a-arc (center p/a) and b-arc (center q/b),")
print("b-arc just to the right, no other center between:")
print("  Delta(a,b,p,q) = q/b - p/a - 1/(14a) - 1/(14b)")
print("This is POSITIVE iff the two 1/14-arcs fail to meet, i.e. the configuration")
print("is loose there. The two competing speeds in the champion are 5 and 36 —")
print("NOT the dropped 12 nor added 36 alone; it's the 5<->36 interface.")

# Now: across ALL single perturbations of the AP, is 1/1260 truly the global min?
# Re-run the exact minimal-positive sweep but with wide w and also report the
# clash decomposition of each champion residual.
def danger(v):
    out=[]; w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(x,y) for x,y in out if y>x]
def union(arcs):
    arcs=sorted((x,y) for x,y in arcs if y>x)
    if not arcs: return []
    res=[]; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=max(ch,hi)
        else: res.append((cl,ch)); cl,ch=lo,hi
    res.append((cl,ch)); return res
def total(a): return sum(y-x for x,y in union(a))
def L_exact(S):
    arcs=[]
    for v in set(S): arcs+=danger(v)
    return F(1)-total(arcs)

print("\n--- Global single-perturbation minimum over AP, w in 14..2000 ---",flush=True)
# fast float screen, exact-confirm only near target
def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=max((14*k-1)*inv,0.0); hi=min((14*k+1)*inv,1.0)
            if lo<hi: arcs.append((lo,hi))
    arcs.sort(); tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=hi if hi>ch else ch
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot
TARGET=1.0/1260.0
C=list(range(1,14)); best=(F(2),None)
ties=[]; cands=[]
for e in C:
    rest=[x for x in C if x!=e]
    for w in range(14,2001):
        if w in rest: continue
        S=rest+[w]
        lf=L_float(S)
        if lf<TARGET*1.5: cands.append((e,w,S))
for e,w,S in cands:
    Le=L_exact(S)
    if Le==0: continue
    if Le<best[0]: best=(Le,(e,w))
    if Le==F(1,1260): ties.append((e,w,tuple(sorted(S))))
print(f"  global min positive L = {best[0]} = {float(best[0]):.6e} at (e,w)={best[1]}")
print(f"  number of (e,w) achieving exactly 1/1260 (w<=2000): {len(ties)}")
for e,w,S in ties[:20]:
    print(f"    drop {e} -> {w}:  S={list(S)}")

# distinct resulting configs achieving 1/1260
cfgs=set(s for _,_,s in ties)
print(f"  distinct configs at L=1/1260: {len(cfgs)}")
for s in sorted(cfgs): print(f"    {list(s)}")
