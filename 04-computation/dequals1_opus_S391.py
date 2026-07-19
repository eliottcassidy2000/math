# opus-2026-07-19-S391 -- HYP-7770: THE D=1 ACTIVE-PAIR QUESTION, RESOLVED STRUCTURALLY.
#
# THM-1205: g(V) = D/(v_i+v_j) at the straddling active pair, so LRC(14) fails
# only if some family has v_i+v_j > 14*D.  At D = 1 that needs an active pair
# summing to >= 15.  WHAT IS D = 1 ACTUALLY?
#
# D = v_i a_j - v_j a_i = 1 forces gcd(v_i,v_j) = 1.  Canonical instance
# v_i = 1, v_j = s-1: then t* = 1/s, both runners sit at distance EXACTLY 1/s,
# and every other runner v is >= 1/s away iff ||v/s|| >= 1/s iff s does not
# divide v.  So the D = 1 pair (1, s-1) IS THE CLASSICAL SIEVE AT MODULUS s,
# and it delivers gap exactly 1/s.
#   s = 14  -> gap 1/14: the tight families ({1..13}, {1..11,13,24}), q0 = 14.
#   s >= 15 -> gap 1/s < 1/14: too weak, so the maximiser must come from
#              ELSEWHERE, i.e. from a pair with D >= 2.
# PREDICTION: families in the HARD STRATUM (q0 > 14, THM-1105 -- every q <= 14
# divides some speed, so no small-modulus sieve exists) must have ACTIVE PAIRS
# WITH D >= 2.  That is the sharp, testable consequence.
from fractions import Fraction as F
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def maximizer(V):
    Dn=set()
    for a,b in combinations(V,2):
        Dn.add(a+b)
        if a!=b: Dn.add(abs(a-b))
    for v in V: Dn.add(2*v)
    best=F(0); arg=None
    for q in sorted(Dn):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,arg=g,F(p,q)
    return best,arg
def active(V,t,g):
    up=[];dn=[]
    for v in V:
        x=v*t; a=int(x); r=x-a
        if r<0: a-=1; r+=1
        if r==g: up.append((v,a))
        if 1-r==g: dn.append((v,a+1))
    return up,dn
def q0(V):
    q=1
    while any(v%q==0 for v in V): q+=1
    return q

print("(1) CONFIRM: the D=1 pair (1,s-1) at t*=1/s IS the classical sieve at s")
for V,nm in [(list(range(1,14)),"{1,...,13}"),
             ([1,2,3,4,5,6,7,8,9,10,11,13,24],"{1..11,13,24}")]:
    g,t=maximizer(V); up,dn=active(V,t,g)
    vi,ai=up[0]; vj,aj=dn[0]; D=vi*aj-vj*ai
    print(f"    {nm:16s} t*={t} g={g}  pair=({vi},{vj}) D={D} s={vi+vj}  q0={q0(V)}"
          f"   {'s == q0 (sieve)' if vi+vj==q0(V) else ''}")

print()
print("(2) THE PREDICTION: hard-stratum families (q0 > 14) must have D >= 2")
random.seed(391)
rows=[]; tries=0
while len([r for r in rows if r[0]>14])<8 and tries<900:
    tries+=1
    V=sorted(random.sample(range(2,60),13))
    z=q0(V)
    if z<=14 and len([r for r in rows if r[0]<=14])>=6: continue
    g,t=maximizer(V); up,dn=active(V,t,g)
    if not up or not dn: continue
    vi,ai=up[0]; vj,aj=dn[0]; D=vi*aj-vj*ai
    rows.append((z,D,vi,vj,g))
print("    q0   D   pair        g          s=v_i+v_j   s<=14D?")
easy=[r for r in rows if r[0]<=14][:6]; hard=[r for r in rows if r[0]>14][:8]
for z,D,vi,vj,g in easy+hard:
    s=vi+vj
    print(f"    {z:3d} {D:3d}  ({vi:3d},{vj:3d})  {str(g):9s}  {s:5d}       {'yes' if s<=14*D else 'NO'}"
          f"   {'<-- HARD stratum' if z>14 else ''}")
d1hard=[r for r in hard if r[1]==1]
print()
print(f"    hard-stratum families with D = 1: {len(d1hard)} of {len(hard)}")
print("    (prediction says 0 -- a D=1 active pair means a sieve at s = q0 > 14,")
print("     whose gap 1/s < 1/14 cannot be the maximum since LRC(14) holds there)")
