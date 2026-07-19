# opus-2026-07-19-S390 -- HYP-7760: THE ACTIVE-PAIR SUM-TO-DETERMINANT RATIO.
#
# Built on boxeph-S120's located-maximizer theorem (HYP-7745).  At the maximizer
# t*, two runners STRADDLE: v_i t* = a_i + g and v_j t* = a_j - g.  Eliminating
# t* gives  g*(v_i+v_j) = v_i a_j - v_j a_i,  i.e.
#       g(V) = D / (v_i + v_j),   D = v_i a_j - v_j a_i  a positive integer.
# THEREFORE
#       LRC(14) for V   <=>   v_i + v_j  <=  14 * D   at the active pair.
# A purely ARITHMETIC condition on ONE pair -- no measures anywhere, so it is on
# the right side of the THM-1185 triage (measure methods are blind to the tight
# families; pointwise ones are not).  This run verifies the formula and measures
# the ratio (v_i+v_j)/D, whose threshold is exactly 14.
from fractions import Fraction as F
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def maximizer(V):
    """exact (g, t*) via the critical-point set (HYP-2059 / THM-1170)."""
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b)
        if a!=b: D.add(abs(a-b))
    for v in V: D.add(2*v)
    best=F(0); arg=None
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,arg=g,F(p,q)
    return best,arg
def active_pair(V,t,g):
    """runners at distance exactly g, split by which side they straddle."""
    up=[];dn=[]
    for v in V:
        x=v*t; a=int(x); r=x-a
        if r<0: a-=1; r+=1
        if r==g: up.append((v,a))          # v t = a + g
        if 1-r==g: dn.append((v,a+1))      # v t = (a+1) - g
    return up,dn

print("(1) VERIFY boxeph's formula g = (v_i a_j - v_j a_i)/(v_i+v_j) at the maximizer")
random.seed(390)
ok=0; n=0; rows=[]
for _ in range(18):
    V=sorted(random.sample(range(1,50),13))
    g,t=maximizer(V)
    up,dn=active_pair(V,t,g)
    if not up or not dn: continue
    vi,ai=up[0]; vj,aj=dn[0]
    D=vi*aj-vj*ai
    n+=1
    pred=F(D,vi+vj)
    if pred==g: ok+=1
    rows.append((g,vi,vj,D,vi+vj))
print(f"    {n} families with a straddling pair; formula exact in {ok}")
print()
print("(2) THE RATIO (v_i+v_j)/D AT THE ACTIVE PAIR -- threshold is exactly 14")
print("    g       v_i  v_j    D   v_i+v_j   ratio   LRC(14)?")
ratios=[]
for g,vi,vj,D,s in rows[:12]:
    r=F(s,D); ratios.append(float(r))
    print(f"    {str(g):9s} {vi:3d}  {vj:3d}  {D:3d}   {s:5d}   {float(r):6.2f}   {'yes' if r<=14 else 'NO'}")
ratios.sort()
print(f"    ratio over {len(ratios)} families: min {ratios[0]:.2f}, median {ratios[len(ratios)//2]:.2f}, max {ratios[-1]:.2f}")

print()
print("(3) THE TIGHT FAMILIES -- ratio must be EXACTLY 14")
for nm,V in [("{1,...,13}",list(range(1,14))),("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
    g,t=maximizer(V); up,dn=active_pair(V,t,g)
    if up and dn:
        vi,ai=up[0]; vj,aj=dn[0]; D=vi*aj-vj*ai
        print(f"    {nm:16s} t*={t}  g={g}  active pair ({vi},{vj})  D={D}  "
              f"(v_i+v_j)/D = {F(vi+vj,D)}  {'== 14 EXACTLY' if F(vi+vj,D)==14 else ''}")
