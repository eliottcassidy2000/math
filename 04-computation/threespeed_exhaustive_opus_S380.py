# opus-2026-07-17-S380 -- THE EXHAUSTIVE THREE-SPEED ENUMERATION.
#
# Chain of bounds (all proved): r <= 3/(4*ell_max(E_ijk)) [3-comb density],
# then s <= 2/(5*ell_max(E')) [2-comb density, S379], then t <= 2*lam/ell_max(E'')
# [single-arc swap bound, S378].  Termination: E' nonempty by LRC(12) (11 speeds
# cannot cover at 1/14), E'' nonempty by LRC(13) (12 speeds cannot).
#
# HARD PRUNE that makes this feasible: with r <= s <= t, the best possible
# two-comb budget for the residue is at s = t = r, giving ell <= 2/(5r).  So if
# ell_max(E') > 2/(5r), NO admissible (s,t) exists and the branch dies at once.
from fractions import Fraction as F
from itertools import combinations
import sys
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
TEETH={}
def T(x):
    if x not in TEETH: TEETH[x]=teeth01(x)
    return TEETH[x]
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def complement(u):
    out=[]; prev=F(0)
    for a,b in u:
        if a>prev: out.append((prev,a))
        prev=b
    if prev<1: out.append((prev,F(1)))
    return out
def subtract(E,D):
    """E \ D, both sorted disjoint interval lists -- LINEAR merge scan."""
    out=[]; j=0; n=len(D)
    for a,b in E:
        cur=a
        while j<n and D[j][1]<=cur: j+=1
        jj=j
        while jj<n and D[jj][0]<b:
            c,d=D[jj]
            if c>cur: out.append((cur,min(c,b)))
            if d>cur: cur=d
            if cur>=b: break
            jj+=1
        if cur<b: out.append((cur,b))
    return out
def circ_lmax(c):
    if not c: return F(0)
    if len(c)>1 and c[0][0]==0 and c[-1][1]==1:
        wrap=(c[0][1]-c[0][0])+(c[-1][1]-c[-1][0])
        inner=max((b-a) for a,b in c[1:-1]) if len(c)>2 else F(0)
        return max(wrap,inner)
    return max(b-a for a,b in c)
def ess(V,drop):
    allv=[]
    for x in V:
        if x not in drop: allv.extend(T(x))
    return complement(union(allv))

base=list(range(1,14))
tight=set(); combos=0; pruned=0; empty1=0; empty2=0; branches=0
for i,j,k in combinations(base,3):
    E=ess(base,{i,j,k}); L=circ_lmax(E)
    if L==0: continue
    rb=int(F(3,4)/L)
    for r in range(1,rb+1):
        Ep=subtract(E,union(T(r)))
        if not Ep: empty1+=1; continue          # would contradict LRC(12)
        Lp=circ_lmax(Ep)
        if Lp > F(2,5*r): pruned+=1; continue   # no (s,t) with s,t >= r can work
        sb=int(F(2,5)/Lp)
        for s in range(r,sb+1):
            branches+=1
            Epp=subtract(Ep,union(T(s)))
            if not Epp: empty2+=1; continue     # would contradict LRC(13)
            Lpp=circ_lmax(Epp)
            tb=int((2*LAM)/Lpp)
            for t in range(s,tb+1):
                combos+=1
                if not subtract(Epp,union(T(t))):
                    V=tuple(sorted(set([x for x in base if x not in (i,j,k)]+[r,s,t])))
                    if len(V)==13: tight.add(V)
print(f"branches explored: {branches}   (triple,r) pruned by the 2/(5r) test: {pruned}")
print(f"empty E'  (would contradict LRC(12)): {empty1}")
print(f"empty E'' (would contradict LRC(13)): {empty2}")
print(f"(i,j,k,r,s,t) combinations fully checked: {combos}")
print(f"DISTINCT TIGHT FAMILIES FOUND: {len(tight)}")
B=tuple(base); T2=(1,2,3,4,5,6,7,8,9,10,11,13,24)
for V in sorted(tight):
    tag = "= {1..13}" if V==B else ("= {1..11,13,24}" if V==T2 else "*** NEW ***")
    print(f"    {list(V)}   [{tag}]")
