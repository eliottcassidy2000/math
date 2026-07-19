# opus-2026-07-17-S381 -- HYP-7650: THE FOUR-SPEED CASE, EXHAUSTIVE.
#
# THE UNIFORM PER-LEVEL BOUND.  At a node with residue E, j combs still to
# place, and all remaining speeds >= r_min, the density bound with all speeds
# equal to r_min gives the most generous budget:
#     ell*(1 - 2*j*lam) <= 2*lam * (j/r_min)   =>   ell <= j / (r_min*(7-j))
# at lam = 1/14.  So if ell_max(E) > j/(r_min*(7-j)) the node is DEAD -- sound,
# because that is the best any admissible assignment could do.  Applying this at
# EVERY level (not just the root) is what makes k=4 feasible.
#   j=1: use the SHARPER single-arc bound ell <= 2*lam/r = 1/(7r)  (THM-1125).
#   j=2: ell <= 2/(5r)   j=3: ell <= 3/(4r)   j=4: ell <= 4/(3r)
# Termination: residues nonempty by LRC(11), LRC(12), LRC(13) -- 10, 11 and 12
# speeds cannot cover at radius 1/14 (gaps >= 1/11, 1/12, 1/13, all > 1/14).
from fractions import Fraction as F
from itertools import combinations
import sys
LAM=F(1,14)
TEETH={}
def T(x):
    if x not in TEETH:
        w=LAM/x; out=[]
        for j in range(0,x+1):
            a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
            if a<b: out.append((a,b))
        TEETH[x]=out
    return TEETH[x]
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
UT={}
def UT_(x):
    if x not in UT: UT[x]=union(T(x))
    return UT[x]
def complement(u):
    out=[]; prev=F(0)
    for a,b in u:
        if a>prev: out.append((prev,a))
        prev=b
    if prev<1: out.append((prev,F(1)))
    return out
def subtract(E,D):
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

def budget(j, r):
    """max coverable component length with j combs all >= r (sharper at j=1)."""
    if j==1: return 2*LAM/r          # single-arc bound, THM-1125
    return F(j, r*(7-j))             # density bound at lam = 1/14

base=list(range(1,14))
stats={'nodes':0,'pruned':0,'empty':0,'checked':0}
tight=set()

def rec(E, j, rmin, chosen, dropped):
    """place j more combs, all >= rmin, to cover E exactly."""
    stats['nodes']+=1
    L=circ_lmax(E)
    if L > budget(j, rmin):
        stats['pruned']+=1; return
    ub=int(budget(j,rmin)/L*rmin) if L>0 else rmin
    # per-level speed bound: w <= (budget coefficient)/L
    if j==1: wb=int((2*LAM)/L)
    else:    wb=int(F(j,(7-j))/L)
    for w in range(rmin, wb+1):
        E2=subtract(E, UT_(w))
        if j==1:
            stats['checked']+=1
            if not E2:
                V=tuple(sorted(set([x for x in base if x not in dropped]+chosen+[w])))
                if len(V)==13: tight.add(V)
        else:
            if not E2:
                stats['empty']+=1     # would contradict LRC(14-j)
                continue
            rec(E2, j-1, w, chosen+[w], dropped)

for quad in combinations(base,4):
    E=ess(base,set(quad))
    if not E: continue
    rec(E, 4, 1, [], set(quad))

print(f"nodes visited: {stats['nodes']}   pruned: {stats['pruned']}")
print(f"empty residues (would contradict LRC(11)/(12)/(13)): {stats['empty']}")
print(f"leaf checks: {stats['checked']}")
print(f"DISTINCT TIGHT FAMILIES: {len(tight)}")
B=tuple(base); T2=(1,2,3,4,5,6,7,8,9,10,11,13,24)
for V in sorted(tight):
    tag = "= {1..13}" if V==B else ("= {1..11,13,24}" if V==T2 else "*** NEW ***")
    print(f"    {list(V)}   [{tag}]")
