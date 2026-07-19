# opus-2026-07-17-S385 -- k=5 REBUILT ON THE SHARPER PER-LEVEL BOUND.
#   j>1 : subtract_local  -- only arcs meeting E  (verified exact, 0/360)
#   j=1 : covered_by_single -- O(1) rounding test per component, NO subtraction
#         (verified exact, 0/360; 360x faster than the old leaf check)
# Bounds unchanged and all proved: per-level ell <= j/(r*(7-j)), sharpened to
# 2*lam/r at j=1 (THM-1125); termination by LRC(10)..LRC(13).
from fractions import Fraction as F
from math import floor, ceil
from itertools import combinations
import time
LAM=F(1,14)
def arcs_meeting(w,a,b): return range(int(floor(w*a-LAM)), int(ceil(w*b+LAM))+1)
def subtract_local(E,w):
    out=[]; hw=LAM/w
    for a,b in E:
        cur=a
        for k in arcs_meeting(w,a,b):
            c=F(k,w)-hw; d=F(k,w)+hw
            if d<=cur: continue
            if c>=b: break
            if c>cur: out.append((cur,min(c,b)))
            if d>cur: cur=d
            if cur>=b: break
        if cur<b: out.append((cur,b))
    return out
def covered_by_single(E,w):
    for a,b in E:
        if ceil(w*b-LAM) > floor(w*a+LAM): return False
    return True
def teeth_full(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
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
        if x not in drop: allv.extend(teeth_full(x))
    return complement(union(allv))

base=list(range(1,14))
st={'nodes':0,'pruned':0,'empty':0,'leaf':0}
tight=set(); t0=time.time()

def rec(E,j,rmin,chosen,dropped):
    st['nodes']+=1
    L=circ_lmax(E)
    B = 2*LAM/rmin if j==1 else F(j, rmin*(7-j))
    if L>B: st['pruned']+=1; return
    wb = int((2*LAM)/L) if j==1 else int(F(j,(7-j))/L)
    for w in range(rmin, wb+1):
        if j==1:
            st['leaf']+=1
            if covered_by_single(E,w):
                V=tuple(sorted(set([x for x in base if x not in dropped]+chosen+[w])))
                if len(V)==13: tight.add(V)
        else:
            E2=subtract_local(E,w)
            if not E2: st['empty']+=1; continue
            rec(E2,j-1,w,chosen+[w],dropped)

quints=list(combinations(base,5))
# expensive quintuples first: dropping LARGE speeds keeps the coarse combs and
# leaves a small residue, which is where the old run stalled.
quints.sort(key=lambda q: -sum(q))
for idx,q in enumerate(quints):
    E=ess(base,set(q))
    if E: rec(E,5,1,[],set(q))
    if idx%100==0:
        print(f"  [{idx}/{len(quints)}] nodes={st['nodes']} pruned={st['pruned']} "
              f"leaf={st['leaf']} t={time.time()-t0:.0f}s", flush=True)
print(f"nodes={st['nodes']} pruned={st['pruned']} empty={st['empty']} leaf={st['leaf']}")
print(f"elapsed {time.time()-t0:.0f}s")
print(f"DISTINCT TIGHT FAMILIES: {len(tight)}")
B0=tuple(base); T2=(1,2,3,4,5,6,7,8,9,10,11,13,24)
for V in sorted(tight):
    tag="= {1..13}" if V==B0 else ("= {1..11,13,24}" if V==T2 else "*** NEW ***")
    print(f"    {list(V)}   [{tag}]")
