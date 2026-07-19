# opus-2026-07-17-S385 -- THE WITNESS-POINT BRANCHING RULE.
#
# Every valid cover of the residue E must cover EVERY point of E.  So fix a
# witness point x in E (take the midpoint of the longest component, the most
# constrained one).  Some placed comb w must satisfy ||w x|| < lam -- a condition
# of density 2*lam = 1/7.  Branching on "which comb covers x" therefore explores
# only ~1/7 of the speeds per level instead of all of them.
#
# COMPLETENESS.  This drops the sorted-order symmetry breaking (the covering comb
# need not be the smallest), so assignments are counted up to j! times; we dedupe
# by canonicalising the final family.  Net factor 7^j / j!  -- at j=5 that is
# 16807/120 = 140x.  Validated below against the KNOWN k=4 answer.
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
def covers_point(w,x):
    r=(w*x) % 1
    return min(r,1-r) < LAM
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
def longest_comp(c):
    best=None
    for a,b in c:
        if best is None or (b-a)>(best[1]-best[0]): best=(a,b)
    return best
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
def run(k, quints_limit=None, report=True):
    st={'nodes':0,'pruned':0,'leaf':0,'empty':0}
    tight=set(); t0=time.time()
    def rec(E,j,chosen,dropped):
        st['nodes']+=1
        L=circ_lmax(E)
        # per-level budget: all remaining speeds >= 1 in the worst case, so use
        # the global bound with r = 1 (weaker but still exact for the ell_max test)
        wb = int((2*LAM)/L) if j==1 else int(F(j,(7-j))/L)
        if wb < 1: st['pruned']+=1; return
        a,b = longest_comp(E)
        x = (a+b)/2                                   # witness point
        for w in range(1, wb+1):
            if not covers_point(w,x): continue          # 1/7 filter, exact
            if j==1:
                st['leaf']+=1
                if covered_by_single(E,w):
                    V=tuple(sorted(set([y for y in base if y not in dropped]+chosen+[w])))
                    if len(V)==13: tight.add(V)
            else:
                E2=subtract_local(E,w)
                if not E2: st['empty']+=1; continue
                rec(E2,j-1,chosen+[w],dropped)
    combs=list(combinations(base,k))
    combs.sort(key=lambda q:-sum(q))
    if quints_limit: combs=combs[:quints_limit]
    for q in combs:
        E=ess(base,set(q))
        if E: rec(E,k,[],set(q))
    el=time.time()-t0
    if report:
        print(f"    k={k}: nodes={st['nodes']} leaf={st['leaf']} pruned={st['pruned']} "
              f"elapsed={el:.1f}s  families={len(tight)}")
        B0=tuple(base); T2=(1,2,3,4,5,6,7,8,9,10,11,13,24)
        for V in sorted(tight):
            tag="{1..13}" if V==B0 else ("{1..11,13,24}" if V==T2 else "*** NEW ***")
            print(f"        {list(V)}  [{tag}]")
    return tight

print("(4) VALIDATE the witness-branching rule against the KNOWN k=4 answer")
print("    (S381: 2242028 nodes, 414499 leaf checks, exactly 2 families)")
t4=run(4)
print()
print("(5) k=3 cross-check (S380 found exactly 2 families)")
t3=run(3)

print()
print("(6) k=5 -- the level that stalled under the old search")
t5=run(5)
