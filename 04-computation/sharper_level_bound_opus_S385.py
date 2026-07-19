# opus-2026-07-17-S385 -- HYP-7700: THE SHARPER PER-LEVEL BOUND.
#
# S382's k=5 enumeration stalled.  Diagnosis: the cost is dominated by nodes
# where the residue E is SMALL, because the bound w <= (j/(7-j))/ell_max then
# admits hundreds of speeds, and each was tested by materialising all w+1 arcs
# of D_w and running a merge scan -- O(w) work per candidate.
#
# TWO SHARPENINGS, both exact.
# (1) LOCALITY.  Only arcs MEETING E matter.  For a component [a,b], the arcs
#     k/w +- lam/w that meet it are exactly  w*a - lam <= k <= w*b + lam,
#     a range of length w*(b-a) + 2*lam.  Summed over components the work is
#     O(w*mu(E) + c) instead of O(w).  Since mu(E) is small at the expensive
#     nodes, this is precisely where the old code was slowest.
# (2) THE j=1 ROUNDING TEST.  At the last level E must lie inside D_w, and by
#     the separation lemma (THM-1125) each component must sit in a SINGLE arc.
#     So containment holds iff for every component [a,b]
#            ceil(w*b - lam)  <=  floor(w*a + lam),
#     an O(1) test per component with NO subtraction at all.  Testing longest
#     components first makes failures immediate.
from fractions import Fraction as F
from math import floor, ceil
from itertools import combinations
import random, time
LAM=F(1,14)

def arcs_meeting(w, a, b):
    """indices k whose arc [k/w - lam/w, k/w + lam/w] can meet [a,b]."""
    return range(int(floor(w*a - LAM)), int(ceil(w*b + LAM)) + 1)

def subtract_local(E, w):
    """E \ badArcs(w), touching only arcs that meet E.  Exact."""
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

def covered_by_single(E, w):
    """True iff E subset badArcs(w): each component inside ONE arc.  O(1)/comp."""
    for a,b in E:
        if ceil(w*b - LAM) > floor(w*a + LAM): return False
    return True

# ---------- reference implementation for the correctness check ----------
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
def subtract_ref(E,D):
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
def complement(u):
    out=[]; prev=F(0)
    for a,b in u:
        if a>prev: out.append((prev,a))
        prev=b
    if prev<1: out.append((prev,F(1)))
    return out
def ess(V,drop):
    allv=[]
    for x in V:
        if x not in drop: allv.extend(teeth_full(x))
    return complement(union(allv))

print("(1) CORRECTNESS: local subtraction vs the reference, on real residues")
random.seed(385)
base=list(range(1,14)); bad=0; n=0
for _ in range(60):
    drop=set(random.sample(base,random.randint(2,5)))
    E=ess(base,drop)
    if not E: continue
    for w in random.sample(range(1,300),6):
        n+=1
        if subtract_local(E,w) != subtract_ref(E,union(teeth_full(w))): bad+=1
print(f"    {n} (residue, w) comparisons; mismatches = {bad}")

print()
print("(2) CORRECTNESS: the j=1 rounding test vs actual containment")
bad2=0; n2=0
for _ in range(60):
    drop=set(random.sample(base,random.randint(2,5)))
    E=ess(base,drop)
    if not E: continue
    for w in random.sample(range(1,300),6):
        n2+=1
        if covered_by_single(E,w) != (len(subtract_local(E,w))==0): bad2+=1
print(f"    {n2} comparisons; mismatches = {bad2}")

print()
print("(3) SPEEDUP on an expensive node (small residue, many candidate speeds)")
drop=set([9,10,11,12,13])           # keep the coarse combs -> small residue
E=ess(base,drop)
mu=sum(b-a for a,b in E)
print(f"    residue: {len(E)} components, measure {float(mu):.6f}")
t0=time.time()
for w in range(1,400): subtract_ref(E,union(teeth_full(w)))
t_ref=time.time()-t0
t0=time.time()
for w in range(1,400): subtract_local(E,w)
t_loc=time.time()-t0
t0=time.time()
for w in range(1,400): covered_by_single(E,w)
t_rnd=time.time()-t0
print(f"    reference subtract : {t_ref:.2f}s")
print(f"    local subtract     : {t_loc:.2f}s   ({t_ref/max(t_loc,1e-9):.1f}x faster)")
print(f"    j=1 rounding test  : {t_rnd:.3f}s  ({t_ref/max(t_rnd,1e-9):.0f}x faster)")
