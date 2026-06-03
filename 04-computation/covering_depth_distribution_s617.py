from fractions import Fraction as F
from itertools import combinations
import math

# ---- circular-arc covering: depth distribution p_k for a speed set V, gap delta=1/(n+1), n=|V| ----
def norm_frac(x):  # distance to nearest integer of a Fraction
    x = x - int(x); 
    if x<0: x+=1
    return min(x, 1-x)

def forbidden_arcs(v, delta):
    # runner v forbidden = {t: ||v t|| < delta} = union_{k=0}^{v-1} (k/v - delta/v, k/v + delta/v) on circle [0,1)
    r = delta/v; arcs=[]
    for k in range(v):
        c = F(k,v); lo, hi = c-r, c+r
        # split wraparound into [0,1)
        lo%=1; hi_full=c+r
        a, b = (c-r)%1, (c+r)
        # represent on [0,1): handle wrap
        L=(c-r); H=(c+r)
        L%=1; H%=1
        if L<=H: arcs.append((L,H))
        else: arcs.append((L,F(1))); arcs.append((F(0),H))
    return arcs

def depth_distribution(V):
    n=len(V); delta=F(1,n+1)
    # events: for each runner, +1 entering any of its arcs (a runner contributes 1 if t in its forbidden set;
    # its arcs are disjoint so a simple +/- per arc works for that runner). depth = # runners in forbidden set.
    events=[]  # (point, +1/-1)
    for v in V:
        for (lo,hi) in forbidden_arcs(v,delta):
            if lo<hi:
                events.append((lo,+1)); events.append((hi,-1))
    pts=sorted(set([p for p,_ in events]+[F(0),F(1)]))
    # sweep
    from collections import defaultdict
    pk=defaultdict(F)
    # build delta at each point
    deltamap=defaultdict(int)
    for p,d in events: deltamap[p]+=d
    depth=0; prev=F(0)
    for p in pts:
        if p>prev:
            pk[depth]+=(p-prev)
        depth+=deltamap.get(p,0)
        prev=p
    return dict(pk), delta

def gap(V, res=200000):
    # sup_t min_i ||v_i t|| via fine grid (for reporting the LRC value)
    n=len(V); best=0.0
    for j in range(res):
        t=j/res
        m=min(abs((v*t)%1 - round((v*t)%1)) for v in V)
        if m>best: best=m
    return best

def p0_exact(V):
    pk,delta=depth_distribution(V)
    return pk.get(0,F(0))

print("=== MASTER OBJECT: covering-depth distribution p_k (delta=1/(n+1)) ===")
tests = {
 "AP n=2 [1,2]":[1,2], "AP n=3 [1,2,3]":[1,2,3], "AP n=4 [1,2,3,4]":[1,2,3,4],
 "AP n=5 [1,2,3,4,5]":[1,2,3,4,5],
 "chain (1,3,4,7)":[1,3,4,7], "chain (1,3,4,5,9)":[1,3,4,5,9],
 "random [1,4,6,9]":[1,4,6,9], "random [2,3,7,8]":[2,3,7,8], "random [1,5,8,11,13]":[1,5,8,11,13],
}
for name,V in tests.items():
    pk,delta=depth_distribution(V); n=len(V)
    p0=pk.get(0,F(0)); Edepth=sum(k*m for k,m in pk.items())
    dist=" ".join(f"p{k}={float(pk[k]):.3f}" for k in sorted(pk))
    print(f"{name:22} n={n} delta={delta}  p0={float(p0):.4f} {'<-- COLLAPSE' if p0==0 else ''}  E[depth]={float(Edepth):.3f} (=2n/(n+1)={2*n/(n+1):.3f})")

print("\n=== first-moment check: E[depth]=2n/(n+1) exactly, p0 union bound 1-2n/(n+1) is VACUOUS for n>=2 ===")
for n in range(2,7):
    print(f"  n={n}: E[depth]={2*n/(n+1):.4f}, naive union bound p0>=1-2n/(n+1)={1-2*n/(n+1):+.4f} (negative=>useless)")

print("\n=== COLLAPSE FAMILY search: all distinct sets, min=1, gcd=1, max<=12, sizes 2..5, p0==0 ===")
collapse=[]
for n in range(2,6):
    for combo in combinations(range(1,13), n):
        if combo[0]!=1: continue
        if math.gcd(*combo)!=1: continue
        if p0_exact(list(combo))==0:
            collapse.append(combo)
print(f"found {len(collapse)} collapse sets (gap exactly delta, lonely set measure 0):")
for c in collapse:
    # additive-relation analysis: which elements = sum of two earlier (or any two others)?
    s=set(c); rels=[]
    for x in c:
        for a,b in combinations([y for y in c if y<x],2):
            if a+b==x: rels.append((a,b,x))
        for y in c:
            if y<x and 2*y==x: rels.append((y,y,x))
    chaindepth = len(set(r[2] for r in rels))
    print(f"  {c}   additive relations a+b=c: {rels if rels else 'NONE'}   ({len(rels)} rels, {chaindepth} tops)")
