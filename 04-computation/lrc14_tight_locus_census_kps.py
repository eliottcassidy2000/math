#!/usr/bin/env python3
"""
lrc14_tight_locus_census_kps  (kind-pasteur 2026-06-16, DISPROVE side, THE CRUX)

THE CRUX: is the PRIMITIVE TIGHT LOCUS finite / bounded-lcm?
A config S (13 distinct pos ints) is TIGHT iff L(S)=0.
If we can exhibit a SEQUENCE of distinct primitive tight configs with lcm -> infinity,
the compactness reduction is in danger and inf L may be 0 (kills the approach).
If the tight configs all sit at bounded lcm, that supports HYP-2561 (finite locus).

Census strategy: tight configs need MASSIVE danger-arc overlap to cover all of [0,1).
The AP {1..13} is the prototype. We enumerate 13-subsets of [1,N] and record all
L=0 (tight) ones, tracking lcm and max element.  We also test the "AP-with-resonant-
replacements" family systematically (replace k of {1..13} by larger resonant values).

Because exhaustive 13-subsets of [1,N] explodes, we use a smart search:
 - START from the AP, do BFS over single-element replacements that PRESERVE L=0.
 - A replacement e->w preserves tightness only if w's danger arcs still cover the
   gap that e covered. We test all w up to W and keep the tight ones, building a graph
   of tight configs and measuring how large lcm/max-element can grow while staying tight.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd
from functools import reduce
from collections import deque

def danger_arcs(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A
def L_exact(S):
    A=[]
    for v in set(S): A.extend(danger_arcs(v))
    A=sorted((a,b) for a,b in A if b>a)
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch: ch=max(ch,b)
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot
def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=(14*k-1)*inv; hi=(14*k+1)*inv
            if lo<0.0: arcs.append((0.0,hi)); arcs.append((1.0+lo,1.0))
            elif hi>1.0: arcs.append((lo,1.0)); arcs.append((0.0,hi-1.0))
            else: arcs.append((lo,hi))
    arcs=[(a,b) for a,b in arcs if b>a]; arcs.sort()
    tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch:
            if hi>ch: ch=hi
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot
def lcm(S): return reduce(lambda a,b:a*b//gcd(a,b),S,1)
def is_prim(S): return reduce(gcd,S)==1

def is_tight(S):
    if L_float(S)>1e-12: return False
    return L_exact(S)==0

BASE=tuple(range(1,14))
W=200   # replacement upper bound

print("="*78)
print("BFS over TIGHT configs from AP {1..13}, single-element replacements, w<=%d"%W)
print("="*78)
visited=set([BASE])
frontier=deque([BASE])
tight_all=[BASE]
max_elt_seen=13
max_lcm_seen=lcm(BASE)
steps=0
while frontier:
    S=frontier.popleft()
    Slist=list(S)
    for idx in range(13):
        e=Slist[idx]
        rest=set(Slist)-{e}
        for w in range(1,W+1):
            if w in rest: continue
            cand=tuple(sorted(rest|{w}))
            if cand in visited: continue
            steps+=1
            if is_tight(cand):
                visited.add(cand)
                frontier.append(cand)
                tight_all.append(cand)
                me=max(cand); lc=lcm(cand)
                if me>max_elt_seen: max_elt_seen=me
                if lc>max_lcm_seen: max_lcm_seen=lc
    if len(tight_all)>4000: break

prim_tight=[S for S in tight_all if is_prim(S)]
print(f"tight configs reachable (w<=%d, capped): {len(tight_all)}"%W)
print(f"  primitive among them: {len(prim_tight)}")
print(f"  max element seen in any tight config: {max_elt_seen}")
print(f"  max lcm seen in any tight config: {max_lcm_seen}")
print(f"  candidates tested: {steps}")

# the distinct lcms among primitive tight configs:
lcms=sorted(set(lcm(S) for S in prim_tight))
print(f"\n  distinct lcms among primitive tight configs: {lcms[:20]}{' ...' if len(lcms)>20 else ''}")
print(f"  number of distinct lcms: {len(lcms)}")

# largest-max-element tight configs (are they just dilations? check primitivity)
print("\n  TIGHT configs with the LARGEST max element (primitive only):")
prim_tight_sorted=sorted(prim_tight, key=lambda S:max(S), reverse=True)
for S in prim_tight_sorted[:12]:
    print(f"     max={max(S):>4} lcm={lcm(S):>8}  S={S}")

# Are ALL tight configs related to AP/sporadic by bounded replacements?
# Report which AP indices admit a tight replacement and to what values.
print("\n  Single tight replacements directly from AP {1..13}  (e -> w):")
rep={}
for idx in range(13):
    e=BASE[idx]; rest=set(BASE)-{e}
    ws=[w for w in range(14,W+1) if w not in rest and is_tight(tuple(sorted(rest|{w})))]
    if ws: rep[e]=ws
for e in sorted(rep):
    print(f"     {e} -> {rep[e]}")
