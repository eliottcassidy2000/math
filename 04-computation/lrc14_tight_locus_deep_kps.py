#!/usr/bin/env python3
"""
lrc14_tight_locus_deep_kps  (kind-pasteur 2026-06-17, DISPROVE side, THE CRUX, deep)

Deeper hunt for TIGHT configs (L=0) with GROWING lcm.  If the primitive tight
locus is finite/bounded-lcm => inf L>0 => LRC(14) via singular series.  If we can
build a sequence of distinct primitive tight configs with lcm -> infinity, the
compactness reduction is threatened.

Searches:
 [1] BFS over tight configs allowing 1- AND 2-element replacements (w<=140) from
     BOTH the AP and the sporadic {1..11,13,24}.  Track max lcm / max element.
 [2] Structured resonant families: {1..13}\{j} ∪ {14m} and ∪{multiples}, exact L.
 [3] "scaled-core" families: does any dilation/partial-dilation stay tight while
     lcm grows?  (Scale invariance L(cS)=L(S) means c*S is tight iff S tight, but
     c*S is NOT primitive; we test PARTIAL scalings that keep gcd=1.)
 [4] exhaustive 13-subsets of [1,18] (smart: only those containing 1 for primitivity,
     pruned by float L screen) -> any tight ones beyond AP/sporadic?
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
SPOR=(1,2,3,4,5,6,7,8,9,10,11,13,24)
W=140

print("="*78)
print("[1] BFS over tight configs: 1- and 2-element replacements (w<=%d) from AP & sporadic"%W)
print("="*78)
seeds=[BASE,SPOR]
visited=set(seeds)
frontier=deque(seeds)
tight_all=list(seeds)
cap=6000
while frontier and len(tight_all)<cap:
    S=frontier.popleft(); Sl=list(S)
    # single replacements
    for idx in range(13):
        e=Sl[idx]; rest=set(Sl)-{e}
        for w in range(1,W+1):
            if w in rest: continue
            cand=tuple(sorted(rest|{w}))
            if cand in visited: continue
            if is_tight(cand):
                visited.add(cand); frontier.append(cand); tight_all.append(cand)
    # double replacements (smaller window to bound cost)
    for i in range(13):
        for j in range(i+1,13):
            e1,e2=Sl[i],Sl[j]; rest=set(Sl)-{e1,e2}
            for w1 in range(1,61):
                if w1 in rest: continue
                for w2 in range(w1+1,61):
                    if w2 in rest: continue
                    cand=tuple(sorted(rest|{w1,w2}))
                    if len(cand)!=13 or cand in visited: continue
                    if is_tight(cand):
                        visited.add(cand); frontier.append(cand); tight_all.append(cand)
            if len(tight_all)>=cap: break
        if len(tight_all)>=cap: break

prim=[S for S in tight_all if is_prim(S)]
print(f"tight configs found: {len(tight_all)} (capped at {cap})")
print(f"primitive tight: {len(prim)}")
print(f"distinct lcms (primitive): {sorted(set(lcm(S) for S in prim))}")
print(f"max lcm (primitive): {max(lcm(S) for S in prim)}")
print(f"max element (primitive): {max(max(S) for S in prim)}")
print("ALL primitive tight configs found:")
for S in sorted(prim):
    print(f"   lcm={lcm(S):>8} max={max(S):>3}  {S}")

print("\n"+"="*78)
print("[2] structured resonant family {1..13}\\{j} ∪ {14m}: any tight?")
print("="*78)
fam_tight=[]
for j in range(1,14):
    core=[x for x in range(1,14) if x!=j]
    for m in range(1,30):
        w=14*m
        if w in core: continue
        S=tuple(sorted(core+[w]))
        if len(S)!=13: continue
        if is_tight(S): fam_tight.append((j,m,S))
print(f"tight in this family: {len(fam_tight)}")
for j,m,S in fam_tight:
    print(f"   drop {j}, +{14*m}: lcm={lcm(S)} prim={is_prim(S)}  {S}")

print("\n"+"="*78)
print("[3] partial-scaling: scale a SUBSET of AP by c, renormalize gcd. tight & growing lcm?")
print("="*78)
hits3=[]
for c in range(2,13):
    for k in range(1,13):  # scale the k largest elements by c
        S0=list(range(1,14))
        for idx in range(13-k,13):
            S0[idx]*=c
        S=tuple(sorted(set(S0)))
        if len(S)!=13: continue
        g=reduce(gcd,S); Sp=tuple(x//g for x in S)
        if is_tight(Sp):
            hits3.append((c,k,Sp))
print(f"tight partial-scalings: {len(hits3)}")
for c,k,S in hits3[:30]:
    print(f"   c={c} scale-top-{k}: lcm={lcm(S)} max={max(S)} {S}")

print("\n"+"="*78)
print("[4] exhaustive 13-subsets of [1,17] containing 1 (float-screened) -> tight?")
print("="*78)
pool=list(range(1,18))  # need 13 of 17, contain 1
cnt=0; tight4=[]
for combo in combinations(range(2,18),12):  # choose 12 from {2..17}, plus 1
    S=(1,)+combo
    cnt+=1
    if L_float(S)>1e-9: continue
    if L_exact(S)==0:
        tight4.append(S)
print(f"13-subsets tested: {cnt}; tight found: {len(tight4)}")
for S in tight4:
    print(f"   lcm={lcm(S)} max={max(S)} prim={is_prim(S)} {S}")
