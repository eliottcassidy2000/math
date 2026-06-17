#!/usr/bin/env python3
"""
lrc14_disprove_counterexample_v2_kps  (kind-pasteur 2026-06-17, DISPROVE side, part B, FAST)

Hunt for GENUINE LRC(14) counterexample: 13-set with max-min < 1/14.
max-min := max_tau min_v ||v tau||.  (This is exactly the Lonely Runner gap for
13 runners with speeds S; LRC(14) <=> max-min >= 1/14 for all 13-sets.)

FAST pipeline:
 1. float screen: evaluate min_v||v tau|| on a fine deterministic grid built from
    ALL kink points {j/v} (the true max sits at a kink or a crossing; the lower
    envelope is concave-ish between kinks so the kink set + midpoints gives a tight
    float estimate). Take the configs with the LOWEST float max-min.
 2. EXACT verify max-min on the survivors using the exact candidate set.
Report the smallest max-min; flag anything < 1/14.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd
from functools import reduce

THRESH=Fr(1,14)
THRESHf=1.0/14.0

def fn_f(x):  # ||x|| float
    x=x-int(x)
    if x<0: x+=1.0
    return x if x<=1-x else 1-x

def maxmin_float(S):
    # candidate tau: all kinks j/v and midpoints between consecutive kinks
    ks=set()
    for v in S:
        for j in range(v+1):
            ks.add(j/v)
    ks=sorted(ks)
    cand=list(ks)
    for i in range(len(ks)-1):
        cand.append(0.5*(ks[i]+ks[i+1]))
    best=0.0; bt=0.0
    for t in cand:
        m=1.0
        for v in S:
            d=fn_f(v*t)
            if d<m: m=d
        if m>best: best=m; bt=t
    return best,bt

def fn(x):  # exact ||x||
    f=x-(x.numerator//x.denominator)
    return f if f<=1-f else 1-f

def maxmin_exact(S):
    S=sorted(set(S))
    cands=set()
    for v in S:
        for j in range(v+1):
            cands.add(Fr(j,v))
    n=len(S)
    for i in range(n):
        a=S[i]
        for jx in range(i+1,n):
            b=S[jx]
            for p in range(a+1):
                for q in range(b+1):
                    if a!=b:
                        t=Fr(p-q,a-b)
                        if 0<=t<=1: cands.add(t)
                    t2=Fr(p+q,a+b)
                    if 0<=t2<=1: cands.add(t2)
    best=Fr(-1); bt=None
    for t in cands:
        m=min(fn(v*t) for v in S)
        if m>best: best=m; bt=t
    return best,bt

def lcm(S): return reduce(lambda a,b:a*b//gcd(a,b),S,1)
def is_prim(S): return reduce(gcd,S)==1

BASE=list(range(1,14))
SPOR=[1,2,3,4,5,6,7,8,9,10,11,13,24]

scored=[]  # (maxmin_float, S, tag)
def add(S,tag):
    Sf=tuple(sorted(set(S)))
    if len(Sf)!=13: return
    mm,_=maxmin_float(Sf)
    scored.append((mm,Sf,tag))

# families (tight-like)
add(BASE,"AP"); add(SPOR,"sporadic")
for e in BASE:
    for w in range(14,201):
        if w in BASE and w!=e: continue
        add([x for x in BASE if x!=e]+[w], f"single {e}->{w}")
for e in SPOR:
    for w in range(14,121):
        if w in SPOR and w!=e: continue
        add([x for x in SPOR if x!=e]+[w], f"spor {e}->{w}")
WR=list(range(14,71))
for e1,e2 in combinations(BASE,2):
    rest=[x for x in BASE if x not in (e1,e2)]
    for w1,w2 in combinations(WR,2):
        if w1 in rest or w2 in rest: continue
        add(rest+[w1,w2], f"double {e1},{e2}->{w1},{w2}")
for w in range(13,201):
    add(list(range(1,13))+[w], f"dense 1..12,{w}")
for a,b in combinations(range(12,45),2):
    add(list(range(1,12))+[a,b], f"1..11,{a},{b}")
# also: small dense sets {1..13} with one element doubled-ish & arithmetic-progression variants
for d in range(2,8):
    for start in range(1,4):
        S=[start+d*i for i in range(13)]
        if reduce(gcd,S)==1:
            add(S, f"AP start {start} diff {d}")

scored.sort(key=lambda x:x[0])
print("="*78)
print("Total configs screened:", len(scored))
print("Lowest FLOAT max-min (candidate counterexamples sit near 1/14=%.7f):"%THRESHf)
print("="*78)
# exact-verify the lowest 120 by float
seen=set(); verified=[]
for mmf,S,tag in scored:
    if S in seen: continue
    seen.add(S)
    mm,t=maxmin_exact(S)
    verified.append((mm,S,tag,t))
    if len(verified)>=120: break

verified.sort(key=lambda x:x[0])
print("\nEXACT max-min, smallest 25 (counterexample iff < 1/14):")
for mm,S,tag,t in verified[:25]:
    flag=""
    if mm<THRESH: flag=" <<<<< GENUINE COUNTEREXAMPLE max-min < 1/14"
    elif mm==THRESH: flag="  (=1/14 exactly: LRC tight)"
    print(f"  maxmin={float(mm):.9g}={mm}{flag}")
    print(f"        tau={t}  S={S} lcm={lcm(S)} prim={is_prim(S)} [{tag}]")

mn=min(v[0] for v in verified)
print("\nSMALLEST exact max-min over verified survivors:", mn, "=", float(mn))
print("1/14 =", float(THRESH))
print("Below 1/14? ", mn<THRESH)
nb=sum(1 for v in verified if v[0]<THRESH)
print("count below 1/14:", nb)
