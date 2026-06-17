#!/usr/bin/env python3
"""
lrc14_disprove_l_to_zero_kps  (kind-pasteur 2026-06-16, DISPROVE side, part A/C)

Try to drive L -> 0 with primitive 13-sets of GROWING lcm, and find the
smallest positive L (does anything beat 1/1260?).  Also stress compactness:
track lcm of low-L configs.

Strategies:
 (A) single-element perturbations e->w of the AP, w huge (resonant 14m vs generic).
 (B) scale a near-tight core by growing factor c, keep one entry primitive (escape).
 (C) double-resonance: two strangers 14m1, 14m2 replacing two AP entries.
 (D) "tight + minimal perturbation" of sporadic tight configs (if any large-lcm).

We record L vs lcm to see whether L floors at 1/1260 or decays.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd
from functools import reduce

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
def lcm(S): return reduce(lambda a,b:a*b//gcd(a,b), S, 1)
def is_prim(S): return reduce(gcd,S)==1

BASE=list(range(1,14))
TARGET=1/1260.0
best=[]   # (Lfloat, lcm, S, tag)

def consider(S, tag):
    Sf=tuple(sorted(set(S)))
    if len(Sf)!=13: return
    if not is_prim(Sf): return
    Lf=L_float(Sf)
    if Lf<=1e-12: return   # tight (L=0) configs are NOT loose; exclude from inf over loose
    best.append((Lf, lcm(Sf), Sf, tag))

print("="*78)
print("[A] SINGLE perturbations e->w of AP, e in 1..13, w in 14..3000")
print("="*78)
for e in BASE:
    for w in range(14,3001):
        if w in BASE and w!=e: continue
        S=[x for x in BASE if x!=e]+[w]
        consider(S, f"single {e}->{w}")

print("[B] scale near-tight core by c, keep entry 1 primitive (escape sequences)")
# Take AP minus j, scale rest by c, keep a '1' to stay primitive. lcm grows ~ c.
for j in BASE:
    core=[x for x in BASE if x!=j]
    for c in range(2,40):
        # scale 11 of the 12 remaining by c, keep '1' fixed, add nothing -> only 12; need 13.
        # Build: keep 1, scale the rest (11 entries) by c -> 12 entries, then add a 13th = j*c? primitive via the kept 1.
        if 1 not in core: continue
        rest=[x for x in core if x!=1]
        scaled=[1]+[c*x for x in rest]   # 12 entries
        for extra in (c*j, c*j+1, c*j-1, j):
            S=scaled+[extra]
            consider(S, f"escape c={c} drop {j} extra {extra}")

print("[C] double-resonance: replace two AP entries by 14*m1,14*m2")
for e1,e2 in combinations(BASE,2):
    rest=[x for x in BASE if x not in (e1,e2)]
    for m1 in range(1,12):
        for m2 in range(m1,12):
            w1,w2=14*m1,14*m2
            if w1==w2: continue
            if w1 in rest or w2 in rest: continue
            S=rest+[w1,w2]
            consider(S, f"double {e1},{e2}->{w1},{w2}")

print("[D] resonant single: e->14*m, m up to 60 (period-resonant strangers)")
for e in BASE:
    for m in range(1,61):
        w=14*m
        if w in BASE and w!=e: continue
        S=[x for x in BASE if x!=e]+[w]
        consider(S, f"resonant {e}->{w}")

best.sort()
print("\n"+"="*78)
print("LOWEST positive-L (loose) configs found (float):")
print("="*78)
seen=set()
shown=0
for Lf,lc,S,tag in best:
    if S in seen: continue
    seen.add(S)
    Le=L_exact(S)
    print(f"  L={float(Le):.9g}={Le}  lcm={lc}  {tag}\n       S={S}")
    shown+=1
    if shown>=25: break

print("\n"+"="*78)
print("COMPACTNESS STRESS: largest-lcm configs among the lowest-L (top 200 by L)")
print("="*78)
low=sorted(set((s for _,_,s,_ in best)), key=lambda s:L_float(s))[:400]
by_lcm=sorted(low, key=lambda s:lcm(s), reverse=True)[:10]
for S in by_lcm:
    Le=L_exact(S)
    print(f"  lcm={lcm(S):>10}  L={float(Le):.9g}={Le}  S={S}")
print("\nlowest L overall (exact):", min(L_exact(s) for _,_,s,_ in best) if best else None)
