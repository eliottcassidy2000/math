#!/usr/bin/env python3
"""
lrc14_tight_locus_deep2_kps  (kind-pasteur 2026-06-17) -- phases [2]-[4] (ASCII out)
Continues the deep tight-locus census after the phase-[1] BFS result
(only 2 primitive tight configs reachable, both lcm=360360).
"""
import sys, io
try: sys.stdout=io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
except Exception: pass
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
def lcm(S): return reduce(lambda a,b:a*b//gcd(a,b),S,1)
def is_prim(S): return reduce(gcd,S)==1
def is_tight(S):
    if L_float(S)>1e-12: return False
    return L_exact(S)==0

print("="*78)
print("[2] structured resonant family AP-drop-j PLUS one 14m: any tight?")
print("="*78)
fam_tight=[]
for j in range(1,14):
    core=[x for x in range(1,14) if x!=j]
    for m in range(1,40):
        w=14*m
        if w in core: continue
        S=tuple(sorted(core+[w]))
        if len(S)!=13: continue
        if is_tight(S): fam_tight.append((j,m,S))
print("tight in this family:", len(fam_tight))
for j,m,S in fam_tight:
    print("   drop %d, +%d: lcm=%d prim=%s  %s"%(j,14*m,lcm(S),is_prim(S),S))

print("\n"+"="*78)
print("[3] partial-scaling: scale top-k AP entries by c, renormalize. tight & growing lcm?")
print("="*78)
hits3=[]
for c in range(2,13):
    for k in range(1,13):
        S0=list(range(1,14))
        for idx in range(13-k,13): S0[idx]*=c
        S=tuple(sorted(set(S0)))
        if len(S)!=13: continue
        g=reduce(gcd,S); Sp=tuple(x//g for x in S)
        if is_tight(Sp): hits3.append((c,k,Sp))
print("tight partial-scalings:", len(hits3))
for c,k,S in hits3[:40]:
    print("   c=%d top-%d: lcm=%d max=%d %s"%(c,k,lcm(S),max(S),S))

print("\n"+"="*78)
print("[4] exhaustive 13-subsets of [1,17] containing 1 (float-screened) -> tight?")
print("="*78)
cnt=0; tight4=[]
for combo in combinations(range(2,18),12):
    S=(1,)+combo; cnt+=1
    if L_float(S)>1e-9: continue
    if L_exact(S)==0: tight4.append(S)
print("13-subsets tested: %d; tight found: %d"%(cnt,len(tight4)))
for S in tight4:
    print("   lcm=%d max=%d prim=%s %s"%(lcm(S),max(S),is_prim(S),S))

print("\n"+"="*78)
print("[5] exhaustive 13-subsets of [1,18] containing 1 (float-screened) -> tight?")
print("="*78)
cnt=0; tight5=[]
for combo in combinations(range(2,19),12):
    S=(1,)+combo; cnt+=1
    if L_float(S)>1e-9: continue
    if L_exact(S)==0: tight5.append(S)
print("13-subsets tested: %d; tight found: %d"%(cnt,len(tight5)))
for S in tight5:
    print("   lcm=%d max=%d prim=%s %s"%(lcm(S),max(S),is_prim(S),S))
