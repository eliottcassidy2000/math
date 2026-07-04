#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
{12,24} (LRC tight far-runners) vs {7,21} (tournament forbidden H): the base + uniquely-exceptional
smallest-special-multiple analogy. opus-2026-07-04-S66 (owner steer).

GROUNDED:
 * {7,21} (THM-115/THM-029): H(T)=prod I(C_i,2) multiplicative over conflict-graph components, monoid
   {1+2^k}={3,5,9,17,..}. 7 = forbidden atom; 21=3*7 forbidden by inheritance. Among k*7 ONLY 21 is a
   PERMANENT gap (35=5*7 transient, filled at n=7). Base 7 + uniquely-persistent multiple 3*7.
 * {12,24} (this script): {1..11,13,X} tight (M=1/14) IFF X in {12,24}, both at q*=14 (14th-root grid).
   12 => AP {1..13}; 24=2*12. Doubling a single AP runner j->2j preserves tightness ONLY for j=12.
   Base 12 + uniquely-tight multiple 2*12.
SHARED SHAPE: base + its uniquely-exceptional smallest-special-multiple (3=1+2 gen for H; 2=2-adic for LRC),
on the domain's special grid (2-power monoid / 14th-root q*=14). RESONANCE: 14=2*7 = (LRC mult)*(tourn base).
Honest: opposite roles (forbidden vs tight), different mechanisms, no bijection.
"""
import sys
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    x=x-int(x)
    if x<0:x+=1
    return min(x,1-x)
def exact_M(S):
    S=sorted(set(S));cands=set()
    for v in S:
        for k in range(v):cands.add(Fr(2*k+1,2*v))
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for den in (S[i]+S[j],abs(S[i]-S[j])):
                if den:
                    for s in range(den):cands.add(Fr(s,den))
    b=Fr(0);arg=None
    for t in cands:
        v=min(norm(x*t) for x in S)
        if v>b:b=v;arg=t
    return b,arg

print("(1) LRC: {1..11,13,X} tight (M=1/14) <=> X in {12,24}, both q*=14:")
T=[]
for X in [12,24,36,48,60,72]:
    S=sorted(set(list(range(1,12))+[13,X]))
    if len(S)!=13: continue
    M,t=exact_M(S); isT=(M==Fr(1,14))
    if isT:T.append(X)
    print("   X=%3d  M=%-8s tight=%-6s q*=%d"%(X,str(M),isT,t.denominator))
print("   tight X =",T)

print("\n(2) Doubling a single AP {1..13} runner j->2j is tight ONLY for j=12 (=> 24):")
T2=[]
for j in range(7,14):
    S=[2*j if x==j else x for x in range(1,14)]
    if len(set(S))!=13: continue
    M,_=exact_M(sorted(S))
    if M==Fr(1,14): T2.append(j)
    print("   double %2d -> %2d : M=%-8s %s"%(j,2*j,str(M),"TIGHT" if M==Fr(1,14) else ""))
print("   tightness-preserving doublings:",T2,"(only j=12)")

print("\n(3) THE ANALOGY (base + uniquely-exceptional smallest-special-multiple):")
print("     {7,21}:  base 7 (forbidden atom) + 3*7 (only permanent multiple; 5*7=35 transient). mult 3=1+2.")
print("     {12,24}: base 12 (AP completion)  + 2*12 (only tight multiple; 3*12=36 loose).       mult 2 (2-adic).")
print("     grid: 2-power monoid {1+2^k} (H) vs 14th-root q*=14 (LRC).  resonance: 14 = 2*7 = mult*base.")
print("     honest: opposite roles (forbidden vs tight), diff mechanisms, no bijection (leads -> backlog).")
print("DONE.")
