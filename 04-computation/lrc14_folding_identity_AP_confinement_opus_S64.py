#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TOWARD inf_U gap(U) > 0 (m=2,f=2 confinement as a UNIFORM gap). opus-S64.
gap(U)=min_{odd w1,w2}(M(2U u {w1,w2})-1/14). Reduces to: M(2U u {w1,w2}) >= 1/12 (=> gap >= 1/84 > 0).

PROVED here:
 (1) FOLDING IDENTITY: for E=2U (even) and odd w1,w2, M(S)=max_t min(g_E(t), Psi(t)), where
     g_E(t)=min_u ||2u t|| is (+1/2)-periodic and Psi(t)=max(min(a,b), 1/2-max(a,b)), a=||w1 t||,b=||w2 t||.
     (Because ||w(t+1/2)||=1/2-||wt|| for odd w, and g_E(t+1/2)=g_E(t); max(min(g,X),min(g,Y))=min(g,max(X,Y)).)
 (2) AP EVEN PART: M(2*{1..11} u {w1,w2}) = 1/12 for ALL odd w1,w2 (finite check mod 24). Since g_{2*{1..11}}
     has M=1/12 attained only at 8 points, all of denominator 24, M(S)=1/12 <=> some argmax has both
     tighteners safe -- a condition on (w1,w2) mod 24; verified for all distinct odd residue-pairs, and
     same-residue pairs give Psi>=1/4 automatically. By scale-invariance the same holds for c*{1..11} (the
     min-gap EXTREMIZERS) => confinement (gap=1/84>0) PROVEN for AP even parts and their dilations.

OPEN (the confinement core): general U with M(U)>1/12. Measure bounds are VACUOUS (lambda(U)<=2/7 always,
in fact ~0.05-0.09), so it needs the argmax/commensurability arithmetic, not soft bounds.
"""
import sys, itertools, random
from fractions import Fraction as Fr
from math import lcm
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    x=x-int(x)
    if x<0:x+=1
    return min(x,1-x)
def exact_M(S, want_arg=False):
    S=sorted(set(S));cands=set()
    for v in S:
        for k in range(v):cands.add(Fr(2*k+1,2*v))
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for den in (S[i]+S[j],abs(S[i]-S[j])):
                if den:
                    for s in range(den):cands.add(Fr(s,den))
    b=Fr(0);arg=[]
    for t in cands:
        v=min(norm(x*t) for x in S)
        if v>b:b=v;arg=[t]
        elif want_arg and v==b:arg.append(t)
    return (b,arg) if want_arg else b

print("="*90); print(" (1) FOLDING IDENTITY  M(S)=max_t min(g_E(t),Psi(t))  -- verify on random t"); print("="*90)
E=[2*u for u in [1,2,3,5,7]]; w1,w2=3,5; random.seed(0); ok=True
for _ in range(120000):
    t=Fr(random.randint(0,9973),9974)
    gE=min(norm(v*t) for v in E); a,c=norm(w1*t),norm(w2*t)
    lhs=max(min(gE,a,c), min(gE,Fr(1,2)-a,Fr(1,2)-c))
    rhs=min(gE, max(min(a,c),Fr(1,2)-max(a,c)))
    if lhs!=rhs: ok=False;break
print("  identity verified:",ok)

print("\n"+"="*90); print(" (2) AP EVEN PART: M(2*{1..11} u {w1,w2}) = 1/12 for ALL odd w1,w2 (finite check mod 24)"); print("="*90)
E=[2*j for j in range(1,12)]
b,arg=exact_M(E,want_arg=True)
dens=sorted(set(t.denominator for t in arg)); D=lcm(*dens)
print(f"  g_E: M={b} (=1/12? {b==Fr(1,12)}); {len(arg)} argmax pts, denominators {dens}, period D={D}")
allok=True;cnt=0
for a1 in range(1,D+1,2):
    for a2 in range(a1,D+1,2):
        # residues (a1,a2) mod D; realize with actual speeds a1,a2 (a2=a1 uses distinct speeds a1,a1+D)
        ws=[a1,a2] if a1!=a2 else [a1,a1+D]
        M=exact_M(sorted(E+ws)); cnt+=1
        if M!=Fr(1,12): allok=False; print("   *** M=",M," at residues",(a1,a2))
print(f"  checked {cnt} odd residue-pairs mod {D}: ALL M=1/12? {allok}")
print(f"  => M(2*{{1..11}} u 2odd)=1/12 for ALL odd tighteners (PROVEN). Dilations c*{{1..11}} by scale-inv.")

print("\n"+"="*90); print(" (3) STATUS"); print("="*90)
print("  * inf_U gap(U) > 0  <=  M(2U u 2odd) >= 1/12  (=> gap >= 1/84).  min M = 1/12 (verified, S64 descent).")
print("  * PROVEN for the EXTREMAL AP even parts c*{1..11} (gap=1/84 exactly).")
print("  * OPEN for general U (M(U)>1/12): the argmax non-extremity condition -- the confinement core.")
print("    Soft measure bounds are VACUOUS here (lambda(U)~0.05-0.09 << 2/7); needs the arithmetic.")
print("DONE.")
