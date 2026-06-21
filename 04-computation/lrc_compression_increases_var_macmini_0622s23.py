#!/usr/bin/env python3
"""POTENTIAL CRACK for the FULL concentration extremality: does COMPRESSION (move the largest runner
to the smallest unused position) ALWAYS increase Var(N) -- even off a WIDE base? If yes, any config
compresses to consec with Var monotone-up => consec is the global Var-max (=> gK8 concentration PROVED).
Test: random configs, repeatedly compress, check Var strictly increases at every step."""
import sys, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(11)
def sector_of(p): return int((p%1)*7)
def varN(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        t=7-len(set(sector_of(e*((x0+x1)/2)) for e in E))
        if t<=6: q[t]+=x1-x0
    EN=sum(t*q[t] for t in range(7)); EN2=sum(t*t*q[t] for t in range(7))
    return EN2-EN*EN
def compress_step(E):
    """move the largest element to the smallest non-present nonneg integer (compress toward {0..k-1})."""
    E=sorted(E); k=len(E)
    target=set(range(k))  # consec target
    missing=sorted(target-set(E))
    extra=sorted(set(E)-target, reverse=True)  # elements above the consec block
    if not extra or not missing: return None
    E2=set(E); E2.discard(extra[0]); E2.add(missing[0])
    return tuple(sorted(E2))
print("Does each compression step strictly INCREASE Var(N)? (testing the monotone path to consec)")
bad=0; tot_paths=0; bad_examples=[]
for trial in range(80):
    k=random.choice([8,9,10])
    E=tuple(sorted([0]+random.sample(range(1,50),k-1)))
    tot_paths+=1
    path_bad=0
    cur=E; steps=0
    while steps<40:
        v0=varN(cur); nxt=compress_step(cur)
        if nxt is None: break
        v1=varN(nxt)
        if v1 < v0 - F(1,10**9):  # Var DECREASED on a compression step
            path_bad+=1
            if len(bad_examples)<5: bad_examples.append((cur,float(v0),nxt,float(v1)))
        cur=nxt; steps+=1
    if path_bad>0: bad+=1
print(f"  paths with a Var-DECREASING compression step: {bad}/{tot_paths}")
if bad_examples:
    print("  counterexamples (compression LOWERED Var -- monotonicity FALSE):")
    for (a,va,b,vb) in bad_examples: print(f"    {a} Var={va:.4f} -> {b} Var={vb:.4f}")
else:
    print("  => COMPRESSION MONOTONICITY HOLDS: every compression step raises Var => consec is global Var-max")
    print("     (would PROVE the gK8 concentration extremality by the compression path)")
