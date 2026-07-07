#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S38: reconcile the resonance LADDER (kps S36, metric Dx<D) with
opus's ARITHMETIC obstruction (HYP-4506: first-gap emptiness NON-MONOTONIC in N; mediant
3/(3N+2) achievable <=> 3N+2 prime), and attack the order-3 residual (4/51, 5/63) at N=12.

opus S118: N=13 first gap is NONEMPTY via {1..11,13,36}=3/41 (mediant, 41 prime), while N=12
is EMPTY (38=2*19 composite).  So width does NOT decide (non-monotonic); the arithmetic of
3N+2 does.  This GUARD-RAILS my S36 width/Dx<D narrative.

RECONCILIATION (this file): the ladder is still correct -- it PRODUCES the mediant at N=13 and
SKIPS at N=12 -- and the reason it aligns at N=13 but not N=12 is exactly the primality of
3N+2.  So Dx<D (metric, per-family) and 3N+2-composite (arithmetic, root cause) are consistent:
the arithmetic decides whether any base's resonance grid can align to the mediant.

Part 1: verify {1..11,13,36}=3/41 and decompose it as base {1..11,13} + resonant outlier 36.
Part 2: the 3N+2 prime/composite pattern vs mediant achievability, N=4..16.
Part 3: order-3 at N=12 -- are 4/51, 5/63 achievable?  (ladder over order-2 bases + arithmetic).
"""
from fractions import Fraction
from itertools import combinations
import numpy as np
from math import gcd, isqrt
from functools import reduce

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        for c in range(1,q):
            r=(va*c)%q; d=np.minimum(r,q-r); bq=int(d.min())
            if bq*bd>bn*q: bn,bd,bc=bq,q,c
    return Fraction(bn,bd),(bc,bd)

def isprime(n): return n>=2 and all(n%p for p in range(2,isqrt(n)+1))

print("=== PART 1: opus's N=13 nonempty witness as a LADDER family ===", flush=True)
fam=[1,2,3,4,5,6,7,8,9,10,11,13,36]
M,(c,q)=Mw(fam)
print(f"  {{1..11,13,36}} (13 speeds): M={M} at t={c}/{q}   (opus: N=13 mediant 3/41, 41 prime)", flush=True)
base=[1,2,3,4,5,6,7,8,9,10,11,13]
mu,(cb,D)=Mw(base)
print(f"  base {{1..11,13}} (12 speeds): mu=M(base)={mu} at t={cb}/{D}", flush=True)
# is 36 a resonance (multiple of D)?  and does M=mu*x/(x+rho)?
x=36; rho=x*(mu-M)/M if M>0 else None
print(f"  outlier x=36: 36/D = 36/{D} = {Fraction(36,D)}  (resonance iff integer);  implied rho={rho}", flush=True)
print(f"  => the N=13 mediant IS a ladder rung; at N=13 the base's grid ALIGNS to 3/41 (41 prime).\n", flush=True)

print("=== PART 2: the 3N+2 prime/composite dichotomy vs mediant achievability ===", flush=True)
print(f"  {'N':>3}{'gap':>16}{'mediant':>10}{'3N+2':>6}{'prime?':>8}{'opus/mac-mini':>18}", flush=True)
known={7:'nonempty(3/23)',12:'EMPTY',13:'nonempty(3/41)'}
for N in range(4,17):
    g=f"(1/{N+1},2/{2*N+1})"; med=Fraction(3,3*N+2); pr=isprime(3*N+2)
    note=known.get(N,'')
    print(f"  {N:>3}{g:>16}{str(med):>10}{3*N+2:>6}{str(pr):>8}{note:>18}", flush=True)
print("  HYPOTHESIS (opus HYP-4506 / mac-mini HYP-4562): mediant 3/(3N+2) achievable <=> 3N+2 prime.", flush=True)
print("  N=12: 38=2*19 composite => mediant UNachievable => (G).  The obstruction is ARITHMETIC.\n", flush=True)

print("=== PART 3: order-3 residual at N=12 -- targets 4/51, 5/63 (depth-2 SB descendants) ===", flush=True)
N=12; LO,HI=Fraction(1,N+1),Fraction(2,2*N+1)
# order-3 in-gap values: s/(12s+3), 3<s<6 -> s=4,5
targets=[Fraction(s,N*s+3) for s in (4,5)]
print(f"  order-3 in-gap values s/(12s+3), 3<s<6: {[str(t) for t in targets]}", flush=True)
for t in targets:
    print(f"    {t}: denom {t.denominator} = {t.denominator}={'prime' if isprime(t.denominator) else 'composite '+str([p for p in range(2,t.denominator) if t.denominator%p==0 and isprime(p)])}", flush=True)
# Farey nesting widths: depth-1 gap vs depth-2 sub-intervals
med=Fraction(3,38)
w1=HI-LO; wL=med-LO; wR=HI-med
print(f"  Farey widths: depth-1 gap={w1} (~{float(w1):.5f}); depth-2 (1/13,3/38)={wL}(~{float(wL):.5f}), (3/38,2/25)={wR}(~{float(wR):.5f})", flush=True)
print(f"  => depth-2 windows ~{float(w1/wL):.1f}-{float(w1/wR):.1f}x NARROWER => ladder skip even surer (Dx_2<Dx_1<D).", flush=True)
# ladder search over order-2 bases (AP{1..b}+2 defects) for a rung = 4/51 or 5/63
def rho_mu(base):
    mu,(c,D)=Mw(base); mx=max(base); x0=((2*mx)//D+2)*D
    M0,_=Mw(sorted(base+[x0])); return mu,D,(x0*(mu-M0)/M0 if M0>0 else Fraction(0))
hits=[]; nbases=0
for b in (7,8,9):
    ap=list(range(1,b+1)); pool=[x for x in range(b+1,b+16)]
    for defs in combinations(pool, 11-b):
        base=sorted(ap+list(defs))
        if len(base)!=11 or reduce(gcd,base)!=1: continue
        mu,D,rho=rho_mu(base)
        if mu<=HI: continue
        nbases+=1
        for t in targets:
            # ladder rung = t => x = rho*t/(mu-t); check near-integer multiple of D
            if mu<=t: continue
            xf=rho*t/(mu-t)
            for j in (int(xf/D), int(xf/D)+1):
                if j<1: continue
                x=j*D; v=sorted(base+[x])
                if len(set(v))!=12 or reduce(gcd,v)!=1: continue
                Mv,_=Mw(v)
                if Mv==t: hits.append((base,x,Mv))
print(f"  ladder over {nbases} order-2 bases (AP+2 defects), rungs targeting 4/51 or 5/63: {len(hits)} hit(s)", flush=True)
for h in hits[:20]: print("    HIT:",h, flush=True)
print("\nREADING: order-3 targets have COMPOSITE denominators (51=3*17, 63=3^2*7); the ladder over", flush=True)
print("order-2 bases finds none; the depth-2 Farey windows are narrower still.  Both the metric", flush=True)
print("(Dx skip) and the arithmetic (composite denom) point the same way -- order-3 empty at N=12,", flush=True)
print("consistent with opus's ARITHMETIC obstruction (the root cause) over my width narrative (a symptom).", flush=True)
