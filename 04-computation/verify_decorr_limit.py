#!/usr/bin/env python3
"""
Confirm P_r(B) is the genuine decorrelated LIMIT of p0(B u F).

A single fixed far set F=[23,41] still has residual arithmetic resonance with
B (its sector pattern as x varies is not perfectly uniform/independent), so
p0(B u F) only approximates P_r(B).  The decorrelated boundary value is the
limit when the far runners' sector phases become equidistributed & independent.

We confirm this two ways:
 (A) AVERAGE p0(B u F) over many far sets F whose elements are large and
     mutually dissociated; the average must -> P_r(B).
 (B) Replace each far runner by a genuinely uniform independent sector offset
     (phase shift) and average -> P_r(B) exactly in expectation.
"""
import random
from math import comb
from fractions import Fraction as Fr

DEN = 7
def sector(e, x):
    return int(((e*x) % 1.0)*DEN) % DEN

def c_t(t, r):
    return sum((-1)**i*comb(t,i)*(1-i/7)**r for i in range(t+1))

def prof_exact(B):
    # exact interval-scan profile
    bps = {Fr(0), Fr(1)}
    for e in B:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(Fr(m,7*e))
    bps=sorted(bps); prof=[Fr(0)]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        hit=len({sector(e,xm) for e in B})
        prof[7-hit]+=x1-x0
    return prof

def P_r_exact(B,r):
    p=prof_exact(B); return sum(p[t]*c_t(t,r) for t in range(7))

def p0_mc(E, N=300000):
    """measure{ x : E covers all 7 sectors } by MC."""
    g=0
    for _ in range(N):
        x=random.random()
        if len({sector(e,x) for e in E})==7: g+=1
    return g/N

random.seed(2024)

for B in [(0,1,2,3,4,5,6,7), (0,1,2,3,4,5,6)]:
    r=2
    target=float(P_r_exact(B,r))
    print(f"\nB={B}  P_r(B) target = {target:.6f}")

    # (A) average p0(B u F) over many random large dissociated F
    Ntrials=200
    acc=0.0
    for _ in range(Ntrials):
        # pick r large primes-ish numbers, far apart, mutually coprime-ish
        F=[]
        while len(F)<r:
            v=random.randrange(500,5000)
            if all(abs(v-f)>50 and v!=f for f in F):
                F.append(v)
        acc+=p0_mc(list(B)+F, N=20000)
    avgA=acc/Ntrials
    print(f"   (A) avg p0(B u F) over {Ntrials} random far sets = {avgA:.6f}   diff={avgA-target:+.4f}")

    # (B) far runners as independent uniform sector phases (theta uniform in [0,1))
    N=2_000_000; g=0
    for _ in range(N):
        x=random.random()
        S={sector(e,x) for e in B}
        for _ in range(r):
            S.add(random.randrange(7))   # uniform independent sector
        if len(S)==7: g+=1
    print(f"   (B) decorrelated-phase MC                  = {g/N:.6f}   diff={g/N-target:+.4f}")
