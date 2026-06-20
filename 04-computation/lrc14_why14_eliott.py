#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify the structural 'why N=14, depth 4' reasoning for general LRC(2p).
Claim: tower = {2^a : 2^a < N} = {1,2,4,8} for N=14 (depth floor(log2 N)=4),
and the mouth hole of the (2p-2)-core avoids all powers of 2 -- so the tower is
the carrier specifically because the mouth-minimizer's hole pattern dodges 2^a."""
from fractions import Fraction
import itertools, math
exec(open('04-computation/lrc14_adversarial_verify_eliott.py').read().split('# ---- DECISIVE')[0].split('if __name__')[0])

def lonely_measure_N(C, N):
    theta = Fraction(1,N)
    bad=[]
    for d in C:
        if d==0: continue
        w=theta/d
        for m in range(0,d+1):
            lo=Fraction(m,d)-w; hi=Fraction(m,d)+w
            for s in (-1,0,1):
                a=max(lo+s,Fraction(0)); b=min(hi+s,Fraction(1))
                if a<b: bad.append((a,b))
    bad.sort(); cur=Fraction(-1); un=[]
    for a,b in bad:
        if a>cur: un.append([a,b]); cur=b
        elif b>cur: un[-1][1]=b; cur=b
    cov=sum(b-a for a,b in un)
    return Fraction(1)-cov

print("N=14: floor(log2 14)=", math.floor(math.log2(14)), " tower {2^a<14}=", [2**a for a in range(4) if 2**a<14])
print("CRT split 1/14 = (1/2)(1/7), apex prime 7 = N/2.")
print()
# For general LRC(2p): full core = {1..2p-1}, size 2p-1. We need (2p-2)-core (drop one),
# find which single-hole core is the global mouth (min lonely measure), and check its
# hole vs powers of two.
print("== For each N=2p, find the single-hole (N-2)-core mouth and its hole ==")
for p in [3,5,7,11,13]:
    N=2*p
    full=list(range(1,N))  # {1..2p-1}, size 2p-2 ... wait core size is N-2=2p-2
    # core has size N-2; full {1..N-1} has size N-2. Single-hole means we drop from a larger?
    # Match the LRC(14) setup: core size = N-2 = 12 for N=14, from {1..13}=size 13, drop 1.
    universe=list(range(1,N))  # size N-1
    best=None
    for hole in universe:
        C=tuple(d for d in universe if d!=hole)
        if reduce(gcd,C)!=1: continue
        L=lonely_measure_N(C,N)
        if best is None or L<best[0]: best=(L,hole,C)
    L,hole,C=best
    pw2 = (hole & (hole-1))==0  # hole is power of 2?
    print(f"  N={N:3d} (p={p:2d}): mouth single-hole = drop-{hole}, meas={float(L):.6f}, "
          f"hole_is_pow2={pw2}, tower={[2**a for a in range(N) if 2**a<N]}")
