#!/usr/bin/env python3
"""Independent exact cross-checks. No asserts; all active under -O."""
from fractions import Fraction as F
from itertools import product
import json

def need(p,label):
    if not p: raise RuntimeError(label)
def valuation(z):
    need(z!=0,'punctured-domain zero')
    return (abs(z)&-abs(z)).bit_length()-1
def step(b,x):
    z=b*x+1
    if isinstance(z,F): a=valuation(z.numerator)-valuation(z.denominator)
    else: a=valuation(z)
    return F(z,2**a),a
def replay(b,x,w):
    for a in w:
        x,actual=step(b,x);need(actual==a,'exact cylinder')
    return x
def source(b,w,terminal):
    for a in reversed(w): terminal=(2**a*terminal-1)/b
    return terminal

def res(x,M): return x.numerator*pow(x.denominator,-1,M)%M
counts={'exhaustive_cylinders':0,'denominator_prehistories':0,'controller_runs':0,'barrier_triggers':0}
# Independently enumerate all odd residues, not just formula-generated sources.
for b in (3,5,9):
    for L in range(1,4):
        for w in product((1,2,3),repeat=L):
            M=2**(sum(w)+1); hits=[]
            for n in range(1,M,2):
                x=F(n); digits=[]
                for _ in w: x,a=step(b,x);digits.append(a)
                if tuple(digits)==w: hits.append(n)
            x=source(b,w,F(1));need(hits==[res(x,M)],'exhaustive unique cylinder')
            # Half-size modulus generally loses the exact odd terminal guard.
            n=hits[0]+M//2;x=F(n);digits=[]
            for _ in w: x,a=step(b,x);digits.append(a)
            need(tuple(digits)!=w,'D rather than D+1 must fail')
            counts['exhaustive_cylinders']+=1
# Include denominator factors shared with b, negative sources, and deep prefixes.
for b in (3,5,9):
    for x in (F(19,101),F(-1,3),F(5,9),F(7,15)):
        for L in range(5):
            for w in product((1,2,7),repeat=L):
                y=source(b,w,x)
                need(y.denominator%x.denominator==0,'denominator persistence')
                need(replay(b,y,w)==x,'inverse exact rational word')
                counts['denominator_prehistories']+=1
# Arbitrary-queue finite experiment using native state repetition only.
for b,r in ((3,3),(5,5)):
    for n0 in range(1,1000,2):
        for q0 in (1,2,3,7,19):
            n=F(n0);q=q0;seen=set();entrance=False
            while q and (n,q) not in seen:
                need(len(seen)<10000,'experiment cap')
                seen.add((n,q));entrance |= n==83
                n,a=step(b,n);q+=2 if a%r==0 else -1
            infinite=bool(q)
            need(infinite==(b==5 and entrance),'arbitrary queue trap classification')
            counts['controller_runs']+=1
counterexample=None
for b in (3,5,7):
    for n0 in range(1,502,2):
        n=F(n0);S=0;descended=False
        for k in range(1,61):
            n,a=step(b,n);S+=a;descended |= n<n0
            if 2**S*n0**k>(b*n0+1)**k:
                need(descended,'barrier prior descent')
                counts['barrier_triggers']+=1
                if n>=n0 and counterexample is None:
                    counterexample={'b':b,'n0':n0,'k':k,'endpoint':int(n),'S':S}
need(counterexample is not None,'endpoint qualification witness')
print(json.dumps({'status':'PASS','counts':counts,'prior_not_endpoint_descent':counterexample,'scope':'Finite checks only; infinite claims separately proof-audited.'},indent=2))
