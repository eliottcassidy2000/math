#!/usr/bin/env python3
"""Finite exact hostile/positive controls; no infinite orbit assertion."""
from fractions import Fraction
import json
from pathlib import Path

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def step(b,n):
    x=b*n+1
    a=(x & -x).bit_length()-1
    return x>>a,a

counts={'affine_checks':0,'barrier_checks':0,'inverse_prefix_checks':0}
for b in (3,5,7):
    for n0 in range(1,2002,2):
        n=n0; S=0; R=Fraction(0); descended=False
        for k in range(1,61):
            R+=Fraction(2**S,b**k)
            n,a=step(b,n); S+=a
            descended |= n<n0
            need(Fraction(b**k,2**S)*(n0+R)==n,'affine identity')
            counts['affine_checks']+=1
            need((n0+R).numerator % (2**S)==0,'2-adic prefix source guard')
            counts['inverse_prefix_checks']+=1
            if 2**S*n0**k>(b*n0+1)**k:
                need(descended,'source barrier')
                counts['barrier_checks']+=1

# Finite words really select one odd source class modulo 2^(S+1).
words=[(1,), (2,), (1,1,5), (1,4), (2,2,2), (1,3,2)]
for b in (3,5,7):
    for word in words:
        S=sum(word); modulus=2**(S+1)
        hits=[]
        for n0 in range(1,modulus,2):
            n=n0; found=[]
            for _ in word:
                n,a=step(b,n); found.append(a)
            if tuple(found)==word:
                hits.append(n0)
        need(len(hits)==1,'unique cylinder')
        R=Fraction(0); s=0
        for j,a in enumerate(word,1):
            R+=Fraction(2**s,b**j); s+=a
        # Exact valuation includes odd residual at depth S, not just divisibility.
        need((hits[0]+R.numerator*pow(R.denominator,-1,modulus))%modulus==2**S,'odd residual guard')

for start,word in [(1,(1,4)),(13,(1,1,5))]:
    n=start; vals=[]
    for _ in word:
        n,a=step(5,n); vals.append(a)
    need(n==start and tuple(vals)==word,'5-cycle hostile')

# Pure clock slope averages alone do not mean endpoints tend to zero.
n=1
for _ in range(100):
    n,a=step(3,n)
    need(n==1 and a==2,'3 fixed-point carry control')

out={'status':'PASS','universe':'b in {3,5,7}; 1001 positive odd sources <=2001; first60 odd steps; 18 finite cylinder words; both stated5 cycles',**counts,'cylinder_checks':18,'limitations':'finite checks supplement proofs; no global termination conclusion'}
Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,sort_keys=True))
