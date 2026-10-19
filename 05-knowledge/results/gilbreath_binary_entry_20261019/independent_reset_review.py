#!/usr/bin/env python3
"""Independent exact review, explicit bounded universes; no Git operations."""
from itertools import product
from fractions import Fraction
from collections import defaultdict
from math import comb
import json

counts=defaultdict(int)
def check(p,kind,detail=None):
    counts[kind]+=1
    if not p: raise RuntimeError((kind,detail))
def evolve(x,u):
    # Imperative independent implementation.
    y=list(x)
    for i in range(len(x)-1): y[i]=abs(x[i]-x[i+1])
    y[-1]=abs(x[-1]-u)
    return tuple(y)
def triangle(row,depth):
    row=list(row)
    for j in range(depth):
        row=[abs(row[k]-row[k+1]) for k in range(len(row)-1)]
    return row

# Rational-grid rigidity: seed values 0,1/4,...,1, implemented by scaling.
for n in range(1,8):
    for x in product(range(5),repeat=n):
        if triangle(x,n-1)[0]==4:
            check(all(a in (0,4) for a in x),'rigidity_endpoint',x)
        counts['rigidity_universe']+=1

# Independent triangle embedding and maximum-drop test, all [0,M]^F.
for M in range(2,6):
    for F in range(1,7):
        if (M+1)**F>50000: continue
        for x in product(range(M+1),repeat=F):
            y=x
            for _ in range(F):y=evolve(y,1)
            z=tuple(triangle(x+(1,)+(0,)*(F-1),F))
            check(y==z,'embedding',(M,F,x,y,z))
            check(max(y)<=M-1,'maximum_drop',(M,F,x,y))

# Full deterministic formula for sparse left defect, including non-dyadic widths.
# Until binary entry, x_0(t)=M - sum_{s<t} [F-1 is a binary submask of s].
for F in range(1,66):
    for M in range(2,26):
        x=(M,)+(0,)*(F-1);hits=0;t=0
        while max(x)>1:
            bit=int((t & (F-1))==F-1)
            check(bit==comb(t,F-1)%2 if t>=F-1 else bit==0,'lucas',(F,t))
            if F>1:
                check(x[1]==bit,'neighbor_trace',(F,M,t,x))
            hits+=bit
            x=evolve(x,1);t+=1
            check(x[0]==M-hits,'sparse_defect_formula',(F,M,t,x,hits))
            check(t<=(M-1)*F,'sparse_bound',(F,M,t))
        if F&(F-1)==0:
            check(t==(M-1)*F,'dyadic_sharpness',(F,M,t))

# Every width-F word beginning with one reduces every state maximum.
# The inverse Pascal tail is reconstructed directly and checked by triangles.
for F in range(1,6):
    for inputs in product((0,1),repeat=F):
        tail=[]
        for t,u in enumerate(inputs):
            v=u
            for j,b in enumerate(tail):
                v^=(comb(t,j)&1)&b
            tail.append(v)
        for t,u in enumerate(inputs):
            check(triangle(tail,t)[0]==u,'inverse_pascal_trace',(F,inputs,tail,t))
        for M in (2,3):
            for x0 in product(range(M+1),repeat=F):
                x=x0
                for u in inputs:x=evolve(x,u)
                y=tuple(triangle(x0+tuple(tail),F))
                check(x==y,'arbitrary_word_embedding',(F,M,x0,inputs,x,y))
                if inputs[0]:
                    check(max(x)<=M-1,'first_bit_drop',(F,M,x0,inputs,x))
            if not inputs[0]:
                x=(M,)+(0,)*(F-1)
                for u in inputs:x=evolve(x,u)
                check(x[0]==M,'first_bit_necessary',(F,M,inputs,x))

# Exact stochastic tail via full-state rational forward evolution.
# Uses independent Bernoulli(q), not a biased spatial sea.
for M in (2,3):
    for F in range(1,4):
        for q in (Fraction(1,3),Fraction(1,2),Fraction(2,3),Fraction(1)):
            r=q
            for x0 in product(range(M+1),repeat=F):
                dist={x0:Fraction(1)}
                for t in range(1,6*F+1):
                    nd=defaultdict(Fraction)
                    for x,w in dist.items():
                        for u,pu in ((0,1-q),(1,q)):
                            if pu: nd[evolve(x,u)]+=w*pu
                    dist=nd
                    if t%F==0:
                        k=t//F
                        tail=sum((w for x,w in dist.items() if max(x)>1),Fraction(0))
                        bound=sum((comb(k,j)*r**j*(1-r)**(k-j)
                                   for j in range(min(k,M-2)+1)),Fraction(0))
                        check(tail<=bound,'negative_binomial_tail',(M,F,q,x0,t,tail,bound))

# Binary F-step channel is bijective in forcing word and start state.
for F in range(1,8):
    states=list(product((0,1),repeat=F))
    for x0 in states:
        destinations=[]
        for inputs in states:
            x=x0
            for u in inputs:x=evolve(x,u)
            destinations.append(x)
        check(len(set(destinations))==2**F,'fair_exact_mixing',(F,x0))
    for inputs in states:
        destinations=[]
        for x0 in states:
            x=x0
            for u in inputs:x=evolve(x,u)
            destinations.append(x)
        check(len(set(destinations))==2**F,'fixed_word_bijection',(F,inputs))

# Essential hostile boundaries.
x=(3,0,0)
for _ in range(100):x=evolve(x,0)
check(x==(3,0,0),'zero_input_nonabsorption')
x=(Fraction(23,10),)
for _ in range(10):x=evolve(x,1)
check(x[0] in (Fraction(3,10),Fraction(7,10)),'noninteger_nonabsorption')
print(json.dumps({'status':'PASS','checks':dict(counts)},indent=2))
