#!/usr/bin/env python3
"""Exact character orthogonality guards via integer polynomial reduction."""
from itertools import product
import json

def check(x,msg):
    if not x: raise RuntimeError(msg)
def zero_at_25_root(coeff):
    a=coeff[:]
    for i in range(len(a)-1,19,-1):
        c=a[i]
        for j in (0,5,10,15,20):a[i-20+j]-=c
    return not any(a)
# Uniform full lifts of C={1,3} under Z/25 -> Z/5.
S=[x for x in range(25) if x%5 in (1,3)]
support=[]
for u in range(25):
    coeff=[0]*25
    for x in S:coeff[u*x%25]+=1
    if not zero_at_25_root(coeff):support.append(u)
check(support==[0,5,10,15,20],'uniform-lift annihilator')
# Affine squareclass code for (-1,2,-2).
colors=[(1,z,1^z) for z in range(2)]
for u in product(range(2),repeat=3):
    actual=sum((-1)**sum(a*b for a,b in zip(u,y)) for y in colors)
    expected=2*(-1)**(u[0]+u[2]) if (u[1]+u[2])%2==0 else 0
    check(actual==expected,'binary affine-code Fourier')
# Exact graph law y=3n+2 modulo7. Exponent counts are balanced precisely off line.
line=[]
for u,v in product(range(7),repeat=2):
    c=[0]*7
    for n in range(7):c[(u*n+v*(3*n+2))%7]+=1
    if (u+3*v)%7:
        check(c==[1]*7,'off-line cancellation')
    else:
        check(c[(2*v)%7]==7 and sum(c)==7,'on-line phase')
        line.append([u,v])
check(len(line)==7,'exact Fourier support size')
print(json.dumps({'status':'PASS','lift_support':support,'affine_graph_support':line,'arithmetic':'exact integer coefficient counts'},indent=2))
