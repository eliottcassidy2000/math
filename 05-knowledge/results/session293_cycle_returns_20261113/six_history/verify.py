#!/usr/bin/env python3
"""Exact, dependency-free six-mode return controls."""
from fractions import Fraction as F
from itertools import permutations

def require(test):
    if not test: raise AssertionError("exact control failed")

def mm(a,b): return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*b)) for row in a)
def eye(n): return tuple(tuple(F(i==j) for j in range(n)) for i in range(n))
def power(a,k):
    r=eye(len(a))
    for _ in range(k): r=mm(a,r)
    return r

def scale(a,c): return tuple(tuple(c*x for x in row) for row in a)
A=((F(0),F(-3)),(F(1),F(3)))
C=((F(0),F(0),F(3)),(F(0),F(-1),F(-3)),(F(1,3),F(2),F(3)))
v=((3,0,0),(0,0,1),(3,-3,3),(9,-6,4),(12,-6,3),(9,-3,1))
D=tuple(zip(*v))
P=tuple(tuple(F(i==(j+1)%6) for j in range(6)) for i in range(6))
require(power(A,6)==scale(eye(2),-27))
require(power(C,6)==eye(3))
require(all(power(C,k)!=eye(3) for k in range(1,6)))
require(mm(C,D)==mm(D,P))
require(mm(((1,3,3),),C)==((F(1),F(3),F(3)),))
require(mm(((1,3,3),),D)==((3,)*6,))
require(all(u*w==b*b and u>=0 and w>=0 for u,b,w in v))
require(power(P,6)==eye(6))
# A 3x3 minor (columns 0,1,2) proves rank3.
a,b,c=v[0],v[1],v[2]
det=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
require(det==9)
# Exhaust all permutations through6: primitive sixth mode iff a cycle length
# is divisible by6; this is deliberately not the same as permutation order6.
for n in range(1,7):
    with_mode=0; order6=0
    from math import lcm
    for p in permutations(range(n)):
        seen=set(); ls=[]
        for i in range(n):
            if i in seen: continue
            k=i; l=0
            while k not in seen: seen.add(k); l+=1; k=p[k]
            ls.append(l)
        with_mode+=any(l%6==0 for l in ls)
        order6+=lcm(*ls)==6
    require(with_mode==(120 if n==6 else 0))
    print('n',n,'order6',order6,'primitive_sixth_mode',with_mode)
# Quadratic orbit supplied by historical theorem.
f=lambda x:x*x-F(29,16)
q=[F(-7,4),F(5,4),F(-1,4)]
require(all(f(q[i])==q[(i+1)%3] for i in range(3)))
print('PASS: inherited quadratic3 versus new Möbius6; exact Sym² mass return; minimal positive mode controls')
require(tuple(v[3][j]+2*v[1][j] for j in range(3))==tuple(2*v[2][j]+v[0][j] for j in range(3)))
require(tuple(v[0][j]+v[3][j] for j in range(3))==tuple(v[1][j]+v[4][j] for j in range(3))==tuple(v[2][j]+v[5][j] for j in range(3)))
print('PASS: balanced nonnegative token exchanges preserve full three-mode observer')
