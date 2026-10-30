#!/usr/bin/env python3
"""Standard-library-only independent regular-lift checks; no producer imports."""
from fractions import Fraction as F
from itertools import product
from math import comb
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def clean(a):return {k:F(v) for k,v in a.items() if v}
def add(a,b):
    d=dict(a)
    for k,v in b.items():d[k]=d.get(k,F(0))+v
    return clean(d)
def scale(a,c):return clean({k:v*c for k,v in a.items()})
def mul(a,b):
    d={}
    for k,v in a.items():
        for l,w in b.items():
            q=tuple(ki+li for ki,li in zip(k,l));d[q]=d.get(q,F(0))+v*w
    return clean(d)
one={(0,0,0):F(1)}
def power(a,n):
    r=one
    for _ in range(n):r=mul(r,a)
    return r
def diff(a,j):
    d={}
    for k,v in a.items():
        if k[j]:
            q=list(k);q[j]-=1;d[tuple(q)]=v*k[j]
    return clean(d)
def bracket(a,b):return add(mul(diff(a,1),diff(b,2)),scale(mul(diff(a,2),diff(b,1)),-1))
identities=0
for m in range(2,11):
    d=m-1
    for s in range(1,d+1):
        if d%s:continue
        e=d//s
        for c,beta in [(2,3),(-3,-1)]:
            U={(1,1,m):F(1),(0,0,s):F(c)};P=add(U,scale(one,beta))
            Q={(k-1,k,(m-s)*k):F(-comb(e,k)*c**(e-k),d) for k in range(1,e+1)}
            for extra in range(3):
                h=power(U,e+extra);q=mul(Q,power(U,extra))
                need(bracket(P,q)==h,'regular generator and target multiple');identities+=1
    P={(1,1,m):F(1),(0,0,0):F(5)};U={(1,1,m):F(1)};Q={(0,1,1):F(-1,d)}
    for extra in range(3):
        need(bracket(P,mul(Q,power(U,extra)))==power(U,extra+1),'constant B generator');identities+=1

def composition_membership_y_power(B,d):
    # Decide whether y^d is a polynomial in B by cancelling leading terms.
    n=max(k[2] for k in B);lc=B[(0,0,n)];remainder={(0,0,d):F(1)}
    while remainder:
        j=max(k[2] for k in remainder)
        if j==0:return True
        if j%n:return False
        e=j//n;c=remainder[(0,0,j)]/lc**e
        remainder=add(remainder,scale(power(B,e),-c))
    return True
composition_cases=0;admitted=0
for n in range(1,5):
    for low in product((-1,0,1),repeat=n):
        for lc in (-1,1):
            B=clean({(0,0,j):v for j,v in enumerate(low+(lc,))})
            nonconstant_support=[k[2] for k in B if k[2]>0]
            for d in range(1,9):
                predicted=len(nonconstant_support)==1 and d%n==0
                got=composition_membership_y_power(B,d)
                need(got==predicted,'composition criterion')
                composition_cases+=1;admitted+=got
print(json.dumps({'status':'PASS','arithmetic':'fractions.Fraction','producer_imports':False,'regular_lift_identities':identities,'composition_cases':composition_cases,'composition_admitted':admitted},indent=2,sort_keys=True))
