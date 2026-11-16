#!/usr/bin/env python3
"""Certify both physical sides of the carry, and a third exact zero step."""
from itertools import product
from check_splice import A,B,coefficients,exact,from_coeff,inflate,lookup,profile,validate

def extend(base, masks):
    m=len(base)
    k=len(masks)
    return tuple(r|sum((not(x>>i&1))<<(m+j) for j,x in enumerate(masks))
                 for i,r in enumerate(base))+tuple(
           x|sum(1<<(m+j) for j in range(i+1,k))
           for i,x in enumerate(masks))

def test(a,goal):
    validate(a)
    if a[2]!=8 or a[3]!=16:raise ValueError("run bounds changed")
    c=coefficients(a)
    for sx,sy,X,Y in product((-1,1),(-1,1),(2,4),(2,4)):
        if from_coeff(c,X,Y,sx,sy)!=goal(X,Y):raise ValueError("polynomial mismatch")
    for s,t in ((1,1),(1,2),(2,1),(2,2),(3,3)):
        if exact(*inflate(a,s,t))!=goal(1<<s,1<<t):raise ValueError("DP mismatch")
    print(a,goal(2,2))

for masks in ((189,189),(191,189),(253,189),(253,253),(253,255)):
    test(extend(A,masks),lambda X,Y:profile(profile(lookup(X,Y),0),0))

for base, choices in ((B, ((124,),(190,),(253,60),(253,188),(253,252),(254,255))),
                       (A, ((60,60),(188,60),(188,188)))):
    for masks in choices:
        test(extend(base,masks),lambda X,Y:profile(profile(lookup(X,Y),2),0))

test(extend(extend(B,(124,)),(509,60)),
     lambda X,Y:profile(profile(profile(lookup(X,Y),2),0),0))
