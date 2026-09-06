#!/usr/bin/env python3
"""Exact all-prime (2,2,1) Smith classification and intrinsic dyadic bit.

No producer imports. Standard-library Fraction DVR Smith path, independent
all-minor gcd path, and a direct inverse-precision path. Reproducible finite
universes are printed, and no assert statement or floating point is used.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import gcd, prod
import json
import random
import sys
sys.stdout.reconfigure(newline='\n')
GATES=0


def require(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def det(v,w):
    return v[0]*w[1]-v[1]*w[0]


def tangent(v):
    a,b=v
    x,y,s,t,r,q=1,0,0,1,a,b
    while q:
        k=r//q
        r,q,x,s,y,t=q,r-k*q,s,x-k*s,t,y-k*t
    require(abs(r)==1,'primitive vector')
    ans=(-y*r,x*r)
    require(det(v,ans)==1,'unimodular tangent')
    return ans


def vp(x,p):
    if x==0:
        return None
    x=F(x)
    a,b,ans=abs(x.numerator),x.denominator,0
    while a%p==0:
        a//=p
        ans+=1
    while b%p==0:
        b//=p
        ans-=1
    return ans


def matrix(vectors,shifts=(0,0,0)):
    rows=[]
    for i,(a,b) in enumerate(vectors):
        w=tangent((a,b))
        c,d=w[0]+shifts[i]*a,w[1]+shifts[i]*b
        rows.append([a**j*b**(4-j) for j in range(5)])
        if i<2:
            rows.append([(j*c*a**(j-1)*b**(4-j) if j else 0)
                         +((4-j)*d*a**j*b**(3-j) if j<4 else 0)
                         for j in range(5)])
    return rows


def smith(a,p):
    a=[[F(x) for x in row] for row in a]
    ans=[]
    for k in range(len(a)):
        i,j=min(((i,j) for i in range(k,len(a)) for j in range(k,len(a)) if a[i][j]),
                key=lambda ij:vp(a[ij[0]][ij[1]],p))
        a[k],a[i]=a[i],a[k]
        for row in a:
            row[k],row[j]=row[j],row[k]
        pivot=a[k][k]
        ans.append(vp(pivot,p))
        for i in range(k+1,len(a)):
            ratio=a[i][k]/pivot
            require(not ratio or vp(ratio,p)>=0,'integral DVR elimination')
            for j in range(k,len(a)):
                a[i][j]-=ratio*a[k][j]
    require(ans==sorted(ans),'Smith order')
    return tuple(ans)


def bareiss(matrix):
    a=[row[:] for row in matrix]
    n=len(a)
    if n==0:
        return 1
    previous,sign=1,1
    for k in range(n-1):
        i=next((i for i in range(k,n) if a[i][k]),None)
        if i is None:
            return 0
        if i!=k:
            a[k],a[i]=a[i],a[k]
            sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=pivot*a[i][j]-a[i][k]*a[k][j]
                require(numerator%previous==0,'Bareiss exact division')
                a[i][j]=numerator//previous
            a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]


def inverse(a):
    n=len(a)
    a=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for k in range(n):
        i=next(i for i in range(k,n) if a[i][k])
        a[k],a[i]=a[i],a[k]
        pivot=a[k][k]
        a[k]=[x/pivot for x in a[k]]
        for i in range(n):
            if i!=k:
                ratio=a[i][k]
                a[i]=[x-ratio*y for x,y in zip(a[i],a[k])]
    return [row[n:] for row in a]


def bit(vectors,reference=None):
    v1,v2,v3=vectors
    w=reference or tangent(v1)
    require(all(det(v,w)%2 for v in vectors),'unit-separated reference')
    tau=F(det(v1,v2)*det(v3,w),2*det(v1,v3)*det(v2,w))
    require(tau.numerator%2 and tau.denominator%2,'unit half cross ratio')
    residue=tau.numerator*pow(tau.denominator,-1,4)%4
    return int(residue==1)


def predicted(vectors,p):
    v1,v2,v3=vectors
    A,B,C=vp(det(v1,v2),p),vp(det(v1,v3),p),vp(det(v2,v3),p)
    e,f=min(A,B,C),max(A,B,C)
    if A==e:
        delta=min(2*e,e+(p==2))
        return (0,0,delta,4*e-delta,2*e+2*f)
    require(A==f and B==C==e,'unique deepest doubled pair')
    if p!=2:
        delta=min(f,2*e)
        return (0,0,delta,f+3*e-delta,3*f+e)
    d=f-e
    if d>=2:
        delta=min(2*e,e+d+1)
        return (0,0,delta,4*e+d+1-delta,4*e+3*d-1)
    if e==0:
        return (0,0,0,2,2)
    if e==1:
        return (0,0,2,5,5)
    epsilon=bit(vectors)
    return (0,0,e+2,3*e+2-epsilon,4*e+epsilon)


def affine_loss(nodes,p):
    # Independent direct specialization of reciprocal inverse coefficients.
    x,y,z=nodes
    a,b=y-x,z-x
    A,B,C=vp(a,p),vp(b,p),vp(b-a,p)
    candidates=[2*A+B,2*A+C,2*B+2*C]
    if a+2*b:
        candidates.append(3*A+2*B-vp(a+2*b,p))
    if 3*a-2*b:
        candidates.append(3*A+2*C-vp(3*a-2*b,p))
    return max(candidates)


def apply(g,v):
    return (g[0][0]*v[0]+g[0][1]*v[1],g[1][0]*v[0]+g[1][1]*v[1])


def main():
    cases=[]
    primes=(2,3,5,7)
    for x,y in combinations(range(-6,7),2):
        for z in range(-6,7):
            if z not in (x,y):
                vectors=((x,1),(y,1),(z,1))
                cases.extend((p,vectors) for p in primes)
    head_count=len(cases)
    for p in primes:
        units=[x for x in (-3,-1,1,2,3,5) if x%p]
        for e,d,u,v in product(range(6),range(1,5),units,units):
            a,b=p**(e+d)*u,p**e*v
            cases.append((p,((0,1),(a,1),(b,1))))
            cases.append((p,((0,1),(b,1),(a,1))))
    deep_count=len(cases)-head_count
    pool=[(a,b) for a in range(-4,5) for b in range(5) if gcd(a,b)==1 and (b>0 or a>0)]
    rng=random.Random(20260906221)
    for _ in range(100):
        vectors=tuple(rng.sample(pool,3))
        cases.extend((p,vectors) for p in primes)
    cases.extend((p,((1,0),(0,1),(1,1))) for p in primes)
    trace=[]
    for p,vectors in cases:
        a=matrix(vectors)
        actual=smith(a,p)
        require(actual==predicted(vectors,p),'complete classified Smith list')
        mass=4*vp(det(vectors[0],vectors[1]),p)+2*vp(det(vectors[0],vectors[2]),p)+2*vp(det(vectors[1],vectors[2]),p)
        require(sum(actual)==mass,'weighted determinant mass')
        if all(v[1]==1 for v in vectors):
            require(actual[-1]==affine_loss([v[0] for v in vectors],p),'independent reciprocal precision')
        trace.append((p,vectors,actual))
    symmetry_rows=0
    for e in range(2,8):
        for u,v in product((-3,-1,1,3,5,7),repeat=2):
            vectors=((0,1),(2**(e+1)*u,1),(2**e*v,1))
            original=bit(vectors)
            for g in (((0,-1),(1,0)),((1,2),(1,3)),((1,0),(0,-1))):
                transformed=tuple(apply(g,v) for v in vectors)
                require(bit(transformed)==original,'projective bit invariance')
                require(bit((transformed[1],transformed[0],transformed[2]))==original,'doubled endpoint symmetry')
                require(smith(matrix(transformed,(3,-2,0)),2)==predicted(vectors,2),'projective full module and tangent shifts')
                reference=tangent(transformed[0])
                for k in (-3,-1,1,4):
                    alternative=tuple(reference[j]+k*transformed[0][j] for j in (0,1))
                    require(bit(transformed,alternative)==original,'reference independence')
                symmetry_rows+=1
    all_minor_count=0
    controls=[((0,1),(8,1),(4,1)),((0,1),(8,1),(-4,1)),
              ((1,0),(1,9),(1,3)),((0,1),(27,1),(3,1)),
              ((1,0),(0,1),(1,1)),((0,1),(4,1),(1,1))]
    for vectors in controls:
        a=matrix(vectors)
        divisors=[1]
        for k in range(1,6):
            common=0
            for rows in combinations(range(5),k):
                for cols in combinations(range(5),k):
                    common=gcd(common,bareiss([[a[i][j] for j in cols] for i in rows]))
                    all_minor_count+=1
            divisors.append(common)
        require(divisors[-1]==abs(det(vectors[0],vectors[1])**4*det(vectors[0],vectors[2])**2*det(vectors[1],vectors[2])**2),'literal determinant identity')
        inv=inverse(a)
        for p in primes:
            actual=tuple(vp(divisors[k],p)-vp(divisors[k-1],p) for k in range(1,6))
            require(actual==predicted(vectors,p),'independent full determinantal ideals')
            require(max(-vp(x,p) for row in inv for x in row if x)==actual[-1],'literal inverse denominator')
    for e in range(2,15):
        one=predicted(((0,1),(2**(e+1),1),(2**e,1)),2)
        zero=predicted(((0,1),(2**(e+1),1),(-2**e,1)),2)
        for N in range(1,4*e+3):
            difference=sum(min(N,x) for x in zero)-sum(min(N,x) for x in one)
            require(difference==int(3*e+2<=N<=4*e),'exact full-kernel ratio interval')
    # At shallow depth the raw half-cross-ratio residue can change chart,
    # and the full Smith list masks it; this is an explicit sharp hostile.
    shallow=((0,1),(4,1),(2,1))
    changed=tuple(apply(((1,0),(1,1)),v) for v in shallow)
    require(bit(shallow)!=bit(changed),'shallow coordinate bit is not intrinsic')
    require(smith(matrix(shallow),2)==smith(matrix(changed),2)==(0,0,2,5,5),'shallow module masks bit')
    print('SCOPE: full (2,2,1) Smith classification at all primes on three primitive projective directions.')
    print('Head universe: doubled x<y in -6..6, simple z in -6..6 distinct, primes 2,3,5,7; rows=',head_count)
    print('Deep universe: e=0..5,d=1..4, unit lifts from (-3,-1,1,2,3,5), doubled or mixed closest pair; rows=',deep_count)
    print('Other controls: 100 seeded primitive-direction triples; every dyadic projective residue class.')
    print('Complete Smith/metric rows:',len(cases),'projective/reference symmetry rows:',symmetry_rows)
    print('Independent all-minor determinants:',all_minor_count)
    print('Hostile e=2 epsilon=1: (0,0,4,7,9); epsilon=0: (0,0,4,8,8).')
    print('Shallow e=1 chart flips raw bit but both lists are (0,0,2,5,5).')
    print('Full-kernel ratio is exactly 2 for epsilon0/epsilon1 at 3e+2<=N<=4e, and 1 otherwise.')
    print('semantic_sha256:',sha256(json.dumps(trace,separators=(',',':')).encode()).hexdigest())
    print('Exact gates:',GATES)
    print('PASS')


if __name__=='__main__':
    main()
