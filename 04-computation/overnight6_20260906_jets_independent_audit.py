#!/usr/bin/env python3
"""Separate Fraction Gauss-Jordan audit of normalized Hasse jet precision.

No primary producer imports, no SymPy, no Smith-form implementation.
The inverse of each literal Hasse matrix supplies all observed denominators.
"""
from fractions import Fraction as F
from math import comb, prod, factorial
import hashlib
import json
import sys

sys.stdout.reconfigure(newline='\n')
CHECKS=0

def check(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:
        raise RuntimeError(label)

def vp(x,p):
    x=F(x)
    if x==0:
        return float('inf')
    n,d=abs(x.numerator),x.denominator
    a=0
    while n%p==0:
        n//=p; a+=1
    while d%p==0:
        d//=p; a-=1
    return a

def inverse(matrix):
    n=len(matrix)
    rows=[[F(x) for x in row]+[F(i==j) for j in range(n)]
          for i,row in enumerate(matrix)]
    for j in range(n):
        k=next(k for k in range(j,n) if rows[k][j])
        rows[j],rows[k]=rows[k],rows[j]
        divisor=rows[j][j]
        rows[j]=[x/divisor for x in rows[j]]
        for i in range(n):
            if i==j or not rows[i][j]:
                continue
            factor=rows[i][j]
            rows[i]=[x-factor*y for x,y in zip(rows[i],rows[j])]
    check(all(rows[i][j]==int(i==j) for i in range(n) for j in range(n)),
          'complete exact Gauss-Jordan identity')
    return [row[n:] for row in rows]

def hasse(nodes,multiplicities):
    M=sum(multiplicities)
    matrix=[[comb(j,r)*x**(j-r) if j>=r else 0 for j in range(M)]
            for x,m in zip(nodes,multiplicities) for r in range(m)]
    return matrix

def multiply(a,b,order):
    return [sum((a[j]*b[k-j] for j in range(k+1)
                 if j<len(a) and k-j<len(b)),F(0)) for k in range(order+1)]

def negbin(alpha,m,order):
    return [(-alpha)**r*comb(m+r-1,r) for r in range(order+1)]

def mod(c,p,k):
    c=F(c)
    check(vp(c,p)>=0,'integral normalized residue')
    return c.numerator*pow(c.denominator,-1,p**k)%(p**k)

def observe(nodes,multiplicities,p):
    M=sum(multiplicities)
    inv=inverse(hasse(nodes,multiplicities))
    full_loss=max(-vp(x,p) for row in inv for x in row)
    data=[]; offset=0
    profile={}
    for i,(x,m) in enumerate(zip(nodes,multiplicities)):
        depth=[vp(x-y,p) for j,y in enumerate(nodes) if j!=i]
        f=max(depth,default=0)
        q0=prod((x-y)**multiplicities[j] for j,y in enumerate(nodes) if j!=i)
        D=vp(q0,p)
        # Top row in the literal inverse is the raw reciprocal coefficient.
        a=[inv[M-1][offset+m-1-l] for l in range(m)]
        b=[a[l]*q0*p**(l*f) for l in range(m)]
        check(b[0]==1 and all(vp(c,p)>=0 for c in b),'inverse gives integral normalized jet')
        shell_product=[F(1)]+[F(0)]*(m-1)
        for h in sorted(set(depth)):
            shell=[F(1)]+[F(0)]*(m-1)
            for j,y in enumerate(nodes):
                if j!=i and vp(x-y,p)==h:
                    shell=multiply(shell,negbin(F(p**h,x-y),multiplicities[j],m-1),m-1)
            attenuated=[c*p**(l*(f-h)) for l,c in enumerate(shell)]
            shell_product=multiply(shell_product,attenuated,m-1)
        check(shell_product==b,'shell product versus full literal inverse')
        powers=[F(0)]+[sum((multiplicities[j]*F(p**f,x-y)**r
                              for j,y in enumerate(nodes) if j!=i),F(0))
                          for r in range(1,m)]
        for l in range(1,m):
            check(l*b[l]==sum((-1)**r*powers[r]*b[l-r] for r in range(1,l+1)),
                  'signed power-moment recursion')
            # Fraction-free Newton recurrence loses at most v_p(l!) digits.
            check(vp(factorial(l)*b[l],p)>=0,'factorial integral moment polynomial target')
        for l,c in enumerate(b):
            if c:
                slope=M-m+l
                intercept=D+l*f-vp(c,p)
                profile[slope]=max(profile.get(slope,-float('inf')),intercept)
        data.append((D,f,b)); offset+=m
    check(len(profile)<=max(multiplicities),'merged slope count')
    Dmax=max(D for D,f,b in data)
    recovered=Dmax
    for D,f,b in data:
        for l,c in enumerate(b):
            k=max(0,D+l*f-Dmax)
            if k:
                residue=mod(c,p,k)
                if residue:
                    recovered=max(recovered,D+l*f-vp(residue,p))
    check(recovered==full_loss,'adaptive finite residue budget versus every inverse entry')
    check(max(profile.values())==full_loss,'profile versus full inverse precision')
    return full_loss,dict(sorted(profile.items())),data

def b2(nodes,x):
    a=[F(2,x-y) for y in nodes if y!=x]
    return F(3,2)*(3*sum(a)**2+sum(t*t for t in a))

def main():
    manifest=[]
    cases=[('uniform-A',(0,1,2),(3,3,3),2),
           ('uniform-B',(0,1,3),(3,3,3),2),
           ('weighted-A',(0,2,1),(2,2,1),2),
           ('weighted-B',(1,3,0),(2,2,1),2),
           ('signed-higher',(-5,1,3),(4,2,3),2),
           ('prime3-heterogeneous',(0,3,12,17),(2,4,1,2),3),
           ('prime5-heterogeneous',(-3,2,7,12),(1,3,2,1),5),
           ('singleton',(-7,),(5,),2)]
    for name,nodes,m,p in cases:
        L,profile,st=observe(nodes,m,p)
        for e in (1,2,4):
            actual,shifted,_=observe(tuple(p**e*x for x in nodes),m,p)
            check(actual==max(intercept+slope*e for slope,intercept in profile.items()),
                  'dilation prediction from base envelope')
            check(shifted=={slope:intercept+slope*e for slope,intercept in profile.items()},
                  'all grouped intercepts under dilation')
        manifest.append([name,L,profile])
    hostile=(0,2,-7,-5,1,7,9)
    L,profile,data=observe(hostile,(3,)*7,2)
    check(L==31,'seven-node full inverse global loss31')
    check([b2(hostile,x) for x in (0,2)]==[F(1452032,33075)]*2,'seven-node exact local coefficient')
    check([vp(b2(hostile,x),2) for x in (0,2)]==[11,11],'seven-node simultaneous local valuations')
    check(max(D+l*f-vp(c,2) for D,f,b in data[:2] for l,c in enumerate(b))==4,
          'seven-node pair-local loss4 versus global31')
    manifest.append(['seven-node',L,profile])
    for h in range(2,8):
        r=2**h-1; spacing=2**(h+2)
        nodes=(0,2)+tuple(1+spacing*j for j in range(r))
        target=6*(r+1)*(3*r+1)
        for x in (0,2):
            actual=b2(nodes,x)
            check(vp(actual,2)==h+2,'exact rational simultaneous family valuation')
            check(vp(actual-target,2)>=h+3,'exact rational family congruence')
            first=-3*sum(F(2,x-y) for y in nodes if y!=x)
            check(vp(first,2)==0,'family first jet remains a unit')
            check(max(3,4-vp(first,2),5-vp(actual,2))==4,'family actual pair-local loss4')
        lower=3*(r-1)*(h+2)
        check(lower>5,'outsider order-zero baseline dominates pair')
        if h==2:
            actual_L,profile,_=observe(nodes,(3,)*len(nodes),2)
            check(actual_L>=lower and actual_L>4,'small family global literal inverse firewall')
            manifest.append(['family-h2',actual_L,profile])
    pa=[F(3)*sum(F(2,-y)**r for y in (2,1)) for r in (1,2)]
    pb=[F(3)*sum(F(2,-y)**r for y in (2,3)) for r in (1,2)]
    check([vp(x,2) for x in pa]==[vp(x,2) for x in pb]==[0,0],
          'power-sum valuations do not retain their sum cancellation')
    check(b2((0,2,1),0)==48 and b2((0,2,3),0)==F(44,3),
          'same moment valuations but different normalized coefficient values')
    check(vp(F(48),2)==4 and vp(F(44,3),2)==2,'local valuation-only jet state fails')
    check(max(6*4+3,7*4+4,8*4+1)==33 and max(6*4+3,7*4+4)==32,
          'pruned-at-base order-two slope returns after dilation')
    # Same reduced numerator residue can hide a digit after nonunit division.
    check(0%2==2%2 and (0//2)%2!=(2//2)%2,'modular Newton division needs extra precision')
    print('INDEPENDENT_LITERAL_INVERSE_PROFILES',json.dumps(manifest,sort_keys=True))
    print('FAMILY h2..7 exact rational pair valuations=h+2; pair-local loss=4; global outsider baseline dominates')
    print('BOUNDARIES current-depth pruning is not permanent; power valuations lose unit cancellation; nonunit division loses precision')
    print('PASS',CHECKS,'explicit gates; standard-library only')
    print('semantic_sha256',hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest())

if __name__=='__main__':
    main()
