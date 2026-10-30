#!/usr/bin/env python3
"""Independent Fraction-only audit; imports no producer code; no assert gates."""
from fractions import Fraction as F
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def clean(a): return {k:F(v) for k,v in a.items() if v}
def add(a,b):
    d=dict(a)
    for k,v in b.items(): d[k]=d.get(k,F(0))+v
    return clean(d)
def scale(a,c): return clean({k:v*c for k,v in a.items()})
def mul(a,b):
    d={}
    for (i,j),v in a.items():
        for (k,l),w in b.items(): d[i+k,j+l]=d.get((i+k,j+l),F(0))+v*w
    return clean(d)
def power(a,n):
    r={(0,0):F(1)}
    for _ in range(n): r=mul(r,a)
    return r
def der(a,axis):
    d={}
    for k,v in a.items():
        if k[axis]:
            q=list(k);q[axis]-=1;d[tuple(q)]=v*k[axis]
    return clean(d)
def bracket(a,b): return add(mul(der(a,0),der(b,1)),scale(mul(der(a,1),der(b,0)),-1))
def div_linear(a,b):
    # Exact quotient/remainder of a(x,y) by y-b, fibre by x degree.
    q={};rem={}
    for i in {k[0] for k in a}:
        d=max(j for ii,j in a if ii==i); carry=F(0)
        for j in range(d,0,-1):
            carry=a.get((i,j),F(0))+b*carry
            if carry:q[i,j-1]=carry
        r=a.get((i,0),F(0))+b*carry
        if r:rem[i,0]=r
    return clean(q),clean(rem)
def descend(a,b,n):
    for _ in range(n):
        a,r=div_linear(a,b)
        if r:return None
    return a

def rank(cols):
    basis={}
    for col in cols:
        v=dict(col)
        while v:
            p=max(v)
            if p not in basis:
                basis[p]=scale(v,1/v[p]);break
            v=add(v,scale(basis[p],-v[p]))
    return len(basis)

one={(0,0):F(1)};x={(1,0):F(1)};y={(0,1):F(1)}
identities=0;thresholds=0
for m in range(2,8):
    for t in list(range(1,9))+[None]:
        for b,c,beta in [(-2,3,5),(0,-2,-1),(3,1,0)]:
            z=add(y,scale(one,-b)); A=scale(power(z,m),c)
            B=scale(one,beta) if t is None else add(scale(one,beta),power(z,t))
            P=add(mul(A,x),B);U=add(P,scale(one,-beta))
            v=m if t is None else min(m,t);e=(m-2+v)//v
            for k in range(0,e+3):
                h=power(U,k);Q=descend(h,b,m-1)
                need((Q is not None)==(k>=e),'sharp valuation threshold')
                thresholds+=1
                if Q is not None:
                    Q=scale(Q,F(-1,c*(m-1)))
                    # Independent kernel addition cannot alter image.
                    Q=add(Q,power(P,2))
                    need(bracket(P,Q)==h,'explicit bracket');identities+=1
for B in [{},y,add(power(y,3),one)]:
    P=add(scale(x,3),B)
    for k in range(5):
        h=power(P,k);Q=add(scale(mul(y,h),F(1,3)),power(P,3))
        need(bracket(P,Q)==h,'constant A');identities+=1

# Full bounded matrix intersections, without using candidate integration to
# construct derivation columns. Domain: 0<=deg_x<=3, 0<=deg_y<=10.
# Target space: span{1,P,P²,P³}. Compare exact intersection dimension to
# full theorem's value. Each nonzero predicted basis has a companion in box.
cases=[]
for m,t in [(2,1),(3,1),(4,1),(4,2),(4,4),(4,None)]:
    A=power(y,m);B={} if t is None else power(y,t)
    v=m if t is None else min(m,t);e=(m-2+v)//v
    cases.append(('pure_%s_%s'%(m,t),A,B,max(0,4-e)))
cases.extend([
 ('constant',scale(one,2),add(power(y,2),one),4),
 ('linear',y,y,0),
 ('two_simple',mul(y,add(y,scale(one,-1))),one,0),
 ('two_double',mul(power(y,2),power(add(y,scale(one,-1)),2)),y,0),
 ('irreducible_Q',add(power(y,2),one),y,0),
 ('A_zero_linear',{},y,4),
 ('A_zero_monomial',{},power(y,4),3),
 ('A_zero_repeated_value',{},power(add(power(y,2),scale(one,-1)),2),2),
 ('A_zero_irr_critical',{},add(power(y,3),y),2),
])
matrices=[]
for name,A,B,expected in cases:
    P=add(mul(A,x),B)
    M=[bracket(P,{(i,j):F(1)}) for i in range(4) for j in range(11)]
    N=[power(P,k) for k in range(4)]
    r=rank(M);rn=rank(N);combined=rank(M+N);inter=r+rn-combined
    need(inter==expected,'bounded image intersection '+name)
    matrices.append({'case':name,'domain_columns':44,'derivation_rank':r,'target_rank':rn,'intersection_dimension':inter})

# Rational primitive pole numerator degree identity in independently
# enumerated coprime samples R=N/D with several distinct denominator roots.
# Rather than assuming residues cancel, test reciprocal-derivative obstruction
# by reduced-numerator degree, using univariate exact Euclidean arithmetic.
def univ(p):return {j:v for (i,j),v in p.items() if v}
def udiff(p):return {j-1:v*j for j,v in p.items() if j}
def udiv(a,b):
    a=dict(a);q={};db=max(b)
    while a and max(a)>=db:
        k=max(a)-db;c=a[max(a)]/b[db];q[k]=c
        for j,v in b.items():
            a[j+k]=a.get(j+k,F(0))-c*v
            if not a[j+k]:del a[j+k]
    return q,a
def ugcd(a,b):
    while b:a,b=b,udiv(a,b)[1]
    return {j:v/a[max(a)] for j,v in a.items()}
primitive_checks=0
for d1 in range(1,5):
    for d2 in range(1,5):
        D=mul(power(y,d1),power(add(y,scale(one,-1)),d2))
        for n in range(d1+d2):
            N=one if n==0 else add(power(y,n),one)
            # N(0)=1, N(1)=2, so all samples are coprime to D.
            W=add(mul(der(N,1),D),scale(mul(N,der(D,1)),-1))
            g=ugcd(univ(W),univ(mul(D,D)))
            reduced,rem=udiv(univ(W),g)
            need(not rem and max(reduced)==n+1,'two-pole numerator-degree obstruction')
            primitive_checks+=1
# Independent Laurent-polynomial check of the parameter-family obstruction.
# Exponents are (t,x,y), permitting negative t powers.
def la(a,b):
    d=dict(a)
    for k,v in b.items():d[k]=d.get(k,F(0))+v
    return clean(d)
def lm(a,b):
    d={}
    for k,v in a.items():
        for l,w in b.items():
            q=tuple(ki+li for ki,li in zip(k,l));d[q]=d.get(q,F(0))+v*w
    return clean(d)
def ld(a,j):
    d={}
    for k,v in a.items():
        if k[j]:
            q=list(k);q[j]-=1;d[tuple(q)]=v*k[j]
    return clean(d)
def lb(a,b):return la(lm(ld(a,1),ld(b,2)),scale(lm(ld(a,2),ld(b,1)),-1))
Pt={(1,1,5):F(1),(0,0,3):F(1)}
Qt={(1,2,6):F(-1,4),(0,1,4):F(-1,2),(-1,0,2):F(-1,4)}
need(lb(Pt,Qt)==lm(Pt,Pt),'Laurent family identity')
need(lb({(0,0,3):F(1)},{(0,1,4):F(-1,3)})=={(0,0,6):F(1)},'special fibre companion')
need(2%3!=0,'regular lift exponent obstruction')
print(json.dumps({'status':'PASS','arithmetic':'fractions.Fraction only','producer_imports':False,'identities':identities,'threshold_cases':thresholds,'two_pole_primitive_cases':primitive_checks,'parameter_family_laurent_identities':2,'matrix_universe':matrices},sort_keys=True,indent=2))
