#!/usr/bin/env python3
"""Exact finite controls; the general proofs are in PROOF.md. No assert statements."""
from fractions import Fraction as Q
from math import factorial
from itertools import combinations
import json

def need(b,msg):
    if not b: raise RuntimeError(msg)
def a(m): return Q((m+1)**(m-1),factorial(m)) if m else Q(1)
N=24
A=[a(m) for m in range(N+1)]
E=[Q(1)]
for n in range(1,N+1):
    E.append(sum(Q(k)*A[k-1]*E[n-k] for k in range(1,n+1))/n)
need(A==E,'F=exp(zF) coefficient control failed')
expected=list(map(Q,['1','3/2','8/3','125/24','54/5','16807/720']))
need(A[1:7]==expected,'user coefficient control failed')
for m in range(N): need(A[m+1]/A[m]==Q(m+2,m+1)**m,'ratio control failed')
forests=[]
for n in range(1,7):
    edges=list(combinations(range(n),2)); unrooted=rooted=0
    for mask in range(1<<len(edges)):
        p=list(range(n)); size=[1]*n
        def find(v):
            while p[v]!=v: v=p[v]
            return v
        for i,(u,v) in enumerate(edges):
            if mask>>i&1:
                u,v=find(u),find(v)
                if u==v: break
                p[u]=v; size[v]+=size[u]
        else:
            unrooted+=1; weight=1
            for v in range(n):
                if find(v)==v: weight*=size[v]
            rooted+=weight
    need(rooted==(n+1)**(n-1),'forest formula finite control failed')
    forests.append({'n':n,'unrooted':unrooted,'rooted':rooted})
need(forests[2]['rooted']==16 and forests[2]['unrooted']==7,'triangle hostile control')
# log(2)=2 sum_{j>=0} (1/3)^(2j+1)/(2j+1), positive rational tail.
K=80
L=2*sum((Q(1,3)**(2*j+1))/Q(2*j+1) for j in range(K))
U=L+2*Q(1,3)**(2*K+1)/(Q(2*K+1)*(1-Q(1,9)))
# e<11/4: sum through k=4 is 65/24; remaining <=(1/120)/(1-1/6)=1/100.
need(Q(65,24)+Q(1,100)<Q(11,4),'rational e upper bound')
def outward(x,upper=False,places=30):
    scale=10**places; num=x.numerator*scale; den=x.denominator
    q=(num+den-1)//den if upper else num//den
    return str(q//scale)+'.'+str(q%scale).zfill(places)
intervals=[]
for n,divisor,M in [(4,2,250),(16,4,60)]:
    zl,zu=L/divisor,U/divisor; q=Q(11,4)*zu
    need(q<1,'tail ratio not contractive')
    lo=sum(a(m)*zl**m for m in range(M+1))
    hi=sum(a(m)*zu**m for m in range(M+1))+a(M+1)*zu**(M+1)/(1-q)
    if n==4: need(lo<Q(2)<hi,'exact r4=2 enclosure')
    intervals.append({'n':n,'M':M,'lower_decimal_illustration':float(lo),'upper_decimal_illustration':float(hi),'width_decimal_illustration':float(hi-lo),'outward_lower_30dp':outward(lo),'outward_upper_30dp':outward(hi,True),'rational_bound_checked':True})
print(json.dumps({'status':'PASS','coefficient_checks_through':N,'coefficients':[str(x) for x in A[:10]],'forests':forests,'rigorous_rational_interval_controls':intervals,'note':'Displayed decimals are illustrations, not outward-rounded endpoints. Exact fractions are evaluated by this verifier.'},indent=2))
