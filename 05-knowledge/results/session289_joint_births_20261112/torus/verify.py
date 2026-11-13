#!/usr/bin/env python3
"""Exact supplied-coset certificates; not a complete exceptional-offset solver."""
from fractions import Fraction as Q
from itertools import product
import json

def chi(b,v):
    z=Q(1)
    for x,y in zip(b,v):z*=Q(x)**y
    return z

def value(t,v):return sum((Q(c)*chi(b,v) for c,b in t),Q(0))
def cert(t,a,gens):
    sums={}
    for c,b in t:
        sig=tuple(chi(b,g) for g in gens)
        sums[sig]=sums.get(sig,Q(0))+Q(c)*chi(b,a)
    return all(x==0 for x in sums.values()),sorted([{'signature':list(map(str,k)),'sum':str(v)} for k,v in sums.items()],key=lambda x:x['signature'])
def run():
    d=[(1,(2,1)),(-1,(1,2))]
    m=[(1,(2,3)),(-3,(2,1)),(-2,(1,3)),(6,(1,1))]
    c=[(1,(2,1,1)),(1,(1,2,1)),(-1,(1,1,2))]
    cases=[('physical_diagonal',d,(1,1),[(1,1)],True),('axis_rejected',d,(1,1),[(1,0)],False),('formal_axis1',m,(1,1),[(1,0)],True),('formal_axis2',m,(1,1),[(0,1)],True),('formal_joint_rejected',m,(1,1),[(1,0),(0,1)],False),('carry',c,(1,1,2),[(1,1,1)],True),('wrong_anchor',c,(1,1,1),[(1,1,1)],False)]
    out=[]
    for name,t,a,g,expected in cases:
        ok,s=cert(t,a,g)
        if ok!=expected:raise RuntimeError(name)
        if ok:
            for coefficients in product(range(-3,4),repeat=len(g)):
                p=tuple(a[i]+sum(z*v[i] for z,v in zip(coefficients,g)) for i in range(len(a)))
                if value(t,p)!=0:raise RuntimeError((name,p))
        out.append({'name':name,'passes':ok,'sums':s})
    counts=[]
    for n in range(2,31):
        x=sum(value(d,(a,n-a))==0 for a in range(1,n));y=sum(value(m,(a,n-a))==0 for a in range(1,n));z=sum(value(c,(a,b,n-a-b))==0 for a in range(1,n-1) for b in range(1,n-a))
        if (x,y,z)!=(int(n%2==0),1 if n==2 else 2,int(n>=4 and n%3==1)):raise RuntimeError(n)
        counts.append([n,x,y,z])
    return {'status':'PASS','cosets':out,'total_length_counts':counts}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
