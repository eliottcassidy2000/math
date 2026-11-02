#!/usr/bin/env python3
"""Exact controls for arbitrary integral residue-coset tensor compensation."""
from fractions import Fraction as F
from itertools import product
from random import Random
import json
rng=Random(139202)
counts={};min_gap={}

def require(ok,where):
    if not ok:raise RuntimeError(where)

def record(key,gap):
    require(gap>=0,(key,str(gap)))
    counts[key]=counts.get(key,0)+1
    min_gap[key]=min(min_gap.get(key,gap),gap)

def prod(vals):
    a=F(1)
    for v in vals:a*=v
    return a

def distance(M,rho):
    a=rho%M
    return min(a,M-a)

def transform(b,x):
    f=list(map(F,b))
    for i,xi in enumerate(x):
        for mask in range(len(f)):
            if mask>>i&1:f[mask]+=xi*f[mask^(1<<i)]
    return f

def direct(b,x):
    out=[]
    for mask in range(len(b)):
        value=F(0);sub=mask
        while True:
            value+=b[mask^sub]*prod(x[i] for i in range(len(x)) if sub>>i&1)
            if not sub:break
            sub=(sub-1)&mask
        out.append(value)
    return out

def test(M,rho,b,x,key):
    require(M>=2 and (b[0]-rho)%M==0 and all(v%M==0 for v in b[1:]),'invalid intended coset')
    require(len(b)==1<<len(x) and all(abs(v)<=1 for v in x),'invalid intended tensor')
    d=distance(M,rho)
    require(d==min(abs(rho+M*k) for k in range(-rho//M-2,-rho//M+3)),'minimum residue distance')
    f=transform(b,x)
    gap=sum(v*v for v in f)-d*d*prod(1+v*v for v in x)
    record(key,gap)
    return gap

xs=[F(-1),F(-2,3),F(0),F(1,3),F(1)]
for M in range(2,10):
    for rho in range(M):
        for s in range(3):
            for k in [-1,0,1]:
                for tail in product([-M,0,M],repeat=(1<<s)-1):
                    b=[rho+k*M,*tail]
                    for x in product(xs,repeat=s):test(M,rho,b,x,'exhaustive_dimension_'+str(s))
# Noncanonical representatives, composite and prime moduli, and large moduli.
moduli=[2,3,4,5,6,7,8,9,10,12,15,16,21,32,49,64,97,121,1000,2**61+1]
for M in moduli:
    for s in range(6):
        for trial in range(6):
            rho=rng.randrange(-4*M,4*M+1)
            b=[rho+M*rng.randrange(-3,4)]+[M*rng.randrange(-3,4) for _ in range((1<<s)-1)]
            x=[rng.choice([F(-1),F(-7,8),F(-1,5),F(0),F(2,7),F(1)]) for _ in range(s)]
            test(M,rho,b,x,'mixed_dimension_'+str(s))
            require(transform(b,x)==direct(b,x),'independent tensor transform formula')
            counts['independent_transform_identity']=counts.get('independent_transform_identity',0)+1
# Equality holds for a minimum representative, at every tested point.
for M in range(2,33):
    for rho in range(-M,2*M):
        a=rho%M;nearest=a if 2*a<=M else a-M
        for s in range(6):
            b=[nearest]+[0]*((1<<s)-1)
            x=[rng.choice(xs) for _ in range(s)]
            gap=test(M,rho,b,x,'sharp_minimum_representative')
            require(gap==0,'minimum-representative equality')
# Negative controls verify that necessary assumptions are not silently dropped.
mutations={}
def reject(name,b,x,false_floor):
    gap=sum(v*v for v in transform(b,x))-false_floor
    require(gap<0,'mutation survived: '+name)
    mutations[name]=str(gap)
reject('outside_cube_M2_rho1',[1,-2],[F(2)],F(5))
reject('nonmultiple_tail_M3_rho1',[1,-1],[F(1)],F(2))
reject('wrong_constant_coset_M5_rho2',[1,0],[F(0)],F(4))
# rho itself cannot replace the minimum representative d.
reject('using_unreduced_residue_M5_rho4',[-1,0],[F(0)],F(16))
print(json.dumps({'status':'PASS','counts':counts,'minimum_gaps':{k:str(v) for k,v in min_gap.items()},'negative_controls':mutations,'all_checks_active_under_optimization':True,'scope':'Arbitrary-modulus tensor lemma only; no new Pfaffian, stability cutoff, or RH claim.'},indent=2,sort_keys=True))
