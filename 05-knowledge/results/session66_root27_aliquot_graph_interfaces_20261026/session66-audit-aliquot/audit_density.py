#!/usr/bin/env python3
"""Independent streaming Chen--Zhao g_M interval audit, exact integers only.

Different summand representation from producer:
positive term = phi(M/d)*(sigma(2d)-4d)/(2*M*(sigma(2d)-2d)).
Each odd prime contributes its complement totient, not a phi/q weight.
The 2-adic factor is handled separately after odd divisor generation.
"""
from fractions import Fraction
from itertools import product
from math import gcd, prod
from pathlib import Path
import importlib.util
import json

ROOT = Path(__file__).resolve().parent
D=10**24
BASE=[(2,6),(3,5),(5,4),(7,3),(11,2),(13,1),(17,1),(19,1),(23,1),(29,1),(31,1),(37,1),(41,1)]
checks=0
def require(ok,msg):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(msg)

def evaluate(fs):
    M=prod(p**e for p,e in fs)
    e2=dict(fs).get(2,0)
    odd=[(p,e) for p,e in fs if p!=2]
    tables=[]
    for p,e in odd:
        tables.append([(p**k,sum(p**j for j in range(k+1)),1 if k==e else (p-1)*p**(e-k-1)) for k in range(e+1)])
    two=[(2**k,2**(k+2)-1,1 if k==e2 else 2**(e2-k-1)) for k in range(e2+1)]
    lower=positive=count=0
    def walk(i,a,sigma_a,phi_complement):
        nonlocal lower,positive,count
        if i<len(tables):
            for q,s,t in tables[i]:
                walk(i+1,a*q,sigma_a*s,phi_complement*t)
            return
        for power,s2,phi2 in two:
            d=a*power
            sigma_2d=sigma_a*s2
            count+=1
            if sigma_2d<=4*d:
                continue
            numerator=phi_complement*phi2*(sigma_2d-4*d)
            denominator=2*M*(sigma_2d-2*d)
            lower+=(D*numerator)//denominator
            positive+=1
    walk(0,1,1,1)
    return {'M':M,'divisor_count':count,'positive_terms':positive,'common_denominator':D,'lower_numerator':lower,'strict_upper_numerator':lower+max(1,positive)}

spec=importlib.util.spec_from_file_location('density_producer', ROOT/'source-snapshot/density_lower_bound.py')
producer=importlib.util.module_from_spec(spec)
spec.loader.exec_module(producer)

# Small hostile cases use literal divisors, totient by gcd counting, and
# Fraction evaluation of the original displayed theorem, not our simplification.
small_cases=[]
for fs in [[],[(2,1)],[(3,1)],[(3,2)],[(2,2),(3,1)],[(2,3),(3,2)],[(2,1),(3,1),(5,1)],[(3,2),(5,1)]]:
    M=prod(p**e for p,e in fs)
    exact=Fraction(0)
    direct_floor=positive=0
    for d in range(1,M+1):
        if M%d:
            continue
        q=M//d
        phi=sum(gcd(a,q)==1 for a in range(1,q+1))
        s=sum(j for j in range(1,2*d) if (2*d)%j==0)
        term=Fraction(phi,q)*max(Fraction(0),Fraction(1,2*d)-Fraction(1,s))
        exact+=term
        if term:
            positive+=1
            direct_floor+=D*term.numerator//term.denominator
    got=evaluate(fs)
    require(got['lower_numerator']==direct_floor,('direct small formula',M))
    require(got['positive_terms']==positive,('positive terms',M))
    require(Fraction(got['lower_numerator'],D)<=exact,('lower rounding',M))
    require(exact<Fraction(got['strict_upper_numerator'],D),('strict upper including exact zero',M))
    p=producer.bound(fs)
    require(all(got[k]==p[k] for k in got),('producer small even/odd/zero-positive',M))
    small_cases.append({'M':M,'exact_g_M':str(exact),'positive_terms':positive})

# Independently compare the original theorem with the gcd-expectation formula,
# and check every divisibility refinement in a named finite universe.
sums=[0]*481
for d in range(1,241):
    for n in range(2*d,481,d):
        sums[n]+=d
H={d:max(Fraction(0),Fraction(1,2)-Fraction(d,sums[2*d])) for d in range(1,241)}
g={M:sum((H[gcd(R,M)] for R in range(M)),Fraction(0))/M for M in range(1,241)}
refinements=0
for M in range(1,241):
    original=sum((Fraction(sum(gcd(a,M//d)==1 for a in range(1,M//d+1)),M//d)*H[d]/d
                  for d in range(1,M+1) if M%d==0),Fraction(0))
    require(g[M]==original,('gcd expectation',M))
    for N in range(M,241,M):
        require(g[M]<=g[N],('divisibility monotonicity',M,N))
        refinements+=1

stored=json.loads((ROOT/'source-snapshot/density_lower_bound.json').read_text())
rows={}
for name,fs in [('base',BASE),('enlarged',BASE+[(43,1)])]:
    got=evaluate(fs)
    require(all(got[k]==stored[name][k] for k in got),('complete independent stored interval',name))
    rows[name]=got
require(rows['enlarged']['lower_numerator']>rows['base']['strict_upper_numerator'],'strict certified improvement')
require(rows['enlarged']['lower_numerator']*10**15>60452936134095*D,'strict displayed lower bound')
print(json.dumps({'status':'PASS','explicit_checks':checks,'method':'streaming complement-totient / fixed denominator 2M; integer floors at 10^24','small_exact_checks':small_cases,'gcd_expectation_M':[1,240],'divisibility_refinements':refinements,'intervals':rows,'derived_bound':'liminf U(X)/X > 0.060452936134095','scope':'CITED Chen-Zhao theorem plus independently certified instantiation; no density existence or numerical record claim'},indent=2,sort_keys=True))
