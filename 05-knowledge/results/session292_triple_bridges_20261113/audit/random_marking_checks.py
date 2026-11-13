#!/usr/bin/env python3
"""Exact probability and formal-log expectation checks on every tournament n<=4."""
from itertools import product
from fractions import Fraction as F
from pathlib import Path
import importlib.util,json
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('colors',ROOT/'triples/colored_fibers.py')
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def vp(m,p):
    e=0
    while m%p==0:e+=1;m//=p
    return e
def expected(t,q,p):
    kind,mask,kids=t
    out=sum((expected(k,q,p) for k in kids),F(0))
    if kind=='prime' and p==2:out+=F(1,q**(mask.bit_count()-1))
    if kind=='linear':
        for i in range(len(kids)):
            size=kids[i][1].bit_count()
            for j in range(i+1,len(kids)):
                size+=kids[j][1].bit_count();ell=j-i+1
                out+=F(vp(ell,p)-vp(ell-1,p),q**(size-1))
    return out
cases=0;expected_checks=0;transitive_checks=0
for n in range(1,5):
    for code in range(1<<(n*(n-1)//2)):
        a=f.rows(n,code);t=f.tree(a,(1<<n)-1)
        transitive=sorted(x.bit_count() for x in a)==list(range(n))
        for q in (1,2,3):
            values=[f.predicted(t,c) for c in product(range(q),repeat=n)]
            pu=F(values.count(1),q**n)
            need(pu>=1-F(n-1,q),'uniform uniqueness bound')
            if transitive:
                need(pu==(1-F(1,q))**(n-1),'transitive exact probability');transitive_checks+=1
            for p in (2,3):
                empirical=F(sum(vp(v,p) for v in values),q**n)
                need(empirical==expected(t,q,p),'formal logarithm expectation')
                expected_checks+=1
            cases+=1
print(json.dumps({'status':'PASS','tournament_palette_probability_cases':cases,
 'transitive_exact_cases':transitive_checks,'prime_valuation_expectation_cases':expected_checks,
 'scope':'Every labelled tournament n<=4; q1,2,3. Exact Fractions; no logarithmic floating point.'},indent=2))
