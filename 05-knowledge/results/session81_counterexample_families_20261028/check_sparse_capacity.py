#!/usr/bin/env python3
"""Exact Burnside counts; no external dependencies. n<=40."""
from fractions import Fraction
from math import factorial, gcd, comb
import json

def partitions(n,lo=1):
    if not n:
        yield ()
    for a in range(lo,n+1):
        for rest in partitions(n-a,a): yield (a,)+rest

def edge_orbits(p):
    out=[]
    for i,a in enumerate(p):
        out += [a]*((a-1)//2)
        if a%2==0: out.append(a//2)
        for b in p[:i]: out += [a*b//gcd(a,b)]*gcd(a,b)
    return out

def z(p):
    ans=1
    for a in set(p): ans*=a**p.count(a)*factorial(p.count(a))
    return ans

def coefficient(lengths,m):
    c=[1]+[0]*m
    for a in lengths:
        for j in range(m,a-1,-1): c[j]+=c[j-a]
    return c[m]

def counts(n):
    g=t=Fraction(0)
    for p in partitions(n):
        ls=edge_orbits(p); den=z(p)
        if sum(ls)!=comb(n,2): raise ValueError('pair orbit partition')
        g+=Fraction(coefficient(ls,n-1),den)
        if all(a%2 for a in p): t+=Fraction(2**len(ls),den)
    if g.denominator!=1 or t.denominator!=1: raise ValueError('integrality')
    return int(g),int(t)

rows=[]
for n in range(1,31):
    g,t=counts(n)
    rows.append(dict(n=n,sparse_graph_types=g,tournament_types=t,
        independent_pair_possible=t<=g*g,independent_triple_possible=t<=g**3,
        labelled_pair_upper_exceeded=t>comb(comb(n,2),n-1)**2,
        labelled_triple_upper_exceeded=t>comb(comb(n,2),n-1)**3))
# Hostile to the false tree implication and canonical small universes.
if [r['sparse_graph_types'] for r in rows[:5]]!=[1,1,1,3,6]: raise ValueError('small sparse count')
if [r['tournament_types'] for r in rows[:6]]!=[1,1,2,4,12,56]: raise ValueError('small tournament count')
print(json.dumps({'method':'exact partition Burnside; ordered independent graph tuples', 'rows':rows},indent=2))
