#!/usr/bin/env python3
"""Independent exact checks: even-run formula, succession sum, tournament DP.
Run: python3 session191-moments/check_moments.py
Universe: all permutations n<=9; all 2^(n choose 2) labeled tournaments n<=6.
"""
from math import factorial, comb
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def run_terms(n):
    return [F(factorial(n-2*j), factorial(n))*sum(2**r*comb(j-1,r-1)*comb(n-2*j,r) for r in range(1,min(j,n-2*j)+1)) for j in range(1,(n-1)//2+1)]

def moment_kernel(n):
    return int(factorial(n)*(1+sum(run_terms(n),F(0))))

def succession_kernel(n):
    ans=0
    for p in permutations(range(n)):
        if any(b==a-1 for a,b in zip(p,p[1:])): continue
        ans += 2**sum(b==a+1 for a,b in zip(p,p[1:]))
    return ans

def odd_composition_kernel(n):
    # Ordered compositions; factorial weights cannot be folded into a scalar
    # recurrence, so retain the number of parts as a state.
    d=[[0]*(n+1) for _ in range(n+1)]
    d[0][0]=1
    for s in range(n+1):
        for k in range(n):
            for a in range(1,n-s+1,2):
                d[s+a][k+1] += d[s][k]*(1 if a==1 else 2)
    return sum(factorial(k)*d[n][k] for k in range(n+1))

def h_dp(n,code):
    incoming=[0]*n
    for bit,(a,b) in enumerate(( (a,b) for a in range(n) for b in range(a+1,n))):
        if code>>bit&1: incoming[b]|=1<<a
        else: incoming[a]|=1<<b
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1:
                prev=mask^(1<<v)
                todo=prev&incoming[v]
                while todo:
                    bit=todo&-todo; todo-=bit
                    dp[mask][v]+=dp[prev][bit.bit_length()-1]
    return sum(dp[-1])

rows=[]
for n in range(1,10):
    s=moment_kernel(n)
    assert s==succession_kernel(n)==odd_composition_kernel(n)
    mu=F(factorial(n),2**(n-1)); eh2=F(factorial(n)*s,4**(n-1))
    row={'n':n,'S_n':s,'mean':str(mu),'second_moment':str(eh2),'normalized_variance':str(sum(run_terms(n),F(0)))}
    if n<=6:
        vals=[h_dp(n,i) for i in range(1<<(n*(n-1)//2))]
        assert F(sum(vals),len(vals))==mu
        assert F(sum(v*v for v in vals),len(vals))==eh2
        row['all_labeled_tournaments']=len(vals)
    rows.append(row)
    print(row,flush=True)

# Tail bound, positivity, and exact one-sided remainder; includes boundary n=1,2.
checks=0
for n in range(1,201):
    terms=run_terms(n)
    if n>=3: assert 0<sum(terms)<F(2,n)
    for J in range(1,len(terms)+2):
        A=F(2**J*factorial(n-J),factorial(J)*factorial(n))
        assert 0<=sum(terms[J-1:])<=A
        checks+=1

# Exact rational expansion (the tail with J=7 is O(n^-7)).
import sympy as sp
x=sp.symbols('x')
n=sp.symbols('n',positive=True)
poly=0
for j in range(1,7):
    tj=sum(2**r*comb(j-1,r-1)*sp.prod(n-2*j-i for i in range(r))/factorial(r) for r in range(1,j+1))/sp.prod(n-i for i in range(2*j))
    poly+=sp.series(tj.subs(n,1/x),x,0,7).removeO()
expansion=str(sp.expand(poly))
assert sp.expand(poly).coeff(x,1)==2
assert sp.expand(poly).coeff(x,2)==0
assert sp.expand(poly).coeff(x,3)==-sp.Rational(14,3)
assert sp.expand(poly).coeff(x,4)==-sp.Rational(58,3)
print('tail_checks',checks)
print('variance expansion x=1/n:',expansion)
report={'status':'PASS','rows':rows,'tail_checks':checks,'tail_check_n_max':200,'variance_expansion_x_equals_1_over_n':expansion}
(ROOT/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
