#!/usr/bin/env python3
"""Standalone complete finite inverse certificate for s(n)=343."""
import json
from math import isqrt, prod

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def fac(n):
    f={}; p=2
    while p*p<=n:
        while n%p==0: f[p]=f.get(p,0)+1; n//=p
        p+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def sigma(n): return prod((p**(a+1)-1)//(p-1) for p,a in fac(n).items())
t=343; bound=(t-1)**2
s=[0]*(bound+1)
for d in range(1,bound//2+1):
    for n in range(2*d,bound+1,d): s[n]+=d
fiber=[n for n in range(1,bound+1) if s[n]==t]
need(all(sigma(n)-n==t for n in fiber),'independent multiplicative reencoding')
three=[n for n in fiber if len(fac(n))==3 and all(a==1 for a in fac(n).values())]
need(three==[1001],'unique squarefree three-prime source')
cases=[]
for p in [3,5,7]:
    rhs=343+p*(p+1)
    pairs=[(u,rhs//u) for u in range(1,isqrt(rhs)+1) if rhs%u==0]
    candidates=[]
    for u,v in pairs:
        q=u-p-1; r=v-p-1
        if p<q<r and fac(q)=={q:1} and fac(r)=={r:1}: candidates.append([p,q,r])
    cases.append({'p':p,'rhs':rhs,'factor_pairs':pairs,'eligible_prime_triples':candidates})
need([tr for case in cases for tr in case['eligible_prime_triples']]==[[7,11,13]],'factor-pair proof')
need(sigma(1001)-1001==343 and sigma(343)==400,'anchor')
# Identity checked by exact coefficient convolution rather than floating point.
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
# s(p(p+4)(p+6)) formal three-distinct-prime divisor expression is 35+23p+3p^2.
left=[35,23,3,-1]
right=[-x for x in mul([-7,1],[5,4,1])]
need(left==right,'constellation polynomial identity')
print(json.dumps({'target':t,'proved_source_bound':bound,'fiber':fiber,'source_factorizations':{str(n):fac(n) for n in fiber},'squarefree_three_prime_sources':three,'factor_pair_cases':cases,'constellation_difference_coefficients':left,'status':'PASS'},indent=2,sort_keys=True))
