#!/usr/bin/env python3
"""Finite exact controls only; no numerical fit supplies the all-order theorem."""
import json
from fractions import Fraction
from pathlib import Path

N=500
charges=[]
for ell in range(5,2001,2):
    c=3*(ell*ell-1)//8-ell
    if c>N:break
    charges.append((ell,c))
# All multisets of actual large-orbit sizes under the relaxed charge bound,
# compressed by maximum number of orbits at each exact charge.
dp=[-10**9]*(N+1);dp[0]=0
for charge in range(N+1):
    if dp[charge]<0:continue
    for ell,c in charges:
        if charge+c<=N:dp[charge+c]=max(dp[charge+c],dp[charge]+1)
rows=[];checks=0
for d in range(N+1):
    feasible=[(dp[D]+2*((d-D)//3)+1,D,dp[D]) for D in range(d+1) if dp[D]>=0 and (d-D)%3==0]
    b=d%3;a=(d-4*b)//3
    expected=None if a<0 else 2*a+b+1
    got=max((v for v,D,R in feasible),default=None)
    if got!=expected:raise RuntimeError(('penalty maximum',d,got,expected))
    checks+=1
    if expected is not None:
        k=(expected+1)//2
        if k!=a+(b+2)//2:raise RuntimeError(('nested palindrome pole',d))
        checks+=1
        rows.append({'slack':d,'ordinary_pole':expected,'self_converse_pole':k,'ordinary_recurrence_order':3*expected,'self_converse_recurrence_order':6*k,'merged_recurrence_order':3*expected+6*k,'attaining_a':a,'attaining_b':b,'safe_ordinary_onset':3*expected if d==0 else 4*d*d+d+3*expected+1,'safe_self_converse_onset':6*k if d==0 else 4*d*d+d+6*k+1,'safe_merged_onset':3*expected+6*k if d==0 else 4*d*d+d+3*expected+6*k+1})
# Polynomial exact arithmetic, coefficients in increasing degree.
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
f2=[1,1,0,1]
rmod2=[x%2 for x in mul(f2,f2)]
if rmod2!=[1,0,1,0,0,0,1]:raise RuntimeError('mod2 factor square')
if any(sum(c*t**i for i,c in enumerate(f2))%2==0 for t in (0,1)):raise RuntimeError('cubic reducible')
checks+=2
# A direct finite coefficient check of the zero, three and four channels.
# Fractions from rational recurrence avoid symbolic-library dependence.
def series(num,den,n):
    v=[Fraction(0) for _ in range(n+1)]
    for i in range(n+1):v[i]=(Fraction(num[i] if i<len(num) else 0)-sum(den[j]*v[i-j] for j in range(1,min(i,len(den)-1)+1)))/den[0]
    return v
p=[1,-1,0,-1];q=[1,0,-1,0,0,0,-1]
A=series([1],p,30)
P=series([1,1,0,1],q,30)
# Palindrome words enumerate directly through vertex degree12.
words=[]
def go(w,n):
    words.append(w)
    for c in (1,3):
        if n+c<=12:go(w+(c,),n+c)
go((),0)
for n in range(13):
    if A[n]!=sum(sum(w)==n for w in words):raise RuntimeError('word coefficients')
    if P[n]!=sum(sum(w)==n and w==w[::-1] for w in words):raise RuntimeError('palindrome coefficients')
    checks+=2
print('Charge-universe upper bound',N)
print('Large-orbit charges',charges)
print('Ordinary, self-converse, merged pole/recurrence rows through24:')
for row in rows:
    if row['slack']<=24:print(json.dumps(row,sort_keys=True))
print('Exact checks',checks)
print('PASS; bounded controls only, all-order justification is SHARP_POLE_ORDER.md')
Path(__file__).with_name('sharp_pole_rows.json').write_text(json.dumps({'bound':N,'charges':charges,'rows':rows,'checks':checks},indent=2)+'\n')
