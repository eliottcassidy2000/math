#!/usr/bin/env python3
"""Exact stdlib tests; failures survive python -O. No primality-distribution claim."""
from math import gcd
from collections import defaultdict
checks=0
def ck(t,msg):
    global checks
    checks+=1
    if not t: raise RuntimeError(msg)
def mul(x,y,m):
    u,v=x;r,s=y
    return ((u*r+v*s)%m,(u*s+v*r+v*s)%m)
def power(x,n,m):
    y=(1%m,0)
    while n:
        if n&1:y=mul(y,x,m)
        x=mul(x,x,m);n//=2
    return y
def tr(x,m):return (2*x[0]+x[1])%m
def norm(x,m):return (x[0]**2+x[0]*x[1]-x[1]**2)%m
def bar(x,m):return ((x[0]+x[1])%m,-x[1]%m)
def fib(n):
    a,b=0,1
    for _ in range(n):a,b=b,a+b
    return a
def factors(n):
    out=[];p=2
    while p*p<=n:
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:out.append(n)
    return out
def order(x,bound,m):
    ck(power(x,bound,m)==(1%m,0),'bad supplied bound')
    d=bound
    for q in factors(bound):
        while d%q==0 and power(x,d//q,m)==(1%m,0):d//=q
    ck(power(x,d,m)==(1%m,0),'order power')
    for q in factors(d):ck(power(x,d//q,m)!=(1%m,0),'order minimality')
    return d
for m in range(1,71):
    for d in range(1,13):
        a=power((1,1),d,m);T=tr(a,m);fibers=defaultdict(int)
        for u in range(m):
            for v in range(m):
                x=(u,v);y=mul(a,x,m);z=mul(a,y,m)
                pair=(tr(x,m),tr(y,m));fibers[pair]+=1
                ck(tr(z,m)==(T*pair[1]-pair[0])%m,'recurrence')
                ck((pair[0]**2+pair[1]**2-T*pair[0]*pair[1]+5*fib(2*d)**2*norm(x,m))%m==0,'trace conic norm')
        g=gcd(m,5*fib(2*d))
        ck(len(fibers)==m*m//g,'image size')
        ck(set(fibers.values())=={g},'uniform kernel size')
# Complete one-trace norm-one fibres at unramified primes.
rows=[]
for p in range(3,101,2):
    if any(p%q==0 for q in range(2,int(p**.5)+1)) or p==5:continue
    leg=1 if pow(5,(p-1)//2,p)==1 else -1
    U=[(u,v) for u in range(p) for v in range(p) if norm((u,v),p)==1]
    ck(len(U)==p-leg,'norm-one group size')
    fs=defaultdict(set)
    for x in U:fs[tr(x,p)].add(x)
    for x in U:ck(fs[tr(x,p)]=={x,bar(x,p)},'one-trace inversion orbit')
    alpha=(0,1);a=(1,1)
    if leg<0:ck(power(alpha,p+1,p)==(p-1,0),'inert norm minus one')
    else:ck(power(alpha,p-1,p)==(1,0),'split clock')
    bound=p-leg;oa=order(a,bound,p)
    rows.append((p,'split' if leg>0 else 'inert',bound,oa))
# Ramification: each whole-history class on norm-one states has five points.
p=5;a=(1,1);fs=defaultdict(list)
for u in range(p):
    for v in range(p):
        x=(u,v)
        if norm(x,p)!=1:continue
        hist=[]
        for k in range(12):hist.append(tr(x,p));x=mul(x,a,p)
        fs[tuple(hist)].append((u,v))
ck(len(fs)==2 and sorted(map(len,fs.values()))==[5,5],'ramified five-fold history')
# Decimation hostile: same every-other trace, different next trace.
x=(1,1);y=bar(x,3);a=(1,1)
ck(x!=y and tr(x,3)==tr(y,3),'decimation distinct states')
ck(tr(mul(a,x,3),3)!=tr(mul(a,y,3),3),'unit step sees sign')
for k in range(12):ck(tr(mul(power(a,2*k,3),x,3),3)==tr(mul(power(a,2*k,3),y,3),3),'stride-two blindness')
# CRT sign branching: no single global inversion sidecar suffices.
m=33;fs=defaultdict(list)
for u in range(m):
    for v in range(m):
        x=(u,v)
        if norm(x,m)==1:fs[tr(x,m)].append(x)
key=next(k for k,v in sorted(fs.items()) if len(v)==4)
ck(len(fs[key])==4,'four CRT branches')
# Two coprime clocks restore all unramified information, including composite m.
for m in range(1,51):
    for ds in [(2,3),(2,4),(3,6),(4,6),(5,7)]:
        images=set()
        for u in range(m):
            for v in range(m):
                x=(u,v)
                images.add((tr(x,m),)+tuple(tr(mul(power((1,1),d,m),x,m),m) for d in ds))
        g=gcd(m,*[5*fib(2*d) for d in ds])
        ck(len(images)==m*m//g,'multi-clock exact quotient')
        ck(g==gcd(m,5*fib(2*gcd(*ds))),'Fibonacci strong divisibility')
print('PASS',checks,'checks; m=1..70, d=1..12; odd primes<=100')
print('p,type,norm-one-size,order(phi^2):',rows)
print('mod33 four-state single-trace fibre:',key,sorted(fs[key]))
print('mod5 full-history fibre sizes:',[5,5])
