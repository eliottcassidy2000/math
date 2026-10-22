#!/usr/bin/env python3
import math,json

def need(b,m):
    if not b: raise RuntimeError(m)

def inv_count(N,q,a):
    t=[pow(a,n*(n-1)//2,q) for n in range(N+1)]; s=[0]*(N+1)
    for n in range(1,N+1): s[n]=(t[n]-sum(math.comb(n,k)*s[k]*t[n-k] for k in range(1,n)))%q
    return t,s
base_checks=connected_checks=0
for p,r in [(3,1),(3,2),(3,3),(5,1),(5,2),(7,2)]:
    q=p**r; P=q//p*(p-1); N=P+2*r+15
    for a in [-5,-2,-1,1,2,3,4,7]:
        if a%p==0: continue
        eps=1 if pow(a,(p-1)//2,p)==1 else -1
        t,s=inv_count(N,q,a)
        for n in range(2*r-1,N-P+1):
            need((s[n+P]-eps*s[n])%q==0,('unit base',a,p,r,n)); base_checks+=1
    t,s=inv_count(N,q,2); b=[0]*(N+1)
    for n in range(1,N+1): b[n]=(t[n]-sum(math.comb(n-1,k-1)*b[k]*t[n-k] for k in range(1,n)))%q
    eps=1 if pow(2,(p-1)//2,p)==1 else -1
    for n in range(2*r,N-P+1):
        need((b[n+P]-eps*b[n])%q==0,('connected',p,r,n)); connected_checks+=1

def factors(n):
    out=[]; d=2
    while d*d<=n:
        if n%d==0:
            r=0
            while n%d==0: n//=d; r+=1
            out.append((d,r))
        d+=1
    if n>1: out.append((n,1))
    return out
crt=[]
for n in range(4,201):
    fs=factors(n)
    if fs==[(n,1)]: continue
    t,s=inv_count(n,n,2)
    def D(k): return (s[k]-(pow(2,(k-1)//2,n) if k%2 else 0))%n
    conditions=[]
    for p,r in fs:
        if p==2: continue
        q=p**r; P=q//p*(p-1); low=2*r-1; z=low+(n-low)%P
        need(low<=z<=q<=n,('representative bound',n,p,r,z))
        conditions.append(D(z)%q==0)
    need((D(n)==0)==all(conditions),('CRT',n,conditions,D(n)))
    crt.append(n)
print(json.dumps(dict(status='PASS',base_checks=base_checks,connected_checks=connected_checks,crt_checks=len(crt),crt_maximum=200),indent=2))
