#!/usr/bin/env python3
"""Independent path DP and slot convolution for the shared-center fan."""
from math import comb

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def graph(n,t):
    a=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            u,v=(j,i) if i==0 and 2<=j<=t+1 else (i,j)
            a[u]|=1<<v
    return a

def wrap(a):
    n=len(a);b=a[:]+[0,0]
    for i in range(n):b[i]|=1<<n;b[n+1]|=1<<i
    b[n]|=1<<(n+1)
    return b

def paths(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,x in enumerate(d[mask]):
            if not x:continue
            avail=a[i]&~mask
            while avail:
                bit=avail&-avail;j=bit.bit_length()-1
                d[mask|bit][j]+=x;avail-=bit
    return [1]+[sum(x) for x in d[1:]]

def profile(a,r):
    h=paths(a);n=len(a);v=[0]*(1<<n);v[0]=1
    for step in range(r):
        z=[0]*(1<<n)
        for mask in range(1<<n):
            sub=mask
            while True:
                z[mask]+=h[sub]*v[mask^sub]
                if not sub:break
                sub=(sub-1)&mask
        v=z
    return v[-1]

def closed(n,t,r):return r**n+2*sum(comb(t,q)*r**(n-q-1) for q in range(1,t+1,2))
for n in range(4,9):
    values=[]
    for t in range(n-1):
        a=graph(n,t)
        for r in range(4):need(profile(a,r)==closed(n,t,r),f'profile n={n} t={t} r={r}')
        h=paths(wrap(a))[-1];need(h==closed(n,t,1)+closed(n,t,2),'wrapped count')
        values.append(h)
    jumps=[v-u for u,v in zip(values,values[1:])]
    for t,j in enumerate(jumps):
        formula=2+2**(n-1) if t==0 else 2**t+2**(n-2-t)*(3**t+1)
        need(j==formula,'jump formula')
    need(jumps[0]==jumps[1] and all(jumps[t]<jumps[t+1] for t in range(1,len(jumps)-1)),'monotonic jumps')
print('PASS: fan profiles, literal wrappers, all jump formulas n=4,...,8')
