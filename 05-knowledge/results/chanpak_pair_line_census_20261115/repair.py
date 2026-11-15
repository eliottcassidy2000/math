#!/usr/bin/env python3
"""Physically grow the marked-reader bank on a prescribed finite gap list."""
import sys
from itertools import combinations
from collections import defaultdict

targets=[32825,34111,35099,35395,35557,36833]

def rows(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1: a[i]|=1<<j
        else: a[j]|=1<<i
    return a

def dp(a):
    n=len(a); f=(1<<n)-1
    h=[1]+[0]*f
    ends=[[0]*n for _ in range(f+1)]
    for i in range(n): ends[1<<i][i]=1
    for s in range(1,f+1):
        h[s]=sum(ends[s])
        for i,x in enumerate(ends[s]):
            more=a[i]&~s
            while x and more:
                z=more&-more; more-=z
                ends[s|z][z.bit_length()-1]+=x
    return h,ends

def coeffs(n,code):
    a=rows(n,code); f=(1<<n)-1
    h,end=dp(a)
    _,start=dp([f^a[i]^(1<<i) for i in range(n)])
    lin=[sum(h[s]*start[f^s][i] for s in range(f)) for i in range(n)]
    pair=[[sum(end[s][i]*start[f^s][j] for s in range(1,f)) for j in range(n)] for i in range(n)]
    for o in range(1,f+1):
        bits=[i for i in range(n) if o>>i&1]
        M=sum(lin[i] for i in bits)
        W=sum(pair[i][j] for i in bits for j in bits)
        yield (h[f]-W,M),o

def close(host,mask,source):
    n=len(host); m=len(source); x=n+m
    bulk=((1<<m)-1)<<n
    out=[a|bulk|((1<<x) if not mask>>i&1 else 0) for i,a in enumerate(host)]
    return out+[(a<<n)|(1<<x) for a in source]+[mask]

def main(path,hosts):
    bank={}
    for n,code in sorted(tuple(map(int,line.split())) for line in open(hosts)):
        for k,mask in coeffs(n,code):
            bank.setdefault(k,(n,code,mask))
    sources=defaultdict(dict)
    lowest_h={}
    best={}
    for line in open(path):
        n,h,p,code,_=line.split()
        n,h,p,code=map(int,(n,h,p,code))
        q=p+2*h
        if q<=max(targets):
            sources[h].setdefault(q,(n,code))
        if h in targets:lowest_h.setdefault(h,(n,code))
    profiles = sorted((q,h,n,c) for h,qs in sources.items() for q,(n,c) in qs.items())
    for (r,M),(n,host,mask) in bank.items():
        if not (M>0 and r+2*M>0 and r>=1-M):
            raise ValueError((r,M,n,host,mask))
        missing=[N for N in targets if N not in best]
        if not missing:break
        for q,h,a,b in profiles:
            if M*q > 2*max(missing):break
            N=r*h+M*q
            if N in missing and N not in best:
                real=close(rows(n,host),mask,rows(a,b))
                assert dp(real)[0][-1]==N
                best[N]=(r,M,n,host,mask,a,b,h,q)
    print("Distinct marked readers:",len(bank))
    print("Repairs: target: (r,M,host_order,host_code,mask,source_order,source_code,H,Q)")
    for N in targets:
        print(N,best.get(N),"direct source",lowest_h.get(N))
if __name__=='__main__':
    main(*sys.argv[1:])
