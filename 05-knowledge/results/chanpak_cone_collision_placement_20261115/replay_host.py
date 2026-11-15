#!/usr/bin/env python3
"""Exact two-apex slot DP and a constructive macroscopic host decoder."""
from collections import defaultdict
from fractions import Fraction
from itertools import product
from random import Random
from sys import argv

# Coverage by marked runs on an increasing spine. States per run:
# 0 = not started, 1/2 = open odd/even, 3/4 = closed odd/even.
def runs(a, b, k, endpoints, accept, empty):
    states={(0,)*len(endpoints):1}
    for i in range(len(a)):
        bit=(a[i],b[i])
        new=defaultdict(int)
        for S,c in states.items():
            new[S]+=k*c
            for j,(u,v) in enumerate(endpoints):
                p=S[j]
                if p>=3 or (not p and not bit[u]):
                    continue
                q=2 if p==1 else 1
                new[S[:j]+(q,)+S[j+1:]]+=c
                if not bit[v]:
                    new[S[:j]+(q+2,)+S[j+1:]]+=c
        states=new
    return sum(c for S,c in states.items()
        if all(x>=3 or (x==0 and y) for x,y in zip(S,empty)) and accept(S))

def profile(a, b, k, arc=True):
    m=len(a)
    odd=lambda x:x==3
    even=lambda S:S==(4,)*len(S)
    single=sum(runs(a,b,k,((j,j),),even,(False,)) for j in (0,1))
    separate=runs(a,b,k,((0,0),(1,1)),even,(False,False))
    cross=runs(a,b,k,((0,1),(1,0)),lambda S:sum(map(odd,S))==1,(arc,not arc))
    return k**(m+2)+2*k*k*single+4*k*k*separate+2*k*cross

def quotient_value(a,b,arc=True):
    return 259*profile(a,b,2,arc)-209*profile(a,b,1,arc)

def out(a,b,arc):
    m=len(a); n=m+2
    A=[sum(1<<j for j in range(i+1,m)) for i in range(m)]+[0,0]
    for l,w in enumerate((a,b),m):
        for i,v in enumerate(w):
            if v: A[l]|=1<<i
            else: A[i]|=1<<l
    if arc:A[m]|=1<<(m+1)
    else:A[m+1]|=1<<m
    return A

def subset_h(A):
    n=len(A); D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for S in range(1,1<<n):
        for i,c in enumerate(D[S]):
            if not c:continue
            t=A[i]&~S
            while t:
                b=t&-t;t-=b
                D[S|b][b.bit_length()-1]+=c
    return [1]+[sum(row) for row in D[1:]]

def L(w,v):
    m=len(w)
    return sum((w[i]*v[j]+v[i]*w[j])*Fraction(1,2)**i*Fraction(3,4)**(j-i-1)
        for i in range(m) for j in range(i+1,m))

def phase(a,b):
    return Fraction(777,128)*L(a,b)*L(tuple(1-x for x in a[::-1]),tuple(1-x for x in b[::-1]))

def choose(q):
    """Choose even m and target-coefficient in [64,1024), then place a box."""
    assert q>=64
    m=0
    while q>=1024*4**m: m+=2
    m=max(m,2)
    left=Fraction(q,1554*4**m)*4 # represent a sum in [0,4]
    assert 0<=left<=4
    x=Fraction(1)
    a=[]
    for _ in range(m//2):
        if left>x*3:
            left-=x;a.append(1)
        else:a.append(0)
        x*=Fraction(3,4)
    a+= [0]*(m-len(a))
    b=[1]*(m//2)+[0]*(m-m//2)
    return a,b,quotient_value(a,b)

def choose_residue(q):
    """Place the macroscopic target and match its residue mod 3 exactly."""
    m=8
    while q>=1024*4**m: m+=4
    left=Fraction(q,1554*4**m)*4
    assert 0<=left<=4
    x=Fraction(1)
    prefix=[]
    for _ in range(m//2-4):
        if left>x*3:
            left-=x;prefix.append(1)
        else:prefix.append(0)
        x*=Fraction(3,4)
    b=[1]*(m//2)+[0]*(m//2)
    candidates=[]
    for repair in product((0,1),repeat=4):
        a=prefix+list(repair)+[0]*(m//2)
        value=quotient_value(a,b)
        if value%3==q%3:
            candidates.append((abs(value-q),a,b,value))
    _,a,b,value=min(candidates)
    return a,b,value

def controls():
    count=0
    for m in range(5):
        for bits in product((0,1),repeat=2*m):
            a,b=bits[:m],bits[m:]
            for arc in (False,True):
                h=subset_h(out(a,b,arc));full=len(h)-1
                assert profile(a,b,1,arc)==h[-1]
                assert profile(a,b,2,arc)==sum(h[S]*h[full^S] for S in range(full+1))
                assert (quotient_value(a,b,arc)%3)==0 if (m%2) else ((quotient_value(a,b,arc)+h[-1])%3==0)
                count+=1
    rng=Random(571)
    for m in range(5,10):
        for _ in range(12):
            a,b=tuple(tuple(rng.randrange(2) for i in range(m)) for j in range(2))
            for arc in (True,False):
                h=subset_h(out(a,b,arc));full=len(h)-1
                assert profile(a,b,1,arc)==h[-1]
                assert profile(a,b,2,arc)==sum(h[S]*h[full^S] for S in range(full+1))
                count+=1
    print(f"Profile transfer matched direct subset path counts on {count} two-apex hosts.")
    for n in (6,10,20,40,80):
        t=123*4**n
        a,b,value=choose(t)
        rel=float(Fraction(value-t,4**n))
        print(f"m={len(a)} q/4^m error at phase 123: {rel:+.10g}")

if __name__=="__main__":
    if len(argv)>1:
        q=int(argv[1]); a,b,value=choose_residue(q)
        print({"target":q,"value":value,"error":value-q,"spine_rows":(a,b),"a_beats_b":True})
    else:
        controls()
