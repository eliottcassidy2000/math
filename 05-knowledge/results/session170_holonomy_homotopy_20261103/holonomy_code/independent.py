#!/usr/bin/env python3
"""Independent dense Cayley-Dickson and subset-path controls. No check.py imports."""
import json
import itertools
from functools import lru_cache

def check(condition, message):
    if not condition:
        raise ValueError(message)

def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def conj(a):return (a[0],)+tuple(-x for x in a[1:])
@lru_cache(None)
def multiply(a,b):
    if len(a)==1:return (a[0]*b[0],)
    m=len(a)//2
    A,B,C,D=a[:m],a[m:],b[:m],b[m:]
    return sub(multiply(A,C),multiply(conj(D),B))+add(multiply(D,A),multiply(B,conj(C)))

def table(n):
    size=1<<n
    units=[tuple(int(i==j) for i in range(size)) for j in range(size)]
    result={}
    for a in range(size):
        for x in range(size):
            value=multiply(units[a],units[x])
            support=[i for i,v in enumerate(value) if v]
            check(len(support)==1 and abs(value[support[0]])==1,'not a signed basis product')
            result[a,x]=(support[0],int(value[support[0]]<0))
    return result

def apply(word,start,T):
    x,s=start,0
    for a,direction in word:
        x,b=T[a,x]
        s^=b^int(direction<0)
    return x,s

def inverse_word(word):return [(a,-d) for a,d in reversed(word)]
def diagonal(word,n,T):
    half=1<<(n-1)
    result=[]
    for x in range(2*half):
        y,s=apply(word,x,T)
        check(y==x,'not a diagonal word')
        result.append(s)
    check(result[:half]==result[half:],'top-bit dependence')
    return result[:half]

def decode(values,d,k):
    coefficients={}
    for A in range(1<<d):
        if A.bit_count()<=k:
            total=0
            B=A
            while True:
                total^=values[B]
                if not B:break
                B=(B-1)&A
            coefficients[A]=total
    return [sum(c for A,c in coefficients.items() if y&A==A)%2 for y in range(1<<d)]

def seed_profile(bits):
    N=5
    edges={pair:bool(bits>>i&1) for i,pair in enumerate(itertools.combinations(range(N),2))}
    def arrow(a,b):return edges[a,b] if a<b else not edges[b,a]
    @lru_cache(None)
    def ending(S,v):
        if S==1<<v:return 1
        return sum(ending(S^(1<<v),u) for u in range(N) if S>>u&1 and u!=v and arrow(u,v))
    def paths(S):return sum(ending(S,v) for v in range(N) if S>>v&1)
    @lru_cache(None)
    def covers(S,k):
        if not S:return int(k==0)
        if not k:return 0
        v=S&-S
        subtotal,Q=0,S
        while Q:
            if Q&v:subtotal+=paths(Q)*covers(S^Q,k-1)
            Q=(Q-1)&S
        return subtotal
    return [covers(31,k) for k in range(1,6)]

def main():
    report={'dense_levels':[]}
    seed=0x170
    def next_int(bound):
        nonlocal seed
        seed=(1664525*seed+1013904223)&0xffffffff
        return seed%bound
    for n in [3,4,5]:
        T=table(n)
        M,d=1<<(n-1),n-1
        W=[(1,1),(M+2,1),(1,-1),(M+2,-1)]
        layer=[]
        for r in range(2,d+1):
            if r>2:
                a=1<<(r-1)
                W=[(a,1)]+W+[(a,-1)]+inverse_word(W)
            h=diagonal(W,n,T)
            check(h==[int(y<1<<r) for y in range(M)],'nested flat witness')
            check(decode(h,d,d-r)==h,'optimal interpolation')
            layer.append({'r':r,'word_length':len(W),'support_size':sum(h)})
        distinct=set()
        for trial in range(256):
            labels=[1+next_int(2*M-1) for _ in range(4)]
            bits=[0]*4
            word=[]
            for _ in range(12):
                i=next_int(4)
                word.append((labels[i],-1 if bits[i] else 1))
                bits[i]^=1
            order=sorted(range(4),key=lambda i:next_int(1<<20))
            for i in order:
                if bits[i]:word.append((labels[i],-1))
            h=diagonal(word,n,T)
            check(sum(h)%2==0,'balanced parity')
            check(all(sum(h[y] for y in range(M) if y>>j&1)%2==0 for j in range(d)),'balanced moment')
            check(decode(h,d,d-2)==h,'balanced interpolation')
            distinct.add(tuple(h))
        report['dense_levels'].append({'n':n,'basis_products':4*M*M,'balanced_paths':256,
                                       'distinct_balanced_returns':len(distinct),'layers':layer})
        if n==4:
            word=[(4,1),(1,1),(10,1),(1,-1),(10,-1),(4,-1)]
            h=diagonal(word,n,T)
            check([i for i,b in enumerate(h) if b]==[4,5,6,7],'shifted hostile')
            report['shifted_hostile']={'e0':apply(word,0,T),'e4':apply(word,4,T),'support':[4,5,6,7]}
    profiles={str(b):seed_profile(b) for b in [8,10]}
    check(all(p==[9,33,31,10,1] for p in profiles.values()),'independent seed profile')
    report['subset_path_profiles']=profiles
    report['status']='PASS'
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
