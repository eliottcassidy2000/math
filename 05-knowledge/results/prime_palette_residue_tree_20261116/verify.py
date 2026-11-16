#!/usr/bin/env python3
import itertools as it
from math import gcd,isqrt
import json
from pathlib import Path

P=(3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61)
B=(11,13,17,19)
def isprime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))
def axes(P):
    A=[(8,2,{1:0,3:1,5:0,7:1}),(8,2,{1:0,3:1,5:1,7:0})]
    for p in P:
        m=p*p
        for g in range(2,m):
            d={pow(g,j,m):j for j in range(p*(p-1))}
            if len(d)==p*(p-1) and 0 not in d and all(k%p for k in d):
                break
        else: raise ValueError(p)
        A.append((m,p*(p-1),d))
    return A
def rank(R,l):
    R=[list(map(lambda x:x%l,r)) for r in R]
    k=0
    for j in range(len(R[0])):
        i=next((i for i in range(k,len(R)) if R[i][j]),None)
        if i is None: continue
        R[k],R[i]=R[i],R[k]
        inv=pow(R[k][j],-1,l)
        R[k]=[(v*inv)%l for v in R[k]]
        for i in range(k+1,len(R)):
            a=R[i][j]
            R[i]=[(x-a*y)%l for x,y in zip(R[i],R[k])]
        k+=1
    return k
def signature(P,B):
    A=axes(P)
    return {l:(rank([[d[b%m] for b in B] for m,n,d in A if n%l==0],l),sum(n%l==0 for m,n,d in A)) for l in range(2,max(P)+1) if isprime(l) and any(n%l==0 for m,n,d in A)}
def complete(P,B):
    return all(r==s for r,s in signature(P,B).values())
W=(67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,157,167)

def H(rows):
    n=len(rows); s=0
    for i,x in enumerate(rows):
        assert x>>n==0 and not (x>>i&1)
        assert all((x>>j&1)+(rows[j]>>i&1)==1 for j in range(n) if i!=j)
    for p in it.permutations(range(n)):
        s+=all(rows[i]>>j&1 for i,j in zip(p,p[1:]))
    return s

def main():
    assert complete((3,5),B) and complete(P,W)
    assert len(W)==len(P)+2
    print('Small bank:',B,signature((3,5),B))
    print('Wide bank:',W,signature(P,W))
    data=json.loads(Path(__file__).with_name('witnesses.json').read_text())
    req=set(B)|set(W)|({x for x in range(1,64,2)}-{7,21})|{343}|{7*p for p in P if p not in (3,7)}
    assert set(map(int,data))==req
    for k,v in data.items(): assert H(v)==int(k),(k,v)
    print(len(data),'literal witnesses independently verified by permutations')
if __name__=='__main__':
    main()
