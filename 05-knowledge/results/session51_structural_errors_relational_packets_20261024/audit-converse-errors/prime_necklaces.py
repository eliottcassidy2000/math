#!/usr/bin/env python3
"""Independent finite controls of the prime-cell nonrigid fixed-class extraction."""
from functools import lru_cache
from itertools import combinations,product
from collections import Counter
from math import comb
import json

def need(x,msg):
    if not x:raise RuntimeError(msg)
@lru_cache(None)
def symmetric(capacities):
    if not capacities:return 1
    first,*rest=capacities;ans=0
    for diagonal in (0,1):
        amount=first-diagonal
        if not 0<=amount<=len(rest):continue
        for chosen in combinations(range(len(rest)),amount):
            nxt=list(rest)
            for j in chosen:nxt[j]-=1
            if min(nxt,default=0)>=0:ans+=symmetric(tuple(sorted(nxt)))
    return ans
@lru_cache(None)
def bipartite(rows,d,capacities):
    if not rows:return int(all(c==0 for c in capacities))
    ans=0
    for chosen in combinations(range(len(capacities)),d):
        nxt=list(capacities)
        for j in chosen:nxt[j]-=1
        if min(nxt)>=0 and max(nxt)<=rows-1:ans+=bipartite(rows-1,d,tuple(sorted(nxt)))
    return ans

def move(M,a,b):
    p=len(M)
    return tuple(tuple(M[(i-a)%p][(j-b)%p] for j in range(p)) for i in range(p))
def converse(M):
    p=len(M)
    return tuple(tuple(M[-j%p][-i%p] for j in range(p)) for i in range(p))
def run():
    records=[]
    for p in (3,5,7):
        for d in range(1,(p+1)//2):
            words={frozenset(c) for c in combinations(range(p),d)};necklaces=[]
            while words:
                D=min(words,key=lambda x:tuple(sorted(x)));orb={frozenset((x+a)%p for x in D) for a in range(p)}
                need(len(orb)==p and orb<=words,'proper prime subset free under rotations');words-=orb;necklaces.append(D)
            counts=Counter();nonrigid=[]
            for slope,D in product(range(1,p),necklaces):
                M=tuple(tuple(int((j-slope*i)%p in D) for j in range(p)) for i in range(p))
                orb={move(M,a,b) for a,b in product(range(p),repeat=2)}
                is_sc=converse(M) in orb;counts[slope]+=is_sc
                need(is_sc==(slope==1 or (slope==p-1 and any(frozenset((-x+a)%p for x in D)==D for a in range(p)))),'full matrix converse iff slope and reflection test')
                need(len(orb)==p,'nonrigid prime matrix orbit size')
                nonrigid.append((slope,is_sc))
            B=comb(p,d);c=comb((p-1)//2,d//2);S=symmetric((d,)*p);A=bipartite(p,d,(d,)*p)
            need(sum(counts.values())==B//p+c,'nonrigid self converse formula')
            need(counts[1]==B//p and counts[p-1]==c,'separate plus and minus slope self converse counts')
            rigid=(A-(p-1)*B)//(p*p);fixed=(S-B-c)//p
            need((A-(p-1)*B)%(p*p)==0 and (S-B-c)%p==0,'rigid counts integral')
            need((rigid+fixed)%2==0,'merged rigid count integral')
            records.append({'p':p,'d':d,'A':A,'S':S,'nonrigid_classes':len(nonrigid),'nonrigid_self_converse':sum(counts.values()),'nonrigid_self_converse_by_slope':dict(sorted(counts.items())),'rigid_ordinary':rigid,'rigid_self_converse':fixed,'rigid_merged':(rigid+fixed)//2})
    return {'status':'PASS','scope':'finite p=3,5,7 coefficient and nonrigid orbit controls; proof of general extraction in audit report','records':records}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
