#!/usr/bin/env python3
"""Exact new-family controls; no unlabeled-poset or conjecture census."""
from functools import lru_cache
from itertools import combinations, product
from math import comb
from pathlib import Path
import json

def require(c,msg):
    if not c: raise AssertionError(msg)

def fib(k):
    if k == -1: return 1
    a,b=0,1
    for _ in range(k):a,b=b,a+b
    return a

def close(n,rel):
    p=[0]*n
    for a,b in rel:
        require(type(a) is int and type(b) is int and 0<=a<n and 0<=b<n and a!=b,'bad relation')
        p[b]|=1<<a
    for k in range(n):
        for j in range(n):
            if p[j]>>k&1:p[j]|=p[k]
    return tuple(p)

def count(p):
    if any(x>>i&1 for i,x in enumerate(p)):return 0
    n=len(p);full=(1<<n)-1
    @lru_cache(None)
    def f(s):
        if s==full:return 1
        return sum(f(s|1<<i) for i in range(n) if not(s>>i&1) and p[i]&s==p[i])
    return f(0)

def extensions(p):
    n=len(p);full=(1<<n)-1
    def rec(s,seq):
        if s==full:
            yield seq;return
        for i in range(n):
            if not(s>>i&1) and p[i]&s==p[i]:yield from rec(s|1<<i,seq+(i,))
    return rec(0,())

def table(n,rel,e,f):
    return [[count(close(n,rel+[e if a else e[::-1],f if b else f[::-1]])) for b in (True,False)] for a in (True,False)]

def deficits(t):
    return [[sum(t[1-i])-t[i][j] for j in range(2)] for i in range(2)]

def gold(t):return min(sum(deficits(t),[]))>=0

def thin(n):return [(i,j) for i in range(n) for j in range(i+2,n)]

def terminal_context(n,r,m,D,core=None):
    rel=thin(n) if core is None else list(core)
    rel += [(n+i,n+i+1) for i in range(r-1)]
    rel += [(d,n+r-1) for d in D]
    rel += [(n+r+i,n+r+i+1) for i in range(m-1)]
    if m:rel += [(v,n+r) for v in range(n)]
    return n+r+m,rel

def boundary_claims(n,r,m):
    T=comb(n+r+m,r);a=fib(n-1);b=fib(n-2);q=comb(r+2,2)
    negative=[[a*(T-1),b*(T-1)],[a*(T-r-1),0]]
    positive=[[a*(T-r-1),b*(T-q)],[a*(T-1),0]]
    A=(T-1)*fib(n-1)-r*fib(n-3)
    B=(T-1)*fib(n-2)-r*fib(n-4)
    right=[[A,B],[A,0]]
    return negative,positive,right

def main():
    checks=0; examples=[]
    for n,r,m in product(range(3,11),range(1,8),range(0,6)):
        neg,pos,right=boundary_claims(n,r,m)
        for mark in (0,1):
            nn,rel=terminal_context(n,r,m,{mark})
            actual=table(nn,rel,(0,1),(1,2))
            expected=neg if mark==0 else pos
            require(actual==expected,('left table',n,r,m,mark,actual,expected));checks+=1
            require(sum(sum(actual,[]))==count(close(nn,rel)),'total');checks+=1
            if mark==0:
                actual_right=table(nn,rel,(n-2,n-1),(n-3,n-2))
                require(actual_right==right and gold(right),'right certificate');checks+=1
                require(min(sum(deficits(actual),[]))==-r*fib(n-1),'exact defect');checks+=1
            else:
                require(gold(actual),'positive certificate');checks+=1
                strict=min(sum(deficits(actual),[]))>0
                require(strict != (r==1 and m==0 and n in (3,4)),'strict boundary');checks+=1
        if (n,r,m) in ((3,1,0),(6,4,0),(6,4,10)):
            examples.append(dict(N=n,r=r,m=m,left_negative=neg,left_positive=pos,right_negative=right))
    # Independent core-fibre check, including non-singleton ideals and arbitrary small cores.
    fibre_cases=0
    for n in range(2,6):
        pairs=list(combinations(range(n),2));cores={close(n,[e for k,e in enumerate(pairs) if mask>>k&1]) for mask in range(1<<len(pairs))}
        for pred in cores:
            rel=[(i,j) for j in range(n) for i in range(n) if pred[j]>>i&1]
            ideals=[s for s in range(1<<n) if all(not(s>>i&1) or pred[i]&s==pred[i] for i in range(n))]
            for Dmask in ideals:
                D={i for i in range(n) if Dmask>>i&1}
                if not D:continue
                for sigma in extensions(pred):
                    d=max(sigma.index(x)+1 for x in D)
                    for r,m in ((1,0),(2,1),(4,2)):
                        nn,rr=terminal_context(n,r,m,D,rel)
                        rr += list(zip(sigma,sigma[1:]))
                        v=count(close(nn,rr));want=comb(n+r+m,r)-comb(r+d-1,r)
                        require(v==want,('fibre',n,D,sigma,r,m,v,want));fibre_cases+=1
    # Genuine five-vertex endpoint-locality obstruction and all bridge lengths.
    bridge=[]
    for ell in range(0,15):
        n=4+ell;rel=[(0,2),(0,3),(1,2),(1,3)]
        if ell:rel += [(0,4)]+[(i,i+1) for i in range(4,n-1)]+[(n-1,2)]
        got=table(n,rel,(0,1),(2,3));want=[[ell+1,comb(ell+2,2)],[1,ell+1]]
        require(got==want and gold(got)==(ell<=1),'bridge');checks+=1
        good=[]
        pairs=list(combinations(range(4),2))
        for e,f in product(pairs,repeat=2):
            if gold(table(n,rel,e,f)):good.append((e,f))
        if ell==1:require(good and all(len(set(e+f))==4 for e,f in good),'endpoint uncrossing')
        if ell>=2:require(not good,'all endpoint strategies fail')
        bridge.append(dict(length=ell,table=got,endpoint_certificates=good))
    # Explicit extreme flatness: normalized max/min fibre distortion is r/(T-r-1).
    flat=[]
    for r in (1,2,5,10):
        m=100;T=comb(m+6+r,r)
        flat.append(dict(r=r,m=m,T=T,unnormalized_defect=r*fib(5),distortion_excess_numerator=r,distortion_excess_denominator=T-r-1))
    out={'status':'exact finite controls for proved parameterized statements','family_parameter_box':{'N':[3,10],'r':[1,7],'m':[0,5]},'family_and_bridge_checks':checks,'core_fibre_cases':fibre_cases,'examples':examples,'bridge':bridge,'flatness_examples':flat}
    # Receipt is printed to stdout; the published copy is in results.
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
