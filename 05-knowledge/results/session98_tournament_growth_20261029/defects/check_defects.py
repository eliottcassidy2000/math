#!/usr/bin/env python3
"""Exact missing-comparison census; standard library, no external data."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import factorial, gcd, lcm, comb
import json

def require(b, msg):
    if not b: raise RuntimeError(msg)
def partitions(n, lo=1):
    if not n: yield ()
    for j in range(lo,n+1):
        for tail in partitions(n-j,j): yield (j,)+tail

def mul(a,b,D):
    c=[0]*(D+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:D+1-i]): c[i+j]+=x*y
    return c

def factor(length,D):
    a=[0]*(D+1); a[0]=2
    if length<=D:a[length]=1
    return a

def closed_profile(lam,D):
    debt=sum(j//2 for j in lam if j%2==0)
    if debt>D:return [0]*(D+1)
    a=[0]*(D+1);a[debt]=1
    for j in lam:
        for _ in range((j-1)//2):a=mul(a,factor(j,D),D)
    for i,j in enumerate(lam):
        for k in lam[i+1:]:
            for _ in range(gcd(j,k)):a=mul(a,factor(lcm(j,k),D),D)
    return a

def census(n,D):
    a=[Fraction(0) for _ in range(D+1)]
    for lam in partitions(n):
        z=1
        for j,m in Counter(lam).items():z*=j**m*factorial(m)
        for d,x in enumerate(closed_profile(lam,D)):a[d]+=Fraction(x,z)
    require(all(x.denominator==1 for x in a),'integrality')
    return [int(x) for x in a]

def signed_profile(p,D):
    """Independent literal pair-orbit walk, tracking orientation reversal."""
    pairs=list(__import__('itertools').combinations(range(len(p)),2)); unseen=set(pairs)
    a=[1]+[0]*D
    while unseen:
        first=min(unseen); u,v=first; length=0;flip=0
        while True:
            unseen.discard(tuple(sorted((u,v))))
            u,v=p[u],p[v];length+=1
            if tuple(sorted((u,v)))==first:
                flip=int((u,v)!=first);break
        if flip:
            b=[0]*(D+1)
            if length<=D:b[length]=1
        else:b=factor(length,D)
        a=mul(a,b,D)
    return a

def cycle_type(p):
    unseen=set(range(len(p)));out=[]
    while unseen:
        v=min(unseen);start=v;k=0
        while True:
            unseen.remove(v);k+=1;v=p[v]
            if v==start:break
        out.append(k)
    return tuple(sorted(out))

def transform(state,p,pairs,index):
    a=[0]*len(state)
    for (u,v),x in zip(pairs,state):
        u,v=p[u],p[v]
        if u>v:u,v=v,u;x=-x
        a[index[u,v]]=x
    return tuple(a)

def brute(n):
    pairs=list(__import__('itertools').combinations(range(n),2)); idx={p:i for i,p in enumerate(pairs)}
    ps=list(permutations(range(n)));reps=set();hist=[0]*(len(pairs)+1);checks=0
    for state in product((-1,0,1),repeat=len(pairs)):
        imgs=[transform(state,p,pairs,idx) for p in ps]
        rep=min(imgs)
        if rep in reps:continue
        reps.add(rep);d=state.count(0);hist[d]+=1
        aut=sum(x==state for x in imgs);v2=(aut&-aut).bit_length()-1
        require(v2<=d,'two-primary symmetry budget');checks+=1
    return hist,checks

def main():
    checks=0
    for n in range(8):
        D=min(4,n*(n-1)//2)
        for p in permutations(range(n)):
            require(signed_profile(p,D)==closed_profile(cycle_type(p),D),'signed orbit profile');checks+=1
    brute_hist={}
    for n in range(5):
        hist,k=brute(n);require(hist==census(n,len(hist)-1),'full brute census');checks+=k;brute_hist[n]=hist
    known=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168]
    require([census(n,0)[0] for n in range(len(known))]==known,'A000568')
    # Odd-cycle-only Burnside fails as soon as one comparison is erased.
    require(census(2,1)==[1,1],'one-hole hostile')
    # All even-cycle mass <=2d, equivalently obligatory debt <=d.
    for n in range(15):
        for lam in partitions(n):
            debt=sum(j//2 for j in lam if j%2==0)
            require(all(x==0 for x in closed_profile(lam,3)[:min(debt,4)]),'debt pruning');checks+=1
    one_hole=[]
    for n in range(2,16):
        edge=Fraction(0);vertex=Fraction(0)
        for lam in partitions(n):
            if any(j%2==0 for j in lam):continue
            z=1
            for j,m in Counter(lam).items():z*=j**m*factorial(m)
            edge+=Fraction(closed_profile(lam,0)[0]*comb(lam.count(1),2),z)
        for lam in partitions(n-1):
            if any(j%2==0 for j in lam):continue
            z=1
            for j,m in Counter(lam).items():z*=j**m*factorial(m)
            vertex+=Fraction(closed_profile(lam,0)[0]*lam.count(1),z)
        require(edge.denominator==vertex.denominator==1,'marked integral')
        require(2*census(n,1)[1]==edge+vertex,'one hole edge and twin recursion')
        one_hole.append([n,int(edge),int(vertex),census(n,1)[1]])
    rows={n:census(n,3) for n in range(21)}
    print(json.dumps({'status':'PASS','checks':checks,'brute_histograms':brute_hist,'defect_0_to_3':rows,'one_hole_edge_vertex':one_hole},sort_keys=True))
if __name__=='__main__':main()
