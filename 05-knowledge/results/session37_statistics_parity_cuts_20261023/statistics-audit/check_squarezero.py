#!/usr/bin/env python3
from itertools import permutations,combinations
from math import factorial
from pathlib import Path
import json
OUT=Path(__file__).resolve().parent

def ck(x,msg):
    if not x:raise AssertionError(msg)
def mul(p,q):
    N=len(p);r=[0]*N
    for i,a in enumerate(p):
        if not a:continue
        for j,b in enumerate(q):
            if b and not(i&j):r[i|j]+=a*b
    return r
def pw(p,k):
    r=[1]+[0]*(len(p)-1)
    for _ in range(k):r=mul(r,p)
    return r
def inv(p):
    ck(p[0]==1,'inverse constant');h=[1-p[0]]+[-a for a in p[1:]];r=[1]+[0]*(len(p)-1);term=r[:]
    for _ in range((len(p)-1).bit_length()):term=mul(term,h);r=[a+b for a,b in zip(r,term)]
    ck(mul(p,r)==[1]+[0]*(len(p)-1),'inverse product');return r

def detper(A,sign=1,det=True):
    n=len(A);r=[0]*(1<<n)
    # Principal-minor expansion makes diagonal 1+Aii xi contributions explicit.
    for mask in range(1<<n):
        vs=[i for i in range(n)if mask>>i&1]
        for p in permutations(vs):
            value=sign**len(vs)
            for i,j in zip(vs,p):value*=A[i][j]
            if det:value*=(-1)**sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))
            r[mask]+=value
    return r

def pathpoly(A):
    n=len(A);r=[0]*(1<<n);r[0]=1
    for mask in range(1,1<<n):
        vs=[i for i in range(n)if mask>>i&1]
        for p in permutations(vs):
            x=1
            for i in range(len(p)-1):x*=A[p[i]][p[i+1]]
            r[mask]+=x
    return r

def packing_value(A,f,b):
    n=len(A);out=0
    for p in permutations(range(n)):
        seen=set();weight=1
        for i in range(n):
            if i in seen:continue
            j=i;k=0;edges=1
            while j not in seen:
                seen.add(j);k+=1
                if j!=p[j]:edges*=A[j][p[j]]
                j=p[j]
            if k>1:weight*=edges*(b+(-1)**(k-1)*f)
        out+=weight
    return out

def exterior_sign(a,b):
    if a&b:return 0
    s=0
    while b:
        low=b&-b;b^=low;s+=(a>>(low.bit_length())).bit_count()
    return (-1)**s

def main():
    cases=0;parameter_checks=0
    for n in range(5):
        edges=list(combinations(range(n),2));N=1<<n;pref=[1]*N
        for mask in range(1<<len(edges)):
            A=[[0]*n for _ in range(n)]
            for k,(i,j)in enumerate(edges):A[i][j]=(mask>>k)&1;A[j][i]=1-A[i][j]
            dplus=detper(A);dminus=detper(A,-1);boson=inv(dminus)
            ck(boson==detper(A,det=False),('MacMahon',n,mask))
            ck(mul(pref,mul(dplus,boson))==pathpoly(A),('path identity',n,mask))
            for f,b in [(0,0),(1,0),(0,1),(1,1),(2,1),(1,2),(2,2)]:
                z=mul(pref,mul(pw(dplus,f),pw(boson,b)))
                ck(z[-1]==packing_value(A,f,b),('f,b',n,mask,f,b));parameter_checks+=1
            cases+=1
    A=[[2,-3,1],[4,5,-2],[3,7,-1]]
    ck(inv(detper(A,-1))==detper(A,det=False),'weighted diagonal MacMahon')
    embedding=0
    for n in range(6):
        emb=lambda s:sum(3<<(2*i)for i in range(n)if s>>i&1)
        for s in range(1<<n):
            for t in range(1<<n):
                ck(exterior_sign(emb(s),emb(t))==(0 if s&t else 1),('even embedding',n,s,t));embedding+=1
    outs=[8,1,3,6];A=[[(outs[i]>>j)&1 for j in range(4)]for i in range(4)]
    paths=[p for p in permutations(range(4))if all(A[p[i]][p[i+1]]for i in range(3))]
    alt=sum((-1)**sum(p[i]>p[j]for i in range(4)for j in range(i+1,4))for p in paths)
    ck(len(paths)==5 and alt==-1,'raw exterior hostile')
    r={'status':'PASS','all_labelled_tournaments_through_n4':cases,'mixed_parameter_checks':parameter_checks,'paired_exterior_product_checks_through_n5':embedding,'arbitrary_signed_weighted_diagonal_MacMahon':True,'raw_odd_hostile':{'outmasks':outs,'H':len(paths),'H_alt':alt,'paths':paths}}
    (OUT/'squarezero.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
