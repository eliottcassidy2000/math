#!/usr/bin/env python3
"""Exact signed binary-height constructors. Outgoing masks use zero-based arcs."""
import argparse
import csv
import json
from itertools import combinations
from math import isqrt
from pathlib import Path

# (D rows, I mask, G, intercept, power shift j), H=G*2^width+intercept.
REC={
 'dyadic+1':((8,1,3,6,15,31,63),3,128,1,7),
 'dyadic+7':((2,132,9,35,15,23,63,125),38,896,7,7),
 'dyadic-1':((32,5,81,7,75,30,171,63),100,1024,-1,10),
 'dyadic-7':((8,1,3,166,15,23,63,119),80,224,-7,5),
 '3+1':((4,1,2),3,6,1,1),
 '3+7':((24,17,3,6,12),5,42,7,1),
 '5+1':((4,1,18,7,11),24,20,1,2),
 '5+7':((32,1,3,7,15,30),48,35,7,0),
 '7+1':((8,1,3,102,15,23,55),48,112,1,4),
 '7+7':((4,17,10,35,13,151,63,95),144,392,7,3),
 '13+1':((16,1,3,7,206,31,47,111),96,208,1,4),
 '13+7':((8,33,3,6,15,29),33,91,7,0),
}

def check(rows):
    n=len(rows)
    assert all(0<=r<1<<n and not r>>i&1 for i,r in enumerate(rows))
    assert all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i,j in combinations(range(n),2))

def dp(rows):
    n=len(rows);E=[[0]*n for _ in range(1<<n)]
    for i in range(n):E[1<<i][i]=1
    for s in range(1,1<<n):
        for i,a in enumerate(E[s]):
            if not a:continue
            Q=rows[i]&~s
            while Q:
                z=Q&-Q;Q^=z;E[s|z][z.bit_length()-1]+=a
    return E
def reads(rows,I):
    check(rows);n=len(rows);f=(1<<n)-1
    E=dp(rows);B=dp([sum(1<<i for i,r in enumerate(rows) if r>>j&1) for j in range(n)])
    G=Q=0
    for s in range(1,f+1):
        e=sum(x for i,x in enumerate(E[s]) if I>>i&1);t=f^s
        G+=e*(sum(B[t]) if t else 1)
        Q+=e*(sum(x for j,x in enumerate(B[t]) if not I>>j&1) if t else 1)
    return sum(E[f]),G,Q
def wrap(D,I,a):
    n=len(D);b=n+a
    return tuple(r|(((I>>i)&1)<<b) for i,r in enumerate(D))+tuple(
        ((1<<n)-1)|(((1<<b)-1)^((1<<(i+1))-1)) for i in range(n,b))+(
        ((1<<b)-1)^I,)
def ordinal_pad(A,k):
    n=len(A);e=n+k
    extra=((1<<k)-1)<<n
    return tuple(r|extra for r in A)+tuple(((1<<e)-1)^((1<<(i+1))-1) for i in range(n,e))

def pair(A,name):
    """Two exact readings, with the same incoming arbitrary module A."""
    key='dyadic' if name.startswith('dyadic') else name
    sign='-' if name.endswith('-') else '+'
    L=REC[key+sign+'1'];R=REC[key+sign+'7']
    j=max(L[4],R[4]);beta=L[3]
    alpha=L[2]<< (j-L[4])
    assert R[3]==7*beta and (R[2] <<(j-R[4]))==7*alpha
    return (alpha,beta,compose(L[0],L[1],ordinal_pad(A,j-L[4])),
            compose(R[0],R[1],ordinal_pad(A,j-R[4])))

def compose(D,I,A):
    n=len(D);b=n+len(A)
    return (tuple(r|(((I>>i)&1)<<b) for i,r in enumerate(D))
            +tuple((r<<n)|((1<<n)-1) for r in A)+(((1<<b)-1)^I,))

def count(rows):
    check(rows);return sum(dp(rows)[-1])

def birthday(k):
    n=k+2;R=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    R[0]^=1<<(n-1);R[-1]|=1
    return tuple(R)

def seeds():
    with (Path(__file__).resolve().parent/'small.tsv').open() as f:
        return {int(a['count']):tuple(map(int,a['outgoing_masks'].split(',')))
                for a in csv.DictReader(f,delimiter='\t')}
def resolve(a,k,sign=1,scale=1):
    if a not in (1,2,3,4,5,7,13) or k<1 or scale not in (1,7) or sign not in (-1,1):
        raise ValueError('outside known address families')
    if sign<0 and a!=1:raise ValueError('negative branch is for dyadic heights')
    N=scale*(a*2**k+sign)
    if N in (7,21):return N,None
    if a in (1,2,4):
        k+=(a.bit_length()-1)
        if sign>0 and scale==1:return N,birthday(k)
        name='dyadic'+('+' if sign>0 else '-')+str(scale)
    else:name=f'{a}+{scale}'
    D,I,G,b,j=REC[name]
    z=k-j
    if z<1:return N,seeds()[N]
    assert G*2**z+b==N
    return N,wrap(D,I,z)

def required():
    N=set()
    for a in (1,2,3,4,5,7,13):
        for k in range(1,13):
            for sign in ((-1,1) if a==1 else (1,)):
                for scale in (1,7):
                    v=scale*(a*2**k+sign)
                    if v in (7,21):continue
                    if a in (1,2,4):
                        if sign>0 and scale==1:continue
                        z=k+a.bit_length()-1-REC['dyadic'+('+' if sign>0 else '-')+str(scale)][4]
                    else:z=k-REC[f'{a}+{scale}'][4]
                    if z<1:N.add(v)
    return N

def verify():
    for name,(D,I,G,b,j) in REC.items():
        d,g,q=reads(D,I)
        assert g==G and d-g+q==b
        for z in (1,2):
            assert count(wrap(D,I,z))==G*2**z+b
        print(name,'h,G,Q',d,g,q,'height shift',j)
    tab=seeds()
    assert tab.keys()==required()
    for N,R in tab.items():assert count(R)==N
    print('checked tail-start exceptions and seeds:',sorted(tab))
    for a,k in [(3,36),(5,39)]:
        p=a*2**k+1
        assert all(p%d for d in range(3,isqrt(p)+1,2))
        print('certified prime pair',p,7*p,'orders',len(resolve(a,k)[1]),len(resolve(a,k,1,7)[1]))
    for name in ('dyadic+','dyadic-','3','5','7','13'):
        A=(2,4,1)
        alpha,beta,L,R=pair(A,name)
        assert count(L)==alpha*12+beta*3
        assert count(R)==7*count(L)
        print('arbitrary core pair', name, alpha,beta)
    A=ordinal_pad((2,4,9,3),32)
    p=192*2**32+5
    assert all(p%d for d in range(3,isqrt(p)+1,2))
    _,_,L,R=pair(A,'3')
    check(L);check(R)
    print('nonunit-offset prime pair',p,7*p,'orders',len(L),len(R))
    k=61;p=2**k-1;r=4
    for _ in range(k-2):r=(r*r-2)%p
    assert r==0
    print('Lucas-Lehmer prime pair',p,7*p,'orders',len(resolve(1,k,-1)[1]),len(resolve(1,k,-1,7)[1]))

if __name__=='__main__':
    cli=argparse.ArgumentParser();cli.add_argument('a',nargs='?',type=int)
    cli.add_argument('k',nargs='?',type=int);cli.add_argument('--sign',type=int,default=1)
    cli.add_argument('--scale',type=int,default=1)
    o=cli.parse_args()
    if o.a is None:verify()
    else:
        n,rows=resolve(o.a,o.k,o.sign,o.scale)
        print(json.dumps(dict(H=n,status='known hole' if rows is None else 'exact witness',rows=rows)))
