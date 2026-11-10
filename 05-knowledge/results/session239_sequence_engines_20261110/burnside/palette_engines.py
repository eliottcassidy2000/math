#!/usr/bin/env python3
"""Integral future-response engines for graphs and reversal palettes."""
from collections import defaultdict, Counter
from math import factorial,gcd
from fractions import Fraction
from itertools import combinations
from engines import require

def palette_series(N,q,t):
    """q pair states, t reversal-fixed states; return coefficients 0..N.
    Graphs (2,2), loopless digraphs (4,2), oriented graphs (3,1),
    tournaments (2,0). Only realizable finite involution palettes accepted.
    """
    require(isinstance(N,int) and N>=0,'nonnegative integer order')
    require(isinstance(q,int) and isinstance(t,int) and q>=t>=0 and (q-t)%2==0,'involution palette')
    fact=[factorial(i) for i in range(N+1)]
    lengths=list(range(1,N+1,2) if t==0 else range(1,N+1))[::-1]
    states={(0,(0,)*len(lengths)):1};stats={'peak_states':1,'transitions':0}
    for i,ell in enumerate(lengths):
        new=defaultdict(int);future=lengths[i+1:]
        for (d,response),w in states.items():
            old=dict(zip((j for j in lengths[i:] if j<=N-d),response))
            for k in range((N-d)//ell+1):
                D=d+ell*k
                du=k*((ell-1)//2)+ell*k*(k-1)//2+k*old.get(ell,0)
                dv=k if ell%2==0 else 0
                C=fact[D]//(fact[d]*ell**k*fact[k])
                R=tuple(old[j]+k*gcd(ell,j) for j in future if j<=N-D)
                new[D,R]+=w*C*q**du*t**dv
                stats['transitions']+=1
        states={key:w for key,w in new.items() if w}
        stats['peak_states']=max(stats['peak_states'],len(states))
    answer=[0]*(N+1)
    for (n,R),w in states.items():
        require(not R and w%fact[n]==0,'palette terminal integrality')
        answer[n]=w//fact[n]
    return answer,stats

def partitions(n,minimum=1):
    if n==0:yield ()
    for k in range(minimum,n+1):
        for rest in partitions(n-k,k):yield (k,)+rest

def direct_palette(N,q,t):
    """Independent uncompressed rational partition sum, all lengths."""
    answer=[]
    for n in range(N+1):
        total=Fraction(0)
        for p in partitions(n):
            u=sum((ell-1)//2 for ell in p)+sum(gcd(a,b) for a,b in combinations(p,2))
            v=sum(ell%2==0 for ell in p)
            z=1
            for ell,k in Counter(p).items():z*=ell**k*factorial(k)
            total+=Fraction(q**u*t**v,z)
        require(total.denominator==1,'direct palette integrality')
        answer.append(total.numerator)
    return answer
