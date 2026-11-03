#!/usr/bin/env python3
"""Check cyclic three-transitive-block run formula against subset DP."""
from math import factorial
from itertools import product

def stirling(n,k):
    row=[1]+[0]*n
    for i in range(1,n+1):
        row=[0]+[k*row[k]+row[k-1] for k in range(1,n+1)]
    return row[k]

def formula(a,b,c):
    dims=(a,b,c)
    return 1+sum(factorial(k)**3*(
        __import__('math').prod(stirling(t+1,k+1) for t in dims)+
        __import__('math').prod(stirling(t,k) for t in dims))
        for k in range(1,min(dims)+1))

def direct(a,b,c):
    n=a+b+c
    if n==0: return 1
    block=[0]*a+[1]*b+[2]*c
    def edge(i,j):
        return i<j if block[i]==block[j] else (block[j]-block[i])%3==1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for j in range(n):
            if dp[mask][j]:
                for k in range(n):
                    if not mask>>k&1 and edge(j,k):
                        dp[mask|1<<k][k]+=dp[mask][j]
    return sum(dp[-1])

if __name__=='__main__':
    cases=0
    for a,b,c in product(range(0,10),repeat=3):
        if a+b+c>9: continue
        actual=direct(a,b,c); predicted=formula(a,b,c)
        if actual != predicted: raise AssertionError((a,b,c,actual,predicted))
        if a>=1 and b==c==1 and actual!=2**a+1: raise AssertionError((a,b,c,actual))
        if a>=1 and b==2 and c==1 and actual!=3*2**a-1: raise AssertionError((a,b,c,actual))
        if a>=1 and b==c==2 and actual!=4*3**a+5*2**a-11: raise AssertionError((a,b,c,actual))
        cases+=1
    print('PASS: %d ordered nonnegative triples with total order <= 9' % cases)
    print('samples:',[(t, direct(*t)) for t in [(1,1,1),(2,2,2),(3,2,2),(3,3,3)]])
