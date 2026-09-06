#!/usr/bin/env python3
"""Independent exact audit; no producer imports or recurrence implementation."""
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement
from math import factorial, prod
from hashlib import sha256
import json

GATES=0
TRACE=[]
def require(ok, label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)

def elementary(values):
    return [sum((prod(c) for c in combinations(values,k)),Q(0)) for k in range(len(values)+1)]

def falling(n,r):
    return factorial(n)//factorial(n-r) if n>=r else 0

def kernel(M,t,d):
    return prod((M-a)**2-d*d for a in range(t))

def solve(matrix, rhs):
    n=len(matrix)
    a=[[Q(x) for x in row]+[Q(y)] for row,y in zip(matrix,rhs)]
    for j in range(n):
        i=next(i for i in range(j,n) if a[i][j])
        a[j],a[i]=a[i],a[j]
        v=a[j][j]
        a[j]=[x/v for x in a[j]]
        for i in range(n):
            if i!=j:
                v=a[i][j]
                a[i]=[x-v*y for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]

carrier_count=window_count=cone_count=0
for n in range(2,8):
    for roots in combinations_with_replacement((-2,-1,1,3),n-1):
        pre=elementary(roots)
        for k in range(1,n):
            if not pre[k] or not pre[k-1]:
                continue
            coeff=elementary((*roots,-pre[k]/pre[k-1]))
            require(coeff[k]==0 and coeff[0]*coeff[-1]!=0,'valid real-rooted carrier')
            require(coeff[k-1]*coeff[k+1]<0,'strict adjacent sign')
            carrier_count+=1
            for r in range(k):
                for s in range(n-k):
                    # Direct diagonal coefficient multiplier, independently of
                    # the producer's derivative/reversal construction.
                    q=[falling(j,r)*falling(n-j,s)*coeff[j] for j in range(r,n-s+1)]
                    target=k-r
                    require(q[target]==0,'zero preserved')
                    require(q[target-1]*q[target+1]<0,'target remains interior')
                    square=sum(q[j]*q[2*target-j] for j in range(len(q)) if 0<=2*target-j<len(q))
                    response=2*sum(kernel(k,r,d)*kernel(n-k,s,d)*coeff[k-d]*coeff[k+d]
                                   for d in range(1,min(k,n-k)+1))
                    require(square==response and square<0,'strict actual quadratic response')
                    window_count+=1
            TRACE.append((n,k,tuple(map(str,coeff))))

for M in range(1,8):
    matrix=[[kernel(M,t,d) for t in range(M)] for d in range(1,M+1)]
    for N in range(M,15):
        for r in range(M):
            for s in range(N):
                values=[kernel(M,r,d)*kernel(N,s,d) for d in range(1,M+1)]
                coeff=solve(matrix,values)
                require(all(x>=0 for x in coeff) and any(coeff),'independent cone interpolation')
                for d,value in enumerate(values,1):
                    require(sum(coeff[t]*kernel(M,t,d) for t in range(M))==value,'interpolation reconstruction')
                cone_count+=1
require(solve([[1,3],[1,0]],[0,1])==[Q(1),Q(-1,3)],'nonnegative kernel lies outside cone')
print('Independent carrier universe: degrees 2..7, prefix multisets from (-2,-1,1,3), tuned final factor.')
print('Carriers:',carrier_count,'strict coupled windows:',window_count)
print('Independent kernel interpolation universe: 1<=M<=7, M<=N<=14, all legal r,s.')
print('Cone rows:',cone_count)
print('Exact gates:',GATES)
print('semantic_sha256:',sha256(json.dumps(TRACE,separators=(',',':')).encode()).hexdigest())
print('PASS')
