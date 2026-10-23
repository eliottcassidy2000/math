#!/usr/bin/env python3
"""Exact killed-kernel and source-tail controls; source weights are explicit."""
from fractions import Fraction as F
from itertools import permutations
import json
from check import partitions, require


def matmul(A, B):
    return [[sum((a*b for a,b in zip(row,col)),F(0)) for col in zip(*B)] for row in A]


def psd(A):
    """Rational Schur elimination, including singular zero-pivot boundaries."""
    A=[list(row) for row in A]
    require(all(A[i][j]==A[j][i] for i in range(len(A)) for j in range(len(A))), 'symmetric PSD input')
    while A:
        pivot=A[0][0]
        if pivot<0:return False
        if pivot==0:
            if any(A[0][1:]):return False
            A=[row[1:] for row in A[1:]]
        else:
            A=[[A[i][j]-A[i][0]*A[0][j]/pivot for j in range(1,len(A))] for i in range(1,len(A))]
    return True


def solve_ones(A):
    """Gauss elimination for the independent expected-time calculation."""
    N=len(A)
    a=[list(row)+[F(1)] for row in A]
    for j in range(N):
        pivot=next(i for i in range(j,N) if a[i][j])
        a[j],a[pivot]=a[pivot],a[j]
        d=a[j][j];a[j]=[x/d for x in a[j]]
        for i in range(N):
            if i!=j:
                d=a[i][j];a[i]=[x-d*y for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]


def audit(Q, pi, gamma):
    N=len(Q); checks=0; subsets=0; tails=0; means=0
    require(sum(pi)==1 and all(sum(row)==1 for row in Q), 'stationary normalization')
    require(all(pi[i]*Q[i][j]==pi[j]*Q[j][i] for i in range(N) for j in range(N)), 'detailed balance')
    require(psd([[pi[i]*Q[i][j] for j in range(N)] for i in range(N)]), 'ambient PSD')
    require(psd([[pi[i]*(F(i==j)-Q[i][j])-gamma*(pi[i]*F(i==j)-pi[i]*pi[j])
                  for j in range(N)] for i in range(N)]), 'ambient gap certificate')
    checks+=4
    for mask in range(1,1<<N):
        S=[i for i in range(N) if mask>>i&1];C=[i for i in range(N) if not mask>>i&1]
        alpha=sum(pi[i] for i in S);rho=1-gamma*alpha
        B=[[Q[i][j] for j in C] for i in C]
        require(psd([[pi[i]*(rho*F(i==j)-Q[i][j]) for j in C] for i in C]), 'killed contraction')
        subsets+=1;checks+=1
        duration=solve_ones([[F(i==j)-Q[i][j] for j in C] for i in C])
        for src in S:
            extra=sum((Q[src][j]*duration[k] for k,j in enumerate(C)),F(0))
            a_squared=(1-alpha)/pi[src]
            require(extra*extra*(gamma*alpha)**2<=a_squared, 'source expectation bound')
            L=0;amplitude=a_squared
            while amplitude>1:
                L+=1;amplitude*=rho*rho
            require(extra<=L+1/(gamma*alpha), 'capped-tail expectation bound')
            means+=1;checks+=2
        # v=B^(t-1)1, formed by actual matrix iteration, not by the bound.
        v=[F(1)]*len(C)
        for t in range(1,13):
            for src in S:
                tail=sum((Q[src][j]*v[k] for k,j in enumerate(C)), F(0))
                require(F(0)<=tail<=1, 'probability range')
                require(tail*tail <= (1-alpha)/pi[src]*rho**(2*(t-1)), 'source-specific squared tail bound')
                tails+=1;checks+=2
            v=[sum((a*b for a,b in zip(row,v)), F(0)) for row in B]
    return dict(subsets=subsets, source_time_checks=tails,
                source_expectation_checks=means, total_checks=checks)


group=list(permutations(range(3)))
gg,go=partitions(3,group);tt,to=partitions(3,group,affine=True)
PG=[[F(gg[i]==gg[j],len(go[gg[i]])) for j in range(8)] for i in range(8)]
PT=[[F(tt[i]==tt[j],len(to[tt[i]])) for j in range(8)] for i in range(8)]
Q=matmul(matmul(PG,PT),PG)
cube=audit(Q,[F(1,8)]*8,F(8,9))
pi=[F(1,6),F(2,6),F(3,6)];stay=F(2,5)
nonuniform=[[stay*F(i==j)+(1-stay)*pi[j] for j in range(3)] for i in range(3)]
weighted=audit(nonuniform,pi,1-stay)
refresh=[[pi[j] for j in range(3)] for i in range(3)]
sharp=audit(refresh,pi,F(1))
for mask in range(1,7):
    C=[i for i in range(3) if not mask>>i&1]
    alpha=sum(pi[i] for i in range(3) if mask>>i&1)
    require(all(sum(refresh[i][j] for j in C)==1-alpha for i in C), 'sharp killed eigenvalue')
print(json.dumps(dict(status='FINITE-EXACT PASS',cube=cube,nonuniform=weighted,
                      refresh=sharp,sharp_subsets=6,maximum_time=12),indent=2,sort_keys=True))
