#!/usr/bin/env python3
"""Independent local-square audit for the reflected skew-double embedding."""
from itertools import combinations,permutations
import json

def check(p,msg):
    if not p:raise RuntimeError(msg)

def double(A):
    n=len(A);B=[[0]*(2*n) for _ in range(2*n)]
    for i in range(n):
        for j in range(n):
            if i==j: B[i][n+j]=1
            else:
                B[i][j]=A[i][j];B[n+i][n+j]=1-A[i][j]
                B[i][n+j]=A[i][j];B[n+i][j]=A[i][j]
    return B

def energy(A,p):return sum(A[p[j]][p[i]] for i in range(len(p)) for j in range(i+1,len(p)))
def swap(p,i):p=list(p);p[i],p[i+1]=p[i+1],p[i];return tuple(p)

total=0;words=0;cores=0
for n in range(1,5):
    pairs=list(combinations(range(n),2))
    for m in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):A[i][j]=(m>>k)&1;A[j][i]=1-A[i][j]
        B=double(A)
        for p in permutations(range(n)):
            F=p+tuple(n+v for v in reversed(p))
            check(energy(B,F)==n*(n-1)//2+2*energy(A,p),'folded energy')
            # Pairwise coordinate equations characterize the affine diagonal section.
            positions={v:i+1 for i,v in enumerate(F)}
            check(all(positions[v]+positions[n+v]==2*n+1 for v in range(n)),'diagonal equations')
            for i in range(n-1):
                pp=swap(p,i);j=2*n-2-i
                q1=swap(F,i);q2=swap(F,j);q12=swap(q1,j);q21=swap(q2,i)
                check(q12==q21==pp+tuple(n+v for v in reversed(pp)),'commuting square')
                check(len({F,q1,q2,q12})==4,'four distinct corners')
                d=energy(A,pp)-energy(A,p)
                check(d in (-1,1),'base edge')
                check(energy(B,q1)-energy(B,F)==d,'first constituent')
                check(energy(B,q2)-energy(B,F)==d,'second constituent')
                check(energy(B,q12)-energy(B,q1)==d,'second after first')
                check(energy(B,q12)-energy(B,q2)==d,'first after second')
                total+=1
            words+=1
        cores+=1
print(json.dumps({'status':'PASS','scope':'all labeled tournament cores n1..4, all permutations, all adjacent interchanges','cores':cores,'words':words,'oriented_base_edges_tested':total,'method':'direct adjacency construction and all-pairs feedback energy; all four square corners tested'},indent=2))
