#!/usr/bin/env python3
"""Independent exact controls for a seven-vertex MRB/star-filter example."""
from itertools import permutations, combinations
from math import comb
N=7

def adjacency(mask):
    A=[[0]*N for _ in range(N)]
    for i in range(5):
        for j in range(i+1,5): A[i][j]=1
    for i in range(5): A[i][5]=1; A[6][i]=1
    A[5][6]=1
    for k,i in enumerate((2,3,4)):
        if mask>>k&1: A[0][i],A[i][0]=A[i][0],A[0][i]
    return A

def dp(A):
    D=[[0]*N for _ in range(1<<N)]
    for i in range(N): D[1<<i][i]=1
    for m in range(1,1<<N):
        for i in range(N):
            c=D[m][i]
            if c:
                for j in range(N):
                    if not (m>>j&1) and A[i][j]: D[m|1<<j][j]+=c
    return sum(D[-1])

def brute(A):
    return sum(all(A[a][b] for a,b in zip(p,p[1:])) for p in permutations(range(N)))

V=[dp(adjacency(m)) for m in range(8)]
W=[brute(adjacency(m)) for m in range(8)]
if V!=W: raise RuntimeError((V,W))
if V!=[33,51,69,69,93,93,111,93]: raise RuntimeError(V)
A=[V[1<<i]-V[0] for i in range(3)]
B=[V[(1<<i)|(1<<j)]-V[1<<i]-V[1<<j]+V[0] for i,j in combinations(range(3),2)]
if A!=[18,36,60] or B!=[-18]*3: raise RuntimeError((A,B))
C=[sum(V[m] for m in range(8) if m.bit_count()==j) for j in range(4)]
if C!=[33,213,273,93]: raise RuntimeError(C)
if sum((-1)**m.bit_count()*V[m] for m in range(8))!=0: raise RuntimeError('degree')
print('subset DP = independent permutation enumeration:',V)
print('singleton coefficients:',A,'pair coefficients:',B)
print('subset-size layer sums:',C)
print('MRB-filtered projection = 33 + 273*cuberoot(3) - 306*sqrt(2)')

def general_dp(n,mask):
    N=n+2
    A=[[0]*N for _ in range(N)]
    for i in range(n):
        for j in range(i+1,n): A[i][j]=1
    for i in range(n): A[i][n]=1;A[n+1][i]=1
    A[n][n+1]=1
    for i in range(1,n):
        if mask>>(i-1)&1:A[0][i],A[i][0]=A[i][0],A[0][i]
    D=[[0]*N for _ in range(1<<N)]
    for i in range(N):D[1<<i][i]=1
    for m in range(1,1<<N):
        for i in range(N):
            for j in range(N):
                if D[m][i] and not (m>>j&1) and A[i][j]:D[m|1<<j][j]+=D[m][i]
    return sum(D[-1])

def formula(n,m):
    def q(d):
        c=1 if d==1 else 2**(d-2)
        return c+2**(n-d-2)*(3**(d-1)+1)
    return 1+2**n+2*sum(q(j-i) for i in range(1,n) for j in range(i+1,n) if not (m>>(i-1)&1) and (m>>(j-1)&1))

k=0
for n in range(2,8):
    for mask in range(1<<(n-1)):
        x=general_dp(n,mask);y=formula(n,mask)
        if x!=y:raise RuntimeError((n,mask,x,y))
        k+=1
print('strong wrapped weighted-cut formula: exact DP checks',k,'over n=2..7')
rows=['0010000','1010011','0000100','1110100','1100000','1011101','1011100']
T=[[int(c) for c in row] for row in rows]
if any(T[i][j]+T[j][i]!=1 for i in range(7) for j in range(i+1,7)):raise RuntimeError('not a tournament')
reach=[[bool(T[i][j] or i==j) for j in range(7)] for i in range(7)]
for k in range(7):
    for i in range(7):
        for j in range(7):reach[i][j]|=reach[i][k] and reach[k][j]
if not all(all(row) for row in reach):raise RuntimeError('not strong')
if dp(T)!=35 or brute(T)!=35:raise RuntimeError('not 35')
for n in range(2,6):
    if 35 in {formula(n,m) for m in range(1<<(n-1))}:raise RuntimeError('35 in star-cut family')
print('35: explicit strong 7-tournament by two enumerations, absent from star-cut family')
