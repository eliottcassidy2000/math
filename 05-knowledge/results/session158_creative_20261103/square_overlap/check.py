#!/usr/bin/env python3
from itertools import permutations
from math import isqrt

def formula(A,B,C):
    return (1+2**(A+B-1)+2**(B+C-1)
       - (2**(A-1) if C==1 else 0)
       - (2**(C-1) if A==1 else 0)
       + 2**(A+C-1)*3**(B-1))

def tournament(A,B,C):
    a=0;c=A;b=A+B;d=A+B+C
    E={(a,b),(c,d)}; n=d+1
    return [[i!=j and ((i<j)^((min(i,j),max(i,j)) in E)) for j in range(n)] for i in range(n)]

def direct_by_pattern(T,A,B,C):
    n=len(T); a,c,b,d=0,A,A+B,A+B+C; sums=[0,0,0,0]
    for p in permutations(range(n)):
        edges=list(zip(p,p[1:]))
        if all(T[x][y] for x,y in edges):
            kind=((b,a) in edges)+2*((d,c) in edges)
            sums[kind]+=1
    return sums

def dp(T):
    n=len(T); table=[[0]*n for i in range(1<<n)]
    for i in range(n): table[1<<i][i]=1
    for S in range(1,1<<n):
        for x in range(n):
            if table[S][x]:
                for y in range(n):
                    if not S>>y&1 and T[x][y]:table[S|1<<y][y]+=table[S][x]
    return sum(table[-1])

def triangle_after_deletion(T,v):
    V=[x for x in range(len(T)) if x!=v]
    for a in V:
        for b in V:
            for c in V:
                if a!=b and b!=c and c!=a and T[a][b] and T[b][c] and T[c][a]:return tuple(x+1 for x in (a,b,c))
    return None

def check():
    ncase=0
    for A in range(1,7):
      for B in range(1,7):
       for C in range(1,7):
        if A+B+C+1>11:continue
        T=tournament(A,B,C)
        want=formula(A,B,C)
        got=dp(T)
        if got!=want:raise ValueError((A,B,C,got,want))
        if len(T)<=8:
            p=direct_by_pattern(T,A,B,C)
            expected=[1,2**(A+B-1)-(2**(A-1) if C==1 else 0),2**(B+C-1)-(2**(C-1) if A==1 else 0),2**(A+C-1)*3**(B-1)]
            if p!=expected:raise ValueError((A,B,C,p,expected))
        ncase+=1
    print('PASS',ncase,'gap triples by endpoint-subset DP, direct permutations by edge pattern for n<=8')
    T=tournament(1,2,1)
    print('H(13 witness)',dp(T),'deletion triangles',[triangle_after_deletion(T,v) for v in range(5)])
    if dp(T)!=13 or not all(triangle_after_deletion(T,v) for v in range(5)):raise ValueError('13 witness')
    T=tournament(5,2,2)
    print('square-sum pair example H',dp(T),'pairs (1,8), (6,10) sums 9,16')
    if dp(T)!=265:raise ValueError('265 witness')
check()
