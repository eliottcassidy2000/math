#!/usr/bin/env python3
"""Exact integral frame and characteristic-five degeneration; stdlib only."""
from itertools import permutations
import json
CHECKS=0

def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)

def mul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]

def det(a):
    n=len(a)
    total=0
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        term=sign
        for i in range(n): term*=a[i][p[i]]
        total+=term
    return total

def rank(a,p):
    a=[[x%p for x in row] for row in a]
    r=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r]
        inv=pow(a[r][j],-1,p)
        a[r]=[(x*inv)%p for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                t=a[i][j]
                a[i]=[(x-t*y)%p for x,y in zip(a[i],a[r])]
        r+=1
    return r

B=[[2*int(i==j)+int((i-j)%5 in (1,4)) for j in range(5)] for i in range(5)]
P=[[1,0,0,1],[-1,0,1,0],[0,-1,-1,0],[0,0,0,-1],[0,1,0,0]]
C=[[1,1,0,0],[1,2,0,0],[0,0,1,1],[0,0,1,2]]
need(all(sum(col)==0 for col in zip(*P)), 'augmentation columns')
need(det(P[:4])==-1, 'integral augmentation basis')
need(mul(B,P)==mul(P,C), 'integral intertwining')
T=[[1]+row for row in P]
need(abs(det(T))==5, 'ambient index')
fullC=[[4,0,0,0,0]]+[[0]+row for row in C]
need(mul(B,T)==mul(T,fullC), 'ambient rational intertwining')
D=[[B[i][j]+int(i==j) for j in range(5)] for i in range(5)]
D2=mul(D,D); D3=mul(D2,D)
ranks=[rank(D,5),rank(D2,5),rank(D3,5)]
need(ranks==[3,1,0], 'ambient nilpotent Jordan sizes 3,2')
need(all(x%5==1 for row in D2 for x in row), 'D squared is all-ones mod5')
E=[[C[i][j]+int(i==j) for j in range(4)] for i in range(4)]
need([rank(E,5),rank(mul(E,E),5)]==[2,0], 'augmentation Jordan sizes 2,2')
coords=[2,1,3,-1]
need([sum(x*y for x,y in zip(row,coords))%5 for row in P]==[1]*5,'constant lies in augmentation mod5')
need(rank(T,5)==4, 'ambient split frame singular mod5')
for p in (2,3,7,11,13): need(rank(T,p)==5, 'split frame invertible away from5')
print(json.dumps({'checks':CHECKS, 'augmentation_minor':det(P[:4]), 'ambient_determinant':det(T), 'ambient_D_ranks_mod5':ranks, 'augmentation_D_ranks_mod5':[2,0], 'constant_coordinates_mod5':coords, 'status':'PASS'},sort_keys=True))
