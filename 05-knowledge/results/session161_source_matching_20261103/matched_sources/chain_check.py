#!/usr/bin/env python3
"""Literal independent controls for two matched-source edge reversals."""
from itertools import combinations,permutations

def need(v,msg):
    if not v: raise RuntimeError(msg)

def graph(n,back):
    a=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            p,q=(j,i) if (i,j) in back else (i,j)
            a[p]|=1<<q
    return a

def pad(a,k):
    n=len(a); b=[0]*(k+n)
    for i in range(k):
        for j in range(i+1,k+n): b[i]|=1<<j
    for i in range(n): b[k+i]=a[i]<<k
    return b

def wrap(a):
    n=len(a); b=a[:]+[0,0]
    for i in range(n): b[i]|=1<<n; b[n+1]|=1<<i
    b[n]|=1<<(n+1)
    return b

def path_subset(a):
    n=len(a); d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,x in enumerate(d[mask]):
            if not x: continue
            possible=a[i]&~mask
            while possible:
                bit=possible&-possible; j=bit.bit_length()-1
                d[mask|bit][j]+=x; possible-=bit
    return [1]+[sum(row) for row in d[1:]]

def slots(a,r):
    h=path_subset(a);n=len(a);v=[0]*(1<<n);v[0]=1
    for unused in range(r):
        z=[0]*(1<<n)
        for mask in range(1<<n):
            sub=mask
            while True:
                z[mask]+=h[sub]*v[mask^sub]
                if not sub:break
                sub=(sub-1)&mask
        v=z
    return v[-1]

def cycles(a):
    n=len(a);r=[]
    for k in range(3,n+1,2):
        for t in combinations(range(n),k):
            for tail in permutations(t[1:]):
                c=(t[0],)+tail
                if all((a[i]>>j)&1 for i,j in zip(c,c[1:]+c[:1])):r.append(c)
    return r

def strong(a):
    n=len(a)
    for root in range(n):
        seen=1<<root
        while True:
            new=seen
            for i in range(n):
                if seen>>i&1:new|=a[i]
            if new==seen:break
            seen=new
        if seen!=(1<<n)-1:return False
    return True

def edge_differences(a,b):
    return [(i,j) for i in range(len(a)) for j in range(i+1,len(a)) if bool(a[i]>>j&1)!=bool(b[i]>>j&1)]

back=[set(),{(0,2)},{(0,2),(0,3)}]
base=[graph(5,x) for x in back]
for j,a in enumerate(base):
    C=cycles(a)
    need(len(C)==j and all(len(c)==3 for c in C),'triangle enumeration')
    for r in range(6):need(slots(a,r)==r**5+2*j*r**3,'slot profile')
for k in (0,3):
    B=[wrap(pad(a,k)) for a in base]
    for j,a in enumerate(B):
        need(strong(a),'wrapper strongness')
        h=path_subset(a)[-1]
        need(h==1+32*2**k+j*(2+16*2**k),'wrapped count')
    need(edge_differences(B[0],B[1])==[(k,k+2)],'first edge')
    need(edge_differences(B[1],B[2])==[(k,k+3)],'second edge')
print('PASS: cycles, profiles, literal strong chains, counts, and named changed edges')
