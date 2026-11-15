#!/usr/bin/env python3
from itertools import combinations
from math import comb, factorial

def tournament(n, code):
    return {(i,j): bool(code >> t & 1)
            for t,(i,j) in enumerate(combinations(range(n),2))}

def reader(n, edges, k):
    m=k-1
    def forward(i,j):
        if j<n: return edges[i,j]
        if i<n: return j<n+m
        if i<n+m<=j: return i-n<=j-n-m
        return True
    return {(i,j):forward(i,j) for i,j in combinations(range(n+2*m),2)}

def adj(n,edges):
    out=[0]*n
    for (i,j),f in edges.items():
        if f: out[i]|=1<<j
        else: out[j]|=1<<i
    return out

def ends(n,edges):
    out=adj(n,edges)
    tab=[[0]*n for _ in range(1<<n)]
    for i in range(n): tab[1<<i][i]=1
    for s in range(1,1<<n):
        for i,h in enumerate(tab[s]):
            if not h: continue
            todo=out[i]&~s
            while todo:
                b=todo&-todo
                tab[s|b][b.bit_length()-1]+=h
                todo-=b
    return list(map(sum,tab))

def H(n,edges): return ends(n,edges)[-1]

def U(n,edges,k):
    h=ends(n,edges)
    f={(0,0):1}
    def g(s,j):
        if (s,j) not in f:
            z=0
            if s and j>0:
                b=s&-s
                r=s-b
                v=r
                while True:
                    z+=h[b+v]*g(r-v,j-1)
                    if not v: break
                    v=(v-1)&r
            f[s,j]=z
        return f[s,j]
    return g((1<<n)-1,k)

def Z(n,e,k): return factorial(k)*factorial(k-1)*U(n,e,k)

rows=[]
for k in range(1,5):
    row=[]
    for n in range(1,k+1):
        e=tournament(n,(1<<comb(n,2))-1)
        rem=H(n+2*(k-1),reader(n,e,k))
        rem-=sum(a*Z(n,e,j) for j,a in enumerate(row,1))
        z=Z(n,e,n)
        assert z and rem%z==0
        row.append(rem//z)
    assert row[-1]==1
    rows.append(row)
checks=0
for n in range(1,5):
    for code in range(1<<comb(n,2)):
        e=tournament(n,code)
        for k,row in enumerate(rows,1):
            assert H(n+2*(k-1),reader(n,e,k)) == sum(
                a*Z(n,e,j) for j,a in enumerate(row,1))
            checks+=1
print('Integral calibrated reader rows:',rows)
print('Independent subset-DP comparisons:',checks)
fixture=tournament(5,8)
assert H(5,fixture)==9
assert Z(5,fixture,2)==66
assert H(7,reader(5,fixture,2))==93
print('Three-adic decrement fixture: H=9 -> 93; Z2=66')
