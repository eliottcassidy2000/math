#!/usr/bin/env python3
from itertools import combinations

def edge(i,j,back):
    return i != j and ((i<j and (j,i) not in back) or (i>j and (i,j) in back))

def triangles(n,back):
    return sum(edge(a,b,back) and edge(b,c,back) and edge(c,a,back) or
               edge(a,c,back) and edge(c,b,back) and edge(b,a,back)
               for a,b,c in combinations(range(n),3))

def strong(n,back):
    for i in range(n):
        seen={i}; todo=[i]
        while todo:
            x=todo.pop()
            for j in range(n):
                if j not in seen and edge(x,j,back): seen.add(j);todo.append(j)
        if len(seen)!=n:return False
    return True

def build(N):
    if N==0: return 1,set()
    back={(2,0)};n=3;x=1
    bits=bin(N)[3:]
    if bits:
        x=2+int(bits[0]); n=x+2; back={(n-1,0)}
    for b in map(int,bits[1:]):
        k=x+b-1; v=n+k; back.add((v,n-1)); n=v+1; x=2*x+b
    if x!=N: raise AssertionError((N,x))
    return n,back

def hamilton_count(n, back):
    full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if not dp[mask][i]: continue
            for j in range(n):
                if not (mask>>j)&1 and edge(i,j,back):
                    dp[mask|1<<j][j]+=dp[mask][i]
    return sum(dp[full])

def predicted_h(N):
    if not N:return 1
    x=1; h=3;bits=bin(N)[3:]
    if bits:
        x=2+int(bits[0]);h=1+2**x
    for b in map(int,bits[1:]):
        h*=1+2**(x+b-1);x=2*x+b
    return h

for N in range(31):
    n,b=build(N)
    if triangles(n,b)!=N or not strong(n,b): raise AssertionError((N,n))
    if N and n!=N+2: raise AssertionError(('size',N,n))
    if N<=7 and hamilton_count(n,b)!=predicted_h(N): raise AssertionError(('path_count',N))
    if N:
        xn=N+2; xb={(xn-1,0)}
        if triangles(xn,xb)!=N or not strong(xn,xb): raise AssertionError(('X',N))
        if {(xn-1-j,xn-1-i) for i,j in xb}!=xb: raise AssertionError(('reflect',N))
print('PASS adaptive and self-converse triangle statistics, N=0..30')
