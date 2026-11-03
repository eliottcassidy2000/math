#!/usr/bin/env python3
"""Independent finite checks of the terminal-carrier/3-cycle identity."""
from itertools import combinations, permutations, product
from math import comb, factorial

def check(condition):
    if not condition: raise AssertionError("check failed")

def tournament(n, mask):
    E=[[False]*n for _ in range(n)]
    for bit,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(j,i) if mask>>bit&1 else (i,j)
        E[u][v]=True
    return E

def profile(E):
    n=len(E); c=[0]*(n+2)
    for p in permutations(range(n)):
        b=sum(not E[p[i]][p[i+1]] for i in range(n-1))
        for k in range(1,n+1):
            # Forced cuts at each backward adjacency; all other cuts free.
            if 0 <= k-1-b <= n-1-b:
                c[k]+=comb(n-1-b,k-1-b)
    return c

def terminal_weight(E, subset):
    vertices=[v for v in range(len(E)) if subset>>v&1]
    if len(vertices)==1: return 1
    total=0
    for a,b in combinations(vertices,2):
        interior=[v for v in vertices if v not in (a,b)]
        for middle in permutations(interior):
            for start,end in ((a,b),(b,a)):
                p=(start,)+middle+(end,)
                total+=all(E[p[i]][p[i+1]] for i in range(len(p)-1))
    return total

def carrier_profile(E):
    # Partition recursion anchors the least remaining label. The block weight
    # explicitly sums its endpoint pairs, endpoint orientations, interior orders.
    n=len(E); full=(1<<n)-1
    w={s:terminal_weight(E,s) for s in range(1,full+1)}
    def rec(s,k):
        if not s: return int(k==0)
        if k==0: return 0
        anchor=s&-s; ans=0; block=s
        while block:
            if block&anchor: ans+=w[block]*rec(s^block,k-1)
            block=(block-1)&s
        return ans
    return [0]+[rec(full,k) for k in range(1,n+1)]

def hamilton_dp(E):
    n=len(E); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            if not dp[s][i]: continue
            for j in range(n):
                if not s>>j&1 and E[i][j]: dp[s|1<<j][j]+=dp[s][i]
    return sum(dp[-1])

def cycle(A,B,C):
    blocks=(A,B,C); n=sum(map(len,blocks)); E=[[False]*n for _ in range(n)]
    off=[]; q=0
    for x in blocks: off.append(q); q+=len(x)
    for z,x in enumerate(blocks):
        for i in range(len(x)):
            for j in range(len(x)): E[off[z]+i][off[z]+j]=x[i][j]
    for z in range(3):
        w=(z+1)%3
        for i in range(len(blocks[z])):
            for j in range(len(blocks[w])): E[off[z]+i][off[w]+j]=True
    return E

def formula(a,b,c):
    get=lambda v,k: v[k] if k<len(v) else 0
    s=0
    for k in range(1,max(map(len,(a,b,c)))+1):
        a0,b0,c0=[get(v,k) for v in (a,b,c)]
        a1,b1,c1=[get(v,k+1) for v in (a,b,c)]
        s+=3*a0*b0*c0+a1*b0*c0+a0*b1*c0+a0*b0*c1+a1*b1*c0+a1*b0*c1+a0*b1*c1
    return s

def carriers(n):
    # Enumerate labeled unordered partitions with each non-singleton block
    # carrying an unordered pair of marked terminals.
    if n==0:return 1
    return carriers(n-1)+sum(comb(n-1,s-1)*comb(s,2)*carriers(n-s) for s in range(2,n+1))

def main():
    check([carriers(n) for n in range(3,6)]==[7,28,126])
    for n in range(1,6):
        for mask in range(1<<(n*(n-1)//2)):
            E=tournament(n,mask); p=profile(E)
            check(p[1]==hamilton_dp(E))
            # C_2 consists of (n-1) splits of each path, plus permutations
            # with exactly one backward adjacent edge.
            if n>=2: check(p[2]>=(n-1)*p[1])
    for n in range(1,5):
        for mask in range(1<<(n*(n-1)//2)):
            E=tournament(n,mask); p=profile(E); cp=carrier_profile(E)
            for k in range(1,n+1): check(cp[k]*factorial(k)==p[k])
    for (n,m,l) in ((1,1,1),(1,2,3),(1,3,3),(2,2,3)):
        choices=[range(1<<(s*(s-1)//2)) for s in (n,m,l)]
        for masks in product(*choices):
            A,B,C=(tournament(s,mask) for s,mask in zip((n,m,l),masks))
            check(hamilton_dp(cycle(A,B,C))==formula(profile(A),profile(B),profile(C)))
    singleton=tournament(1,0)
    for mask, expected_c2, expected_h in [(40,78,123),(76,90,135)]:
        B=tournament(5,mask); p=profile(B)
        check(p[1:3]==[15,expected_c2])
        check(hamilton_dp(cycle(singleton,B,singleton))==expected_h)
    print('PASS: carrier 7,28,126 and weighted profiles through n=4; all tournaments through n=5; 3-cycle DP; flexibility witnesses 123,135')
if __name__=='__main__': main()
