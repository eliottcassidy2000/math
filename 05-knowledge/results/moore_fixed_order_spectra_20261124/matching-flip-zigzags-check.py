from itertools import permutations
from math import factorial, comb
from fractions import Fraction

def paths(A, start=None):
    n=len(A)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        if start is None or i==start:
            dp[1<<i][i]=1
    for m in range(1,1<<n):
        for i,c in enumerate(dp[m]):
            if c:
                for j,v in enumerate(A[i]):
                    if v and not m>>j&1:
                        dp[m|1<<j][j]+=c
    return dp[-1]

def tournament(n):
    return [[i<j for j in range(n)] for i in range(n)]

def flip(A,i,j):
    A[i][j],A[j][i]=A[j][i],A[i][j]

def cube(k,n):
    p=[]
    q=[]
    for s in range(1<<k):
        A=tournament(n)
        flip(A,0,n-1)
        for i in range(k):
            if s>>i&1:
                flip(A,i,i+k)
        p.append(sum(paths(A)))
        q.append(sum(x for i,x in enumerate(paths(A,0)) if A[i][0]))
    return p,q

def mixed(p,k):
    return sum((-1)**(k-s.bit_count())*v for s,v in enumerate(p))

def contracted(k):
    n=2*k+1
    A=tournament(n)
    flip(A,0,n-1)
    dp=[[0]*(n-1) for _ in range(1<<k)]
    dp[1][k]=-1
    for m in range(1,1<<k):
        for i,v in enumerate(dp[m]):
            if v:
                for j in range(k):
                    if not m>>j&1:
                        if A[i][j]:
                            dp[m|1<<j][j+k]-=v
                        if A[i][j+k]:
                            dp[m|1<<j][j]+=v
    return sum(v for i,v in enumerate(dp[-1]) if A[i][n-1])

def gadgets(n,mask):
    A=tournament(n)
    pairs=[(i,j) for i in range(n) for j in range(i+2,n)]
    for b,(i,j) in enumerate(pairs):
        if mask>>b&1:
            flip(A,i,j)
    p=paths(A,0)[-1]
    for i in range(n-1):
        flip(A,i,i+1)
    return p,paths(A,n-1)[0]

def bounds():
    assert gadgets(6,163)==(3,3) and gadgets(8,1172854)==(24,23)
    assert 43206+157<=comb(295,2)
    assert comb(7,2)*7+comb(5,2)==157
    assert 3*24**7>2**33 and 3*23**7>2**33
    n=296
    scale=10**7
    lo,hi=27182818,27182819
    partial=sum(Fraction(1,factorial(i)) for i in range(15))
    assert Fraction(lo,scale)<partial
    assert partial+Fraction(16,15*factorial(15))<Fraction(hi,scale)
    assert all(factorial(d)*hi**d<=((d+5)*scale)**d for d in range(1,n+1))
    a=(n*(n+3)+10*(n+1))*scale
    b=2*(n+1)*lo
    assert a**(n+1)<(b**(n+1))<<1730

if __name__=='__main__':
    p,q=cube(3,7)
    assert p==[33,33,69,65,69,71,117,113]
    assert q==[1,1,5,4,5,5,13,11]
    for k in range(1,5):
        p,q=cube(k,2*k+1)
        assert mixed(p,k)==2*mixed(q,k)==2*contracted(k)
        if k>=2:
            p,q=cube(k,2*k)
            assert mixed(p,k)==0
    for m in range(9):
        e=sum(all((a>b)==(i%2==0) for i,(a,b) in enumerate(zip(p,p[1:])))
              for p in permutations(range(m)))
        assert -contracted(m+1)==e
    bounds()
    print('exact flip cubes, Euler zigzags, and marked-code bounds verified')
