from itertools import combinations
from functools import cache
from math import isqrt
import sys

def prime(n):return n>1 and all(n%k for k in range(2,isqrt(n)+1))
def irreducible(n):
    return prime(n) or n%7==0 and prime(n//7) and n//7 not in (3,7)
def count_tables(a):
    n=len(a);F=(1<<n)-1
    f=[[0]*n for _ in range(1<<n)]
    H=[1]+[0]*F
    for S in range(1,1<<n):
        if S&(S-1)==0:
            H[S]=f[S][S.bit_length()-1]=1;continue
        for i in range(n):
            if S>>i&1:
                U=S^(1<<i)
                f[S][i]=sum(f[U][j] for j in range(n) if (U&a[i])>>j&1)
        H[S]=sum(f[S])
    return H
def coeff(a,w,S):
    # p[k] counts original-vertex runs of w; all other sites are singletons.
    n=len(a)
    # Place the singleton vertices and marked runs one by one, never two adjacent marked runs.
    out=[0]*((n+1)//2+1)
    def rec(L,last,r):
        if not L:
            if r:out[r]+=1
            if last!=w and (a[last]>>w&1):out[r+1]+=1
            return
        if last!=w and (last<0 or (a[last]>>w&1)):
            rec(L,w,r+1)
        choices=L if last<0 else L&a[last]
        while choices:
            b=choices&-choices;choices-=b
            rec(L^b,b.bit_length()-1,r)
    rec(S^(1<<w),-1,0)
    return out
def source(file,n):
    if n==7:
        for line in open(file):
            yield tuple(map(int,line.split()))
    else:
        with open(file) as f:
            next(f)
            for line in f:
                code=int(line.split(',')[0])
                yield tuple((code>>(8*i))&255 for i in range(n))
def evaluate(c,r):
    out=[0]*(len(c))
    out[1]=1
    for _ in range(r-1):out=[0]+[j*(out[j]+out[j-1]) for j in range(1,len(c))]
    return sum(v*w for v,w in zip(c,out))
n=int(sys.argv[1]);file=sys.argv[2]
N=J=P=0;good={}
for a in source(file,n):
    N+=1;F=(1<<n)-1;H=count_tables(a);h=H[F]
    for i,j in combinations(range(n),2):
        A,B=F^(1<<i),F^(1<<j)
        p,q=H[A],H[B]
        if min(p,q)<11 or p+q-1!=h or not (prime(p) and prime(q)):continue
        J+=1
        for w in range(n):
            if w in (i,j):continue
            s,x,y=coeff(a,w,F),coeff(a,w,A),coeff(a,w,B)
            if any(s[t]!=x[t]+y[t]-(t==1) for t in range(1,len(s))):continue
            P+=1
            if not any(s[2:]):continue
            for r in range(2,13):
                nh,np,nq=evaluate(s,r),evaluate(x,r),evaluate(y,r)
                if irreducible(nh) and prime(np) and prime(nq):
                    key=(a,i,j,w,tuple(s),tuple(x),tuple(y))
                    good.setdefault(key,[]).append((r,nh,np,nq))
print(n,N,J,P,'ports with further irreducible prime summand instances',len(good))
for k,t in list(good.items())[:8]:print(k,t)
