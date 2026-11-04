#!/usr/bin/env python3
"""Exact standard-library audit; no assertions removed by -O."""
import itertools,json

def require(b,msg):
    if not b: raise RuntimeError(msg)

def stream(w):
    P=L=0; h=1
    for c in w:
        h+=(1-c)*(P+L)
        P=2*P+c; L=c
    return h

def poly(w):
    n=len(w)
    return 1+sum((1<<(n-i-1))*w[i] for i in range(n-1))-sum((1<<(j-i-1))*w[i]*w[j] for i in range(n) for j in range(i+1,n))-sum(w[i]*w[i+1] for i in range(n-1))

def literal(w):
    n=len(w); N=n+1
    edge=lambda i,j: (i<j if i<n and j<n else (bool(w[j]) if i==n else not w[i]))
    dp=[[0]*N for _ in range(1<<N)]
    for i in range(N):dp[1<<i][i]=1
    for S in range(1,1<<N):
        for i in range(N):
            if not dp[S][i]:continue
            for j in range(N):
                if not S>>j&1 and edge(i,j):dp[S|1<<j][j]+=dp[S][i]
    return sum(dp[-1])

def coef(i,j):return -(1<<(j-i-1))-(j==i+1)

count=0
for n in range(9):
    for w in itertools.product((0,1),repeat=n):
        require(stream(w)==poly(w)==literal(w),('Hamilton',w));count+=1
# Every derivative is constant; independently recover on all small base assignments.
derivatives=0
for n in range(2,8):
    for i in range(n):
        for j in range(i+1,n):
            others=[r for r in range(n) if r not in (i,j)]
            for bits in itertools.product((0,1),repeat=len(others)):
                w=[0]*n
                for r,b in zip(others,bits):w[r]=b
                vals={}
                for a,b in itertools.product((0,1),repeat=2):
                    w[i]=a;w[j]=b;vals[a,b]=stream(w)
                require(vals[1,1]-vals[1,0]-vals[0,1]+vals[0,0]==coef(i,j),('derivative',n,i,j));derivatives+=1
cuts=0
for n in range(2,101):
    for k in range(1,n):
        for i in range(k):
            for j in range(k,n):
                require(coef(i,j)==-(1<<(k-i-1))*(1<<(j-k))-(i==k-1 and j==k),('factor',n,k,i,j))
        if k>=2 and n-k>=2:
            a,b=k-2,k-1;c,d=k,k+1
            require(coef(a,c)*coef(b,d)-coef(a,d)*coef(b,c)!=0,('rank',n,k))
        cuts+=1
prefixes=0
for k in range(13):
    signatures=set()
    for w in itertools.product((0,1),repeat=k):
        P=sum(b<<(k-i-1) for i,b in enumerate(w))
        sig=(stream(w+(0,0)),stream(w+(0,1)))
        require(sig[0]-sig[1]==2*P,('residual',w))
        signatures.add(sig);prefixes+=1
    require(len(signatures)==1<<k,('collision',k))
print(json.dumps(dict(status='PASS',literal_tournaments=count,derivative_base_assignments=derivatives,cut_factorizations=cuts,prefixes=prefixes,scope='finite audits supporting separately proved all-length statements'),indent=2))
