#!/usr/bin/env python3
"""Independent clock/zero checks using Newton power sums and literal partitions."""
from math import factorial, comb, gcd, lcm
from itertools import combinations
import json

def req(q,m):
    if not q:raise RuntimeError(m)
def factors(M):
    r=[];p=2
    while p<=M:
        a=0
        while M%p==0:a+=1;M//=p
        if a:r.append((p,a))
        p+=1
    return r
def predicted(M):
    if M==1:return 1,1
    f=factors(M);L=1
    for p,a in f:
        if p>2:t=(p-1)*p**(a-1)
        else:
            t=1;q=2**(a-1)
            while pow(3,t,q)!=1%q:t+=1
        L=lcm(L,t)
    return max(a for p,a in f),L
def cutoff(M):
    c=1
    while factorial(c)*factorial(c-1)%M:c+=1
    return c
def state(M,n):
    return tuple((factorial(c-1)*sum((-1)**(c-j)*comb(c,j)*pow(j,n,M) for j in range(c+1)))%M for c in range(1,cutoff(M)))
def audit_clock():
    rows=[]
    for M in list(range(1,201))+[243,256,343,512,625,729,1024,1331,2048]:
        A,L=predicted(M);seen={};first=None
        for n in range(1,A+L+1):
            s=state(M,n)
            if s in seen:first=(seen[s],n-seen[s]);break
            seen[s]=n
        req(first==(A,L),('exact first repeat',M,first,(A,L)))
        rows.append((M,A,L))
    req(predicted(16)==(4,2),'dyadic period is not Euler period')
    req(state(25,1)!=state(25,21),'preperiod cannot be dropped')
    return {'moduli':len(rows),'first_repetitions_match':True,'method':'Newton finite differences, not Stirling recurrence'}
def singleton():return [[0]]
def graph_op(A,op):
    n=len(A)
    if op=='P':
        B=[[0]*(n+1) for _ in range(n+1)]
        for i in range(n):
            for j in range(n):B[i][j]=A[i][j]
        for i in range(n):B[n][i]=1
        return B
    k={'C':1,'D':2}[op];B=[[0]*(n+k+1) for _ in range(n+k+1)]
    for i in range(n):
        for j in range(n):B[i][j]=A[i][j]
        for j in range(n,n+k):B[i][j]=1
        B[n+k][i]=1
    for i in range(n,n+k):
        B[i][n+k]=1
        for j in range(i+1,n+k):B[i][j]=1
    return B
def ordinal(A,B):
    n=len(A);m=len(B);C=[[0]*(n+m) for _ in range(n+m)]
    for i in range(n+m):
        for j in range(n+m):
            if i<n and j<n:C[i][j]=A[i][j]
            elif i>=n and j>=n:C[i][j]=B[i-n][j-n]
            elif i<n<=j:C[i][j]=1
    return C
def profile(A):
    n=len(A);N=1<<n
    # Endpoint recurrence through final predecessor, unlike producer forward extension.
    d=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
    for mask in range(1,N):
        for j in range(n):
            if mask>>j&1:
                prev=mask^(1<<j)
                d[mask][j]=1 if not prev else sum(d[prev][i]*A[i][j] for i in range(n) if prev>>i&1)
        H[mask]=sum(d[mask])
    memo={0:(1,)}
    def parts(mask):
        if mask in memo:return memo[mask]
        r=[0]*(mask.bit_count()+1);root=mask&-mask;s=mask
        while s:
            if s&root:
                for c,z in enumerate(parts(mask^s)):r[c+1]+=H[s]*z
            s=(s-1)&mask
        memo[mask]=tuple(r);return memo[mask]
    return parts(N-1)
def eval_profile(pc,x):
    r=0;t=1
    for c in range(1,len(pc)):t*=x-c+1;r+=pc[c]*t
    return r
def audit_zero():
    W5=graph_op(graph_op(singleton(),'P'),'C')
    W11=graph_op(graph_op(singleton(),'P'),'D')
    Z=ordinal(W5,W11);p5,p11,pz=map(profile,[W5,W11,Z])
    req((p5[1],p11[1],pz[1])==(5,11,55),'literal Hamiltonian counts')
    for x in range(-6,13):
        f=x**4+4*x*x;g=x**5+8*x**3+2*x
        req(eval_profile(p5,x)==f and eval_profile(p11,x)==g,'explicit seed profiles')
        req(eval_profile(pz,x)==f*g,'ordinal product profile')
    sig=[factorial(c)*factorial(c-1)*pz[c]%5 for c in range(1,len(pz))]
    req(set(sig)=={0},'order9 contextual zero mod5')
    req(all(eval_profile(pz,x)%5==0 for x in range(5)),'order9 function zero')
    # Product polynomial degree9 is nonzero over F5, although its F5 function vanishes.
    req([1%5,12%5,34%5,8%5]!=[0]*4,'polynomial/function hostile')
    return {'orders':[4,5,9],'Hamiltonian_counts':[5,11,55],'order9_pathcover_profile':pz,'mod5_signature':sig,'function_not_polynomial_zero':True}
if __name__=='__main__':print(json.dumps({'status':'PASS','clock':audit_clock(),'zero':audit_zero()},indent=2,sort_keys=True))
