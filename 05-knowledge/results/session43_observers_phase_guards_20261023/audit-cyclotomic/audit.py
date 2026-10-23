#!/usr/bin/env python3
"""Independent audit: direct permutations, cyclotomic norms, tangent ODE."""
from itertools import combinations,permutations
from math import factorial,comb
import json,hashlib, pathlib,subprocess

def check(x,msg):
    if not x:raise RuntimeError(msg)
def vprime(n,p):
    if not n:return None
    n=abs(n);v=0
    while n%p==0:n//=p;v+=1
    return v
def ds(n,p):
    return sum_digits(n,p)
def sum_digits(n,p):
    a=[]
    while n:a.append(n%p);n//=p
    return sum(a)
def bound(n):return min(n-vprime(n+1,2), n+1-((n-1).bit_length()-1))
def matrix(n,code):
    A=[[False]*n for i in range(n)]
    for j,(a,b) in enumerate(combinations(range(n),2)):
        A[a][b]=bool(code&(1<<j));A[b][a]=not A[a][b]
    return A

def direct_f(A):
    n=len(A);out=[0]*n
    for seq in permutations(range(n)):
        k=sum(A[a][b] for a,b in zip(seq,seq[1:]));out[k]+=1
    return out

def det(M):
    A=[row[:] for row in M];n=len(A);old=1;sign=1
    if n==1:return A[0][0]
    for k in range(n-1):
        if not A[k][k]:
            pivot=next((j for j in range(k+1,n) if A[j][k]),None)
            if pivot is None:return 0
            A[k],A[pivot]=A[pivot],A[k];sign=-sign
        t=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                val=A[i][j]*t-A[i][k]*A[k][j]
                check(val%old==0,'Bareiss exact division');A[i][j]=val//old
        for i in range(k+1,n):A[i][k]=0
        old=t
    return sign*A[-1][-1]
def norm(f,p):
    # Multiplication by F(zeta) on integral power basis; determinant is field norm.
    M=[[0]*(p-1) for _ in range(p-1)]
    for col in range(p-1):
        totals=[0]*p
        for i,c in enumerate(f):totals[(i+col)%p]+=c
        for row in range(p-1):M[row][col]=totals[row]-totals[p-1]
    return det(M)

def odd_cycles(A):
    n=len(A);out=[]
    for k in range(3,n+1,2):
        for V in combinations(range(n),k):
            for tail in permutations(V[1:]):
                seq=V[:1]+tail
                if all(A[seq[i]][seq[(i+1)%k]] for i in range(k)):
                    out.append(sum(1<<v for v in V))
    return out

def eulerian(n):
    a=[1]
    for m in range(2,n+1):
        b=[0]*m
        for j,c in enumerate(a):b[j]+=(j+1)*c;b[j+1]+=(m-j-1)*c
        a=b
    return a

def expansion(A):
    n=len(A);cs=odd_cycles(A);out=[0]*n
    def walk(start,mask,r):
        d=n-mask.bit_count()+r;gap=n-d
        for j,c in enumerate(eulerian(d)):
            for k in range(gap+1):out[j+k]+=2**r*c*comb(gap,k)*(-1)**(gap-k)
        for i in range(start,len(cs)):
            if not (cs[i]&mask):walk(i+1,mask|cs[i],r+1)
    walk(0,0,0)
    return out

counts=[];total=0;exp_cases=0;norm_cases=0
for n in range(1,6):
    low=None;vals=[]
    for code in range(1<<(n*(n-1)//2)):
        A=matrix(n,code);f=direct_f(A);check(f==expansion(A),('cycle identity',n,code));exp_cases+=1
        for p in (2,3,5,7):
            z=norm(f,p);vp=vprime(z,p);norm_cases+=1
            E=n-max(ds(m,p) for m in range(1,n+1))
            if vp is not None:
                check(vp>=E,('digit bound',n,code,p))
                if p!=2:check(vp%2==0 and vp>=2*((E+1)//2),'odd prime even bound')
            if p==2:
                if n%2==0:check(z==0,'even antipodal zero')
                elif vp is not None:
                    low=vp if low is None else min(low,vp)
                    if n>=3:check(vp>=bound(n),'binary improved bound')
        total+=1
    if n>=3 and n%2:check(low==bound(n),'binary sharp exhaustive')
    counts.append({'n':n,'tournaments':1<<(n*(n-1)//2),'minimum_binary_valuation':low})
# Deterministic independent actual larger tournaments, including extremizers.
for n in (6,7,8,9):
    a=(n-1).bit_length()-1
    q=n-(1<<a)+2 if n%2 and n-vprime(n+1,2)>n+1-a else 0
    N=n*(n-1)//2;codes={0,(1<<N)-1}
    if q:codes.add((1<<N)-1 ^ (1<<list(combinations(range(n),2)).index((0,q-1))))
    for code in sorted(codes):
        A=matrix(n,code);f=direct_f(A)
        if n<=7:check(f==expansion(A),'larger exact identity');exp_cases+=1
        val=sum((-1)**j*c for j,c in enumerate(f))
        if n%2==0:check(val==0,'larger even zero')
        else:
            witness=((1<<N)-1 ^ (1<<list(combinations(range(n),2)).index((0,q-1)))) if q else (1<<N)-1
            check(vprime(val,2)>=bound(n),'larger bound')
            if code==witness:check(vprime(val,2)==bound(n),'actual sharp witness')
        for p in (2,3,5,7,11):
            vp=vprime(norm(f,p),p);norm_cases+=1
            E=n-max(ds(m,p) for m in range(1,n+1))
            if vp is not None:check(vp>=E and (p==2 or vp%2==0),'larger norm bound')
        total+=1
# Independent tangent construction from differential equation tan'=1+tan^2.
T=[0]*402;T[1]=1
for n in range(1,401):T[n+1]=sum(comb(n,k)*T[k]*T[n-k] for k in range(n+1))
sharp=[]
for n in range(3,402,2):
    check(vprime(T[n],2)==n-vprime(n+1,2),'tangent valuation')
    empty=n-vprime(n+1,2);a=(n-1).bit_length()-1
    if empty<=n+1-a:q=0;value=(-1)**((n-1)//2)*T[n]
    else:
        q=n-(1<<a)+2
        value=(-1)**((n-1)//2)*T[n]
        for ell in range(3,q+1,2):
            d=n-ell+1;value+=comb(q-2,ell-2)*2**ell*(-1)**((d-1)//2)*T[d]
    check(vprime(value,2)==bound(n),('allsize scalar sharp',n,q))
    z=n+1
    repunit=(z & (z-1))==0
    double_triple=(z%3==0 and (z//3 & (z//3-1))==0 and z//3>=2)
    check((q==0)==(repunit or double_triple),'transitive classification')
    sharp.append({'n':n,'q':q,'valuation':vprime(value,2)})
print(json.dumps({'status':'PASS','direct_permutation_tournaments':total,'cycle_expansion_checks':exp_cases,'independent_norm_checks':norm_cases,'exhaustive_orders':counts,'tangent_ODE_sharp_checks':sharp},sort_keys=True,indent=2))
