#!/usr/bin/env python3
"""Actual ordinal-profile bank, exact finite certificates, and witness constructor."""
from itertools import combinations
from math import gcd,isqrt
import argparse

def counts(rows):
    n=len(rows);N=1<<n
    end=[[0]*n for _ in range(N)]
    for i in range(n):end[1<<i][i]=1
    for s in range(1,N):
        for u,c in enumerate(end[s]):
            if c:
                z=rows[u]&~s
                while z:
                    b=z&-z
                    z-=b
                    end[s|b][b.bit_length()-1]+=c
    h=[1]+list(map(sum,end[1:]));F=N-1
    return h[F],sum(h[s]*h[F^s] for s in range(F+1))

def bank():
    rows_for={}
    for n in range(1,7):
        E=[(i,j) for j in range(n) for i in range(j) if j>i+1]
        for k in range(1<<len(E)):
            rows=[(2<<i) if i<n-1 else 0 for i in range(n)]
            for b,(i,j) in enumerate(E):
                if k>>b&1:rows[i]|=1<<j
                else:rows[j]|=1<<i
            rows_for.setdefault(counts(rows),tuple(rows))
    assert len(rows_for)==39
    return rows_for

def prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))

def det(a,b):
    return a[0]*b[1]-a[1]*b[0]

def logtable(p):
    for g in range(2,p):
        table={};x=1
        while x not in table:
            table[x]=len(table)
            x=x*g%p
        if len(table)==p-1:return table
    raise ValueError('p is not an odd prime')

def setup(p,R):
    units=sorted((d,q) for d,q in R if d*q%p)
    m=p-1; L=logtable(p)
    cols=[(m,0),(0,m)]+[(L[d%p],L[q%p]) for d,q in units]
    F=[((pow(d,m,p*p)-1)//p,(pow(q,m,p*p)-1)//p) for d,q in units]
    return units,L,cols,F

def egcd(a,b):
    if not b:return abs(a),1 if a>0 else -1,0
    g,x,y=egcd(b,a%b)
    return g,y,x-(a//b)*y

def exponent_vector(p,e,d,q,R):
    """Return (profile, rows, exponent) list whose ordinal join gives (d,q)."""
    mod=p**e;m=p-1
    if p%2==0 or not prime(p) or e<1 or gcd(d*q,p)!=1:
        raise ValueError('need odd prime p, positive exponent e, and unit target')
    units,L,C,F=setup(p,R)
    d%=mod;q%=mod
    target=(L[d%p],L[q%p])
    G=0; coeff=[0]*len(C)
    for i,j in combinations(range(len(C)),2):
        delta=det(C[i],C[j])
        if not delta:continue
        G,old,new=egcd(G,delta)
        coeff=[old*x for x in coeff]
        coeff[i]+=new*det(target,C[j])
        coeff[j]+=new*det(C[i],target)
        if G==1:break
    if G!=1:raise ValueError('bank is not certified over this prime')
    E=coeff[2:]
    A,B=1,1
    for (a,b),k in zip(units,E):
        A=A*pow(a,k,mod)%mod
        B=B*pow(b,k,mod)%mod
    if e>1:
        pair=next(((i,j) for i,j in combinations(range(len(F)),2) if det(F[i],F[j])%p),None)
        if pair is None:raise ValueError('bank has no two-direction lift over this prime')
        i,j=pair
        inv=pow(det(F[i],F[j]),-1,p)
        for level in range(1,e):
            scale=p**level; small=scale*p
            t=((d*pow(A,-1,small)%small-1)//scale,(q*pow(B,-1,small)%small-1)//scale)
            for index,digit in ((i,det(t,F[j])*inv%p),(j,det(F[i],t)*inv%p)):
                step=digit*m*p**(level-1)
                a,b=units[index]
                E[index]+=step
                A=A*pow(a,step,mod)%mod
                B=B*pow(b,step,mod)%mod
    assert (A,B)==(d,q)
    period=m*p**(e-1)
    out=[(profile,R[profile],k%period) for profile,k in zip(units,E) if k%period]
    if not out:
        profile=units[0]
        out=[(profile,R[profile],period)]
    return out

def vp(x,p):
    v=0
    while x%p==0:x//=p;v+=1
    return v

def closure_size(q,unit_profiles):
    reached={(1,1)};todo=[(1,1)]
    G=[(a%q,b%q) for a,b in unit_profiles]
    for d,k in todo:
        for a,b in G:
            z=(d*a%q,k*b%q)
            if z not in reached:reached.add(z);todo.append(z)
    return len(todo)

def joint_certificate(S,R):
    B=1
    for p in S:
        B*=p**(1+max((vp(r-1,p) for r in S if r!=p),default=0))
    U=[(a,b) for a,b in R if gcd(a*b,B)==1]
    phi=sum(gcd(x,B)==1 for x in range(B))
    assert closure_size(B,U)==phi*phi
    for p in S:
        F=[((pow(a,p-1,p*p)-1)//p,(pow(b,p-1,p*p)-1)//p) for a,b in U]
        assert any(det(a,b)%p for a,b in combinations(F,2))
    return B,phi*phi

def check(R):
    primes=[p for p in range(3,1000,2) if prime(p)]
    for p in primes:
        _,_,C,F=setup(p,R)
        assert gcd(*(det(a,b) for a,b in combinations(C,2)))==1
        assert any(det(a,b)%p for a,b in combinations(F,2))
    print(len(primes),'odd primes below 1000: full prime-power unit pairs')
    for q in range(3,256,2):
        U=[(d,k) for d,k in R if gcd(d*k,q)==1]
        phi=sum(gcd(a,q)==1 for a in range(q))
        assert closure_size(q,U)==phi*phi
    print('127 odd moduli 3 to 255: full unit-pair groups, checked directly')
    for p in primes:
        e=4;mod=p**e
        d=(1-p)%mod;q=(1+2*p+p*p)%mod
        witness=exponent_vector(p,e,d,q,R)
        x=y=1
        for (a,b),rows,k in witness:
            assert counts(rows)==(a,b)
            x=x*pow(a,k,mod)%mod;y=y*pow(b,k,mod)%mod
        assert (x,y)==(d,q)
    print('167 explicit constructor checks modulo p^4')
    for S in ((3,5),(3,7),(5,7),(3,5,7),(3,11),(5,11),(7,11),(3,5,11),(3,7,11)):
        B,size=joint_certificate(S,R)
        print('all joint powers over',S,'via base',B,'with',size,'unit pairs')

if __name__=='__main__':
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--target',nargs=4,type=int,metavar=('P','E','D','Q'))
    args=arg.parse_args()
    R=bank()
    if args.target:
        for profile,rows,k in exponent_vector(*args.target,R):
            print(k, 'copies of',rows,'with profile',profile)
    else:check(R)
