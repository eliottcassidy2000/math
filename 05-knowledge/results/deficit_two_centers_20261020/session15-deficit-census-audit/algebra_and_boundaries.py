#!/usr/bin/env python3
import itertools,json
from pathlib import Path
import sympy as sp
from reconstruct_spaces import affine,rows,require
ROOT=Path(__file__).resolve().parent
z=sp.symbols('z');P=1-z-z**3;Q=1-z**2-z**6;N=2*z**4-z**5-z**7
require(sp.cancel(sp.Rational(1,2)*(1/P-(1+z+z**3)/Q)+z**4/Q-N/(P*Q))==0,'leaf rational mismatch')
require(sp.gcd(N,P*Q)==1,'denominator not reduced')
require(sp.expand(P*Q)==1-z-z**2+z**5-z**6+z**7+z**9,'denominator expansion')
require(sp.rem(z**7+2*z**2-1,P)==0,'even cancellation polynomial')
A=[0]*501;L=[0]*501;M=[0]*501;A[0]=1
at=lambda seq,j:seq[j] if j>=0 else 0
for n in range(501):
    if n:A[n]=at(A,n-1)+at(A,n-3)
    M[n]=at(M,n-2)+at(M,n-6)+(1 if n==4 else 0)
    L[n]=at(L,n-1)+at(L,n-2)-at(L,n-5)+at(L,n-6)-at(L,n-7)-at(L,n-9)+(2 if n==4 else -1 if n in (5,7) else 0)
for n in range(501):
    r=n//2
    alt=(A[n]-A[r])//2+at(A,r-2) if n%2==0 else (A[n]-A[r]-at(A,r-1))//2
    require((A[n]-A[r]-(at(A,r-1) if n%2 else 0))%2==0,'nonintegral parity formula')
    require(L[n]==alt,'parity formula mismatch')
    require(M[n]==(at(A,r-2) if n%2==0 else 0),'mixed parity mismatch')
    if n%2==0 and r>=5:require(2*L[n]==A[n]-at(A,r-7),'even cancellation mismatch')
# Explicit q=2 boundary: the enumerated family really contains non-reflection frames.
q2=[]
for n,q in [(10,2),(14,2)]:
    pairs,eqs,basis,dirs,free=affine(n,q,None);anti=0;nonanti=0;witness=None
    for mask in range(1<<len(dirs)):
        x=basis
        for k,d in enumerate(dirs):
            if mask>>k&1:x^=d
        a=rows(n,pairs,x)
        isanti=all(((a[i]>>j)&1)!=((a[n-1-i]>>(n-1-j))&1) for i,j in pairs)
        if isanti:anti+=1
        else:
            nonanti+=1
            if witness is None:witness=dict(binary_assignment=mask,outneighbors=a)
    require(nonanti>0,'dangerous hidden reflection normalization')
    q2.append(dict(n=n,total=anti+nonanti,standard_reflection_anti=anti,standard_reflection_not_anti=nonanti,nonanti_control=witness))
# Positive atom: independently count all permutations on explicit strong S4.
a=[2,4,9,3];perms=list(itertools.permutations(range(4)))
arc=lambda u,v: bool(a[u]>>v&1)
autos=[p for p in perms if all(arc(u,v)==arc(p[u],p[v]) for u in range(4) for v in range(4) if u!=v)]
paths=[p for p in perms if all(arc(p[i],p[i+1]) for i in range(3))]
blue=[]
for p in paths:
    rev={p[i]:p[3-i] for i in range(4)}
    if all(arc(u,v)!=arc(rev[u],rev[v]) for u in range(4) for v in range(4) if u!=v):blue.append(p)
require((len(autos),len(paths),len(blue))==(1,5,3),'S4 control fails')
result=dict(algebra_checked_through=500,total_denominator_degree=sp.degree(P*Q),total_gcd=1,mixed_denominator_degree=sp.degree(Q),homogeneous_total_onset=8,homogeneous_mixed_onset=5,coefficients_0_through_30=L[:31],q2_boundaries=q2,S4=dict(automorphisms=len(autos),hamilton_paths=len(paths),blue_paths=len(blue),deficit=2))
(ROOT/'algebra_and_boundaries.json').write_text(json.dumps(result,indent=2,default=int)+'\n')
print(json.dumps(result,indent=2,default=int))
