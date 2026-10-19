#!/usr/bin/env python3
"""Independent canonical-tournament and integral divided-power review."""
from itertools import permutations, combinations
from functools import lru_cache
from collections import defaultdict
from math import comb
import json
from sympy import Matrix
counts=defaultdict(int)
def check(p,kind,detail=None):
    counts[kind]+=1
    if not p:raise RuntimeError((kind,detail))
@lru_cache(None)
def pairs(n):return tuple(combinations(range(n),2))
def arcs(n,code):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs(n)):
        a[i][j]=(code>>k)&1;a[j][i]=1-a[i][j]
    return a
def encode(a):return sum(a[i][j]<<k for k,(i,j) in enumerate(pairs(len(a))))
@lru_cache(None)
def canon(n,code):
    a=arcs(n,code)
    return min((sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs(n))) for p in permutations(range(n))),default=0)
def induced(n,code,vs):
    a=arcs(n,code);return canon(len(vs),encode([[a[i][j] for j in vs]for i in vs]))
def source(n,code):
    a=arcs(n,code)
    for row in a:row.append(0)
    a.append([1]*n+[0]);return canon(n+1,encode(a))
N=5
classes=[[0]]
for n in range(1,N+1):
    vals=set()
    for code in classes[-1]:
        a=arcs(n-1,code)
        for mask in range(1<<(n-1)):
            b=[row[:] + [(mask>>i)&1]for i,row in enumerate(a)]
            b.append([1-((mask>>i)&1)for i in range(n-1)]+[0])
            vals.add(canon(n,encode(b)))
    classes.append(sorted(vals))
check([len(c)for c in classes]==[1,1,1,2,4,12],'class_counts')
ids=[{c:i for i,c in enumerate(cs)}for cs in classes]
@lru_cache(None)
def D(n,m):
    out=[]
    for code in classes[n+m]:
        row=[0]*len(classes[n])
        for vs in combinations(range(n+m),n):row[ids[n][induced(n+m,code,vs)]]+=1
        out.append(row)
    return Matrix(out)
@lru_cache(None)
def H(n):
    if n==0:return Matrix.zeros(0,1)
    out=[]
    for code in classes[n-1]:
        row=[0]*len(classes[n]);row[ids[n][source(n-1,code)]]=1;out.append(row)
    return Matrix(out)
def hpow(n,k):
    out=Matrix.eye(len(classes[n]));degree=n
    for _ in range(k):
        if degree<0:return Matrix.zeros(0,len(classes[n]))
        out=H(degree)*out;degree-=1
    return out
# Integral identities on all class bases within the available graded range.
for n in range(N+1):
    for a in range(N-n+1):
        for b in range(N-n-a+1):
            check(D(n+b,a)*D(n,b)==comb(a+b,a)*D(n,a+b),'divided_power',(n,a,b))
    for m in range(1,N-n+1):
        right=D(n-1,m)*H(n) if n else Matrix.zeros(len(classes[m-1]),1)
        check(H(n+m)*D(n,m)-right==D(n,m-1),'source_commutator',(n,m))
# Integral P, source-free basis, inverse extraction, and diagonalization.
bases=[];core_ranks=[];labels=[];dets=[]
for n in range(N+1):
    cols=[];labs=[]
    for r in range(n+1):
        for k,code in enumerate(classes[r]):
            a=arcs(r,code)
            if any(sum(row)==r-1 for row in a):continue
            v=Matrix.zeros(len(classes[r]),1);v[k]=1
            cols.append(D(r,n-r)*v);labs.append((r,k))
    B=Matrix.hstack(*cols);bases.append(B);labels.append(labs)
    determinant=int(B.det());dets.append(determinant)
    check(abs(determinant)==1,'integral_unimodularity',n)
    P=Matrix.zeros(len(classes[n]))
    for j in range(n+1):P+=(-1)**j*D(n-j,j)*hpow(n,j)
    check(H(n)*P==Matrix.zeros(max(0,len(classes[n-1]) if n else 0),len(classes[n])),'hP_zero',n)
    check(P*P==P,'P_idempotent',n)
    core_ranks.append(int(P.rank()))
for n in range(N+1):
    for m in range(N-n+1):
        expected=Matrix.zeros(len(classes[n+m]),len(classes[n]))
        for j,(r,k)in enumerate(labels[n]):expected[j,j]=comb(n-r+m,m)
        check(bases[n+m].inv()*D(n,m)*bases[n]==expected,'integral_diagonalization',(n,m))
# Full-function contraction cases, not just constant-function tests.
for p in (2,3,5):
    for s in (1,p,p*p):
        if (p-1)*s>N:continue
        for n in range(N-(p-1)*s+1):
            total=Matrix.zeros(len(classes[n]))
            for j in range(p):
                degree=n;v=Matrix.eye(len(classes[n]))
                for _ in range(j):v=D(degree,s)*v;degree+=s
                v=hpow(degree,(p-1)*s)*v;degree-=(p-1)*s
                if degree<0:continue
                for _ in range(p-1-j):v=D(degree,s)*v;degree+=s
                total-=v
            check(all(int(x)%p==int(y)%p for x,y in zip(total,Matrix.eye(len(classes[n])))),'p_contraction',(p,s,n))
# Scalar Frobenius/carry and borrow checks independently via exact binomials.
for p in (2,3,5,7):
    for a in range(65):
        for b in range(65):
            aa=a;bb=b;carry=False
            while aa or bb:
                carry|=(aa%p+bb%p>=p);aa//=p;bb//=p
            check((comb(a+b,a)%p!=0)==(not carry),'carry_scalar',(p,a,b))
    for e in range(3):
        s=p**e
        for f in range(3):
            t=p**f
            for n in range(32):
                # h^s D_t - D_t h^s, on the degree-n constant function.
                lhs=(comb(n+t,t) if n+t>=s else 0)-(comb(n-s+t,t) if n>=s else 0)
                rhs=0 if e>f else (1 if e==f else comb(n+t-s,t-s))
                check((lhs-rhs)%p==0,'borrow_commutator',(p,e,f,n))
# Aut divisibility for every oriented graph on the cyclic three-vertex host.
cyclic=[[0,1,0],[0,0,1],[1,0,0]]
from itertools import product
for options in product((0,1,2),repeat=3):
    constraints=[]
    for (i,j),v in zip(pairs(3),options):
        if v:constraints.append((i,j)if v==1 else(j,i))
    embeddings=sum(all(cyclic[perm[i]][perm[j]]for i,j in constraints)for perm in permutations(range(3)))
    check(embeddings%3==0,'cyclic_embedding_blindness',options)
print(json.dumps({'status':'PASS','checks':dict(counts),'class_counts':[len(c)for c in classes],'source_free_ranks':core_ranks,'basis_determinants':dets},indent=2))
