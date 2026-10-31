#!/usr/bin/env python3
"""Clean-room exact rational checks; imports no primary checker/results."""
from fractions import Fraction as Q
from itertools import product
import json

def fail(label): raise RuntimeError(label)
def eye(n): return [[Q(i==j) for j in range(n)] for i in range(n)]
def add(A,B): return [[x+y for x,y in zip(ar,br)] for ar,br in zip(A,B)]
def neg(A): return [[-x for x in r] for r in A]
def scale(q,A): return [[q*x for x in r] for r in A]
def mul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def inverse(A):
    n=len(A); T=[list(map(Q,row))+r for row,r in zip(A,eye(n))]
    for k in range(n):
        pivot=next((i for i in range(k,n) if T[i][k]),None)
        if pivot is None: return None
        T[pivot],T[k]=T[k],T[pivot]
        d=T[k][k]; T[k]=[x/d for x in T[k]]
        for i in range(n):
            if i!=k:
                d=T[i][k]; T[i]=[x-d*y for x,y in zip(T[i],T[k])]
    return [r[n:] for r in T]
def solve(A,b):
    R=inverse(A)
    if R is None: fail('unexpected singular periodic system')
    return [row[0] for row in mul(R,[[x] for x in b])]
def owner(x):
    if x<=Q(-1,6): return 0
    if x<=Q(1,6): return 1
    return 2

words=valid=fixed=0
for c in product((Q(-1,3),Q(0),Q(1,3)),repeat=3):
    for l in (Q(1,4),Q(1,2),Q(3,4),Q(1)):
        fixed+=sum(owner(v)==i for i,v in enumerate(c))
        a=1-l
        for p in range(1,5):
            for code in product(range(3),repeat=p):
                # Solve all simultaneous cyclic affine equations, not formula (6).
                A=eye(p); rhs=[]
                for j in range(p):
                    A[j][(j-1)%p]-=a
                    rhs.append(l*c[code[(j-1)%p]])
                states=solve(A,rhs)
                admissible=all(owner(v)==i for v,i in zip(states,code))
                if admissible:
                    valid+=1
                    for j in range(p):
                        if states[(j+1)%p]!=a*states[j]+l*c[owner(states[j])]:
                            fail('cyclic equation not realized')
                words+=1
if words!=12960:fail('wrong declared word universe')

# Independent sign-sum lemma controls on all nonzero integer triples summing 0.
triples=0
for u,v in product(range(-9,10),repeat=2):
    t=-u-v
    if not u or not v or not t:continue
    signs=sum(1 if a>0 else -1 for a in (u,v,t))
    if abs(signs)!=1:fail('zero-sum sign triple')
    triples+=1

# Clean-room direct Gaussian inverse vs rank-two assembled update at 81 gains.
L=[[Q(-n*n if i==j else 0) for j in range(6)] for i,n in enumerate((1,2,3,6,7,8))]
U=list(map(lambda r:list(map(Q,r)),[[1,0],[0,1],[1,1],[2,-1],[-1,2],[1,3]]))
C=list(map(lambda r:list(map(Q,r)),[[2,1],[-3,4]]))
V=list(map(lambda r:list(map(Q,r)),[[1,2,0,-1,0,1],[0,1,1,0,2,-1]]))
R=inverse(L); B=mul(mul(U,C),V)
count=singular=0
for n in range(-40,41):
    gamma=Q(n,3)
    left=add(L,scale(-gamma,B))
    Sinv=inverse(add(eye(2),scale(-gamma,mul(mul(mul(V,R),U),C))))
    direct=inverse(left)
    if (Sinv is None)!=(direct is None):fail('singularity mismatch')
    if Sinv is None:
        singular+=1; continue
    correction=mul(mul(mul(mul(mul(mul(R,U),C),Sinv),V),R),eye(6))
    candidate=add(R,scale(gamma,correction))
    if candidate!=direct:fail('independent inverse mismatch')
    count+=1

# Excluded endpoint limit and two-cycle controls, without symbolic algebra.
for l in (Q(1,4),Q(1,2),Q(3,4),Q(1)):
    d=l/(3*(2-l))
    T=lambda b:(1-l)*b+l*(Q(1,3) if b<=0 else Q(-1,3))
    if T(-d)!=d or T(d)!=-d:fail('two-cycle')
b=Q(1,4)
for k in range(1,31):
    b=b/2
    if b!=Q(1,4*2**k) or b<=0:fail('excluded endpoint trajectory')
if Q(1,2)*Q(1,4)!=Q(1,8):fail('nonfixed limit')

print(json.dumps({'status':'PASS','arithmetic':'stdlib Fraction only',
 'abstract_templates':27,'gains':4,'periodic_words_checked':words,
 'admissible_words':valid,'branch_fixed_instances':fixed,
 'zero_sum_sign_triples':triples,'independent_gaussian_resolvent_gains':count,
 'singular_resolvent_gains':singular},sort_keys=True,indent=2))
