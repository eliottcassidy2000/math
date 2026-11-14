"""Literal profile routines adapted from the S326 replay."""
from itertools import permutations,product
from collections import defaultdict
def extend(A,cut):
    n=len(A)
    return tuple(r|(((cut>>i&1)^1)<<n) for i,r in enumerate(A))+(cut,)

def canon(A):
    groups=defaultdict(list)
    for i,r in enumerate(A): groups[r.bit_count()].append(i)
    best=None
    for runs in product(*(tuple(permutations(groups[d])) for d in sorted(groups))):
        p=tuple(x for run in runs for x in run)
        row=tuple(sum((A[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
        if best is None or row<best: best=row
    return best

def strong(A):
    for i in range(len(A)):
        rb=1<<i
        while True:
            b=rb
            for j,r in enumerate(A):
                if rb>>j&1: b|=r
            if b==rb: break
            rb=b
        if rb!=(1<<len(A))-1: return False
    return True

def hp_all(A):
    n=len(A); N=1<<n
    E=[[0]*n for _ in range(N)]; H=[1]+[0]*(N-1)
    for i in range(n): E[1<<i][i]=1
    for s in range(1,N):
        H[s]=sum(E[s])
        for i,e in enumerate(E[s]):
            if not e: continue
            W=A[i]&~s
            while W:
                b=W&-W; W-=b
                E[s|b][b.bit_length()-1]+=e
    return H

def profile(A):
    H=hp_all(A)
    return H[-1],sum(a*b for a,b in zip(H,H[::-1]))

