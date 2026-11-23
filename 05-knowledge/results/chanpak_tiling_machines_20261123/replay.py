"""Finite controls for the tiling-language, counter, and complement claims."""
from itertools import combinations, permutations
from functools import cache

def from_arcs(n,E):
    A=[0]*n
    for i,j in E: A[i]|=1<<j
    return tuple(A)

def tt(n):
    return from_arcs(n,combinations(range(n),2))

def join(U,V):
    n=len(U); m=len(V)
    return tuple(x|(((1<<m)-1)<<n) for x in U)+tuple(y<<n for y in V)

def hp(A):
    return {P for P in permutations(range(len(A)))
            if all(A[i]>>j&1 for i,j in zip(P,P[1:]))}

@cache
def count(A):
    if not A: return 1
    return len(hp(A))

def residuals(language):
    rows={}
    for w in language:
        for n in range(len(w)+1):
            rows.setdefault(w[:n],set()).add(w[n:])
    return {frozenset(s) for s in rows.values()}|{frozenset()}

for m in range(1,10):
    W=tuple(range(m)); a,b=m,m+1
    X={tuple(v for v in W if S>>v&1)+(b,a)+tuple(v for v in W if not S>>v&1)
       for S in range(1<<m)}
    X.add((a,)+W+(b,))
    assert len(X)==(1<<m)+1
    assert len(residuals(X))==3*(1<<m)+m+2
    if m<=4:
        A=from_arcs(m+2,list(combinations(W,2))+[(b,a)]+[(a,v) for v in W]+[(v,b) for v in W])
        assert hp(A)==X
print('exact streaming-state count, m=1 through 9; literal tournaments through m=4')

def sub(A,S):
    return tuple(sum((A[i]>>j&1)<<k for k,j in enumerate(S)) for i in S)

C=(2,4,1)
def z2(A):
    n=len(A)
    return sum(count(sub(A,[i for i in range(n) if s>>i&1]))*
               count(sub(A,[i for i in range(n) if not s>>i&1]))
               for s in range(1<<n))
for a in range(4):
    for b in range(4):
        A=join(join(tt(a),C),tt(b))
        assert count(A)==3 and z2(A)==12*2**(a+b)
print('literal two-gap counts all H=3; complete two-slot reads see only total padding')

def phi(A,P):
    n=len(A); f=(1<<n)-1
    B=[f^x^(1<<i) for i,x in enumerate(A)]
    for i,j in zip(P,P[1:]):
        assert A[i]>>j&1
        B[i]^=1<<j;B[j]^=1<<i
    return tuple(B)

def scc_width(A,drop):
    n=len(A); M=((1<<n)-1)^sum(1<<v for v in drop)
    R={}
    for i in range(n):
        if not M>>i&1: continue
        s=1<<i
        for _ in range(n):
            for j in range(n):
                if s>>j&1:s|=A[j]&M
        R[i]=s
    return max(sum((s>>j&1) and (t>>i&1) for j,t in R.items())
               for i,s in R.items())
def w(A,r):
    return min(scc_width(A,X) for X in combinations(range(len(A)),r))
def degrees(A):return [x.bit_count() for x in A]
def dv(P):
    return [int(P[0]==i)-int(P[-1]==i) for i in range(len(P))]

# Exact two-line controls on all 1024 five-vertex tournaments.
E=list(combinations(range(5),2))
for x in range(1<<len(E)):
    A=from_arcs(5,[(u,v) if x>>i&1 else (v,u) for i,(u,v) in enumerate(E)])
    P=min(hp(A)); B=phi(A,P)
    assert degrees(B)==[4-u+v for u,v in zip(degrees(A),dv(P))]
    for Q in hp(B):
        F=phi(B,Q)
        assert degrees(F)==[u+v-r for u,v,r in zip(degrees(A),dv(Q),dv(P))]
        assert sum((u^v).bit_count() for u,v in zip(A,F))%4==0
# Width-one controls on a single coherent 13-letter frame.
W=()
for _ in range(4): W=join(W,C)
n=len(W); z=n
# Apex beats first block and loses to the other three.
S=7
A=tuple(x|((1<<n) if not S>>i&1 else 0) for i,x in enumerate(W))+(S,)
assert w(A,1)<=3
for i in range(n+1):
    Q=tuple(range(n+1)); Q=Q[i:]+Q[:i]
    assert w(phi(A,Q),1)>=4
print('two-line parity and degree identities; one-deletion complement width on all cyclic frames')
