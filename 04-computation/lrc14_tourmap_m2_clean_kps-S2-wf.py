"""
Clean M2 test: a TRUE binding-pair crossing for difference speed d=v_i-v_j is
snrm(d*tau)=0 (the runners coincide) OR snrm(d*tau)=+/-1/2 (antipodal). BOTH are
degenerate for an orientation rule. A clean M2 tournament requires snrm(d*tau)
NOT in {0, 1/2} for every pair. Recheck the SC property excluding ALL degenerate
(0 and 1/2) pairs. Does M2 stay universally SC at n=6 when fully clean?
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def snrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else r-1
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def score_seq(A):
    n=len(A); return tuple(sorted(sum(A[i]) for i in range(n)))
def count_3cycles(A):
    n=len(A); c=0
    for i in range(n):
        for j in range(n):
            if i==j or not A[i][j]: continue
            for k in range(n):
                if k==i or k==j: continue
                if A[j][k] and A[k][i]: c+=1
    return c//3
def canon(A):
    n=len(A); best=None
    for p in permutations(range(n)):
        flat=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or flat<best: best=flat
    return best
def complement(A):
    n=len(A); B=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j: B[i][j]=A[j][i]
    return B
def is_SC(A): return canon(A)==canon(complement(A))
def is_primitive(S):
    g0=0
    for v in S: g0=gcd(g0,v)
    return g0==1
def gen_speedsets(n,vmax):
    return [c for c in combinations(range(1,vmax+1),n) if is_primitive(c)]
def adj_m2_clean(S,tau):
    """Valid only if snrm(d*tau) not in {0, 1/2} for every pair."""
    S=sorted(set(S)); n=len(S); A=[[0]*n for _ in range(n)]; valid=True
    for a in range(n):
        for b in range(a+1,n):
            i,j=S[a],S[b]
            s=snrm((i-j)*tau)
            if s==0 or s==F(1,2) or s==F(-1,2): valid=False
            elif s>0: A[a][b]=1
            else: A[b][a]=1
    return A,valid

print("CLEAN M2 (exclude snrm in {0,1/2}): is it universally SC?")
for n_base,vmax in [(5,11),(6,10),(6,11)]:
    tot=0; sc=0; first=None
    for S in gen_speedsets(n_base,vmax):
        for tau in cand(S):
            A,valid=adj_m2_clean(S,tau)
            if not valid: continue
            tot+=1
            if is_SC(A): sc+=1
            elif first is None: first=(S,tau,score_seq(A),count_3cycles(A))
    print(f"  n={n_base} vmax={vmax}: clean valid={tot}, SC={sc}, ALL SC? {sc==tot}")
    if first: print(f"    first non-SC: {first}")

print()
print("CLEAN M2 at OPTIMAL tau only, n=6:")
for vmax in [10,11,12]:
    tot=0; sc=0; realized={}; first=None
    for S in gen_speedsets(6,vmax):
        _,tau=M(S)
        A,valid=adj_m2_clean(S,tau)
        if not valid: continue
        tot+=1
        c=canon(A); realized[c]=is_SC(A)
        if not is_SC(A) and first is None: first=(S,tau,score_seq(A),count_3cycles(A))
    nonsc=sum(1 for v in realized.values() if not v)
    print(f"  n=6 vmax={vmax}: clean valid at optimal tau={tot}, distinct classes={len(realized)}, "
          f"non-SC classes among realized={nonsc}")
    if first: print(f"    first non-SC at optimal: {first}")

print("\nDONE.")
