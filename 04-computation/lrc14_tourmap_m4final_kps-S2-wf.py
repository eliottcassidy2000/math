"""M4 final: confirm the loneliness lever. The regular class (rotational T5,
c3=5) is the UNIQUE class M4 cannot reach at optimal lonely tau but CAN reach
at non-optimal crossings. This is the cleanest 'loneliness forbids a class'
signal. Quantify across vmax and report exact counts."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
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
def c3(A):
    n=len(A); c=0
    for i in range(n):
        for j in range(n):
            if i==j or not A[i][j]: continue
            for k in range(n):
                if k==i or k==j: continue
                if A[j][k] and A[k][i]: c+=1
    return c//3
def is_primitive(S):
    g0=0
    for v in S: g0=gcd(g0,v)
    return g0==1
def gen(n,vmax): return [c for c in combinations(range(1,vmax+1),n) if is_primitive(c)]
def m4(S,tau):
    S=sorted(set(S)); n=len(S); A=[[0]*n for _ in range(n)]
    dd={(a,b):nrm((a-b)*tau) for a in S for b in S}
    for a in range(n):
        for b in range(a+1,n):
            i,j=S[a],S[b]; wi=wj=0
            for k in S:
                if k==i or k==j: continue
                if dd[(i,k)]<dd[(j,k)]: wi+=1
                elif dd[(j,k)]<dd[(i,k)]: wj+=1
            if wi>wj or (wi==wj and i<j): A[a][b]=1
            else: A[b][a]=1
    return A
# regular class is the only one with score (2,2,2,2,2)
print("M4: count speed sets whose optimal-tau tournament is REGULAR (2,2,2,2,2):")
for vmax in [11,14,17,20,23]:
    sets=gen(5,vmax); reg_opt=0; reg_any=0
    for S in sets:
        _,tau=M(S); A=m4(S,tau)
        if score_seq(A)==(2,2,2,2,2): reg_opt+=1
    print(f"  vmax={vmax}: {len(sets)} primitive 5-speed sets; #regular at OPTIMAL lonely tau = {reg_opt}")
