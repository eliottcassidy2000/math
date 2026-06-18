#!/usr/bin/env python3
"""
Final: M5 mechanism + N-dependence, and a combined statement.
M5 arc i->j (i<j) iff #{a in {1,3,5}: (v_i a)%14 < (v_j a)%14} is ODD.
This is a GF(2) parity aggregation of 3 half-grid section orders.
Check:
  1. Is M5's regular-T_5 forbiddance N-dependent? (control N=15,20 with half-units)
  2. M5 nontriviality confirmed (it has 3-cycles), and forbids exactly the
     regular T_5 over primitive sets vmax=40.
  3. Combined: at modulus 14 BOTH pinning maps forbid the H-maximizing
     self-complementary regular tournament; at non-14 moduli they don't.
"""
from math import gcd
from itertools import combinations
from functools import reduce
def is_primitive(S): return reduce(gcd,S)==1
def score(adj,k): return tuple(sorted(sum(adj[i][j] for j in range(k) if j!=i) for i in range(k)))
def is_tour(adj,k): return all(adj[i][j]!=adj[j][i] for i in range(k) for j in range(i+1,k))
def num3(adj,k):
    c=0
    for i,j,l in combinations(range(k),3):
        if adj[i][j] and adj[j][l] and adj[l][i]:c+=1
        if adj[i][l] and adj[l][j] and adj[j][i]:c+=1
    return c

def m5(S,N,half):
    k=len(S); adj=[[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1,k):
            c2=sum(1 for a in half if (S[i]*a)%N < (S[j]*a)%N)
            if c2%2==1: adj[i][j]=True
            else: adj[j][i]=True
    return adj if is_tour(adj,k) else None

# half = first half of units(N)
def half_units(N):
    Uu=[a for a in range(1,N) if gcd(a,N)==1]
    return Uu[:len(Uu)//2]

print("### M5 regular-T_5 reachability over primitive 5-sets, various N ###")
for N in (14,15,16,20,21):
    half=half_units(N)
    found=None
    nt=False
    for S in combinations(range(1,30),5):
        if not is_primitive(S): continue
        adj=m5(S,N,half)
        if adj is None: continue
        if num3(adj,5)>0: nt=True
        if score(adj,5)==(2,2,2,2,2):
            found=S; break
    print(f"  N={N} half-units={half}: nontrivial={nt}  regular-T5 found={found}")

print("\n### M3 vs M5 forbidden-class agreement at N=14 (vmax=40 primitive) ###")
# already established: both forbid regular T_5. confirm jointly here briefly at vmax=20.
N=14
U=[a for a in range(1,N) if gcd(a,N)==1]
def depth(x): r=x%N;return min(r,N-r)
def dvec(v): return tuple(depth((v*a)%N) for a in U)
def margin(ri,rj):
    Di,Dj=dvec(ri),dvec(rj)
    return sum((Di[c]>Dj[c])-(Di[c]<Dj[c]) for c in range(len(U)))
def m3(S):
    k=len(S);adj=[[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i==j:continue
            m=margin(S[i]%N,S[j]%N)
            if m>0:adj[i][j]=True
            elif m<0:adj[i][j]=False
            else:adj[i][j]=S[i]<S[j]
    return adj if is_tour(adj,k) else None
m3reg=False;m5reg=False
half=(1,3,5)
for S in combinations(range(1,21),5):
    if not is_primitive(S):continue
    a3=m3(S)
    if a3 and score(a3,5)==(2,2,2,2,2): m3reg=True
    a5=m5(S,14,half)
    if a5 and score(a5,5)==(2,2,2,2,2): m5reg=True
print(f"  vmax=20: M3 realizes regular T_5={m3reg}; M5 realizes regular T_5={m5reg}")
print("  (both False => both pinning maps forbid the H-maximizing regular tournament at N=14)")
