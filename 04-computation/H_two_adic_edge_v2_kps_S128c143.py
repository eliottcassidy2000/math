#!/usr/bin/env python3
"""
THE 2-ADIC EDGE OF H, clean version — how many low bits of H does the SPECTRUM pin?
kind-pasteur-2026-07-21-S128c143 (v2, after the arb-root bug: arb-rooted-at-0 is NOT iso-invariant;
here arb_inv = sorted tuple of per-root arborescence counts, a genuine invariant).

Bulletproof route to the user's "H sits at the formula/no-formula edge": group iso classes by an
EXACT poly-time invariant; within each group H varies; the gcd of within-group H-differences = M;
v2(M) = number of low 2-adic bits of H pinned by that invariant. Redei => H odd => v2 >= 1 always.

Partitions tested (all EXACT, all poly-time): specA (spectrum = THM-1780 ladder), (specA,score),
(specA,specS,score,disc,arb_inv) [full poly-tuple]. The edge = how fast v2 -> 1 as n grows.
Everything exact integer (Faddeev-LeVerrier, Bareiss); no floats.
"""
import itertools, math
from functools import reduce
from fractions import Fraction
from collections import defaultdict

def setup(n):
    pairs=list(itertools.combinations(range(n),2)); pos={pr:k for k,pr in enumerate(pairs)}
    maps=[[ (pos[(p[i],p[j])],0) if p[i]<p[j] else (pos[(p[j],p[i])],1) for (i,j) in pairs]
          for p in itertools.permutations(range(n))]
    return pairs,maps
def canon_orbit(bits,maps,npairs):
    seen=set()
    for m in maps:
        v=0
        for k in range(npairs):
            b=(bits>>k)&1; tb,inv=m[k]
            if inv:b^=1
            v|=b<<tb
        seen.add(v)
    return min(seen),len(seen)
def adj(bits,pairs,n):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs):
        if (bits>>k)&1:A[i][j]=1
        else:A[j][i]=1
    return A
def matmul(A,B,n): return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def charpoly_int(M,n):
    c=[1]; Mk=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for k in range(1,n+1):
        Mk=matmul(M,Mk,n); tr=sum(Mk[i][i] for i in range(n)); ck=-tr//k; c.append(ck)
        for i in range(n): Mk[i][i]+=ck
    return tuple(c)
def bareiss(M,n):
    if n==0: return 1
    M=[r[:] for r in M]; sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            sw=next((r for r in range(k+1,n) if M[r][k]!=0),None)
            if sw is None: return 0
            M[k],M[sw]=M[sw],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]
def score(A,n): return tuple(sorted(sum(A[i]) for i in range(n)))
def Kmat(A,n): return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
def disc(A,n):
    cK=charpoly_int(Kmat(A,n),n); val=sum(cK[k]*((-1)**(n-k)) for k in range(n+1))
    return Fraction(abs(((-1)**n)*val),2**(n-1))
def arb_root(A,n,r):
    indeg=[sum(A[x][c] for x in range(n)) for c in range(n)]
    Lm=[[ (indeg[i] if i==j else 0)-A[j][i] for j in range(n)] for i in range(n)]
    minor=[[Lm[i][j] for j in range(n) if j!=r] for i in range(n) if i!=r]
    return bareiss(minor,n-1)
def arb_inv(A,n): return tuple(sorted(arb_root(A,n,r) for r in range(n)))  # ISO-INVARIANT
def Hcount(A,n):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            if not (mask>>v)&1 or dp[mask][v]==0: continue
            for w in range(n):
                if (mask>>w)&1 or not A[v][w]: continue
                dp[mask|(1<<w)][w]+=dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))

def census(n):
    pairs,maps=setup(n); npairs=len(pairs); reps={}
    for bits in range(2**npairs):
        c,osz=canon_orbit(bits,maps,npairs)
        if c not in reps: reps[c]=bits
    data=[]
    for c,bits in reps.items():
        A=adj(bits,pairs,n)
        data.append(dict(score=score(A,n),specA=charpoly_int(A,n),specS=charpoly_int(Kmat(A,n),n),
                         disc=disc(A,n),arb=arb_inv(A,n),H=Hcount(A,n)))
    return data
def gcd_list(xs):
    xs=[x for x in xs if x]; return reduce(math.gcd,xs) if xs else 0
def v2(m):
    if m==0: return None
    k=0
    while m%2==0: m//=2; k+=1
    return k

def edge(data,keys,label,n):
    cells=defaultdict(list)
    for d in data: cells[tuple(d[k] for k in keys)].append(d["H"])
    diffs=[]; split=0; wit=[]
    for Hs in cells.values():
        u=sorted(set(Hs))
        if len(u)>1:
            split+=1; wit.append(u)
            for i in range(len(u)):
                for j in range(i+1,len(u)): diffs.append(u[j]-u[i])
    M=gcd_list(diffs)
    if split==0:
        print(f"  n={n} [{label}]: {len(cells)} cells, H-split=0  => this invariant DETERMINES H")
    else:
        vv=v2(M)
        print(f"  n={n} [{label}]: cells={len(cells)} H-split={split}  M=gcd(H-diffs)={M} v2={vv}"
              f"  => spectrum pins H mod {2**vv}; H mod {2**(vv+1)} FREE (#P). witnesses={wit[:5]}")
    return M

print("="*78)
print("THE 2-ADIC EDGE OF H — how many low bits does a poly-time invariant pin?  (exact)")
print("="*78)
for n in range(4,7):
    data=census(n)
    edge(data,("specA",),"spectrum",n)
    edge(data,("specA","score"),"spectrum+score",n)
    edge(data,("specA","specS","score","disc","arb"),"full poly-tuple",n)
    print()
print("Redei = 'H mod 2 = 1', the parity bit, is the ONE universally-formula bit.")
print("The edge question: does the poly-determined depth v2 -> 1 (only parity) as n grows?")
print("DONE.")
