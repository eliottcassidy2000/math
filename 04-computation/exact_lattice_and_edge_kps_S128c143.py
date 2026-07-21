#!/usr/bin/env python3
"""
EXACT-ARITHMETIC recompute of the invariant lattice + the 2-adic edge of H.
kind-pasteur-2026-07-21-S128c143.

Motivation: the S128c142 census computed disc, arb via numpy round(det) — a FLOATING-POINT hazard
(the sampler S128c143 found exact poly-tuple collisions at n=6 that numpy had spuriously split).
This redoes the WHOLE exhaustive census n=3..6 with EXACT integer arithmetic (Faddeev-LeVerrier char
polys, Bareiss determinant, exact disc via charpoly_K(-1)) and:
  (A) prints |disc|,|arb| class-counts EXACT vs the numpy values (correction check for THM-1965),
  (B) recomputes the disc/arb refinement relations exactly (fixing any float artifact),
  (C) computes the EXACT 2-adic edge: group iso classes by the exact poly-tuple P=(score,specA,
      specS,disc,arb); within each cell M_cell = gcd of H-differences; M = gcd over cells;
      v2(M) = poly-determined low bits of H.  (Redei => v2>=1.)
"""
import itertools, math
from functools import reduce
from fractions import Fraction
from collections import defaultdict

def setup(n):
    pairs = list(itertools.combinations(range(n), 2)); pos = {pr: k for k, pr in enumerate(pairs)}
    maps = []
    for p in itertools.permutations(range(n)):
        m = [ (pos[(p[i],p[j])],0) if p[i]<p[j] else (pos[(p[j],p[i])],1) for (i,j) in pairs ]
        maps.append(m)
    return pairs, maps

def canon_orbit(bits, maps, npairs):
    seen=set()
    for m in maps:
        v=0
        for k in range(npairs):
            b=(bits>>k)&1; tb,inv=m[k]
            if inv: b^=1
            v|=b<<tb
        seen.add(v)
    return min(seen), len(seen)

def adj(bits, pairs, n):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs):
        if (bits>>k)&1: A[i][j]=1
        else: A[j][i]=1
    return A

def matmul(A,B,n): return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def charpoly_int(M,n):
    c=[1]; Mk=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for k in range(1,n+1):
        Mk=matmul(M,Mk,n); tr=sum(Mk[i][i] for i in range(n)); ck=-tr//k
        c.append(ck)
        for i in range(n): Mk[i][i]+=ck
    return tuple(c)
def bareiss(M,n):
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
def arb(A,n):
    indeg=[sum(A[r][c] for r in range(n)) for c in range(n)]
    Lm=[[ (indeg[i] if i==j else 0)-A[j][i] for j in range(n)] for i in range(n)]
    return bareiss([r[1:] for r in Lm[1:]],n-1)
def cyc(A,n):
    out=[]
    for k in range(3,n+1):
        cnt=0
        for S in itertools.combinations(range(n),k):
            s0=S[0]
            for perm in itertools.permutations(S[1:]):
                seq=(s0,)+perm
                if all(A[seq[i]][seq[(i+1)%k]] for i in range(k)): cnt+=1
        out.append(cnt)
    return tuple(out)
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
def Rsigned(A,n):
    dp=[defaultdict(int) for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        if not dp[mask]: continue
        for v,val in list(dp[mask].items()):
            if val==0: continue
            for w in range(n):
                if (mask>>w)&1 or not A[v][w]: continue
                inv=bin(mask>>(w+1)).count("1"); sgn=-1 if inv&1 else 1
                dp[mask|(1<<w)][w]+=val*sgn
    return abs(sum(dp[(1<<n)-1].values()))

def census(n):
    pairs,maps=setup(n); npairs=len(pairs); reps={}
    for bits in range(2**npairs):
        c,osz=canon_orbit(bits,maps,npairs)
        if c not in reps: reps[c]=(bits,osz)
    data={}
    for c,(bits,osz) in reps.items():
        A=adj(bits,pairs,n)
        data[c]=dict(score=score(A,n),specA=charpoly_int(A,n),specS=charpoly_int(Kmat(A,n),n),
                     cyc=cyc(A,n),H=Hcount(A,n),R=Rsigned(A,n),disc=disc(A,n),arb=arb(A,n),
                     aut=math.factorial(n)//osz)
    return data

def ncls(data,f): return len(set(d[f] for d in data.values()))
def refines(data,f,g):
    b=defaultdict(set)
    for d in data.values(): b[d[f]].add(d[g])
    return all(len(s)==1 for s in b.values())
def gcd_list(xs):
    xs=[x for x in xs if x]; return reduce(math.gcd,xs) if xs else 0
def v2(m):
    if m==0: return math.inf
    k=0
    while m%2==0: m//=2; k+=1
    return k

INVS=["score","specA","specS","cyc","H","R","disc","arb","aut"]
print("="*76); print("EXACT-ARITHMETIC recompute (fixes numpy disc/arb float artifact)"); print("="*76)
CEN={}
NUMPY={4:dict(disc=2,arb=3),5:dict(disc=2,arb=7),6:dict(disc=5,arb=32)}  # S128c142 numpy values
for n in range(3,7):
    CEN[n]=census(n)
    line=f"n={n}: iso={len(CEN[n])}  "+"  ".join(f"|{f}|={ncls(CEN[n],f)}" for f in INVS)
    print(line)
    if n in NUMPY:
        for f in ("disc","arb"):
            ex=ncls(CEN[n],f); nv=NUMPY[n][f]
            flag="  <-- numpy said "+str(nv)+" (FLOAT ARTIFACT)" if ex!=nv else "  (matches numpy)"
            print(f"      exact |{f}|={ex}{flag}")

print()
print("="*76); print("(A) EXACT refinement relations for disc/arb (correcting THM-1965)"); print("="*76)
data6=CEN[6]
for f,g in [("disc","H"),("H","disc"),("disc","specS"),("specS","disc"),("arb","specA"),
            ("specA","arb"),("arb","cyc"),("cyc","arb"),("arb","H"),("H","arb"),
            ("disc","specA"),("specA","disc"),("disc","arb"),("arb","disc")]:
    fails=None
    for n in range(3,7):
        if not refines(CEN[n],f,g): fails=n; break
    print(f"  {f:>6} refines {g:<6}? "+("through n=6 (YES)" if fails is None else f"FAILS at n={fails}"))

print()
print("="*76); print("(C) THE EXACT 2-ADIC EDGE  (poly-tuple P=(score,specA,specS,disc,arb))"); print("="*76)
POLY=("score","specA","specS","disc","arb")
for n in range(4,7):
    data=CEN[n]; cells=defaultdict(list)
    for c,d in data.items(): cells[tuple(d[f] for f in POLY)].append(d)
    diffs=[]; split=0; wit=[]
    for key,mem in cells.items():
        Hs=sorted(set(m["H"] for m in mem))
        if len(Hs)>1:
            split+=1; wit.append(Hs)
            for i in range(len(Hs)):
                for j in range(i+1,len(Hs)): diffs.append(Hs[j]-Hs[i])
    M=gcd_list(diffs)
    print(f" n={n}: iso={len(data)}  EXACT poly-cells={len(cells)}  H-split cells={split}")
    if split:
        print(f"        M=gcd(H-diffs)={M}  v2(M)={v2(M)}  odd={M>>v2(M)}  "
              f"=> H mod {2**v2(M)} poly-time, H mod {2**(v2(M)+1)} NOT")
        print(f"        witnesses (H-values per split cell): {wit[:6]}")
    else:
        print(f"        poly-tuple DETERMINES H (and iso) exactly at n={n}")
print("DONE.")
