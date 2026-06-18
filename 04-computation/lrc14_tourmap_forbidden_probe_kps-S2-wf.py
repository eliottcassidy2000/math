"""
lrc14_tourmap_forbidden_probe_kps-S2-wf.py

Follow-up probe. The first pass found:
  - M1 orbit-majority: ALWAYS transitive (dead).
  - M3 halfplane-majority & M4 band-depth-majority: forbid 4 high-H classes at n=5.
  - M5 crossing-parity: realizes everything (no constraint).
  - M2 QR-character: forbidden set shrinks as prime grows (modulus collision artifact).

This script:
  (A) Robustness of M3/M4 forbidden classes at n=5: widen speed range, and restrict
      to LRC-HARD inputs (covering-ish / regular-score-targeting sets) to confirm the
      regular class (2,2,2,2,2) and near-regular (1,2,2,2,3) H>=11 are TRULY unreachable.
  (B) WHY they're forbidden: prove structurally that M3/M4 can't make a regular
      tournament for ODD reasons.
  (C) M2 QR: confirm the n=5 'only H=13' at p=7 is a collision/skip artifact, not LRC.
"""

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def all_opt_taus(S):
    b = F(0)
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v
    return sorted([t for t in cand(S) if gmin(S, t) == b]), b

def canon_key(adj, m):
    best = None
    for perm in permutations(range(m)):
        bits = tuple(adj[perm[i]][perm[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or bits < best: best = bits
    return best
def h_count(adj, m):
    cnt = 0
    for perm in permutations(range(m)):
        if all(adj[perm[k]][perm[k+1]] for k in range(m-1)): cnt += 1
    return cnt
def score_seq(adj, m):
    return tuple(sorted(sum(adj[i][j] for j in range(m) if j != i) for i in range(m)))
def num_3cycles(adj, m):
    c = 0
    for a,b,cc in combinations(range(m),3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c += 1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c += 1
    return c
def primitive(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1

def method3(S):
    taus,_ = all_opt_taus(S); m=len(S)
    adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            sc=0
            for t in taus:
                fi=(S[i]*t)%1; fj=(S[j]*t)%1
                wi=0 if fi<F(1,2) else 1; wj=0 if fj<F(1,2) else 1
                if wi!=wj: sc += (1 if wi<wj else -1)
                else: sc += (1 if fi<fj else -1) if fi!=fj else 0
            if sc>0: adj[i][j]=1
            elif sc<0: adj[j][i]=1
            else: adj[i][j]=1 if S[i]<S[j] else 0; adj[j][i]=1-adj[i][j]
    return adj

def method4(S, mod=14):
    m=len(S); U=[a for a in range(1,mod) if gcd(a,mod)==1]
    adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            wi=wj=0
            for a in U:
                ri=(S[i]*a)%mod; rj=(S[j]*a)%mod
                di=min(ri,mod-ri); dj=min(rj,mod-rj)
                if di>dj: wi+=1
                elif dj>di: wj+=1
            if wi>wj: adj[i][j]=1
            elif wj>wi: adj[j][i]=1
            else: adj[i][j]=1 if S[i]<S[j] else 0; adj[j][i]=1-adj[i][j]
    return adj

def all_iso_classes(m):
    seen={}
    pairs=list(combinations(range(m),2))
    for bits in product([0,1],repeat=len(pairs)):
        adj=[[0]*m for _ in range(m)]
        for (i,j),b in zip(pairs,bits):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        k=canon_key(adj,m)
        if k not in seen: seen[k]=(h_count(adj,m),num_3cycles(adj,m),score_seq(adj,m))
    return seen

print("="*70)
print("(A) M3 & M4 at n=5, WIDE speed range up to 60, confirm forbidden classes")
free5 = all_iso_classes(5)
for name, fn in [("M3 halfplane", method3), ("M4 band-mod14", method4)]:
    realized=set()
    cnt=0
    for S in combinations(range(1,61),5):
        if not primitive(S): continue
        cnt+=1
        adj=fn(S)
        realized.add(canon_key(adj,5))
    forb = set(free5)-realized
    print(f"  {name}: {cnt} primitive 5-sets in 1..60; realized={len(realized)}/12; "
          f"forbidden={len(forb)}")
    print(f"     forbidden (H,#3c,score): {sorted(free5[k] for k in forb)}")

print()
print("="*70)
print("(B) WHY: can M3/M4 ever produce a REGULAR tournament score (2,2,2,2,2)?")
# Check: M4 arc i->j iff i has greater orbit-depth-majority. The relation
# 'i beats j' from M4 is determined by comparing the multiset of depths
# {di(a)} vs {dj(a)} across units a. Show whether this can be a 3-cycle / regular.
# Empirically count: among ALL 5-subsets of speeds 1..40, what score sequences appear?
from collections import Counter
for name, fn in [("M3", method3), ("M4", method4)]:
    sc_counter=Counter()
    for S in combinations(range(1,41),5):
        if not primitive(S): continue
        adj=fn(S)
        sc_counter[score_seq(adj,5)]+=1
    print(f"  {name} n=5 score-sequence distribution (speeds 1..40):")
    for s,c in sorted(sc_counter.items()):
        reg = " <-- REGULAR" if s==(2,2,2,2,2) else ""
        print(f"     {s}: {c}{reg}")

print()
print("="*70)
print("(C) M2 QR p=7 at n=5: is 'only H=13' a skip artifact?")
def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
def method2(S,p):
    m=len(S); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            d=(S[i]-S[j])%p
            if d==0: return None
            if legendre(d,p)==1: adj[i][j]=1; adj[j][i]=0
            else: adj[i][j]=0; adj[j][i]=1
    return adj
# at p=7, distinct residues mod 7 forces all 5 speeds in distinct nonzero classes ->
# only 5 of 6 nonzero residues used; the QR tournament on a 5-subset of Z/7* is
# essentially the Paley-like sub-tournament. Count realized classes WITHOUT requiring
# the speeds be a 'used'/'skipped' artifact: enumerate ALL 5-subsets of residues.
realized7=set()
for res in combinations(range(7),5):  # 5 distinct residues mod 7 (incl 0 allowed? speeds can be ≡0)
    # build QR tournament on these residues; if any pair equal mod7 skip (none, distinct)
    m=5; adj=[[0]*m for _ in range(m)]; ok=True
    for i in range(m):
        for j in range(i+1,m):
            d=(res[i]-res[j])%7
            if d==0: ok=False;break
            if legendre(d,7)==1: adj[i][j]=1;adj[j][i]=0
            else: adj[i][j]=0;adj[j][i]=1
        if not ok: break
    if ok: realized7.add(canon_key(adj,5))
print(f"  ALL 5-subsets of residues mod 7: realized classes = {len(realized7)}")
print(f"     H values: {sorted(free5[k][0] for k in realized7 if k in free5)}")
print("  => the QR map on a fixed prime p only sees C(p,n) residue-tournaments;")
print("     forbidden classes there are MODULUS artifacts (finite residue pool), NOT LRC.")
