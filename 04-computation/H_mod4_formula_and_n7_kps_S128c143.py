#!/usr/bin/env python3
"""
(a) find the explicit spectral quantity giving H mod 4; (b) test at n=7 whether H mod4 stays
spectral and H mod8 fails. kind-pasteur-2026-07-21-S128c143.
"""
import sys, itertools, math, random
from collections import defaultdict
sys.path.insert(0,'04-computation')
from H_two_adic_edge_v2_kps_S128c143 import (setup, canon_orbit, adj, charpoly_int, Kmat,
                                             Hcount, score)

def cyc_k(A,n,k):
    cnt=0
    for S in itertools.combinations(range(n),k):
        s0=S[0]
        for perm in itertools.permutations(S[1:]):
            seq=(s0,)+perm
            if all(A[seq[i]][seq[(i+1)%k]] for i in range(k)): cnt+=1
    return cnt

print("="*72); print("(a) EXPLICIT spectral quantity for (H-1)/2 mod 2  (n<=6, exhaustive)"); print("="*72)
for n in [4,5,6]:
    pairs,maps=setup(n); npairs=len(pairs); reps={}
    for bits in range(2**npairs):
        c,_=canon_orbit(bits,maps,npairs)
        if c not in reps: reps[c]=bits
    rows=[]
    for c,bits in reps.items():
        A=adj(bits,pairs,n)
        H=Hcount(A,n); spA=charpoly_int(A,n)
        cyc={k:cyc_k(A,n,k) for k in range(3,n+1)}
        rows.append((H,spA,cyc))
    g=lambda H:(H-1)//2 % 2                       # the 2nd 2-adic bit of H
    # candidates for g: c3,c4,c5,c6 mod2, sums, char-coeff mod2
    cand={}
    for k in range(3,n+1):
        cand[f'c{k} mod2']=lambda H,cyc,k=k:(cyc.get(k,0)%2)
    cand['c3+c4 mod2']=lambda H,cyc:((cyc.get(3,0)+cyc.get(4,0))%2)
    for k in range(1,n+1):
        cand[f'charA[{k}] mod2']=lambda H,spA,k=k:(spA[k]%2)
    print(f" n={n}:")
    for name,f in cand.items():
        try:
            ok=all(g(H)==(f(H,cyc) if 'char' not in name else f(H,spA)) for H,spA,cyc in rows)
        except Exception:
            ok=False
        if ok: print(f"    (H-1)/2 == {name}  for ALL classes  <== MATCH")
    # also just report g vs each cycle count correlation table absent; note if none matched
print()
print("="*72); print("(b) n=7 SAMPLING: is H mod4 spectral? is H mod8 spectral?"); print("="*72)
def rand_A(n,rng):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if rng.getrandbits(1): A[i][j]=1
            else: A[j][i]=1
    return A
n=7; rng=random.Random(7)
bucket=defaultdict(lambda:[set(),set()])   # specA -> [set(H%4), set(H%8)]
N=120000
for _ in range(N):
    A=rand_A(n,rng); spA=charpoly_int(A,n); H=Hcount(A,n)
    bucket[spA][0].add(H%4); bucket[spA][1].add(H%8)
mod4_spectral=all(len(v[0])==1 for v in bucket.values())
mod8_spectral=all(len(v[1])==1 for v in bucket.values())
split4=[k for k,v in bucket.items() if len(v[0])>1]
split8=[k for k,v in bucket.items() if len(v[1])>1]
print(f" n=7: samples={N}  distinct specA buckets={len(bucket)}")
print(f"      H mod4 spectral (constant on cospectral)? {mod4_spectral}  ({len(split4)} buckets violate)")
print(f"      H mod8 spectral? {mod8_spectral}  ({len(split8)} buckets violate)")
# global mod4 distribution at n=7
from collections import Counter
allres=Counter()
rng2=random.Random(70)
for _ in range(20000):
    A=rand_A(n,rng2); allres[Hcount(A,n)%4]+=1
print(f"      H mod4 global distribution (20k samples): {dict(allres)}  (both 1 and 3 => nontrivial)")
print("DONE.")
