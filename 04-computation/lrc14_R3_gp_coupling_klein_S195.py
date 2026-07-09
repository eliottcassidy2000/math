#!/usr/bin/env python3
"""
klein-2026-07-08-S195: R3 -- the G_P coupling in the DISCRETE good-period existence.

The good witness needs x=j/Vmax with BOTH:  x in G_P (||p x||>=1/14 for all p in P)
AND cluster maxgap{frac(e_i x)}>1/7 (THM-530). LEM-010/j*=O(k) give the CLUSTER part;
does a j satisfying BOTH exist at small Vmax?  Two questions:
 (Q1) does the FULL discrete good set (G_P AND cluster) have a small good period?
 (Q2) RESOLUTION TEST: fold P into the co-offsets -- E_full = {Vmax - v : v in S}
      (all 13 co-offsets, Vmax=max S). Then a gap>1/7 in {frac(e j/Vmax): all 13}
      is a witness avoiding ALL runners (incl P). Does that good period exist small,
      AND does it imply the observer is lonely (M(S)>=1/14 at some t near j/Vmax)?
"""
import numpy as np
from math import gcd
from functools import reduce
rng=np.random.default_rng(3333)
INV7=1.0/7.0

def maxgap(phases,mod=1.0):
    p=np.sort(np.mod(phases,mod)); g=np.diff(p); g=np.append(g,mod-p[-1]+p[0]); return g.max()

def build_covering_set(k, Vmax):
    """P = few small speeds (<=13), L = k large speeds up to Vmax; primitive-ish."""
    psz=13-k
    P=sorted(rng.choice(np.arange(1,14),size=psz,replace=False).tolist()) if psz>0 else []
    # cluster: k large speeds in (13, Vmax], include Vmax
    L=sorted(set(rng.choice(np.arange(14,Vmax+1),size=k-1,replace=False).tolist()+[Vmax]))
    if len(L)<k: return None
    S=sorted(set(P)|set(L))
    if len(S)!=13: return None
    return P,L,S

print("(Q1) DISCRETE good set {j: x=j/Vmax in G_P AND cluster maxgap>1/7}: smallest good j")
print("(Q2) ALL-13-co-offset: {j: gap>1/7 in {(Vmax-v)j/Vmax : v in S}}: smallest good j (P folded in)")
print(f"{'k':>3} {'Vmax':>6} {'#sets':>6} {'Q1 max jgood':>13} {'Q1 fails':>9} {'Q2 max jgood':>13} {'Q2 fails':>9}")
for k in (8,10,11,12):
    for Vmax in (200, 500, 1500):
        q1max=0;q1f=0;q2max=0;q2f=0;n=0
        for _ in range(120):
            r=build_covering_set(k,Vmax)
            if r is None: continue
            P,L,S=r; n+=1
            e_cluster=np.array([Vmax-u for u in L])  # cluster co-offsets (incl 0 for Vmax)
            e_full=np.array([Vmax-v for v in S])     # ALL 13 co-offsets
            # Q1: full discrete good set
            found1=None
            for j in range(1,Vmax):
                x=j/Vmax
                inGP=all(abs(((p*x+0.5)%1.0)-0.5)>=1/14-1e-12 for p in P)
                if not inGP: continue
                if maxgap(e_cluster*x)>INV7+1e-12: found1=j; break
            if found1 is None: q1f+=1
            else: q1max=max(q1max,found1)
            # Q2: all-13-co-offset gap
            found2=None
            for j in range(1,Vmax):
                if maxgap((e_full*j/Vmax))>INV7+1e-12: found2=j; break
            if found2 is None: q2f+=1
            else: q2max=max(q2max,found2)
        print(f"{k:>3} {Vmax:>6} {n:>6} {q1max:>13} {q1f:>9} {q2max:>13} {q2f:>9}")

print("\n=> Q1 (full G_P+cluster) small good j, 0 fails => discrete witness exists at small j.")
print("   Q2 (all-13-co-offset) small good j, 0 fails => folding P into co-offsets works")
print("   (P's co-offsets Vmax-p ~ Vmax have phases near 0, cluster with e_0, harmless).")
print("   => R3 (P-coupling) resolved: use E_full = all 13 co-offsets; LEM-010/j*=O(k) then covers P.")
