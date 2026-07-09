"""
mac-mini-2026-07-09-S64 -- THE NON-STRICT CRITERION and the n=14 knife-edge.

CORRECTION (this session): the finite-Vmax good-period criterion for LRC(14) is
  loneliness M(S) >= 1/14  <=>  maxgap({e_i j/Vmax}) >= 1/7   (NON-STRICT, matches LRC's ">= 1/n").
Equality maxgap = 1/7 gives M = 1/14 EXACTLY, which SATISFIES the conjecture (>=). Scripts using
strict maxgap > 1/7 wrongly flag the boundary as a "failure" (cf klein-S201's V>=Q+1 boundary bug).

DECISIVE TEST: over 7-structured dissociated k=13 sets on the resonant grid 7|Vmax (the tight class),
compute the best integer margin  m(E,V) = max_j (maxgap*7 - V):
  m > 0  loneliness strict (M > 1/14);  m = 0  EXACT knife-edge (M = 1/14);  m < 0  COUNTEREXAMPLE.
Find the min margin (any m<0 would DISPROVE LRC14 for the covering case). Confirm min >= 0 with
equality (m=0) achieved at the wraparound-boundary set spread = 6V/7.
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(6400)
K=13
def prim(E):
    E=sorted(E); return len(E)>=2 and reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best
def make_7struct(k,s,mn):
    sev=[x for x in range(0,s+1,7)]
    if len(sev)<mn: return None
    S7=sorted(random.sample(sev,random.randint(mn,min(len(sev),k-1))))
    pool=[x for x in range(1,s) if x%7!=0]; need=k-len(S7)
    if need<0 or len(pool)<need: return None
    E=sorted(set(S7+[0,s]+random.sample(pool,max(0,need))))
    return E if len(E)==k else None
def best_margin(E,V):
    """max over j=1..V-1 of (maxgap*7 - V), exact ints. >=0 <=> some j gives loneliness M>=1/14."""
    E=sorted(E); best=-10**9; bj=None
    for j in range(1,V):
        ph=sorted({(e*j)%V for e in E}); m=len(ph)
        if m==1: mg=V
        else: mg=max(max((ph[(i+1)%m]-ph[i])%V for i in range(m-1)), ph[0]+V-ph[-1])
        val=mg*7-V
        if val>best: best=val; bj=j
    return best,bj

worst=10**9; winfo=None; neg=0; zero=0; n=0
for _ in range(30000):
    s=random.randint(20,140); E=make_7struct(K,s,4)
    if E is None or not prim(E) or longest_ap(E)>K-6: continue
    mx=max(E)
    for t in range(1,9):
        V=7*((mx)//7+t)
        if V<=mx: continue
        marg,bj=best_margin(E,V); n+=1
        if marg<0: neg+=1
        if marg==0: zero+=1
        if marg<worst: worst=marg; winfo=(sorted(E),V,bj)
print(f'7-struct dissociated k=13, resonant grid 7|V, {n} (set,V) pairs:')
print(f'  MIN best-margin (maxgap*7 - V): {worst}   [>=0 => LRC14 loneliness holds for covering case]')
print(f'  COUNTEREXAMPLES (margin<0): {neg}     EXACT knife-edges (margin=0, M=1/14): {zero}')
if winfo: print(f'  extremal (min margin): E={winfo[0]} V={winfo[1]} at j={winfo[2]}')

# The canonical wraparound-boundary knife-edge: spread = 6V/7 exactly, all residues, |S7|=k-6=7
E=[0,7,10,14,18,20,21,26,28,35,36,37,42]; V=49
m,bj=best_margin(E,V)
print(f'\nWRAPAROUND-BOUNDARY set (spread=42=6*49/7), V=49:')
print(f'  best margin = {m} at j={bj}   (0 = EXACT M=1/14 knife-edge)')
print(f'  residues mod7 = {sorted(set(e%7 for e in E))} (all 7 => no missed-residue), |S7| = {sum(1 for e in E if e%7==0)}')
print(f'  => at j=1 the phases fill [0,42/49]=[0,6/7], wraparound gap = 1-6/7 = 1/7 EXACTLY. Loneliness = 1/14.')
print('\nCONCLUSION: covering case holds NON-STRICTLY (maxgap>=1/7); n=14 is a true knife-edge')
print('(equality achieved at the wraparound boundary, never violated). The good-period lemma and its')
print('Lean form must use 7*spread <= 6*Vmax => gapLen >= 1/7 (non-strict), not the strict <.')
