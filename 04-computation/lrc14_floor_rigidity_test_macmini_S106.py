#!/usr/bin/env python3
"""mac-mini-S106: SKEPTICAL test of the S105 multi-killer floor conjecture (M=1/13 iff near-dilate).
Concern: the binding subset is a dilated block c*{1..12} (M=1/13); ANY coprime killer w safe at the
block's tight point t=1/(13c) keeps M=1/13 -- not just the antipode w=13c+1. Construct counterexamples."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
def M_exact(S):
    best=F(0); denoms=set()
    for a,b in combinations(sorted(set(S)),2):
        denoms.add(a+b)
        if b>a: denoms.add(b-a)
    for q in denoms:
        for m in range(1,q):
            num=q
            for v in S:
                r=(v*m)%q; d=min(r,q-r)
                if d<num: num=d
                if num*13<q: break
            if F(num,q)>best: best=F(num,q)
    return best
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def primitive(S): return reduce(gcd,S)==1
def is_near_dilate(S):
    S=sorted(S)
    for L in range(2,max(S)+1):
        if sorted([i*L for i in range(1,13)]+[13*L+1])==S: return True
    return False
print("The dilated covering block c*{1..12} + a coprime killer w safe at t=1/(13c): M=1/13?")
c=26  # smallest c with 13|c and (2|c or 7|c) so the block covers {2..14}
block=[c*i for i in range(1,13)]
print(f"c={c}: block={block}; block covering? {covering(block)}; M(block)=M({{1..12}})=1/13")
print(f"  block tight at t=1/(13c)=1/{13*c}; killer w is SAFE there iff ||w/{13*c}||>=1/13 iff (w mod {13*c}) in [{13*c//13},{13*c-13*c//13}]")
found=[]
for w in range(15, 13*c+50):
    if w in block: continue
    S=sorted(block+[w])
    if len(set(S))!=13 or not covering(S) or not primitive(S): continue
    M=M_exact(S)
    if M==F(1,13):
        found.append((w, is_near_dilate(S)))
print(f"\nkillers w giving EXACT M=1/13 (primitive+covering): {len(found)}")
print(f"  {[w for w,_ in found][:20]}")
nd=[w for w,isd in found if isd]; nnd=[w for w,isd in found if not isd]
print(f"  near-dilate killers (w=13c+1={13*c+1}): {nd}")
print(f"  NON-near-dilate killers with M=1/13: {nnd[:20]} ({len(nnd)} total)")
if nnd:
    w=nnd[0]; S=sorted(block+[w]); print(f"\n  COUNTEREXAMPLE: {S}")
    print(f"    M={M_exact(S)} EXACT, primitive={primitive(S)}, covering={covering(S)}, near-dilate={is_near_dilate(S)}")
    print("  => S105 CONJECTURE 'M=1/13 iff near-dilate' is FALSE. The minimizers are a whole FAMILY:")
    print("     {dilated covering 12-block c*{1..12}} + {any coprime killer w safe at 1/(13c)}. Antipode is one.")
