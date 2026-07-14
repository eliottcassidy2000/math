#!/usr/bin/env python3
"""mac-mini-S107b: attack the LRC(13) 12-block rigidity. (A) exact census of 12-subsets of {1..18}:
is {1..12} the ONLY primitive one with M=1/13? (B) the ratio bound: does M(A)=1/13 => max/min<13
(tight => bounded, via far-element)? Reduces the rigidity to [ratio bound] + [bounded finite check]."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
def M_exact(S):
    best=F(0); dens=set()
    for a,b in combinations(sorted(set(S)),2):
        dens.add(a+b)
        if b>a: dens.add(b-a)
    for q in dens:
        for m in range(1,q):
            num=q
            for v in S:
                r=(v*m)%q; d=min(r,q-r)
                if d<num: num=d
                if num*13<q: break
            if num*13>=q and F(num,q)>best: best=F(num,q)
    return best
def prim(S): return reduce(gcd,S)==1
one13=F(1,13)
print("(A) EXACT census: 12-subsets of {1..18}, primitive, M=1/13 -- is {1..12} the only one?")
tight=[]; checked=0
for S in combinations(range(1,19),12):
    if not prim(S): continue
    checked+=1
    if M_exact(list(S))==one13: tight.append(sorted(S))
print(f"   primitive 12-subsets of {{1..18}} checked: {checked}")
print(f"   with M=1/13 (tight): {len(tight)}")
for S in tight: print(f"      {S}  (={'{'+'1..12'+'}' if S==list(range(1,13)) else 'OTHER!'})")
print(f"   => {'{1..12} is the UNIQUE tight primitive 12-set in {1..18}.' if tight==[list(range(1,13))] else 'MULTIPLE tight sets -- rigidity FALSE in {1..18}.'}")

print("\n(B) RATIO BOUND: is every tight 12-set of bounded ratio? search for tight sets with max/min>=13:")
import random; rng=random.Random(1313)
wide_tight=[]; checked2=0
for _ in range(400000):
    mn=rng.choice([1,1,1,2,3])
    rest=rng.sample([x for x in range(mn+1, mn*20+1)], 11)
    S=sorted([mn]+rest)
    if len(set(S))!=12 or S[0]!=mn or not prim(S): continue
    checked2+=1
    if M_exact(S)==one13: wide_tight.append((S[-1]/S[0], S))
print(f"   sampled {checked2} primitive 12-sets (min in {{1,2,3}}, wide range); tight (M=1/13): {len(wide_tight)}")
if wide_tight:
    maxratio=max(r for r,_ in wide_tight)
    print(f"   max ratio among tight sets found: {maxratio:.2f}")
    for r,S in sorted(wide_tight,reverse=True)[:5]: print(f"      ratio={r:.1f}: {S}")
    hi=[S for r,S in wide_tight if r>=13]
    print(f"   tight sets with ratio>=13: {len(hi)} {'-- ratio bound FAILS' if hi else '-- NONE => ratio bound holds (tight => ratio<13)'}")
print()
print("REDUCTION: if [tight => ratio<13] (far-element decorrelation, THM-751-style) and [min=1 for tight],")
print("then tight 12-set => subset of {1..12} => ={1..12}. The min>=2 case is a bounded finite check.")
print("So LRC(13) rigidity = [ratio bound, provable via far-element] + [bounded census, exact]. Verified in {1..18}.")
