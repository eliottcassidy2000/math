#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S16 (HYP-4432) -- the WITNESS-DENOMINATOR LEMMA (verification).

LEMMA. If M(S)=c/q in lowest terms, then q divides (v_i+v_j) or (v_i-v_j) for some
pair, or 2*v_i for some i.
PROOF. f(t)=min_i||v_i t|| is piecewise-linear; its max is at a crossing
||v_i t||=||v_j t|| (=> (v_i-+v_j)t in Z => t=k/(v_i-+v_j)) or a single-runner peak
t=(2k+1)/(2 v_i).  Either way q | v_i+-v_j or 2 v_i.
CONSEQUENCE: q <= 2*max(v_i) -- bounding height bounds the witness denominator =>
(G) is a finite check (exact M computable).  The additive realization of kps's
'bound the off-13-grid denominators' guidance.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys
sys.path.insert(0,'04-computation')
from lonely_profile import profile

def M_exact(S):
    S=sorted(set(abs(x) for x in S))
    for cap in (14,11,8,6,4,3,2):
        p=profile(S,F(1,cap)); m=p.M()
        if m is not None: return m

print("CLAIM: witness denom q (of M=c/q) DIVIDES some (v_i+v_j),(v_i-v_j), or 2*v_i.\n")
fams={
 "AP {1..12}": list(range(1,13)),
 "doubled-apex {1..11,24}": list(range(1,12))+[24],
 "block {..17,19}": [1,2,3,5,7,8,9,10,11,12,17,19],
 "deep well {1..11,168}": list(range(1,12))+[168],
 "{1..11,23}": list(range(1,12))+[23],
 "{2..13}": list(range(2,14)),
 "{1..10,12,25}": list(range(1,11))+[12,25],
}
for name,S in fams.items():
    M=M_exact(S); q=M.denominator
    hits=[]
    for i in range(len(S)):
        if 2*S[i]%q==0: hits.append(f"2*{S[i]}")
    for a,b in combinations(S,2):
        if (a+b)%q==0: hits.append(f"{a}+{b}")
        if (a-b)%q==0: hits.append(f"{a}-{b}")
    print(f"{name:26s} M={str(M):>7} q={q:>4}: q|sum/diff/2v? {'YES' if hits else 'NO'}   e.g. {hits[:4]}")
print("\n=> q | pairwise sum/difference.  For a gap member (q>=38) a pair must sum to a")
print("   multiple of q>=38 => a LIFTED runner is forced; the covering-multiplicative")
print("   structure constrains q (parity of the sum, 13|q via a pair summing to a mult of 13).")
