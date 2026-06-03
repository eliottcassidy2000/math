#!/usr/bin/env python3
"""lrc_single_swap_classification_s612b.py -- complete classification of tight
AP single-swaps (exploration alongside HYP-2177).

RESULT [VERIFIED n=4..26]: the tight single-swaps AP[a->a'] (a in AP={1..n-1},
a' not in AP) are EXACTLY:
  * FAMILY II (doubling V*): a = n-2, a' = 2(n-2), tight  <=>  3 | (2n-1)
    -- the infinite family n in {8,14,20,26,...} (n == 2 mod 3); PROVED HYP-2177.
  * FAMILY I (reflection): a = 2, a' = 2n-3 (== -2 mod 2n-1), tight ONLY at
    n = 5 and n = 6 -- exactly the classic sporadics (1,3,4,7) [n=5] and
    (1,3,4,5,9) [n=6]. Exceptional (small-n), does NOT continue.
  * NOTHING ELSE.
So the famous sporadics (1,3,4,7),(1,3,4,5,9) are AP reflection-swaps a=2->-2,
and V* is the AP doubling-swap a=n-2->2a; both are single swaps governed by the
residue of the swapped element mod 2n-1. (Multi-swap sporadics, e.g. n=8s
(1,4,5,6,7,11,13), exist separately.)

Session: claude-2026-06-03-S612b (lrc-single-swap-classification).
"""
import sys; sys.stdout.reconfigure(line_buffering=True)
from math import gcd
from functools import reduce
def lcm(a,b): return a*b//gcd(a,b)
def is_tight(V):
    n=len(V); L=reduce(lcm,V); D=(n+1)*L; ivs=[]
    for v in V:
        Lv=L//v
        for j in range(v+1):
            lo=((j*(n+1)-1)*Lv)%D; hi=lo+2*Lv
            if hi<=D: ivs.append((lo,hi))
            else: ivs.append((lo,D)); ivs.append((0,hi-D))
    ivs.sort(); pos=0
    for s,e in ivs:
        if s>pos: return False
        if e>pos: pos=e
    return pos>=D
fam1=[]; fam2=[]; other=[]
for n in range(4,27):
    AP=set(range(1,n)); Q=2*n-1
    for a in range(1,n):
        for ap in range(n,5*n+1):
            if ap in AP: continue
            V=tuple(sorted((AP-{a})|{ap}))
            if reduce(gcd,V)==1 and is_tight(V):
                if ap==2*a: fam2.append((n,a,ap))
                elif ap%Q==(Q-a)%Q: fam1.append((n,a,ap))
                else: other.append((n,a,ap))
print("Complete classification of tight AP single-swaps, n=4..26 (a' up to 5n):")
print(f"  Family II (doubling a=n-2 -> 2(n-2), needs 3|2n-1): {[(n,a,ap) for n,a,ap in fam2]}")
print(f"    all have a=n-2? {all(a==n-2 for n,a,ap in fam2)}; all 3|2n-1? {all((2*n-1)%3==0 for n,a,ap in fam2)}")
print(f"  Family I (reflection a=2 -> 2n-3 == -2): {[(n,a,ap) for n,a,ap in fam1]}")
print(f"  Other (neither): {other}")
