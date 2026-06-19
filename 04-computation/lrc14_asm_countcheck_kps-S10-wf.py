#!/usr/bin/env python3
"""Determine the enumeration that yields the claimed primitive-shape counts
11432/6435/715 for k=8,9,10 at B=16,15,13.  Two readings of 'span':
  (a) span = max element = B  (E={0}u(k-1 subset of {1..B}), max==B not required)
  (b) span = max element <= B
  (c) primitivity = gcd(E)==1
Print counts under each so we can pin the exact enumeration the claim audits.
"""
from itertools import combinations
from math import gcd
from functools import reduce

def prim(E):
    return reduce(gcd, [e for e in E if e>0])==1

for k,B in [(8,16),(9,15),(10,13)]:
    # (b) max element <= B, all subsets, primitive
    cnt_le=0
    cnt_le_all=0
    for rest in combinations(range(1,B+1), k-1):
        cnt_le_all+=1
        E=[0]+list(rest)
        if prim(E): cnt_le+=1
    # (a) max element == B exactly (span exactly B), primitive
    cnt_eq=0
    for rest in combinations(range(1,B+1), k-1):
        if rest[-1]!=B: continue
        E=[0]+list(rest)
        if prim(E): cnt_eq+=1
    # non-primitive too, max<=B
    print(f"k={k} B={B}: maxElt<=B all={cnt_le_all}  primitive={cnt_le}  "
          f"| span==B primitive={cnt_eq}  | claim={ {8:11432,9:6435,10:715}[k] }")
