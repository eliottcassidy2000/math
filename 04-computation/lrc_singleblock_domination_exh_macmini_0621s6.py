#!/usr/bin/env python3
"""Exhaustive-ish single-block domination: for many bounded cores B and far sizes, does the
single consec block maximize p0 over all same-size far-shapes (at matched offset window)?"""
import itertools, sys, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(3)
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot
# For far window [m, m+W], among all s-subsets, is the consec block max p0? (small W exhaustive)
viol=0; tot=0; worst=None
for B in [(0,),(0,1),(0,1,2),(0,1,2,3),(0,1,2,3,4),(0,2,4),(0,1,3)]:
    for s in [2,3]:
        for m in [15,25,40]:
            W=s+3
            block=tuple(range(m,m+s))  # consec block
            pblock=measS7(B+block)
            for far in itertools.combinations(range(m,m+W+1), s):
                if far==block: continue
                pf=measS7(B+far); tot+=1
                if pf>pblock+F(1,10**7):
                    viol+=1
                    if worst is None or pf-pblock>worst[0]: worst=(pf-pblock,B,block,far)
print(f"single consec block NOT max in {viol}/{tot} (block dominated by a same-size shape in window)")
if worst: print(f"  worst beater: B={worst[1]} block={worst[2]} beaten by {worst[3]} (by {float(worst[0]):.5f})")
else: print("  => consec block is the per-window same-size MAXIMIZER (0 violations) -- domination holds")
