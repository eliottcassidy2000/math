#!/usr/bin/env python3
"""
lrc_splitting_domination_macmini_0621s6.py  (mac-mini-2026-06-21-S6)
TEST the "splitting lowers cover" / single-block domination (Route E, the SOLE OPEN LEMMA route).
Claim: for far runners, a SINGLE coherent block maximizes p0; splitting into >=2 separated blocks
LOWERS p0. If true, multi-far <= single-block (closed form), closing OPEN-Q-108.
Test: fix a bounded core B; compare p0(B u single-block) vs p0(B u split-blocks) at matched size/spread.
"""
import itertools, sys, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(11)
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
B=(0,1,2,3,4)  # bounded core
print("Single coherent block vs split blocks (core B=(0,1,2,3,4)), p0:")
print(f"{'config':<46}{'p0':>9}{'note':>14}")
# single block of 4 consecutive far runners at offset m; vs 2 blocks of 2 at m and m+G
tests=[]
for m in [15,30,60]:
    single=(m,m+1,m+2,m+3)
    for G in [10,30,100]:
        split2=(m,m+1,m+G,m+G+1)
        tests.append((m,G,single,split2))
worst_split_exceeds=0; total=0
for (m,G,single,split2) in tests:
    p_single=measS7(B+single); p_split=measS7(B+split2); total+=1
    exc = p_split>p_single
    if exc: worst_split_exceeds+=1
    print(f"m={m} block(m..m+3)            -> {float(p_single):.5f}")
    print(f"m={m},G={G} split (m,m+1)+(m+G,m+G+1)  -> {float(p_split):.5f}  {'SPLIT>SINGLE!' if exc else 'single>=split OK':>14}")
# random stress: single block of size s vs random split of same s, matched count
print("\nRandom stress (single consec block vs scattered, same far-count):")
viol=0; N=0
for _ in range(60):
    s=random.choice([3,4,5]); m=random.randint(15,40)
    block=tuple(range(m,m+s))
    scattered=tuple(sorted(random.sample(range(15,200),s)))
    pb=measS7(B+block); ps=measS7(B+scattered); N+=1
    if ps>pb+F(1,10**6): viol+=1
print(f"  scattered p0 > block p0 in {viol}/{N} (if 0: single block dominates scattered)")
print(f"\nSplit>single in {worst_split_exceeds}/{total} structured tests.")
print("If single block dominates => multi-far reduces to single-block closed form (OPEN-Q-108 route).")
