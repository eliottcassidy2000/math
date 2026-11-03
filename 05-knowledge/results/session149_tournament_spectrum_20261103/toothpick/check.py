#!/usr/bin/env python3
"""Exact integer-coordinate check; see CORRECTION.md for conventions."""
from collections import Counter
from itertools import combinations, permutations

def require(condition, message):
    if not condition:
        raise AssertionError(message)

def seg(a,b): return tuple(sorted((a,b)))
def rot(s): return seg(*[(-x,-y) for x,y in s])
def quarter(s): return seg(*[(-y,x) for x,y in s])
def refl(s): return seg(*[(-x,y) for x,y in s])
def mid(s): return ((s[0][0]+s[1][0])//2, (s[0][1]+s[1][1])//2)
def run(N):
    seed=seg((0,-1),(0,1)); allseg={seed}; latest={seed}; totals=[]
    for t in range(1,N+1):
        totals.append(len(allseg))
        require({rot(s) for s in allseg} == allseg, "half-turn invariance")
        require({refl(s) for s in allseg} == allseg, "reflection invariance")
        require([s for s in allseg if rot(s)==s] == [seed], "half-turn fixed segments")
        require(sum(refl(s)==s for s in allseg)%2 == 1, "reflection fixed count")
        require(quarter(seed) not in allseg, "quarter-turn seed")
        endpoints=Counter(p for s in allseg for p in s)
        middles={mid(s) for s in allseg}
        fresh=set()
        for s in latest:
            for p in s:
                if endpoints[p]!=1 or p in middles: continue
                x,y=p
                if s[0][0]==s[1][0]: z=seg((x-1,y),(x+1,y))
                else: z=seg((x,y-1),(x,y+1))
                require(z not in allseg, "duplicate segment")
                fresh.add(z)
        allseg.update(fresh); latest=fresh
    return totals

q={(0,1),(0,2),(1,2),(1,3),(2,3),(3,0)}
p=(3,2,1,0)
require(all((p[b],p[a]) in q for a,b in q), "Q4 converse")
require({(0,1),(1,3),(3,0)} <= q, "Q4 triangle")
actual=run(17)
expected=[1,3,7,11,15,23,35,43,47,55,67,79,95,123,155,171,175]
require(actual == expected, f"stage totals: {actual}")
print('Q4 reversal and directed triangle: PASS')
print('17 stages: totals match recorded A139250 values; half-turn has exactly one fixed segment; reflections have odd fixed count; quarter-turn image of seed absent: PASS')
print('totals:', actual)
