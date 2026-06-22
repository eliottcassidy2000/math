#!/usr/bin/env python3
"""CHECK 4 fixed: cap-bound counterexample search. kind-mendel-2026-06-22-S1."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
def lcm(a,b): return a*b//gcd(a,b)
def gl(xs): return reduce(lcm,[x for x in xs if x],1)
def gall(xs): return reduce(gcd,[x for x in xs if x],0)   # gcd of nonzero

def p0_fast(E):
    pos=[e for e in E if e!=0]
    if not pos: return F(0)
    D=7*gl(pos)
    bset=set([0,D])
    for e in pos:
        step=D//(7*e); x=0
        while x<=D: bset.add(x); x+=step
    bps=sorted(bset); num=0
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid2=a+b; sec=set()
        for e in E:
            sec.add((7*((e*mid2)%(2*D)))//(2*D))
            if len(sec)==7: break
        if len(sec)==7: num+=b-a
    return F(num,D)

capvals={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
print("=== CHECK 4 (fixed): max_E p0(E) over primitive E=0+(k-1)-subset of {1..S} ===")
for k,S in [(8,16),(9,16),(10,16),(8,20)]:
    cap=capvals[k]; best=(F(-1),None); over=0; cnt=0
    for sub in combinations(range(1,S+1),k-1):
        E=[0]+list(sub)
        if gall(E)!=1: continue
        cnt+=1; v=p0_fast(E)
        if v>best[0]: best=(v,E)
        if v>=cap: over+=1
    print(f"k={k} spread<={S}: prim sets={cnt:6d}  MAX p0={float(best[0]):.5f} at {best[1]}  "
          f"cap={float(cap):.5f}  margin={float(cap-best[0]):+.5f}  #(p0>=cap)={over}  "
          f"argmax_is_consec:{best[1]==list(range(k))}")
