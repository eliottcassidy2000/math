from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def M(V):
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b)
        if a!=b: D.add(abs(a-b))
    for v in V: D.add(2*v)
    best=F(0)
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best=g
    return best
print("(4) AP UNIQUENESS ACROSS n: is {1,...,n-1} the ONLY set with M = 1/n?")
print("    search: single substitutions of {1,...,n-1}, replacement r <= 60")
print("     n  speeds  floor   AP M     other families at the floor")
for n in range(10,16):
    base=list(range(1,n))
    fl=F(1,n)
    mb=M(base)
    ties=[]
    for i in base:
        for r in range(n,61):
            V=sorted([x for x in base if x!=i]+[r])
            if len(set(V))!=len(base): continue
            if M(V)==fl: ties.append((i,r,V))
    par = "even" if n%2==0 else "odd "
    print(f"    {n:3d} {len(base):5d}   1/{n:<3d}  {str(mb):7s}  "
          f"{len(ties)} tie(s) {par}  {[(i,r) for i,r,_ in ties][:3]}")
print()
print("(5) THE TIES IN FULL (single substitutions achieving the floor exactly)")
for n in range(10,16):
    base=list(range(1,n)); fl=F(1,n)
    for i in base:
        for r in range(n,61):
            V=sorted([x for x in base if x!=i]+[r])
            if len(set(V))!=len(base): continue
            if M(V)==fl:
                print(f"    n={n:3d}: {i} -> {r}   V = {V}")
