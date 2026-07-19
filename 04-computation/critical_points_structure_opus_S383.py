from fractions import Fraction as F
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def cands(V):
    D={}
    for v in V: D.setdefault(2*v,set()).add(('peak',v))
    for a,b in combinations(V,2):
        D.setdefault(a+b,set()).add(('sum',a,b))
        if a!=b: D.setdefault(abs(a-b),set()).add(('diff',a,b))
    D.pop(0,None)
    return D
def best_point(V):
    D=cands(V); best=F(0); arg=None; kinds=None
    for q,ks in D.items():
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,arg,kinds=g,F(p,q),ks
    return best,arg,kinds

print("(3) THE TIGHT FAMILIES -- does the critical-point picture reproduce g = 1/14?")
for nm,V in [("{1,...,13}",list(range(1,14))),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("2*{1,...,13}",[2*i for i in range(1,14)])]:
    g,arg,ks=best_point(V)
    print(f"    {nm:16s} g = {g} ({float(g):.8f})  at t = {arg}   1/14 = {float(F(1,14)):.8f}")
    print(f"        critical-point type(s) at that denominator: {sorted(ks)[:3]}")

print()
print("(4) WHICH KIND OF CRITICAL POINT WINS?  (sum v_i+v_j / diff v_i-v_j / peak 2v)")
random.seed(3832)
from collections import Counter
C=Counter(); slack=[]
for _ in range(40):
    V=sorted(random.sample(range(1,80),13))
    D=cands(V); best=F(0); bk=None
    nge=0; tot=0
    for q,ks in D.items():
        for p in range(1,q):
            g=gap_at(V,F(p,q)); tot+=1
            if g>=F(1,14): nge+=1
            if g>best: best,bk=g,ks
    for k in bk: C[k[0]]+=1
    slack.append((nge,tot))
print(f"    winning critical-point kind over 40 families: {dict(C)}")
med=sorted(slack)[len(slack)//2]
print(f"    candidate points achieving gap >= 1/14: median {med[0]} of {med[1]} tested")
print(f"    (so the conjecture has substantial slack per family -- many witnesses, not one)")
