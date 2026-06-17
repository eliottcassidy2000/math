from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations
import sys

def danger_arcs(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A

# cache danger arcs
_cache={}
def darcs(v):
    if v not in _cache: _cache[v]=danger_arcs(v)
    return _cache[v]

def L_exact(S):
    A=[]
    for v in S: A.extend(darcs(v))
    A=sorted((a,b) for a,b in A if b>a)
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch:
            if b>ch: ch=b
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

def lcm(a,b): return a*b//gcd(a,b)
def lcmset(S): return reduce(lcm,S,1)
def prim(S): g=reduce(gcd,S); return tuple(x//g for x in S)

# Exhaustive census: all 13-subsets of [1..18] (must contain 1 for primitivity-ish; check all)
tight=set()
cnt=0
for S in combinations(range(1,19),13):
    cnt+=1
    if L_exact(S)==0:
        tight.add(prim(S))
print(f"[1..18] checked {cnt} subsets; tight primitive configs:")
for t in sorted(tight): print("  ", t)
sys.stdout.flush()
