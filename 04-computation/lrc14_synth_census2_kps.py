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

# Targeted: 12-core {1..11,13} fixed, vary W in 14..3000, record L and lcm
core=[1,2,3,4,5,6,7,8,9,10,11,13]
print("Family {1..11,13,W}: W where L=0, and L trend vs lcm")
tightW=[]
samples={}
for W in range(14,3001):
    S=core+[W]
    L=L_exact(S)
    if L==0: tightW.append(W)
    if W in (24,36,72,144,150,300,600,1200,2400,3000):
        samples[W]=(L, lcmset(S))
print("  tight W (L=0):", tightW)
for W in sorted(samples):
    L,lc=samples[W]
    print(f"  W={W:5d}: L={str(L):>16s} = {float(L):.6f}  lcm={lc}")
sys.stdout.flush()

# k=1 replacement census of AP: drop one of {1..13}, add w<=300
print()
print("k=1 replacements of AP (drop e in 1..13, add w in 14..300): tight configs")
APset=set(range(1,14))
found=set()
for e in range(1,14):
    base=sorted(APset-{e})
    for w in range(14,301):
        if w in APset: continue
        S=base+[w]
        if L_exact(S)==0:
            found.add(prim(S))
for t in sorted(found): print("  ", t, " lcm=", lcmset(t))
sys.stdout.flush()
