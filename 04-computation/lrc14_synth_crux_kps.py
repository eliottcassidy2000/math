from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations
import sys, random

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

# meas(G_C) = lonely measure of a 12-set C  (this is L_exact applied to a 12-set)
# PROVE-side claim: min over 12-subsets of {1..13} of meas(G_C); decoupling floors >= 1/143.
print("Lonely measure meas(G_C) of 12-element cores {1..13}\{e}:")
floors=[]
for e in range(1,14):
    C=sorted(set(range(1,14))-{e})
    m=L_exact(C)
    floor=Fr(6,7)*m   # decoupling limit (6/7)*meas(G_C)
    floors.append((floor,e,m))
    print(f"  drop e={e:2d}: meas(G_C)={str(m):>16s}={float(m):.6f}  (6/7)*meas={float(floor):.6f}")
mn=min(floors)
print(f"  MIN decoupling floor = {mn[0]} = {float(mn[0]):.6f} at e={mn[1]}  (claim: 1/143)")
print(f"  1/143 = {float(Fr(1,143)):.6f};  min floor == 1/143? {mn[0]==Fr(1,143)}")
sys.stdout.flush()

# DISPROVE-side attack: minimize meas(G_C) over ALL 12-sets (large entries allowed).
# Claim: min is 7/858 at {1,2,3,4,5,7,8,9,10,11,12,13}. And random large-entry starts flow back up.
print()
test12=[1,2,3,4,5,7,8,9,10,11,12,13]
print(f"  meas(G_C) of {{1..5,7..13}} = {L_exact(test12)} = {float(L_exact(test12)):.6f}; ==7/858? {L_exact(test12)==Fr(7,858)}")
# Search all 12-subsets of [1..16] for the minimum
best=None
for C in combinations(range(1,17),12):
    m=L_exact(C)
    if best is None or m<best[0]:
        best=(m,C)
print(f"  min meas(G_C) over 12-subsets of [1..16]: {best[0]}={float(best[0]):.6f} at {best[1]}")
# any tight 12-set (meas=0) in [1..16]?
tight12=[C for C in combinations(range(1,17),12) if L_exact(C)==0]
print(f"  tight 12-sets (meas=0) in [1..16]: {len(tight12)}")
sys.stdout.flush()

# Large-lcm 12-set: does meas stay above the small-lcm floor?
random.seed(7)
print()
print("  Random large-entry 12-sets meas(G_C) (should be >= small-lcm floor ~0.0082, never ->0):")
mins=[]
for _ in range(2000):
    C=sorted(random.sample(range(1,121),12))
    mins.append(L_exact(C))
print(f"    over 2000 random 12-subsets of [1..120]: min meas = {min(mins)}={float(min(mins)):.6f}")
