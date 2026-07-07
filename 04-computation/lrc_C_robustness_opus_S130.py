from fractions import Fraction as F
from math import gcd, lcm, ceil, floor
from functools import reduce
from itertools import combinations, product
import random
from sympy import nextprime

def dist_int(num, den):
    r = num % den
    return F(min(r, den-r), den)

def reach_M(V):
    V=[abs(v) for v in V if v!=0]; dens=set(); n=len(V)
    for i in range(n):
        dens.add(2*V[i])
        for j in range(i+1,n):
            dens.add(V[i]+V[j])
            if V[i]!=V[j]: dens.add(abs(V[i]-V[j]))
    best=F(0)
    for d in dens:
        if d==0: continue
        for m in range(1,d):
            mn=None
            for v in V:
                dd=dist_int(v*m,d)
                if mn is None or dd<mn: mn=dd
                if mn==0: break
            if mn is not None and mn>best: best=mn
    return best

def clears_at(V,q):
    lo=ceil(2*q/25); hi=floor(23*q/25)
    if lo>hi: return False
    for c in range(1,q):
        if gcd(c,q)!=1: continue
        if all(lo<=(v*c)%q<=hi for v in V): return True
    return False

lo,hi=F(1,13),F(2,25)
def in_gap(M): return lo<M<hi
random.seed(12345)

print("=== (C) ROBUSTNESS -- exact-M search for gap members (opus-S130) ===")
found=[]; tested=0
def consider(V):
    global tested
    V=[v for v in V]
    if any(v==0 for v in V) or len(set(map(abs,V)))<len(V): return
    g=reduce(gcd,[abs(v) for v in V]); 
    if g!=1: return
    M=reach_M(V); tested+=1
    if in_gap(M): found.append((tuple(V),float(M)))

AP=list(range(1,13))
# A) {1..11,x}
for x in range(13,400): consider(list(range(1,12))+[x])
# B) {1..10,x,y}
for x in range(11,60):
    for y in range(x+1,80): consider(list(range(1,11))+[x,y])
# C) dilated AP k*{1..12} with up to 2 defects
for k in [2,3,5,7]:
    base=[k*i for i in range(1,13)]
    for i in range(12):
        for d in range(-6,7):
            V=base.copy(); V[i]+=d; consider(V)
    for (i,j) in combinations(range(12),2):
        for di,dj in product([-2,-1,1,2],repeat=2):
            V=base.copy(); V[i]+=di; V[j]+=dj; consider(V)
# D) random small primitive families (broad net), max<=45
for _ in range(4000):
    V=random.sample(range(1,46),12); consider(V)
# E) random near-AP: AP with 1-4 random entries bumped by multiples of 13 (divisibility-preserving)
for _ in range(3000):
    V=AP.copy()
    for _ in range(random.randint(1,4)):
        i=random.randrange(12); V[i]=AP[i]+13*random.randint(1,8)
    consider(V)

print(f"  exact-M tested: {tested}")
print(f"  GAP MEMBERS (would refute (C)): {len(found)}")
for V,M in found[:20]: print("    ", V, M)

print("\n=== ESCAPE-FAMILY probe: do compressed varying-k families ALWAYS clear at some q? ===")
# a family escaping ALL q up to a large bound would be a gap-member candidate at huge height
maxq=120; noclear=[]
probes=0
for Q0 in [12,20,25,30,37]:
    L=lcm(*range(2,Q0+1))
    for _ in range(60):
        k=[random.randint(1,4) for _ in range(12)]
        if len(set(k))==1: k[0]+=1   # force varying (non-translate)
        V=[(i+1)+L*k[i] for i in range(12)]
        probes+=1
        # find any clearing q in [6,maxq]
        cq=next((q for q in range(6,maxq+1) if clears_at(V,q)), None)
        if cq is None: noclear.append((Q0,k))
print(f"  varying-k compressed escape families probed: {probes}")
print(f"  families NOT clearing at any q<=120 (gap-member candidates): {len(noclear)}")
for Q0,k in noclear[:10]: print("    Q0=",Q0,"k=",k)
print("\n  => if 0 gap members AND 0 non-clearing escape families: (C) strongly confirmed TRUE,")
print("     escape families all LOOSE (clear at some q, just unbounded). Only the finite-COVERING fails.")
