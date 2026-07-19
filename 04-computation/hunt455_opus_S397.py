from fractions import Fraction as F
from itertools import combinations
import random
TARGET=F(4,55); LO=F(1,14); HI=F(3,41)
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def beats(V, thr):
    """True if SOME critical point gives gap > thr -- early exit, no full M."""
    Dn=set()
    for a,b in combinations(V,2):
        Dn.add(a+b)
        if a!=b: Dn.add(abs(a-b))
    for v in V: Dn.add(2*v)
    for q in sorted(Dn):
        if q<2: continue
        for p in range(1,q):
            t=F(p,q); g=min(fn(t*v) for v in V)
            if g>thr: return True
    return False
def M_in_interval(V):
    """cheap test: is M strictly inside (1/14, 3/41)?"""
    if beats(V, HI - F(1,10**9)): return None      # M >= 3/41, reject fast
    # M <= 3/41; now confirm M > 1/14 by finding a point above the floor
    Dn=set()
    for a,b in combinations(V,2):
        Dn.add(a+b)
        if a!=b: Dn.add(abs(a-b))
    for v in V: Dn.add(2*v)
    best=F(0)
    for q in sorted(Dn):
        if q<2: continue
        for p in range(1,q):
            g=min(fn(F(p,q)*v) for v in V)
            if g>best: best=g
    return best if LO<best<HI else None

print("HUNTING the mediant 4/55 and anything else inside (1/14, 3/41)")
print("  Farey: |1*41 - 3*14| = 1, so 1/14 and 3/41 are NEIGHBOURS and")
print("  4/55 = (1+3)/(14+41) is the unique least-denominator fraction inside.")
print()
random.seed(3971)
band=[u for u in range(4,52)]          # residues mod 55 at distance >= 4
found=[]; tried=0
# (a) residue-band construction at t* = p/55
for _ in range(500):
    p=random.choice([x for x in range(1,55) if __import__('math').gcd(x,55)==1])
    inv=pow(p,-1,55)
    V=set()
    while len(V)<13:
        u=random.choice(band)
        v=(inv*u)%55 + 55*random.randint(0,3)
        if v>0: V.add(v)
    V=sorted(V)
    if len(V)!=13: continue
    tried+=1
    m=M_in_interval(V)
    if m: found.append((m,V))
# (b) near-extremal shapes with two free slots
for x in range(14,140):
    for y in range(x+1,160):
        V=sorted(set(list(range(1,12))+[x,y]))
        if len(V)!=13: continue
        tried+=1
        m=M_in_interval(V)
        if m: found.append((m,V))
print(f"  {tried} families tested (residue-band at s=55, plus {{1..11,x,y}} shapes)")
print(f"  families with M strictly inside (1/14, 3/41): {len(found)}")
for m,V in sorted(found)[:6]:
    print(f"    M = {m} = {float(m):.8f}   {V}")
if not found:
    print("  NONE -- including no realisation of the mediant 4/55")
