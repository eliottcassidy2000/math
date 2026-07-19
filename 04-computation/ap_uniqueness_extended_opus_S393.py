from fractions import Fraction as F
from itertools import combinations
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
print("(6) EXTEND THE RANGE: is n=14 still the only exceptional case?")
print("     n   ties (single substitution, r <= 3n)")
for n in [16,17,18]:
    base=list(range(1,n)); fl=F(1,n); ties=[]
    for i in base:
        for r in range(n,3*n+1):
            V=sorted([x for x in base if x!=i]+[r])
            if len(set(V))!=len(base): continue
            if M(V)==fl: ties.append((i,r))
    print(f"    {n:3d}   {len(ties)} {ties[:3]}")
print()
print("(7) WHY 12 -> 24 AT n=14: the essential-region criterion (THM-1125)")
print("    the swap i -> r preserves the floor iff E_i is contained in D_r,")
print("    where E_i is the part of the circle only speed i keeps clear.")
print("    Check which i admit ANY r at each n (r <= 3n):")
for n in [12,13,14,15]:
    base=list(range(1,n)); fl=F(1,n); who=set()
    for i in base:
        for r in range(n,3*n+1):
            V=sorted([x for x in base if x!=i]+[r])
            if len(set(V))!=len(base): continue
            if M(V)==fl: who.add(i)
    print(f"    n={n:3d}: swappable speeds = {sorted(who) if who else 'NONE (AP unique)'}")
