from fractions import Fraction as F
from math import gcd
from itertools import combinations
def fd(x):
    r=x-(x.numerator//x.denominator); return min(r,1-r)
def fmin(sp,t): return min(fd(F(v)*t) for v in sp)
def covering(sp): return all(any(v%q==0 for v in sp) for q in range(2,13))
def primitive(sp):
    g=0
    for v in sp: g=gcd(g,v)
    return g==1
def M_pinch(sp):
    sums={sp[i]+sp[j] for i in range(12) for j in range(i+1,12)}
    best=F(0)
    for q in sums:
        for m in range(1,q):
            v=fmin(sp,F(m,q))
            if v>best: best=v
            if best>F(1,13): return best  # early exit: not tight
    return best
THIRT=F(1,13); N=20; tight=[]; ncov=0
for combo in combinations(range(1,N+1),12):
    if not primitive(combo) or not covering(combo): continue
    ncov+=1
    M=M_pinch(list(combo))
    if M==THIRT: tight.append(combo)
print(f"N={N}: primitive covering 12-subsets={ncov}; with M==1/13: {len(tight)}")
for t in tight: print("  ",t," == {1..12}:",t==tuple(range(1,13)))
print("DONE")
