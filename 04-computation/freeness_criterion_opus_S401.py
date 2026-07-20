from fractions import Fraction as F
from itertools import combinations
from math import gcd
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def maximizers(V):
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b)
        if a!=b: D.add(abs(a-b))
    for v in V: D.add(2*v); D.add(2)
    best=F(0); args=[]
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,args=g,[F(p,q)]
            elif g==best and F(p,q) not in args: args.append(F(p,q))
    return best, sorted(set(args))
print("(3) WHEN IS THE Z/2 ACTION FREE?  t=1/2 is the unique fixed point of t -> 1-t.")
print("    At t=1/2: ||v/2|| = 0 if v EVEN, 1/2 if v ODD.  So t=1/2 is a maximizer")
print("    exactly when ALL speeds are odd (then M = 1/2, the global maximum).")
print()
print("    family                     all odd?  M          1/2 a maximizer?  free?")
for nm,V in [("{1,...,13}",list(range(1,14))),
             ("{1,3,5,...,25} all odd",[2*i+1 for i in range(13)]),
             ("{1,3,...,23,26}",[2*i+1 for i in range(12)]+[26]),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("{1..11,13,36}",[1,2,3,4,5,6,7,8,9,10,11,13,36])]:
    allodd=all(v%2==1 for v in V)
    M,args=maximizers(V)
    half = F(1,2) in args
    free = not half
    print(f"    {nm:26s} {str(allodd):8s} {str(M):10s} {str(half):16s} {free}")
print()
print("  => the free-Z/2 hypothesis of Borsuk-Ulam holds EXACTLY when the family")
print("     contains an even speed.  All-odd families have the reflection-FIXED")
print("     witness t = 1/2 with M = 1/2 -- Brouwer's regime, not Borsuk-Ulam's.")
print("     (And all-odd families are trivial for LRC: 1/2 >> 1/14.)")
