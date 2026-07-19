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
def dbl(V): return F(len({a+b for a in V for b in V}), len(V))
def q0(V):
    q=1
    while any(v%q==0 for v in V): q+=1
    return q
print("(4) THE MECHANISM: translation preserves doubling but destroys the floor")
print("    family        doubling   M         q0   (q0 is NOT translation-invariant)")
base=list(range(1,14))
for k in range(0,5):
    V=[v+k for v in base]
    print(f"    {{{1+k},...,{13+k}}}".ljust(17)
          + f"{float(dbl(V)):.3f}     {str(M(V)):9s} {q0(V):3d}")
print()
print("    => all translates share doubling 1.923 (minimal), but M and q0 vary wildly.")
print("       Any TRANSLATION-INVARIANT invariant is therefore blind to the floor,")
print("       and doubling -- the Freiman invariant -- is translation-invariant.")
print()
print("(5) SAME POINT UNDER DILATION (which LRC DOES respect, THM-1050)")
for k in [1,2,3]:
    V=[k*v for v in base]
    print(f"    {k}*{{1..13}}: doubling {float(dbl(V)):.3f}  M = {M(V)}  q0 = {q0(V)}")
print("    => dilation preserves BOTH doubling and M, consistent with THM-1050.")
print()
print("(6) THE DIVISIBILITY INVARIANT SEPARATES WHAT DOUBLING CANNOT")
for nm,V in [("{1,...,13}",base),("{2,...,14}",[v+1 for v in base]),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("AP d=3",[1+3*i for i in range(13)])]:
    print(f"    {nm:16s} doubling {float(dbl(V)):.3f}  q0 = {q0(V):3d}  M = {M(V)}")
