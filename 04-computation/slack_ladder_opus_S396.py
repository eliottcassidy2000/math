# opus-2026-07-19-S396 -- HYP-7820: THE SLACK LADDER, AND WHY SLACK 1 IS OPTIMAL.
#
# slack := 14D - s at the active pair, an INTEGER (D, s integers).  Then
#     slack >= 0  <=>  LRC(14) holds for the family
#     slack =  0  <=>  M = D/s = D/(14D) = 1/14  <=>  EXTREMAL
#     slack =  1  <=>  M = D/(14D - 1)   -- the tightest NON-extremal value
# So "beating slack 1" is impossible except by being extremal: 1 is the minimum
# over non-extremal families, and 3/41 (D=3) already attains it.
#
# THE REAL QUESTION IS THE SLACK-1 LADDER  M = D/(14D-1):
#     D=1 -> 1/13,  D=2 -> 2/27,  D=3 -> 3/41,  D=4 -> 4/55, ... -> 1/14.
# If arbitrarily large D are REALISED by 13-families, then 1/14 is an
# ACCUMULATION POINT of non-extremal values from above, and there is no gap
# above the floor at all -- a sharp structural fact about LRC(14).
# NOTE slack is NOT dilation-invariant: D->kD, s->ks gives slack->k*slack, so
# slack is meaningful only for PRIMITIVE families (dilates inflate it).
from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def Mfull(V):
    Dn=set()
    for a,b in combinations(V,2):
        Dn.add(a+b)
        if a!=b: Dn.add(abs(a-b))
    for v in V: Dn.add(2*v)
    best=F(0); arg=None
    for q in sorted(Dn):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,arg=g,F(p,q)
    return best,arg
def act(V,t,g):
    up=[];dn=[]
    for v in V:
        x=v*t; a=int(x); r=x-a
        if r<0: a-=1; r+=1
        if r==g: up.append((v,a))
        if 1-r==g: dn.append((v,a+1))
    return up,dn
def info(V):
    M,t=Mfull(V); up,dn=act(V,t,M)
    if not up or not dn: return M,None,None,None
    vi,ai=up[0]; vj,aj=dn[0]; D=vi*aj-vj*ai; s=vi+vj
    return M,D,s,14*D-s

print("(1) slack = 0 IFF extremal (M = 1/14); slack SCALES under dilation")
for nm,V in [("{1,...,13}",list(range(1,14))),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("2*{1,...,13}",[2*i for i in range(1,14)]),
             ("3*{1,...,13}",[3*i for i in range(1,14)]),
             ("{1..11,13,36} (3/41)",[1,2,3,4,5,6,7,8,9,10,11,13,36]),
             ("2*{1..11,13,36}",[2*x for x in [1,2,3,4,5,6,7,8,9,10,11,13,36]])]:
    M,D,s,sl=info(V)
    print(f"    {nm:24s} M={str(M):9s} D={str(D):4s} s={str(s):5s} slack={str(sl):4s}"
          f"  {'EXTREMAL' if M==F(1,14) else ''}")

print()
print("(2) THE SLACK-1 LADDER  M = D/(14D-1) -- which rungs are REALISED?")
print("    D    M = D/(14D-1)    value      realised by")
targets={D:F(D,14*D-1) for D in range(1,9)}
found={}
# known witnesses + a targeted family sweep
cands=[list(range(1,13))+[14], [1,2,3,4,5,6,7,8,9,10,11,13,36],
       list(range(1,13))+[15], [1,2,3,4,5,6,7,8,9,10,11,13,48],
       [1,2,3,4,5,6,7,8,9,10,11,13,60], list(range(1,14))]
random.seed(396)
for _ in range(220):
    cands.append(sorted(random.sample(range(1,50),13)))
for V in cands:
    if len(set(V))!=13: continue
    M,D,s,sl=info(V)
    for d,tv in targets.items():
        if M==tv and d not in found: found[d]=(list(V),D,s,sl)
for d in sorted(targets):
    if d in found:
        V,D,s,sl=found[d]
        print(f"    {d}  {str(targets[d]):9s} {float(targets[d]):.6f}   {V}  (D={D}, s={s}, slack={sl})")
    else:
        print(f"    {d}  {str(targets[d]):9s} {float(targets[d]):.6f}   -- not found in this sweep")
print()
print(f"    ladder limit: D/(14D-1) -> 1/14 = {float(F(1,14)):.8f} as D -> infinity")
print(f"    D=8 already gives {float(F(8,111)):.8f}, within {float(F(8,111)-F(1,14)):.2e} of the floor")
print()
print("(3) CAN SLACK BE 0 WITH D >= 2 (a 'new' extremal)?  slack=0 => M=1/14, so any")
print("    such family is extremal -- and dilation is how D>=2 extremals arise:")
for k in [1,2,3]:
    V=[k*i for i in range(1,14)]
    M,D,s,sl=info(V)
    print(f"      {k}*{{1..13}}: M={M} D={D} s={s} slack={sl}  (primitive: {reduce(gcd,V)==1})")
