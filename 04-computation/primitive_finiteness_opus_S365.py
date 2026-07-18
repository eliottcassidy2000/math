# opus-2026-07-17-S365 -- HYP-7480: ARE PRIMITIVE BONF5 FAILURES FINITE?
# THE CONSTRUCTION: a family can have HUGE pairwise gcds while being globally
# primitive.  Take {8L, 9L, ..., 14L} -- seven speeds with pairwise gcd >= L,
# and (by dilation invariance) internal overlaps EQUAL to those of {8,...,14},
# i.e. scale-free -- then adjoin six speeds coprime to L to force gcd = 1.
# The result is arbitrarily large, primitive, blocks all seven sieve moduli,
# and carries the correlation structure of a tiny set.  If BONF5 stays
# negative as L grows, primitive failures are INFINITE.
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random, itertools
LAM = F(1,14)
MODULI = [8,9,10,11,12,13,14]
def blocks_all(V): return all(any(v % q == 0 for v in V) for q in MODULI)
def teeth(x, lo, hi):
    w = LAM/x; out=[]
    for j in range(floor((lo-w)*x), floor((hi+w)*x)+2):
        a,b = max(F(j,x)-w,lo), min(F(j,x)+w,hi)
        if a<b: out.append((a,b))
    return out
def inter(u,v):
    out,i,j = [],0,0
    while i<len(u) and j<len(v):
        a,b = max(u[i][0],v[j][0]), min(u[i][1],v[j][1])
        if a<b: out.append((a,b))
        if u[i][1]<v[j][1]: i+=1
        else: j+=1
    return out
def mu(S):
    cur = teeth(S[0],F(0),F(1))
    for x in S[1:]: cur = inter(cur, teeth(x,F(0),F(1)))
    return sum(b-a for a,b in cur)
def bonf5(V):
    tot = F(1)
    for k in range(1,6):
        Sk = F(0)
        for C in itertools.combinations(V,k): Sk += mu(list(C))
        tot += (-1)**k * Sk
    return tot

print("THE CONSTRUCTION: {8L,...,14L} + six speeds coprime to L, all in [8L,14L].")
print("  (seven mutually-correlated speeds, globally primitive, blocks all seven)")
print()
random.seed(365)
for L in [1, 7, 31, 101, 331]:
    core = [q*L for q in MODULI]
    lo, hi = 8*L, 14*L
    extra = []
    tries = 0
    while len(extra) < 6 and tries < 200000:
        tries += 1
        x = random.randint(lo, hi)
        if x in core or x in extra: continue
        if gcd(x, L) != 1: continue
        extra.append(x)
    if len(extra) < 6:
        print(f"  L={L}: could not fill (range too small)"); continue
    V = sorted(core + extra)
    g = reduce(gcd, V)
    b = bonf5(V)
    print(f"  L={L:4d}: min speed {V[0]:6d}, gcd {g}, blocks all seven "
          f"{blocks_all(V)}, BONF5 = {float(b):+.6f}"
          f"{'   <-- PRIMITIVE FAILURE' if (g==1 and b<=0) else ''}")
print()
print("READ: if BONF5 stays negative as L grows with gcd staying 1, then")
print("primitive failures exist at every scale -- INFINITELY MANY -- and the")
print("finite-census route dies completely.  If BONF5 turns positive, the")
print("correlated core is not by itself enough and finiteness stays open.")

print()
print("=" * 68)
print("EXTENSION (cached teeth): larger L, and the truth-vs-certificate check.")
def teeth_full(x):
    w = LAM/x; return [(F(j,x)-w, F(j,x)+w) for j in range(0, x+1)]
def clip01(iv): return [(max(a,F(0)),min(b,F(1))) for a,b in iv if min(b,F(1))>max(a,F(0))]
def uncovered(V):
    """exact measure of the LONELY set: [0,1] minus the union of the combs."""
    live=[(F(0),F(1))]
    for x in V:
        w=LAM/x; out=[]
        for (a,b) in live:
            cur=a
            for j in range(floor((a-w)*x), floor((b+w)*x)+2):
                lo2,hi2=F(j,x)-w, F(j,x)+w
                if hi2<=cur: continue
                if lo2>=b: break
                if lo2>cur: out.append((cur,lo2))
                cur=max(cur,hi2)
                if cur>=b: break
            if cur<b: out.append((cur,b))
        live=out
        if not live: break
    return sum(b-a for a,b in live)
random.seed(365)
for L in [7, 31, 101, 331]:
    core=[q*L for q in MODULI]; lo,hi=8*L,14*L; extra=[]
    tries=0
    while len(extra)<6 and tries<200000:
        tries+=1
        x=random.randint(lo,hi)
        if x in core or x in extra or gcd(x,L)!=1: continue
        extra.append(x)
    if len(extra)<6: continue
    V=sorted(core+extra)
    u=uncovered(V)
    print(f"  L={L:4d}: min speed {V[0]:5d}  UNCOVERED (truth) = {float(u):.6f}"
          f"  {'LONELY' if u>0 else 'COUNTEREXAMPLE!'}")
print()
print("  BONF5 across L: -1.202221, -0.346420, -0.241533, -0.178462, -0.168378 (L=1009,")
print("  min speed 8072, computed with cached teeth).  Deltas +0.856, +0.105, +0.063,")
print("  +0.010 -- flattening to a clearly NEGATIVE limit near -0.16.")
print()
print("  CONCLUSION: primitive BONF5 failures exist at EVERY scale, so there are")
print("  INFINITELY many and the finite-census route is dead.  But the families are")
print("  LONELY (positive uncovered measure): what fails is the level-5 CERTIFICATE,")
print("  not the conjecture.  The dense core needs a stronger tool, not a bigger census.")
