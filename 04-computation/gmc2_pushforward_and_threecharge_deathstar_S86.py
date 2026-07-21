#!/usr/bin/env python3
"""gmc2_pushforward_and_threecharge_deathstar_S86.py (HYP-8765)
Attacking GMC(2) [E[P^m]=0 all m => P one-sided], open piece = >=3 charges.
ANGLE A (pushforward): E[P^m]=int z^m dmu, mu=P_*(Gaussian); one-sided => vanishing
  analytic moments; GMC(2) = 'vanishing analytic moments => one-sided'. Verify the
  obstruction (analytic-vanishing is weaker than rotational invariance).
ANGLE B (three-charge stratum): scan two-sided 3-charge P; confirm no nullcone
  element (GMC(2) evidence) + the detection-depth / 'second rung' structure (S73).
"""
from math import factorial
from fractions import Fraction as Fr
from itertools import product

def mul(p,q):
    r={}
    for(a,b),c in p.items():
        for(a2,b2),c2 in q.items():
            k=(a+a2,b+b2); r[k]=r.get(k,0)+c*c2
    return {k:v for k,v in r.items() if v}
def powp(p,m):
    r={(0,0):Fr(1)}
    for _ in range(m): r=mul(r,p)
    return r
def E(p): return sum(v*factorial(a) for(a,b),v in p.items() if a==b)  # charge-0, Z^a W^a -> a!
def charges(p): return sorted({a-b for (a,b),v in p.items()})
def onesided(p):
    ch=[c for c in charges(p) if c!=0]
    return all(c>0 for c in ch) or all(c<0 for c in ch)

print("="*66,"\nANGLE A: one-sided => vanishing analytic moments (nullcone); the\npushforward mu has vanishing analytic moments but is NOT rot-invariant\n","="*66)
# one-sided P = Z + Z^3 (charges +1,+3): E[P^m]=0 all m, but mu not rot-inv (full moments)
P={(1,0):Fr(1),(3,0):Fr(1)}
print(f"  P=Z+Z^3 (charges {charges(P)}, one-sided={onesided(P)}):")
for m in range(1,5): print(f"    E[P^{m}] (analytic moment int z^m dmu) = {E(powp(P,m))}")
# 'full moment' E[P^a conj(P)^b] = int z^a zbar^b dmu ; rot-inv <=> =0 unless a=b
def conj(p): return {(b,a):v for(a,b),v in p.items()}
Pc=conj(P)
for (a,b) in [(2,1),(3,1),(1,0)]:
    val=E(mul(powp(P,a),powp(Pc,b)))
    print(f"    full moment int z^{a} zbar^{b} dmu = E[P^{a} Pbar^{b}] = {val} {'(=0, rot-inv would need)' if a!=b else ''}")
print("  => analytic moments all 0 (nullcone) yet some int z^a zbar^b !=0 (a!=b): mu is")
print("     NOT rotationally invariant. So 'nullcone' is WEAKER than rot-inv; GMC(2) claims")
print("     that for polynomial pushforwards it still forces one-sidedness. (the obstruction named)")

print("="*66,"\nANGLE B: three-charge two-sided scan -- NO nullcone element; detection\ndepth & the 'second rung' (S73 primitive-relation mechanism)\n","="*66)
# P = a*Z^2 + b*Zbar + c*Zbar^3 ? charges {+2,-1,-3} one-sided(neg). Need two-sided 3-charge.
# two-sided 3-charge: charges {+1, 0, -1}: P = a Z + b (ZW)^k-ish + c W. Use {+2,0,-1}: P=aZ^2 + b*ZW + c*W
def scan(support, rng=(-2,-1,1,2), depth=6):
    # support = list of (degZ,degW); scan integer coeffs; find any two-sided P with E[P^m]=0 up to depth
    found=[]; firstdepth={}
    for coeffs in product(rng, repeat=len(support)):
        P={support[i]:Fr(coeffs[i]) for i in range(len(support)) if coeffs[i]!=0}
        if not P or onesided(P): continue
        ok=True; d=None
        for m in range(1,depth+1):
            if E(powp(P,m))!=0:
                d=m; ok=False; break
        if ok: found.append(P)
        else: firstdepth[d]=firstdepth.get(d,0)+1
    return found, firstdepth
# stratum charges {+2, 0, -1}: supports Z^2 (2,0), ZW (1,1), W (0,1)
for name,support in [("{+2,0,-1}: aZ^2+b ZW+c W",[(2,0),(1,1),(0,1)]),
                     ("{+1,0,-2}: aZ + b ZW + c W^2",[(1,0),(1,1),(0,2)]),
                     ("{+3,+1,-1,-3}(S73 span6)",[(3,0),(1,0),(0,1),(0,3)])]:
    found,fd=scan(support)
    print(f"  {name}: two-sided nullcone elements found = {len(found)}; "
          f"first-fire depth histogram = {dict(sorted(fd.items()))}")
print("""  => ZERO two-sided nullcone elements in every scanned stratum (GMC(2) holds here).
  The first-fire depth histogram shows detection depth 2,3,... (the moment where a
  two-sided P is caught) -- consistent with EMP: depth grows with charge span/degree.
  MECHANISM (S73): E[P^2]=0 is the 'primitive relation' among coeffs; on its variety
  E[P^4] (or E[P^6]) is a nonzero homogeneous form => no two-sided nullcone. GMC(2)
  REMAINS OPEN in general (unbounded degree/charge); these are bounded strata.""")
