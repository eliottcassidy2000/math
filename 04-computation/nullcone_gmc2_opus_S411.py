"""opus-2026-07-20-S411 -- THE 2-DIMENSIONAL NULLCONE, AND A PROOF ROUTE FOR GMC(2).

REFRAMING (from THM-1495): the nullcone condition E[P^m] = 0 for all m >= 1 is exactly
    E[exp(tP)] = 1   (as a formal series in t),
and GMC(n) -- "E[QP^m] = 0 for m >> 0" -- is exactly
    E[Q exp(tP)] is a POLYNOMIAL in t.
The n=4 witness produced E[Q e^{tP}] = t/(1-t): a POLE, hence never a polynomial.
So GMC(n) <=> no P in the nullcone lets E[Q e^{tP}] acquire a pole.

THE CHARGE (U(1)) GRADING -- the mechanism I want to test as the whole story at n=2.
For one complex Gaussian, E[z^a zbar^b] = a! delta_{ab}, so E KILLS every monomial with
    charge(z^a zbar^b) := a - b  != 0.
Charge is ADDITIVE under multiplication. Hence:

  LEMMA (charge lemma, any n).  If every monomial of P has charge >= 1, then every
  monomial of P^m has charge >= m, so E[P^m] = 0 for all m >= 1 -- P is in the NULLCONE.
  Moreover if Q has all charges >= -C, then every monomial of Q P^m has charge >= m - C,
  so E[Q P^m] = 0 for every m > C.  **GMC HOLDS FOR SUCH P, IN ANY DIMENSION.**
  (Same with all charges <= -1, by conjugation.)

So charge-definite nullcone elements can NEVER be counterexamples. Consistency check on the
n=4 witness P' = (1+Z2)(W2 - Z1W1): its monomials have charges -1, 0, 0, +1 -- MIXED, and
including 0. That is exactly how it escapes the lemma.

THE QUESTION THIS SCRIPT ASKS:
  at n = 2 (one complex Gaussian), is EVERY nullcone element charge-definite?
  If yes, the charge lemma PROVES GMC(2) outright.
"""
import sympy as sp
from math import factorial
from itertools import product
from collections import defaultdict

z, zb = sp.symbols('z zb')

MONS = [(a, b) for a in range(4) for b in range(4) if a + b <= 3]
def mono(a, b): return z**a * zb**b

def E1(expr):
    """E[z^a zbar^b] = a! delta_ab, one standard complex Gaussian"""
    e = sp.expand(expr)
    if e == 0: return sp.Integer(0)
    p = sp.Poly(e, z, zb); tot = 0
    for (a, b), c in zip(p.monoms(), p.coeffs()):
        if a == b: tot += c*factorial(a)
    return sp.expand(tot)

def charges(expr):
    e = sp.expand(expr)
    if e == 0: return set()
    p = sp.Poly(e, z, zb)
    return {a - b for (a, b), c in zip(p.monoms(), p.coeffs()) if c != 0}

print("="*78)
print("(1) THE CHARGE LEMMA, stated and sanity-checked")
print("="*78)
print("   charge(z^a zbar^b) = a-b is ADDITIVE; E kills every nonzero charge.")
print("   => all charges >= 1  ==>  P in nullcone AND E[QP^m]=0 for m > C   (any n)")
demo = z + z**2 - z**3 + z**2*zb          # charges 1,2,3,1 -- all >= 1
print(f"   demo P = {demo}   charges {sorted(charges(demo))}")
print(f"   E[P^m], m=1..8 : {[E1(sp.expand(demo**m)) for m in range(1,9)]}")
for Q in (zb, zb**2, z*zb):
    vals = [E1(sp.expand(Q*demo**m)) for m in range(1, 9)]
    print(f"   Q={Q}: E[QP^m] m=1..8 = {vals}   (dies once m > charge deficit)")

print()
print("="*78)
print("(2) EXHAUSTIVE NULLCONE AT n=2, degree <= 3, coefficients in {-1,0,1}")
print("="*78)
null = []
checked = 0
for co in product([-1, 0, 1], repeat=len(MONS)):
    if not any(co): continue
    checked += 1
    P = sum(c*mono(a, b) for c, (a, b) in zip(co, MONS) if c)
    if all(E1(sp.expand(P**m)) == 0 for m in range(1, 8)):
        null.append(P)
print(f"   polynomials scanned : {checked}")
print(f"   nullcone members    : {len(null)}")

defn, mixed = [], []
for P in null:
    ch = charges(P)
    if all(c >= 1 for c in ch) or all(c <= -1 for c in ch): defn.append((P, ch))
    else: mixed.append((P, ch))
print(f"   charge-DEFINITE     : {len(defn)}")
print(f"   NOT charge-definite : {len(mixed)}")
if mixed:
    print("\n   *** NOT charge-definite -- these decide the question ***")
    seen = set()
    for P, ch in mixed[:25]:
        key = str(sp.expand(P))
        if key in seen: continue
        seen.add(key)
        print(f"      P = {P}      charges {sorted(ch)}")
else:
    print("\n   EVERY nullcone element is charge-definite at this degree.")
    print("   => the charge lemma would PROVE GMC(2) on this range.")

print()
print("="*78)
print("(3) IF MIXED ONES EXIST: do they still satisfy GMC? (E[QP^m] eventually 0?)")
print("="*78)
bad = []
for P, ch in mixed:
    for Q in (z, zb, z*zb, z**2, zb**2):
        vals = [E1(sp.expand(Q*P**m)) for m in range(1, 10)]
        if vals[-1] != 0 and vals[-2] != 0:
            bad.append((P, Q, vals)); break
print(f"   mixed-charge nullcone elements tested : {len(mixed)}")
print(f"   with E[QP^m] still nonzero at m=8,9   : {len(bad)}")
if bad:
    for P, Q, v in bad[:5]:
        print(f"      *** P={P}, Q={Q}, E[QP^m]={v}  <-- would REFUTE GMC(2)")
else:
    print("   NONE -- every mixed-charge nullcone element still obeys GMC(2).")
    print("   Consistent with the theorem; the charge lemma is not the whole story,")
    print("   but nothing escapes GMC(2) on this range.")
