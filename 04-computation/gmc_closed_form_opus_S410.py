"""opus-2026-07-20-S410 -- THE GMC(4) COUNTEREXAMPLE: CLOSED FORM, MINIMALITY, AND THE n=3 GAP.

Owner supplied an outside counterexample to the Gaussian Moment Conjecture at n >= 4:
    P = (1+Z2)(W1(1-Z1) + W2),  Q = Z2,   with W_j = conj(Z_j)
    (2 complex Gaussians = 4 REAL Gaussians; P is cubic with 6 terms)
    E[P^m] = 0 for all m >= 1,  but  E[Q P^m] = m! != 0.

REPO STATUS (checked): canon records Zhao's Image Conjecture / Mathieu subspaces as
"COROLLARY-false, no witness" -- false as consequences of the JC counterexample (Alpoge),
but with NO explicit witness. A VC-witness transport was attempted and only PARTLY executed
(dimension ~76 -> ~20 -> contingent, MISTAKE-201). So we had the CONCLUSION but not the
OBJECT, and this witness is n=4 / degree 3 / 6 terms and UNCONDITIONAL -- it needs no JC
input at all. That is stronger than what we had in the ways that matter: smaller, explicit,
and independent of the deep theorem.

WHAT THIS SCRIPT ADDS.
 (A) A CLOSED-FORM PROOF, replacing moment-by-moment verification:
        E[exp(tP)] = 1        exactly, and
        E[Z2 exp(tP)] = t/(1-t) = sum_{m>=1} t^m,  i.e. E[Z2 P^m] = m! for ALL m.
     Derivation (verified below): integrate Z1 first. With c = 1+Z2, d = conj(Z2),
        P = c*conj(Z1)*(1-Z1) + c*d
        E_{Z1}[exp(tP)] = exp(t c d) / (1 + t c)          <- a RESOLVENT appears
     then integrate Z2 using  int w^k e^{a wbar} e^{-L|w|^2} d2w/pi = (1/L)(a/L)^k :
        E[exp(tP)] = (1/(1-t^2)) * sum_k [-t^2/(1-t^2)]^k = (1/(1-t^2))*(1-t^2) = 1
        E[Z2 exp(tP)] = (1/(1-t^2))*(t/(1-t))*(1-t^2) = t/(1-t)
     THE MECHANISM IS A TWO-STAGE CASCADE: Z1 manufactures a pole 1/(1+tc); the Z2 integral
     resums the geometric series to exactly 1 for E[e^{tP}], but the extra factor of w
     shifts the series by one term and leaves t/(1-t), whose pole at t=1 is why NO moment
     vanishes. You need TWO complex Gaussians to run this. That is a structural reason the
     boundary sits between n=2 and n=4.
 (B) MINIMALITY probes: is 6 terms / degree 3 needed?
 (C) The n = 2 case (GMC(2) is a THEOREM, Derksen-van den Essen-Zhao): search for any
     counterexample in ONE complex Gaussian -- expect none, as a consistency check.
"""
import sympy as sp
from math import factorial
from itertools import combinations, product

z1, z1b, z2, z2b, t = sp.symbols('z1 z1b z2 z2b t')

def E(expr, vars2=((z1, z1b), (z2, z2b))):
    """exact expectation, Z_j iid standard complex: E[z^a zb^b] = a! delta_ab"""
    gens = [g for pair in vars2 for g in pair]
    e = sp.expand(expr)
    if e == 0: return sp.Integer(0)
    p = sp.Poly(e, *gens)
    tot = 0
    for mono, c in zip(p.monoms(), p.coeffs()):
        ok = True; w = 1
        for i in range(0, len(gens), 2):
            a, b = mono[i], mono[i+1]
            if a != b: ok = False; break
            w *= factorial(a)
        if ok: tot += c*w
    return sp.expand(tot)

TERMS = [-z1*z1b*z2, -z1*z1b, z1b*z2, z1b, z2*z2b, z2b]
P = sum(TERMS)

print("="*78)
print("(A) CLOSED FORM -- check E[P^m]=0 and E[Z2 P^m]=m! far beyond the original range")
print("="*78)
allok = True
for m in range(1, 13):
    Pm = sp.expand(P**m)
    e1 = E(Pm); e2 = E(sp.expand(z2*Pm))
    good = (e1 == 0 and e2 == factorial(m))
    allok &= good
    print(f"   m={m:2d}   E[P^m]={e1}   E[Z2 P^m]={e2}   m!={factorial(m)}   {'OK' if good else 'FAIL'}")
print(f"   => {'closed form confirmed through m=12' if allok else 'DISCREPANCY'}")

print()
print("   intermediate step of the derivation, checked symbolically:")
print("   E_{Z1}[exp(tP)] should equal exp(t*c*d)/(1+t*c) with c=1+Z2, d=conj(Z2).")
c, d = 1 + z2, z2b
lhs_series = 0
for m in range(0, 9):
    lhs_series += t**m/factorial(m)*sp.expand(P**m)
# integrate out Z1 only
def E1(expr):
    e = sp.expand(expr)
    if e == 0: return sp.Integer(0)
    p = sp.Poly(e, z1, z1b)
    tot = 0
    for (a, b), co in zip(p.monoms(), p.coeffs()):
        if a == b: tot += co*factorial(a)
    return sp.expand(tot)
lhs = sp.expand(E1(lhs_series))
rhs = sp.series(sp.exp(t*c*d)/(1 + t*c), t, 0, 9).removeO()
diff = sp.simplify(sp.expand(lhs - sp.expand(rhs)))
print(f"   difference of the two t-series through t^8: {diff}")
print(f"   => {'RESOLVENT FORM CONFIRMED' if diff == 0 else 'mismatch'}")

print()
print("="*78)
print("(B) MINIMALITY -- can any PROPER SUBSET of the 6 terms still work?")
print("="*78)
found = []
for r in range(1, 6):
    for sub in combinations(range(6), r):
        Ps = sum(TERMS[i] for i in sub)
        if Ps == 0: continue
        if all(E(sp.expand(Ps**m)) == 0 for m in range(1, 6)):
            nz = [m for m in range(1, 6) if E(sp.expand(z2*Ps**m)) != 0]
            if nz: found.append((r, sub, nz))
print(f"   proper subsets that still give a counterexample: {len(found)}")
for r, sub, nz in found[:6]:
    print(f"      {r} terms {sub}: E[Q P^m] nonzero at m={nz}")
if not found:
    print("   NONE -- all 6 terms are needed on this support. The example is term-minimal here.")

print()
print("="*78)
print("(C) n = 2 (ONE complex Gaussian): GMC(2) is a THEOREM (Derksen-van den Essen-Zhao)")
print("="*78)
print("   consistency check: brute-force search for ANY counterexample in one complex")
print("   Gaussian, coefficients in {-1,0,1}, degree <= 3.")
mon = [z1, z1b, z1**2, z1*z1b, z1b**2, z1**3, z1**2*z1b, z1*z1b**2, z1b**3]
def E_one(expr):
    e = sp.expand(expr)
    if e == 0: return sp.Integer(0)
    p = sp.Poly(e, z1, z1b); tot = 0
    for (a, b), co in zip(p.monoms(), p.coeffs()):
        if a == b: tot += co*factorial(a)
    return sp.expand(tot)
hits = 0; checked = 0
for coeffs in product([-1, 0, 1], repeat=len(mon)):
    if not any(coeffs): continue
    checked += 1
    Ps = sum(cc*mm for cc, mm in zip(coeffs, mon))
    if any(E_one(sp.expand(Ps**m)) != 0 for m in (1, 2, 3, 4)): continue
    for Q in (z1, z1b, z1*z1b):
        if any(E_one(sp.expand(Q*Ps**m)) != 0 for m in (1, 2, 3, 4)):
            hits += 1
            if hits <= 3: print(f"      CANDIDATE: P={Ps}, Q={Q}")
            break
print(f"   searched {checked} polynomials; counterexample candidates found: {hits}")
print(f"   => {'consistent with GMC(2) TRUE' if hits == 0 else 'CANDIDATES -- inspect (likely higher m needed)'}")
